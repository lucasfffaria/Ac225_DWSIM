"""Build, solve, reload and sweep the connected DWSIM purification subflowsheet."""
import os,sys,json,hashlib,time,argparse
from pathlib import Path
import numpy as np
import pandas as pd
from process_model import ROOT,BASE,reference,ISO,LAM,lp
sys.stdout.reconfigure(encoding='utf8')
parser=argparse.ArgumentParser();parser.add_argument('--quick',action='store_true');args=parser.parse_args()
BIN=Path(os.environ.get('DWSIM_BIN',Path(os.environ['LOCALAPPDATA'])/'DWSIM'))
sys.path.insert(0,str(BIN));os.chdir(BIN)
import clr
from System import AppDomain,ResolveEventHandler,Array,Double
from System.Reflection import Assembly
def resolve(sender,args):
    for folder in [BIN,BIN/'unitops',BIN/'ppacks']:
        p=folder/(str(args.Name).split(',')[0]+'.dll')
        if p.exists():return Assembly.LoadFrom(str(p))
AppDomain.CurrentDomain.AssemblyResolve+=ResolveEventHandler(resolve)
for n in ['DWSIM.Interfaces','DWSIM.GlobalSettings','DWSIM.Automation']:assembly=clr.AddReference(str(BIN/(n+'.dll')))
from DWSIM.Automation import Automation3
from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType
from DWSIM.Interfaces.Enums import ForcedPhase
auto=Automation3();fs=auto.CreateFlowsheet().__implementation__
fs.Options.SimulationName='Ac225 connected DGA-LN purification: batch material balances, hydraulics and isotope clock'
for n in lp.COMPONENTS:fs.AddCompound(n)
pp=fs.CreateAndAddPropertyPackage("Raoult's Law").__implementation__
pp.Tag='Isothermal liquid carrier; empirical Laliberte density override'
prop_code='PARAMS = '+repr(lp.PARAMS)+'\n'+(ROOT/'code/liquid_properties.py').read_text(encoding='utf8')
for phase in ['Mixture','OverallLiquid','Liquid1']:
    pp.PropertyOverrides[phase+'/density']=prop_code+"\npropval=liquid_density([float(matstr.Phases[0].Compounds[n].MassFraction) for n in COMPONENTS],float(matstr.GetTemperature()))\n"
kernel=prop_code+'\n'+(ROOT/'code/dwsim_process_unit.py').read_text()
units={};streams={};edges=[]
eta=(.636/.70)**.5
def unit(tag,kind,x,y,**values):
    u=fs.AddObject(ObjectType.CustomUO,x,y,tag).__implementation__
    u.ScriptText=kernel
    for k,v in {'kind':kind,**values}.items():u.InputVariables[k]=float(v)
    u.GraphicObject.Width=90;u.GraphicObject.Height=65
    units[tag]=u;return u
def stream(tag,x,y,m=None):
    s=fs.AddObject(ObjectType.MaterialStream,int(x),int(y),tag).__implementation__
    s._ppid=pp.UniqueID;s.ForcePhase=ForcedPhase.Liquid
    s.OverrideSingleCompoundFlashBehavior=True
    s.SetTemperature(298.15);s.SetPressure(101325)
    if m is None:m=lp.solution(10)
    s.SetMassFlow(sum(m)/3600)
    s.SetOverallMassComposition(Array[Double]([a/sum(m) for a in m]))
    for k in ISO:s.AddExtraProperty(k,0.)
    s.AddExtraProperty('age_d',0.)
    s.AddExtraProperty('batch_basis_s',3600.)
    s.AddExtraProperty('batch_volume_mL',lp.liquid_volume(m))
    streams[tag]=s;return s
def connect(a,b,outidx=0,inidx=0):
    try:fs.ConnectObjects(a.GraphicObject,b.GraphicObject,outidx,inidx)
    except Exception:
        print('Failed connection',a.GraphicObject.Tag,b.GraphicObject.Tag,outidx,inidx)
        print('Inputs',[(i,c.ConnectorName,str(c.Type),c.IsAttached) for i,c in enumerate(b.GraphicObject.InputConnectors)])
        raise
    edges.append({'from':a.GraphicObject.Tag,'to':b.GraphicObject.Tag,'outlet':outidx,'inlet':inidx})
u0=unit('P-100 PRIMARY BOUNDARY',0,30,330,cooling_d=12.6472106068,r226=10.8,r227=.0014,eta_primary=.70)
u1=unit('M-101 DILUTION',1,300,330,duration_min=10)
u2=unit('C-201 DGA2',2,650,330,flow_mL_min=3,eta_Ac=eta,bed_mL=5,diameter_mm=11,particle_um=100,voidage=.4,carry_HCl_M=.1)
u3=unit('R-201 ACID ADJUSTMENT',3,1000,330,duration_min=10)
u4=unit('C-202 LN',4,1350,330,flow_mL_min=3,eta_Ac=eta,bed_mL=2,diameter_mm=8,particle_um=100,voidage=.4)
u5=unit('T-301 PROTOCOL AND QUEUE',5,1680,330,protocol_allowance_min=BASE['protocol_allowance_min'],extra_queue_h=0)
up=unit('EVIDENCE: CITRATE AG50 pH',6,1050,650,load_pH=2-.828571428571*.2,pH_halfwidth=.2)
path=[('PRIMARY-FRACTION',u0,u1),('DILUTED-FEED',u1,u2),('DGA-CUT',u2,u3),('LN-FEED',u3,u4),('LN-CUT',u4,u5)]
for tag,ua,ub in path:
    s=stream(tag,(ua.GraphicObject.X+ub.GraphicObject.X)/2+20,345)
    connect(ua,s);connect(s,ub)
final=stream('FINAL-PRODUCT',1910,345);connect(u5,final)
for tag,u,ix,x,y in [('DI-WATER',u1,1,210,150),('DGA-PRE',u2,1,440,65),('DGA-RINSE',u2,2,465,140),
 ('DGA-ELUENT',u2,4,490,215),('NH3-ADDITION',u3,1,900,150),('LN-PRE',u4,1,1140,65),
 ('LN-RINSE',u4,2,1165,140),('LN-ELUENT',u4,4,1190,215)]:
    s=stream(tag,x,y,BASE['streams'][tag]);connect(s,u,0,ix)
for tag,u,x in [('DGA-WASTE',u2,760),('LN-WASTE',u4,1460)]:
    s=stream(tag,x,525);connect(u,s,1,0)

def solve():
    errs=auto.CalculateFlowsheet4(fs)
    if errs and errs.Count:raise RuntimeError('\n'.join(str(e) for e in errs))
def extract(s):
    m=[float(s.Phases[0].Compounds[n].MassFraction)*s.GetMassFlow()*3600 for n in lp.COMPONENTS]
    return {'stream':s.GraphicObject.Tag,'mass_kg_batch':sum(m),
      **{n+'_kg':m[i] for i,n in enumerate(lp.COMPONENTS)},
      'volume_mL_batch':s.GetVolumetricFlow()*3600*1e6,'density_kg_m3':float(s.Phases[0].Properties.density),
      'age_d':float(s.GetExtraPropertyValue('age_d')),**{k:float(s.GetExtraPropertyValue(k)) for k in ISO}}
solve()
rows=[extract(s) for s in streams.values()]
outdir=ROOT/'results';outdir.mkdir(exist_ok=True)
pd.DataFrame(rows).to_csv(outdir/'connected_stream_table.csv',index=False)
outputs={tag:{str(k):float(u.OutputVariables[k]) for k in u.OutputVariables.Keys} for tag,u in units.items()}
(outdir/'connected_unit_outputs.json').write_text(json.dumps(outputs,indent=2))
(outdir/'connected_edges.json').write_text(json.dumps(edges,indent=2))
ref=reference();checks=[]
for row in rows:
    key=row['stream'];expected=ref['recipe']['streams'].get(key)
    if expected is not None:
        for i,n in enumerate(lp.COMPONENTS):
            assert abs(row[n+'_kg']-expected[i])<1e-10,(key,n,row[n+'_kg'],expected[i])
        assert abs(row['volume_mL_batch']-lp.liquid_volume(expected))<1e-6,(key,row['volume_mL_batch'])
for tag,label in [('PRIMARY-FRACTION','PRIMARY-BOUNDARY'),('DILUTED-FEED','DILUTION'),('DGA-CUT','DGA2'),
 ('LN-FEED','CONDITIONING'),('LN-CUT','LN'),('FINAL-PRODUCT','PROTOCOL-HOLD')]:
    native=extract(streams[tag]);oracle=next(x for x in ref['stage_inventory'] if x['stage']==label)
    for k in ISO:
        err=abs(native[k]-oracle[k])/max(abs(oracle[k]),1)
        assert err<1e-10,(tag,k,native[k],oracle[k]);checks.append(err)
    assert abs(native['age_d']-oracle['age_d'])<1e-10
# Verify balances at every reactive/non-reactive operation and the global subflowsheet boundary.
assert max(abs(o.get('total_mass_residual_kg',0)) for o in outputs.values())<1e-12
assert min(o.get('heavy_descendant_sink_atoms',0) for o in outputs.values())>=-1e-3
feed_names=['PRIMARY-FRACTION']+BASE['fresh_tags'];exit_names=['DGA-WASTE','LN-WASTE','FINAL-PRODUCT']
index={r['stream']:r for r in rows}
mass_residual=sum(index[t]['mass_kg_batch'] for t in feed_names)-sum(index[t]['mass_kg_batch'] for t in exit_names)
assert abs(mass_residual)<1e-12
dest=ROOT/'dwsim/Ac225_Connected_Purification.dwxml'
auto.SaveFlowsheet(fs,str(dest),False)
loaded=auto.LoadFlowsheet(str(dest))
errs=auto.CalculateFlowsheet4(loaded)
if errs and errs.Count:raise RuntimeError('\n'.join(str(e) for e in errs))
reload_max=0.
for obj in loaded.SimulationObjects.Values:
    obj=obj.__implementation__
    if obj.GraphicObject.Tag in streams:
        row=extract(obj);old=index[row['stream']]
        for k in ['mass_kg_batch','volume_mL_batch','density_kg_m3','age_d']+ISO:
            err=abs(row[k]-old[k])/max(abs(old[k]),1);reload_max=max(reload_max,err)
assert reload_max<1e-10
print('Baseline solved and reloaded:',len(units),'custom operations,',len(streams),'material streams,',len(edges),'connections',flush=True)
if not args.quick:
    sweep=[]
    for q in [1.5,3.,4.5,6.,9.]:
        for h in [0.,6.,12.,24.,48.]:
            u2.InputVariables['flow_mL_min']=q;u4.InputVariables['flow_mL_min']=q;u5.InputVariables['extra_queue_h']=h
            solve();r=extract(final);oracle=reference(q,h)
            relative=abs(r['A225_Bq']-oracle['activity_Bq'][0])/oracle['activity_Bq'][0]
            assert relative<1e-10;checks.append(relative)
            for unit_tag,key in [('C-201 DGA2','DGA2'),('C-202 LN','LN')]:
                native_p=float(units[unit_tag].OutputVariables['max_pressure_drop_Pa'])
                assert abs(native_p-oracle['pressure_Pa'][key])<1e-6
            sweep.append({'flow_mL_min':q,'extra_queue_h':h,'process_h':(r['age_d']-float(u0.InputVariables['cooling_d']))*24,
              'A225_Bq':r['A225_Bq'],'A225_Bq_mL':r['A225_Bq']/r['volume_mL_batch'],
              'DGA2_pressure_Pa':float(u2.OutputVariables['max_pressure_drop_Pa']),
              'LN_pressure_Pa':float(u4.OutputVariables['max_pressure_drop_Pa']),
              'DGA2_pressure_upper_Pa':float(u2.OutputVariables['max_pressure_drop_upper_Pa']),
              'LN_pressure_upper_Pa':float(u4.OutputVariables['max_pressure_drop_upper_Pa'])})
        print('Native flow/queue sweep complete at',q,'mL/min',flush=True)
    pd.DataFrame(sweep).to_csv(outdir/'dwsim_process_sweep.csv',index=False)
    # Recalculate the measured-matrix unit for a control-setpoint sweep without re-running unrelated streams.
    phrows=[]
    up.InputVariables['pH_halfwidth']=0.
    for ph in np.linspace(1.5,2.5,101):
        up.InputVariables['load_pH']=float(ph);up.Calculate(None)
        eta_native=float(up.OutputVariables['Ac_recovery'])
        assert abs(eta_native-lp.ph_recovery(float(ph)))<1e-12
        phrows.append({'load_pH':ph,'recovery_fraction':eta_native})
    pd.DataFrame(phrows).to_csv(outdir/'dwsim_pH_sweep.csv',index=False)
report={'dwsim_version':str(assembly.GetName().Version),'custom_operations':len(units),'material_streams':len(streams),
 'connections':len(edges),'native_flowsheet_errors':0,'saved_reload_verified':True,'max_reload_relative_error':reload_max,
 'max_isotope_matrix_oracle_relative_error':max(checks),'global_mass_residual_kg':mass_residual,
 'batch_basis_s':3600,'actual_pumping_rate_is_separate':True,'flow_queue_scenarios':0 if args.quick else 25,
 'native_pH_scenarios':0 if args.quick else 101,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),
 'scope':'Isothermal liquid analytical components, embedded density correlation, empirical collection cuts, conditional ammonia addition, hydraulic screening; no VLE/enthalpy/retention prediction.'}
(outdir/'connected_dwsim_verification.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2),flush=True)
# DWSIM may keep CLR worker threads alive after the headless automation finishes.
# End only this Python/CLR worker process after all outputs have been written.
from System import Environment
Environment.Exit(0)
