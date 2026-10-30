"""Execute native DWSIM custom operations, save, reload, and compare to CPython.

This model intentionally has no surrogate radionuclide thermodynamics.
It is a reusable DWSIM inventory workspace, not a column-retention simulation.
"""
import sys,os,json,csv,hashlib,platform
from pathlib import Path
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[1]
BIN=Path(os.environ.get('DWSIM_BIN',Path(os.environ['LOCALAPPDATA'])/'DWSIM'))
if not (BIN/'DWSIM.Automation.dll').exists(): raise FileNotFoundError(BIN)
sys.path.insert(0,str(BIN))
os.chdir(BIN)
import clr
from System import AppDomain,ResolveEventHandler
from System.Reflection import Assembly
def resolve(sender,args):
    for folder in [BIN,BIN/'unitops',BIN/'ppacks']:
        p=folder/(str(args.Name).split(',')[0]+'.dll')
        if p.exists(): return Assembly.LoadFrom(str(p))
    return None
AppDomain.CurrentDomain.AssemblyResolve+=ResolveEventHandler(resolve)
for name in ['DWSIM.Interfaces','DWSIM.GlobalSettings','DWSIM.Automation']:
    assembly=clr.AddReference(str(BIN/(name+'.dll')))
from DWSIM.Automation import Automation3
from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType
from model import Decay, released_inventory, purity, best_age, total_window, run_recipe

cfg=json.loads((ROOT/'data/source_inputs.json').read_text())
recipe=list(csv.DictReader((ROOT/'data/reagent_recipe.csv').open(encoding='utf8')))
auto=Automation3();fs=auto.CreateFlowsheet().__implementation__
fs.Options.SimulationName='Actinium 225 isotope timing and documented reagent inventory'
kernel=(ROOT/'code/dwsim_inventory_unit.py').read_text()
base={'half225_d':9.92,'half226_d':1.224,'half227_d':21.7*365.25,'halfRa225_d':14.9,
 'A225_reference_Bq':1e6,'delivery_age_d':14.,'parent_removal_age_d':10.,'eta_Ac':.636,
 'parent_ratio':0.,'R226_limit':.01,'R227_limit':.01,'P_Ac_limit':.99,
 'delay_min_d':3.,'delay_max_d':5.,'eta_lower':.6042}
units=[]
for i,c in enumerate(cfg['isotope_cases']):
    for j,parent in enumerate([0.,.6]):
        tag=c['label']+(' - parent free' if parent==0 else ' - Ra225 scenario')
        u=fs.AddObject(ObjectType.CustomUO,80+350*j,80+180*i,tag).__implementation__
        u.ScriptText=kernel
        vals={**base,'r226':c['r226'],'r227':c['r227'],'parent_ratio':parent}
        for k,v in vals.items():u.InputVariables[k]=v
        u.Calculate(None)
        units.append((u,vals))
c=cfg['isotope_cases'][0]
u=fs.AddObject(ObjectType.CustomUO,780,80,'Proton robust upper input corner').__implementation__
u.ScriptText=kernel
vals={**base,'r226':c['r226']+2*c['r226_sigma'],'r227':c['r227']+2*c['r227_sigma']}
for k,v in vals.items():u.InputVariables[k]=v
u.Calculate(None);units.append((u,vals))
u=fs.AddObject(ObjectType.CustomUO,780,260,'Documented liquid and time ledger').__implementation__
u.ScriptText=(ROOT/'code/dwsim_reagent_unit.py').read_text()
for k,v in {'recipe_rows':len(recipe),'internal_transfer_mL':90,'flow_mL_min':3,'reported_process_min':540,'final_collection_mL':10}.items():u.InputVariables[k]=float(v)
for i,row in enumerate(recipe):
    for key in ['mL','HNO3_M','HCl_M']:
        u.InputVariables[f'row_{i}_{key}']=float(row['volume_mL'] if key=='mL' else row[key])
u.Calculate(None)
units.append((u,None))

# Run the native flowsheet solver, not merely a hand-written file export.
errors=auto.CalculateFlowsheet4(fs)
if errors and errors.Count: raise RuntimeError('\n'.join(str(x) for x in errors))
modeldir=ROOT/'dwsim';modeldir.mkdir(exist_ok=True)
dest=modeldir/'Ac225_Isotope_Timing.dwxml'
auto.SaveFlowsheet(fs,str(dest),False)
output={obj.GraphicObject.Tag:{str(k):float(obj.OutputVariables[k]) for k in obj.OutputVariables.Keys} for obj,_ in units}
checks=[]
for obj,vals in units:
    got=output[obj.GraphicObject.Tag]
    if vals:
        expected=released_inventory(vals['delivery_age_d'],vals['parent_removal_age_d'],vals['r226'],vals['r227'],vals['parent_ratio'],vals['eta_Ac'])*1e6
        for key,want in zip(['A225_delivery_Bq','A226_delivery_Bq','A227_delivery_Bq'],expected):
            err=abs(got[key]-float(want))/max(abs(float(want)),1e-14)
            assert err<1e-10
            checks.append({'unit':obj.GraphicObject.Tag,'field':key,'relative_error':err})
        win=total_window(vals['r226'],vals['r227'])
        if win:
            for key,want in zip(['parent_free_P_Ac_lower_age_d','parent_free_P_Ac_upper_age_d'],win):
                assert abs(got[key]-want)<1e-9
    else:
        want=run_recipe(recipe)
        assert abs(got['analytical_HNO3_mol']-want['HNO3_mol'])<1e-12
        assert abs(got['documented_pumped_volume_mL']-want['known_pumped_mL'])<1e-9

# Persistence check: a saved .dwxml must contain executable scripts and inputs.
loaded=auto.LoadFlowsheet(str(dest))
errs=auto.CalculateFlowsheet4(loaded)
if errs and errs.Count:raise RuntimeError('\n'.join(str(x) for x in errs))
for obj in loaded.SimulationObjects.Values:
    unit=obj.__implementation__
    for k in unit.OutputVariables.Keys:
        val=float(unit.OutputVariables[k]);ref=output[unit.GraphicObject.Tag][str(k)]
        assert abs(val-ref)<=1e-12*max(1,abs(ref))
report={'dwsim_version':str(assembly.GetName().Version),'runtime':'Automation3 and native IronPython CustomUO',
        'objects':len(units),'flowsheet_solver_errors':0,'saved_reload_verified':True,
        'max_CPython_relative_error':max(x['relative_error'] for x in checks),
        'dwxml_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'checks':checks,
        'thermodynamic_property_package':'none; inventory and analytical liquid-accounting model only',
        'outputs':output}
(ROOT/'results/dwsim_verification.json').write_text(json.dumps(report,indent=2))
print(json.dumps({k:v for k,v in report.items() if k not in ['outputs','checks']},indent=2))
