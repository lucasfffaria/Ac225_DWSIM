# Native IronPython CustomUO body. liquid_properties.py and PARAMS precede this code.
from System import Array,Double
from math import exp,log

_basis=3600.0
_keys=['A225_Bq','A226_Bq','A227_Bq','Ra225_Bq']
_lam=[log(2)/9.92,log(2)/1.224,log(2)/(21.7*365.25),log(2)/14.9]
_p=Me.InputVariables
_result={}
def get_mass(s):
    rate=s.GetMassFlow()
    return [float(s.Phases[0].Compounds[n].MassFraction)*rate*_basis for n in COMPONENTS]
def get_iso(s):
    return [float(s.GetExtraPropertyValue(k) or 0.0) for k in _keys]
def get_age(s):return float(s.GetExtraPropertyValue('age_d') or 0.0)
def put(s,m,A,age):
    s.SetTemperature(298.15);s.SetPressure(101325.0)
    s.SetOverallMassComposition(Array[Double]([v/sum(m) for v in m]))
    s.SetMassFlow(sum(m)/_basis)
    for k,v in zip(_keys,A):s.SetExtraPropertyValue(k,float(v))
    s.SetExtraPropertyValue('age_d',float(age))
    s.SetExtraPropertyValue('batch_basis_s',_basis)
    s.SetExtraPropertyValue('batch_volume_mL',liquid_volume(m))
def decay(A,minutes):
    td=minutes/1440.
    out=[A[i]*exp(-_lam[i]*td) for i in range(4)]
    out[0]+=A[3]*_lam[0]/(_lam[0]-_lam[3])*(exp(-_lam[3]*td)-exp(-_lam[0]*td))
    return out
def close(inmass,outmass,Ain,Aout):
    _result['total_mass_residual_kg']=sum(inmass)-sum(outmass)
    _result['heavy_descendant_sink_atoms']=sum(Ain[i]/(_lam[i]/86400.)-Aout[i]/(_lam[i]/86400.) for i in range(4))

_kind=int(_p['kind'])
if _kind==0: # source boundary: primary fraction is an input to this subflowsheet
    _age=_p['cooling_d']+5./24
    _A=[1e6,1e6*_p['r226'],1e6*_p['r227'],0.]
    _A=decay(_A,_age*1440.)
    _A=[_A[i]*_p['eta_primary'] if i<3 else 0. for i in range(4)]
    put(oms1,solution(80.,{'Nitric acid':10.}),_A,_age)
    _result['primary_duration_h']=5.
elif _kind==1: # isothermal conditioning mixer
    _m1=get_mass(ims1);_m2=get_mass(ims2);_m=mix_mass(_m1,_m2)
    _A0=get_iso(ims1);_dt=_p['duration_min'];_A=decay(_A0,_dt)
    put(oms1,_m,_A,get_age(ims1)+_dt/1440.)
    _result['mixed_volume_mL']=liquid_volume(_m)
    _result['analytical_HNO3_M']=_m[1]/(MW[1]/1000)/(liquid_volume(_m)/1000.)
    close(mix_mass(_m1,_m2),_m,_A0,_A)
elif _kind in [2,4]: # sequential batch-column collection/raffinate ledger
    _ss=[ims1,ims2,ims3,ims4];_mm=[get_mass(s) for s in _ss]
    _mi=mix_mass(*_mm)
    _pumped=sum(liquid_volume(m) for m in _mm)
    _dt=_pumped/_p['flow_mL_min']
    _cut=solution(10.,{'Hydrogen chloride':_p['carry_HCl_M']}) if _kind==2 else solution(10.,{'Nitric acid':.7})
    _waste=subtract_mass(_mi,_cut)
    _A0=get_iso(ims1);_aged=decay(_A0,_dt)
    _prod=[v*_p['eta_Ac'] if i<3 else 0. for i,v in enumerate(_aged)]
    _lost=[_aged[i]-_prod[i] for i in range(4)]
    _age=get_age(ims1)+_dt/1440.
    put(oms1,_cut,_prod,_age);put(oms2,_waste,_lost,_age)
    _result['pumped_volume_mL']=_pumped
    _result['max_pressure_drop_Pa']=max(ergun(m,_p['flow_mL_min'],_p['bed_mL'],_p['diameter_mm'],_p['particle_um'],_p['voidage']) for m in _mm)
    _result['max_pressure_drop_upper_Pa']=max(ergun(m,_p['flow_mL_min'],_p['bed_mL'],_p['diameter_mm'],50.,.35) for m in _mm)
    _result['component_mass_max_residual_kg']=max(abs(_mi[i]-_cut[i]-_waste[i]) for i in range(5))
    close(_mi,mix_mass(_cut,_waste),_A0,_aged)
elif _kind==3: # specified NH3 addition, stoichiometric analytical neutralization
    _f=get_mass(ims1);_b=get_mass(ims2);_m=mix_mass(_f,_b)
    _n=min(_m[2]/(MW[2]/1000),_m[3]/(MW[3]/1000))
    _m[2]-=_n*MW[2]/1000;_m[3]-=_n*MW[3]/1000;_m[4]+=_n*(MW[2]+MW[3])/1000
    _dt=_p['duration_min'];_A0=get_iso(ims1);_A=decay(_A0,_dt)
    put(oms1,_m,_A,get_age(ims1)+_dt/1440.)
    _result['analytical_free_HCl_M']=_m[2]/(MW[2]/1000)/(liquid_volume(_m)/1000)
    _result['NH3_consumed_mol']=_n
    close(mix_mass(_f,_b),_m,_A0,_A)
elif _kind==5: # explicit protocol allowance and additional queue
    _m=get_mass(ims1);_A0=get_iso(ims1);_dt=_p['protocol_allowance_min']+60*_p['extra_queue_h']
    _A=decay(_A0,_dt);put(oms1,_m,_A,get_age(ims1)+_dt/1440.)
    _result['activity_concentration_Bq_mL']=_A[0]/liquid_volume(_m)
    close(_m,_m,_A0,_A)
elif _kind==6: # separate literature-matrix evidence, not a fitted Tan column
    _result['Ac_recovery']=ph_recovery(_p['load_pH'])
    _result['worst_recovery_in_control_interval']=min(ph_recovery(_p['load_pH']-_p['pH_halfwidth']),ph_recovery(_p['load_pH']+_p['pH_halfwidth']))

if _kind not in [0,6]:_result['duration_min']=_dt
for _name,_value in _result.items():globals()[str(_name)]=float(_value)
