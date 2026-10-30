"""Connected purification extension: documented recipe plus labelled scenarios."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import expm
import liquid_properties as lp
ROOT=Path(__file__).resolve().parents[1]
lp.PARAMS=json.loads((ROOT/'data/liquid_property_parameters.json').read_text())['parameters']
LAM=np.log(2)/np.array([9.92,1.224,21.7*365.25,14.9])
ISO=['A225_Bq','A226_Bq','A227_Bq','Ra225_Bq']

def isotope_step(A,minutes):
    # Independent matrix form; native custom operations use scalar Bateman equations.
    mat=np.diag(-LAM);mat[0,3]=LAM[0]
    return expm(mat*minutes/1440.)@np.asarray(A)

def recipe(carry_M=.10,dilution_mL=80.):
    f=lp.solution(80,{'Nitric acid':10})
    water=lp.solution(dilution_mL)
    mixed=lp.mix_mass(f,water)
    cut=lp.solution(10,{'Hydrogen chloride':carry_M})
    base,conditioned,c=lp.condition(cut)
    feeds={'PRIMARY-FRACTION':f,'DI-WATER':water,
      'DGA-PRE':lp.solution(25,{'Nitric acid':5}),
      'DGA-RINSE':lp.solution(15,{'Hydrogen chloride':6}),
      'DGA-ELUENT':lp.solution(10),
      'NH3-ADDITION':base,
      'LN-PRE':lp.solution(10,{'Hydrogen chloride':.05}),
      'LN-RINSE':lp.solution(14,{'Nitric acid':.05}),
      'LN-ELUENT':lp.solution(10,{'Nitric acid':.7})}
    dga_in=lp.mix_mass(mixed,*[feeds[k] for k in ['DGA-PRE','DGA-RINSE','DGA-ELUENT']])
    ln_in=lp.mix_mass(conditioned,*[feeds[k] for k in ['LN-PRE','LN-RINSE','LN-ELUENT']])
    prod=feeds['LN-ELUENT'][:]
    streams={**feeds,'DILUTED-FEED':mixed,'DGA-CUT':cut,'LN-FEED':conditioned,
      'DGA-WASTE':lp.subtract_mass(dga_in,cut),'LN-WASTE':lp.subtract_mass(ln_in,prod),'FINAL-PRODUCT':prod}
    dga_pumped=lp.liquid_volume(mixed)+50
    ln_pumped=lp.liquid_volume(conditioned)+34
    # The published 4 h is reconciled at 3 mL/min; this fixed allowance is not claimed to be idle time.
    allowance=240.-20.-(dga_pumped+ln_pumped)/3.
    return {'streams':streams,'fresh_tags':[k for k in feeds if k!='PRIMARY-FRACTION'],
      'dga_pumped_mL':dga_pumped,'ln_pumped_mL':ln_pumped,
      'protocol_allowance_min':allowance,'NH3_volume_mL':lp.liquid_volume(base),
      'conditioned_acid_M':c}

BASE=recipe()

def reference(flow=3.,queue_h=0.,cooling_d=12.6472106068,eta_dga=None,eta_ln=None,carry_M=.10,dilution_mL=80.,rho_parent=0.):
    R=recipe(carry_M,dilution_mL)
    if eta_dga is None:eta_dga=(.636/.70)**.5
    if eta_ln is None:eta_ln=(.636/.70)**.5
    source_A=np.array([1.,10.8,.0014,rho_parent])*1e6
    # Source boundary conditions apply the empirical primary transfer to all Ac isotopes only.
    A=isotope_step(source_A,(cooling_d*24+5)*60)
    A[:3]*=.70;A[3]=0. # main connected module starts after primary Ra removal
    t=(cooling_d*24+5)*60
    stage=[('DILUTION',10.,1.),('DGA2',R['dga_pumped_mL']/flow,eta_dga),
      ('CONDITIONING',10.,1.),('LN',R['ln_pumped_mL']/flow,eta_ln),
      ('PROTOCOL-HOLD',BASE['protocol_allowance_min']+queue_h*60,1.)]
    rows=[{'stage':'PRIMARY-BOUNDARY','age_d':t/1440.,**dict(zip(ISO,A))}]
    for name,dur,eta in stage:
        A=isotope_step(A,dur);A[:3]*=eta;t+=dur
        rows.append({'stage':name,'duration_min':dur,'age_d':t/1440.,**dict(zip(ISO,A))})
    process_h=5+sum(s[1] for s in stage)/60
    fluids_dga=[R['streams'][x] for x in ['DILUTED-FEED','DGA-PRE','DGA-RINSE','DGA-ELUENT']]
    fluids_ln=[R['streams'][x] for x in ['LN-FEED','LN-PRE','LN-RINSE','LN-ELUENT']]
    pressure={}
    for label,fluids,bed,diam in [('DGA2',fluids_dga,5,11),('LN',fluids_ln,2,8)]:
        pressure[label]=max(lp.ergun(m,flow,bed,diam) for m in fluids)
        pressure[label+'_upper']=max(lp.ergun(m,flow,bed,diam,50,.35) for m in fluids)
        pressure[label+'_lower']=max(lp.ergun(m,flow,bed,diam,150,.45) for m in fluids)
    return {'recipe':R,'stage_inventory':rows,'process_h':process_h,'pressure_Pa':pressure,'activity_Bq':A.tolist()}
