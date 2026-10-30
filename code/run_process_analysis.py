"""Engineering extension: property check, material balances, hydraulics, robust controls."""
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar,brentq
from scipy.interpolate import PchipInterpolator
from thermo import electrochem
from process_model import ROOT,BASE,reference,lp,LAM,ISO
from model import total_window
R=reference();out=ROOT/'results'
# Independent experimental density differences, published after the 2009 property correlation.
den=pd.read_csv(ROOT/'data/nitric_density_check.csv')
rho_w=lp.liquid_density([1,0,0,0,0])
den['acid_mass_fraction']=den.molality_mol_kg_water*.063013/(1+den.molality_mol_kg_water*.063013)
den['predicted_delta_density_kg_m3']=[lp.liquid_density([1-w,w,0,0,0])-rho_w for w in den.acid_mass_fraction]
den['residual_kg_m3']=den.predicted_delta_density_kg_m3-den.delta_density_kg_m3
den['relative_density_error_percent']=100*den.residual_kg_m3/(rho_w+den.delta_density_kg_m3)
den.to_csv(out/'density_independent_check.csv',index=False)
# Compare independently coded formulas with the established open-source implementation.
checks=[]
for name,cas in [('Nitric acid','7697-37-2'),('Hydrogen chloride','7647-01-0'),('Ammonia','7664-41-7'),('Ammonium Chloride','12125-02-9')]:
    for w in [.001,.01,.10,.20]:
        m=[1-w,0,0,0,0];m[lp.COMPONENTS.index(name)]=w
        for method,oracle in [(lp.liquid_density,electrochem.Laliberte_density),(lp.liquid_viscosity,electrochem.Laliberte_viscosity)]:
            actual=method(m);expected=oracle(298.15,[w],[cas]);checks.append(abs(actual-expected)/expected)
assert max(checks)<1e-12
# Material balance, reaction stoichiometry and acid-equivalent balance.
src=lp.mix_mass(BASE['streams']['PRIMARY-FRACTION'],*[BASE['streams'][t] for t in BASE['fresh_tags']])
dst=lp.mix_mass(*[BASE['streams'][t] for t in ['DGA-WASTE','LN-WASTE','FINAL-PRODUCT']])
extent=BASE['streams']['NH3-ADDITION'][3]/(lp.MW[3]/1000)
stoich=np.array([0,0,-lp.MW[2]/1000,-lp.MW[3]/1000,(lp.MW[2]+lp.MW[3])/1000])
res=np.array(dst)-np.array(src)-extent*stoich
assert max(abs(res))<1e-12
pd.DataFrame({'component':lp.COMPONENTS,'in_kg':src,'reaction_change_kg':extent*stoich,'out_kg':dst,'residual_kg':res}).to_csv(out/'component_balance.csv',index=False)
pd.DataFrame(R['stage_inventory']).to_csv(out/'process_stage_inventory.csv',index=False)
# Flow changes are unvalidated hydraulic/scheduling scenarios. No eta(flow) is fitted.
flowrows=[]
for q in np.linspace(1.5,9,76):
    r=reference(flow=float(q));saved_h=R['process_h']-r['process_h']
    flowrows.append({'flow_mL_min':q,'process_h':r['process_h'],'saved_h':saved_h,
      'same_start_activity_gain_percent':100*np.expm1(LAM[0]*saved_h/24),
      'break_even_recovery_loss_percent':100*(-np.expm1(-LAM[0]*saved_h/24)),
      **r['pressure_Pa']})
pd.DataFrame(flowrows).to_csv(out/'hydraulic_flow_screen.csv',index=False)
# Tie deterministic DWSIM processing times to the original delivery uncertainty model.
window=total_window(10.8+2*.83,.0014+2*.00011)
ext_min=3-9/24;ext_max=5-9/24
schedule=[]
for q in [1.5,3,6,9]:
    p=reference(flow=q)['process_h']/24
    for width_h in [48,24,12,0]:
        dmin=ext_min+p;dmax=dmin+width_h/24
        c=window[0]-dmin
        guarantee=.6042*np.exp(-LAM[0]*(c+dmax))
        schedule.append({'flow_mL_min':q,'delay_width_h':width_h,'process_h':p*24,'robust_start_d':c,
          'worst_delivery_age_d':c+dmax,'guaranteed_A225_over_A0':guarantee})
sch=pd.DataFrame(schedule);sch.to_csv(out/'process_robust_schedule.csv',index=False)
assert sch.groupby('delay_width_h').guaranteed_A225_over_A0.std().max()<1e-14
# Sparse pH data: interpolation/control demonstration is restricted to this different citrate matrix.
phdata=pd.DataFrame({'load_pH':[1.5,2.,2.5],'Ac_recovery_fraction':[.95,.98,.66],
 'source':['Fitzsimmons et al 2019, Molecules 24,1921, section 2.7.2; 1 M citrate load']*3})
phdata.to_csv(ROOT/'data/pH_recovery_literature.csv',index=False)
control=[];pchip=PchipInterpolator(phdata.load_pH,phdata.Ac_recovery_fraction)
for delta in [0,.05,.10,.15,.20,.25]:
    for kind,fun in [('piecewise_linear',lp.ph_recovery),('PCHIP',lambda x:float(pchip(x)))]:
        def worst(c):return min(fun(c-delta),fun(c+delta))
        opt=minimize_scalar(lambda c:-worst(c),bounds=(1.5+delta,2.5-delta),method='bounded',options={'xatol':1e-12})
        control.append({'interpolation':kind,'pH_halfwidth':delta,'optimal_setpoint_pH':opt.x,
          'worst_recovery_optimal':worst(opt.x),'worst_recovery_at_pH2':worst(2.)})
pd.DataFrame(control).to_csv(out/'pH_control_sensitivity.csv',index=False)
# Boundary assumptions: vary unreported DGA cut acidity and inter-column recovery allocation.
boundary=[]
for c in [.05,.1,.2,.5]:
    r=reference(carry_M=c)
    boundary.append({'DGA_cut_HCl_M':c,'NH3_solution_mL':r['recipe']['NH3_volume_mL'],
      'LN_feed_volume_mL':lp.liquid_volume(r['recipe']['streams']['LN-FEED']),
      'process_h':r['process_h']})
pd.DataFrame(boundary).to_csv(out/'unreported_acid_carryover_sensitivity.csv',index=False)
allocation=[]
for eta1 in [.92,.95,.98,1.]:
    eta2=(.636/.7)/eta1;r=reference(eta_dga=eta1,eta_ln=eta2)
    allocation.append({'eta_DGA2':eta1,'eta_LN':eta2,'A225_final_Bq':r['activity_Bq'][0]})
assert np.std([r['A225_final_Bq'] for r in allocation])<1e-8
pd.DataFrame(allocation).to_csv(out/'recovery_factorization_check.csv',index=False)
full=reference(flow=6)
summary={'density_check_points':len(den),'density_max_absolute_residual_kg_m3':float(abs(den.residual_kg_m3).max()),
 'density_max_relative_error_percent':float(abs(den.relative_density_error_percent).max()),
 'density_RMSE_kg_m3':float(np.sqrt(np.mean(den.residual_kg_m3**2))),
 'property_formula_max_relative_error':max(checks),'component_balance_max_residual_kg':float(max(abs(res))),
 'mixed_feed_volume_mL':lp.liquid_volume(BASE['streams']['DILUTED-FEED']),
 'mixed_feed_HNO3_M':.8/(lp.liquid_volume(BASE['streams']['DILUTED-FEED'])/1000),
 'mixing_contraction_percent':100*(1-lp.liquid_volume(BASE['streams']['DILUTED-FEED'])/160),
 'fresh_solution_volume_mL':sum(lp.liquid_volume(BASE['streams'][t]) for t in BASE['fresh_tags']),
 'fresh_HNO3_mol':sum(BASE['streams'][t][1] for t in BASE['fresh_tags'])/(lp.MW[1]/1000),
 'fresh_HCl_mol':sum(BASE['streams'][t][2] for t in BASE['fresh_tags'])/(lp.MW[2]/1000),
 'partial_waste_mass_g':1000*sum(sum(BASE['streams'][t]) for t in ['DGA-WASTE','LN-WASTE']),
 'nominal_conditioner_NH3_solution_mL':BASE['NH3_volume_mL'],
 'purification_pumped_volume_mL':BASE['dga_pumped_mL']+BASE['ln_pumped_mL'],
 'unassigned_protocol_allowance_min':BASE['protocol_allowance_min'],
 'baseline_process_h':R['process_h'],'q6_process_h':full['process_h'],
 'q6_same_start_gain_percent':100*np.expm1(LAM[0]*(R['process_h']-full['process_h'])/24),
 'q6_break_even_relative_recovery_loss_percent':100*(-np.expm1(-LAM[0]*(R['process_h']-full['process_h'])/24)),
 'q6_pressure_upper_kPa':max(full['pressure_Pa']['DGA2_upper'],full['pressure_Pa']['LN_upper'])/1000,
 'illustrative_50kPa_upper_packing_flow_limit_mL_min':brentq(lambda q:max(reference(flow=q)['pressure_Pa'][k] for k in ['DGA2_upper','LN_upper'])-50000,3,9),
 'queue_width_48_to24_gain_percent':100*np.expm1(LAM[0]),'pH_control_delta02':[r for r in control if r['pH_halfwidth']==.2],
 'assumptions':'Steady batch-equivalent carrier balances; instantaneous empirical collection at operation end; no measured column-specific efficiencies or DGA carryover acidity; pH data belong to an independent citrate/AG50 experiment.'}
(out/'process_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
