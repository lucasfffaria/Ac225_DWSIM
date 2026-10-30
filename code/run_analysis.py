from pathlib import Path
import csv,json,sys,platform
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from model import *

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
cfg=json.loads((ROOT/'data/source_inputs.json').read_text())
d=Decay(); rng=np.random.default_rng(cfg['random_seed'])
cases=cfg['isotope_cases']; n=cfg['monte_carlo_samples']
def quant(x):
    return {k:float(v) for k,v in zip(['q025','median','q975'],np.quantile(x,[.025,.5,.975]))}
def savecsv(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary={'model_version':'1.0.0','seed':cfg['random_seed'],'n':n,'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},'cases':{}}
curves=[]; mc={}; corr=[]
for c in cases:
    r6,r7=c['r226'],c['r227'];tstar=float(best_age(r6,r7,d));lo,hi,feas=individual_window(r6,r7,d=d)
    x14,y14=ratios(14,r6,r7,d)
    out={'t_best_purity_d':tstar,'pmax':float(purity(tstar,r6,r7,d)),
         'individual_lower_d':float(lo),'individual_upper_d':float(hi),'individual_feasible':bool(feas),
         'total_99_window_d':total_window(r6,r7,.99,d),'total_995_window_d':total_window(r6,r7,.995,d),
         'r226_day14':float(x14),'r227_day14':float(y14),'p_day14':float(purity(14,r6,r7,d)),
         'critical_r227_individual':float(critical_r227(r6,d=d)),
         'critical_parent_total_99':critical_parent_ratio(r6,r7,4.,d=d),
         'critical_parent_individual':critical_parent_ratio(r6,r7,4.,mode='individual',d=d)}
    for rho in [-.8,0.,.8]:
        cov=np.array([[1,rho],[rho,1]])
        z=rng.multivariate_normal([0,0],cov,size=n)
        rs=np.column_stack([r6+c['r226_sigma']*z[:,0],r7+c['r227_sigma']*z[:,1]])
        if np.any(rs<=0): raise RuntimeError('Nonpositive Gaussian isotope draw; implement explicit positive sampling')
        ts=best_age(rs[:,0],rs[:,1],d); ps=purity(ts,rs[:,0],rs[:,1],d)
        ls,hs,fs=individual_window(rs[:,0],rs[:,1],d=d)
        corr.append({'case':c['id'],'rho':rho,'p_individual_feasible':float(fs.mean()),'p_total99_feasible':float((ps>=.99).mean()),'pmax_q025':float(np.quantile(ps,.025)),'pmax_q975':float(np.quantile(ps,.975))})
        if rho==0:
            out['MC']={'pmax':quant(ps),'best_age':quant(ts),'lower_age':quant(ls),'upper_age':quant(hs),'p_individual_feasible':float(fs.mean()),'p_total99_feasible':float((ps>=.99).mean())}
            mc[c['id']]=rs[:3000]
    win=out['total_99_window_d']
    if win:
        out['release_opt_99_d']=win[0]
        out['start_opt_99_d']=win[0]-4
        out['A225_opt_99_per_A0']=float(.636*np.exp(-d.lam[0]*win[0]))
        out['A225_day14_per_A0']=float(.636*np.exp(-d.lam[0]*14))
        out['relative_gain_vs_day14']=float(np.exp(d.lam[0]*(14-win[0]))-1)
    out['new_half_life_sensitivity']={'t_best_d':float(best_age(r6,r7,Decay(t225=9.9172))),'pmax':float(purity(float(best_age(r6,r7,Decay(t225=9.9172))),r6,r7,Decay(t225=9.9172)))}
    summary['cases'][c['id']]=out
    for t in np.linspace(0,45,901):
        x,y=ratios(t,r6,r7,d)
        curves.append({'case':c['id'],'age_d':t,'r226':float(x),'r227':float(y),'p_actinium':float(1/(1+x+y)),'A225_per_A0':float(.636*np.exp(-d.lam[0]*t))})
np.savez_compressed(OUT/'mc_isotope_ratio_draws.npz',**mc)
savecsv('time_curves.csv',curves);savecsv('correlation_sensitivity.csv',corr)

train=[x for x in cfg['recovery_cases'] if x['role']=='calibration']
eta_cal=float(np.mean([x['reported_yield_pct']/100 for x in train]))
validation=[]
for row in cfg['recovery_cases']:
    pred=eta_cal*row['feed_kBq']
    validation.append({**row,'predicted_kBq':pred,'relative_residual_pct':100*(pred-row['product_kBq'])/row['product_kBq']})
summary['recovery']={'eta_calibration':eta_cal,'eta_scenarios':cfg['nominal_recovery'],'check_rows':validation,'max_double_counting_effect_9h_relative':float(np.exp(d.lam[0]*9/24)-1)}
savecsv('source_recovery_check.csv',validation)

recipe=list(csv.DictReader((ROOT/'data/reagent_recipe.csv').open(encoding='utf8')))
ledger=run_recipe(recipe)
ledger['known_pumped_minutes_at_3mL_min']=ledger['known_pumped_mL']/3
ledger['reported_chromatography_minutes']=540.
ledger['unallocated_minutes']=540-ledger['known_pumped_mL']/3
ledger['fraction_time_unallocated']=ledger['unallocated_minutes']/540
summary['recipe_ledger']=ledger
rows=[]
for group in ['primary','purification']:
    selected=[x for x in recipe if x['stage']==group]
    fresh=sum(float(x['volume_mL']) for x in selected)
    transfer=0. if group=='primary' else 90.
    rows.append({'stage':group,'known_fresh_mL':fresh,'known_internal_transfer_mL':transfer,'known_pumped_min':(fresh+transfer)/3,'reported_total_min':300 if group=='primary' else 240})
savecsv('stage_time_ledger.csv',rows)

effects=[]
for hrs in [1.8/3600,1/60,30/60,1,4,8,9,24,48,72,96]:
    effects.append({'time_saved_h':hrs,'gain_at_equal_recovery_pct':100*np.expm1(d.lam[0]*hrs/24)})
summary['intensification']={'effects':effects,'eta_99_vs_96_break_even_hours':24*float(recovery_time_equivalent(.99,.96,d)),
   'eta_636_vs_92_break_even_days':float(recovery_time_equivalent(.92,.636,d)),
   'fixed_release_age_gain_at_equal_recovery':0.,'eliminate_entire_9h_gain_pct':100*float(np.expm1(d.lam[0]*9/24)),
   'double_flow_max_known_pumping_gain_pct':100*float(np.expm1(d.lam[0]*(ledger['known_pumped_mL']/6)/1440))}
savecsv('time_savings.csv',effects)

parent_curves=[]
for c in cases:
    for p in [0,.1,.25,.5,1.]:
        for t in np.linspace(4.01,45,411):
            act=released_inventory(t,t-4,c['r226'],c['r227'],p,d=d)
            parent_curves.append({'case':c['id'],'parent_ratio':p,'delivery_age_d':t,'separation_age_d':t-4,'A225_per_A0':float(act[0]),'p_actinium':float(act[0]/sum(act))})
savecsv('parent_scenarios.csv',parent_curves)

# Additional prior-art benchmark: Fitzsimmons et al. 2019, Table 2 and Fig. 1.
# These are published calculated values, NOT new experimental measurements.
db=Decay(t225=10.,t226=29.37/24,t227=21.772*365.25)
before=np.array([206.91,433.01,.42]);reported=np.array([127.37,8.21,.42])
after=before*np.exp(-db.lam[:3]*7)
tb=float(best_age(before[1]/before[0],before[2]/before[0],db))
summary['prior_art_benchmark']={'source':'Fitzsimmons et al. 2019 Table 2, PDF p3; calculated EOB and final-separation activities',
 'EOB_mCi':before.tolist(),'separation_day':7,'reported_final_mCi':reported.tolist(),'predicted_final_mCi':after.tolist(),
 'relative_difference_pct':(100*(after-reported)/reported).tolist(),'best_age_after_separation_d':tb-7,
 'maximum_radioisotopic_purity':float(purity(tb,before[1]/before[0],before[2]/before[0],db)),
 'source_rounded_optimum':'about 9 days after separation; 99.3%'}

# Robust scheduling contribution: isotope interval + variable delivery delay.
robust_rows=[];schedule_curves=[]
for c in cases:
    robust=robust_start_window(c['r226'],c['r227'],c['r226_sigma'],c['r227_sigma'])
    meanwin=total_window(c['r226'],c['r227'])
    if not robust:
        summary['cases'][c['id']]['robust_schedule']={'feasible':False}
        continue
    eta_min=.636*.95
    nominal_start=meanwin[0]-4
    robust_start=robust[0]
    worst_release=robust_start+5
    guaranteed=eta_min*np.exp(-d.lam[0]*worst_release)
    draws=mc[c['id']]
    # Separate 100000 fresh Gaussian draws and uniform 3-5 d scenario delays.
    z=rng.normal(size=(n,2))
    rr=np.array([c['r226'],c['r227']])+z*np.array([c['r226_sigma'],c['r227_sigma']])
    delay=rng.uniform(3,5,n)
    def passprob(start):
        return float((purity(start+delay,rr[:,0],rr[:,1],d)>=.99).mean())
    rob={'feasible':True,'input_box_halfwidth_sigma':2,'delay_bounds_d':[3,5],
         'start_window_d':robust,'earliest_start_d':robust_start,'worst_delivery_d':worst_release,
         'eta_lower_scenario':eta_min,'guaranteed_A225_per_A0':float(guaranteed),
         'nominal_earliest_start_d':nominal_start,
         'nominal_start_scenario_pass_fraction':passprob(nominal_start),
         'robust_start_scenario_pass_fraction':passprob(robust_start),
         'relative_activity_cost_vs_nominal':float(1-guaranteed/(.636*np.exp(-d.lam[0]*meanwin[0])))}
    robust_age=total_window(c['r226']+2*c['r226_sigma'],c['r227']+2*c['r227_sigma'])
    rob['robust_delivery_age_window_d']=robust_age
    rob['penalty_factors']={'chemical_recovery_lower_bound':.95,'isotope_uncertainty':float(np.exp(-d.lam[0]*(robust_age[0]-meanwin[0]))),'delivery_interval_width':float(np.exp(-d.lam[0]*2))}
    summary['cases'][c['id']]['robust_schedule']=rob
    for start in np.linspace(9,25,161):
        t1,t2=start+3,start+5
        worstp=min(purity(t1,c['r226']+2*c['r226_sigma'],c['r227']+2*c['r227_sigma']),purity(t2,c['r226']+2*c['r226_sigma'],c['r227']+2*c['r227_sigma']))
        schedule_curves.append({'start_age_d':start,'nominal_release_P_Ac':float(purity(start+4,c['r226'],c['r227'])),'box_worst_P_Ac':float(worstp),'worst_A225_per_A0':float(eta_min*np.exp(-d.lam[0]*(start+5)))})
    for delay_width in np.linspace(0,7,29):
        for zbox in [0,1,2,3]:
            w=robust_start_window(c['r226'],c['r227'],c['r226_sigma'],c['r227_sigma'],zbox,4-delay_width/2,4+delay_width/2)
            robust_rows.append({'delay_width_d':delay_width,'box_sigmas':zbox,'start_lower_d':w[0] if w else None,'start_upper_d':w[1] if w else None,'window_width_d':w[1]-w[0] if w else 0})
savecsv('robust_schedule_curves.csv',schedule_curves);savecsv('robust_schedule_sensitivity.csv',robust_rows)

# Numerical verification is distinct from comparison with experimental literature.
errs=[]; ode_errs=[]; atom_errs=[]; semigroup=[]; cancel=[]
for _ in range(200):
    t=float(rng.uniform(.1,80));s=float(rng.uniform(0,t));r6=float(10**rng.uniform(-1,2));r7=float(10**rng.uniform(-4,-1));p=float(rng.uniform(0,2));eta=float(rng.uniform(.1,1))
    closed=released_inventory(t,s,r6,r7,p,eta,d)
    matrix=matrix_release(t,s,r6,r7,p,eta,d)
    errs.append(float(np.max(np.abs(closed-matrix)/np.maximum(closed,1e-14))))
    v=np.array([1.,r6,r7,p]);M=activity_matrix(d)
    solution=solve_ivp(lambda _,x:M@x,[0,s],v,rtol=1e-10,atol=1e-13).y[:,-1]
    solution[:3]*=eta;solution[3]=0
    after=solve_ivp(lambda _,x:M@x,[s,t],solution,rtol=1e-10,atol=1e-13).y[:3,-1]
    ode_errs.append(float(np.max(np.abs(closed-after)/np.maximum(closed,1e-7))))
    semigroup.append(float(np.max(np.abs(expm(M*(t-s))@expm(M*s)-expm(M*t)))))
    n0=v/d.lam
    K=-np.diag(np.append(d.lam,0.));K[0,3]=d.lam[3];K[4,:3]=d.lam[:3]
    nf=expm(K*t)@np.append(n0,0.)
    atom_errs.append(abs(sum(nf)-sum(n0))/sum(n0))
    a=released_inventory(t,s,r6,r7,0,eta,d);x,y=ratios(t,r6,r7,d)
    cancel.append(float(max(abs(a[1]/a[0]-x),abs(a[2]/a[0]-y))))
verification={'closed_vs_matrix_max_relative_error':max(errs),'closed_vs_ode_max_scaled_error':max(ode_errs),'atom_ledger_closure_max_relative_error':max(atom_errs),'matrix_semigroup_max_abs_error':max(semigroup),'eta_cancellation_max_abs_error':max(cancel),'random_cases':200}
assert max(errs)<1e-10
assert max(ode_errs)<1e-5
assert max(atom_errs)<1e-12
assert max(cancel)<1e-10
assert individual_window(10.8,.0014)[2]
assert not individual_window(6.62,.0045)[2]
assert total_window(10.8,.0014,.99) is not None
assert total_window(6.62,.0045,.99) is None
assert all(x['p_total99_feasible']==0 for x in corr if x['case']=='nitrogen')
pr=cases[0];rob=summary['cases']['proton']['robust_schedule']
# Verify the robust guarantee at all ratio-box corners and over the entire
# delay interval; the mathematical proof is monotonicity and convexity.
corner_min=1.
for sign6 in [-1,1]:
    for sign7 in [-1,1]:
        for delay in np.linspace(3,5,101):
            pp=purity(rob['earliest_start_d']+delay,pr['r226']+sign6*2*pr['r226_sigma'],pr['r227']+sign7*2*pr['r227_sigma'])
            corner_min=min(corner_min,float(pp))
assert corner_min>=.99-1e-12
wshift=robust_start_window(pr['r226'],pr['r227'],pr['r226_sigma'],pr['r227_sigma'],delay_min=2.625,delay_max=4.625)
assert abs((wshift[0]+4.625)-rob['worst_delivery_d'])<1e-10
verification['robust_corner_and_delay_min_purity']=corner_min
verification['robust_time_translation_invariance_verified']=True
(OUT/'verification.json').write_text(json.dumps(verification,indent=2))
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({'cases':summary['cases'],'ledger':ledger,'recovery':summary['recovery'],'intensification':summary['intensification'],'verification':verification},indent=2))
