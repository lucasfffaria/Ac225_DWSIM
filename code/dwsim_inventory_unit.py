# DWSIM CustomUO kernel. IronPython 2.7 compatible and self-contained.
# No external modules, property packages, compound surrogates, or paid APIs.
# Activities are in Bq, times in days, and ratios are dimensionless.
import math

v = Me.InputVariables
o = Me.OutputVariables
o.Clear()
la = math.log(2.0) / v['half225_d']
l6 = math.log(2.0) / v['half226_d']
l7 = math.log(2.0) / v['half227_d']
lr = math.log(2.0) / v['halfRa225_d']
a, b = l6-la, la-l7
r6, r7 = v['r226'], v['r227']
t, s = v['delivery_age_d'], v['parent_removal_age_d']
eta, p = v['eta_Ac'], v['parent_ratio']
if t < s or s < 0 or eta <= 0 or eta > 1:
    raise ValueError('Require 0 <= parent removal <= delivery and 0 < eta <= 1')
gain = 1.0 + p * la/(la-lr) * (math.exp((la-lr)*s)-1.0)
A0 = v['A225_reference_Bq']
A225 = A0 * eta * math.exp(-la*t) * gain
A226 = A0 * eta * r6 * math.exp(-l6*t)
A227 = A0 * eta * r7 * math.exp(-l7*t)
o['A225_delivery_Bq'] = A225
o['A226_delivery_Bq'] = A226
o['A227_delivery_Bq'] = A227
o['A225_delivery_per_reference'] = A225/A0
o['R226_delivery'] = A226/A225
o['R227_delivery'] = A227/A225
o['P_Ac_delivery'] = A225/(A225+A226+A227)
o['parent_Ac_gain_factor'] = gain
# All following timing windows refer to the parent-free baseline.
tstar = max(0.0,math.log(a*r6/(b*r7))/(a+b))
o['parent_free_best_purity_age_d'] = tstar
o['parent_free_max_P_Ac'] = 1.0/(1.0+r6*math.exp(-a*tstar)+r7*math.exp(b*tstar))
lo = max(0.0,math.log(r6/v['R226_limit'])/a)
hi = math.log(v['R227_limit']/r7)/b
o['parent_free_individual_lower_age_d'] = lo
o['parent_free_individual_upper_age_d'] = hi
o['parent_free_individual_window_exists'] = 1.0 if hi>=lo else 0.0
limit = 1.0/v['P_Ac_limit']-1.0
def fn(age):
    return r6*math.exp(-a*age)+r7*math.exp(b*age)-limit
def bisect(left,right):
    fl=fn(left)
    for k in range(100):
        mid=(left+right)/2.0
        fm=fn(mid)
        if fl*fm<=0: right=mid
        else: left=mid;fl=fm
    return (left+right)/2.0
if fn(tstar) <= 0:
    tl=0.0 if fn(0)<=0 else bisect(0,tstar)
    right=max(tstar+1,2*tstar)
    while fn(right)<0: right*=2
    o['parent_free_P_Ac_window_exists']=1.0
    o['parent_free_P_Ac_lower_age_d']=tl
    o['parent_free_P_Ac_upper_age_d']=bisect(tstar,right)
    o['parent_free_start_window_lower_d']=max(0.0,tl-v['delay_min_d'])
    o['parent_free_start_window_upper_d']=o['parent_free_P_Ac_upper_age_d']-v['delay_max_d']
    o['parent_free_start_window_exists']=1.0 if o['parent_free_start_window_upper_d']>=o['parent_free_start_window_lower_d'] else 0.0
    o['parent_free_earliest_start_worst_A225_per_reference']=v['eta_lower']*math.exp(-la*(o['parent_free_start_window_lower_d']+v['delay_max_d']))
else:
    o['parent_free_P_Ac_window_exists']=0.0
    # -1 explicitly denotes no solution, not zero days.
    o['parent_free_P_Ac_lower_age_d']=-1.0
    o['parent_free_P_Ac_upper_age_d']=-1.0
    o['parent_free_start_window_exists']=0.0
    o['parent_free_start_window_lower_d']=-1.0
    o['parent_free_start_window_upper_d']=-1.0
o['relative_gain_if_9h_saved']=math.exp(la*9.0/24.0)-1.0
o['parent_free_critical_r227_for_individual_limits']=v['R227_limit']*math.exp(-b*lo)
# CustomUO exports numeric scope variables after execution.
for _name in list(o.Keys):
    globals()[str(_name)] = o[_name]
