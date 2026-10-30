"""Transparent inventory model. Time: days. Activities may use any common unit.

This is not a thermodynamic or chromatographic-retention predictor.
The parent activity ratio is an explicit scenario input, not a fitted datum.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.linalg import expm

@dataclass(frozen=True)
class Decay:
    t225: float=9.92
    t226: float=1.224
    t227: float=21.7*365.25
    tra: float=14.9
    @property
    def lam(self): return np.log(2)/np.array([self.t225,self.t226,self.t227,self.tra])
    @property
    def a(self): return self.lam[1]-self.lam[0]
    @property
    def b(self): return self.lam[0]-self.lam[2]

def ratios(t,r226,r227,d=Decay()):
    return r226*np.exp(-d.a*np.asarray(t)), r227*np.exp(d.b*np.asarray(t))

def purity(t,r226,r227,d=Decay()):
    x,y=ratios(t,r226,r227,d)
    return 1/(1+x+y)

def best_age(r226,r227,d=Decay()):
    return np.maximum(0,np.log(d.a*np.asarray(r226)/(d.b*np.asarray(r227)))/(d.a+d.b))

def individual_window(r226,r227,u226=.01,u227=.01,d=Decay()):
    lo=np.maximum(0,np.log(np.asarray(r226)/u226)/d.a)
    hi=np.log(u227/np.asarray(r227))/d.b
    return lo,hi,hi>=lo

def total_window(r226,r227,pmin=.99,d=Decay()):
    t=float(best_age(r226,r227,d))
    limit=1/pmin-1
    f=lambda t: sum(ratios(t,r226,r227,d))-limit
    if f(t)>1e-13: return None
    if abs(f(t))<1e-13: return (t,t)
    lo=0.0 if f(0)<=0 else brentq(f,0,t,xtol=1e-12)
    hi=max(t+1,2*t)
    while f(hi)<0: hi*=2
    return float(lo),float(brentq(f,t,hi,xtol=1e-12))

def critical_r227(r226,u226=.01,u227=.01,d=Decay()):
    """Largest reference r227 allowing both individual ratio limits."""
    tlo=np.maximum(0,np.log(np.asarray(r226)/u226)/d.a)
    return u227*np.exp(-d.b*tlo)

def parent_gain(s,parent_ratio,d=Decay()):
    la,_,_,lr=d.lam
    delta=la-lr
    return 1+parent_ratio*la*np.expm1(delta*np.asarray(s))/delta

def released_inventory(t,s,r226,r227,parent_ratio=0.,eta=.636,d=Decay()):
    """Ac only after instantaneous parent removal at s, common Ac recovery eta.

    All parent-derived Ac is formed before s. Later Ra-derived Ac stays outside
    the delivered product. No Ra carryover or second generator harvest is added.
    """
    if np.any(np.asarray(s)<0) or np.any(np.asarray(t)<np.asarray(s)):
        raise ValueError('Require 0 <= separation age <= delivery age')
    gain=parent_gain(s,parent_ratio,d)
    a225=eta*np.exp(-d.lam[0]*np.asarray(t))*gain
    a226=eta*r226*np.exp(-d.lam[1]*np.asarray(t))
    a227=eta*r227*np.exp(-d.lam[2]*np.asarray(t))
    return np.array([a225,a226,a227])

def critical_parent_ratio(r226,r227,post_separation_days=4.,mode='total',pmin=.99,d=Decay()):
    """Minimum EOB-equivalent A_Ra225/A_Ac225 opening any delivery window.

    Fixed time from parent removal to release, latest permitted separation.
    Threshold is conditional on complete retention of Ac and removal of Ra.
    """
    def required(t):
        x,y=ratios(t,r226,r227,d)
        required_gain=(x+y)/(1/pmin-1) if mode=='total' else max(x/.01,y/.01)
        unit_gain=parent_gain(t-post_separation_days,1,d)-1
        return max(0,(required_gain-1)/unit_gain)
    opt=minimize_scalar(required,bounds=(post_separation_days+1e-5,80),method='bounded',options={'xatol':1e-11})
    return float(opt.fun),float(opt.x)

def activity_matrix(d=Decay()):
    """Independent activity-space matrix for 225Ac, 226Ac, 227Ac, 225Ra."""
    la=d.lam
    M=-np.diag(la)
    M[0,3]=la[0]
    return M

def matrix_release(t,s,r226,r227,parent_ratio=0.,eta=.636,d=Decay()):
    M=activity_matrix(d)
    before=expm(M*s)@np.array([1.,r226,r227,parent_ratio])
    transfer=np.diag([eta,eta,eta,0.])
    return (expm(M*(t-s))@transfer@before)[:3]

def recovery_time_equivalent(eta_slow,eta_fast,d=Decay()):
    return np.log(eta_slow/eta_fast)/d.lam[0]

def run_recipe(rows):
    volumes=np.array([float(r['volume_mL']) for r in rows])
    hno3=sum(float(r['volume_mL'])*float(r['HNO3_M'])/1000 for r in rows)
    hcl=sum(float(r['volume_mL'])*float(r['HCl_M'])/1000 for r in rows)
    fresh=float(volumes.sum())
    # Primary Ac eluate and concentrated DGA2 eluate are passed once more
    # through a downstream column; these are transfers, not fresh additions.
    internal=80.+10.
    return {'known_fresh_mL':fresh,'known_internal_transfer_mL':internal,
            'known_pumped_mL':fresh+internal,'HNO3_mol':hno3,'HCl_mol':hcl,
            'known_nonproduct_collection_mL':fresh-10.,
            'boundary_note':'Partial documented liquid additions only; excludes target feed, sulfate, initial water, dead-volume rinses and NH4OH adjustment. Collection volume assumes additive volumes and no unrecorded losses.'}

def robust_start_window(r226,r227,sigma226,sigma227,halfwidth_sigmas=2.,delay_min=3.,delay_max=5.,pmin=.99,d=Decay()):
    """Box-robust scheduling, not a 95% joint confidence region.

    Monotonicity makes simultaneous upper isotope ratios the worst case.
    The full delay interval must fit inside the resulting delivery-age window.
    """
    window=total_window(r226+halfwidth_sigmas*sigma226,r227+halfwidth_sigmas*sigma227,pmin,d)
    if window is None:return None
    lo=max(0.,window[0]-delay_min);hi=window[1]-delay_max
    return (lo,hi) if lo<=hi else None
