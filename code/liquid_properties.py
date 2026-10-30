"""Isothermal aqueous analytical-component properties (Python 2/3 compatible).

Laliberte 2009, doi:10.1021/je8008123. PARAMS is injected into native scripts.
No acid speciation, VLE, retention, or mixing-enthalpy model is implied.
"""
from math import exp,log,pi,sqrt
COMPONENTS=['Water','Nitric acid','Hydrogen chloride','Ammonia','Ammonium Chloride']
MW=[18.01528,63.01284,36.46094,17.03052,53.49146] # g/mol; last = HCl + NH3

def liquid_density(masses,T=298.15):
    total=sum(masses)
    if total<=0: raise ValueError('A liquid stream needs positive mass')
    w=[x/total for x in masses]; ww=w[0]; t=T-273.15
    rho_w=(((((-2.8054253e-10*t+1.0556302e-7)*t-4.6170461e-5)*t-.0079870401)*t+16.945176)*t+999.83952)/(1+.01687985*t)
    invrho=ww/rho_w
    for i,name in enumerate(COMPONENTS[1:],1):
        if w[i]<=0:continue
        p=PARAMS[name]; low,high,wmax=p['density_range']
        if t<low-1e-8 or t>high+1e-8 or w[i]>wmax+1e-8:raise ValueError('Density correlation outside domain: '+name)
        c0,c1,c2,c3,c4=p['density']
        rho=((c0*(1-ww)+c1)*exp(1e-6*(t+c4)**2))/((1-ww)+c2+c3*t)
        invrho+=w[i]/rho
    return 1/invrho

def liquid_viscosity(masses,T=298.15):
    total=sum(masses);w=[x/total for x in masses];ww=w[0];t=T-273.15
    mu_w=(t+246)/((.05594*t+5.2842)*t+137.37) # mPa s
    log_mu=ww*log(mu_w)
    for i,name in enumerate(COMPONENTS[1:],1):
        if w[i]<=0:continue
        p=PARAMS[name];low,high,wmax=p['viscosity_range']
        if t<low-1e-8 or t>high+1e-8 or w[i]>wmax+1e-8:raise ValueError('Viscosity correlation outside domain: '+name)
        v1,v2,v3,v4,v5,v6=p['viscosity']
        # Evaluate logarithm directly; pure-water branch avoids negative powers at zero.
        log_mu+=w[i]*((v1*(1-ww)**v2+v3)/(v4*t+1)-log(v5*(1-ww)**v6+1))
    return exp(log_mu)*1e-3

def liquid_volume(masses):
    return sum(masses)/liquid_density(masses)*1e6 # kg / kg m-3 -> mL

def solution(volume_mL, concentrations=None):
    """Masses kg from final solution volume and analytical molarities at 25 C."""
    if concentrations is None: concentrations={}
    masses=[0.0]*5
    for name,C in concentrations.items():
        i=COMPONENTS.index(name);masses[i]=C*volume_mL/1000*MW[i]/1000
    solute=sum(masses)
    # Bracket water mass inside each binary correlation's domain.
    min_water=1e-12
    for i,name in enumerate(COMPONENTS[1:],1):
        if masses[i]>0:min_water=max(min_water,masses[i]/PARAMS[name]['density_range'][2]-solute)
    lo=max(min_water,1e-12);hi=max(volume_mL*.002,lo*2)
    masses[0]=lo
    if liquid_volume(masses)>volume_mL+1e-8:raise ValueError('Requested concentration exceeds validated domain')
    for k in range(80):
        masses[0]=(lo+hi)/2
        if liquid_volume(masses)>volume_mL:hi=masses[0]
        else:lo=masses[0]
    masses[0]=(lo+hi)/2
    return masses

def mix_mass(*vectors):
    return [sum(v[i] for v in vectors) for i in range(5)]

def subtract_mass(a,b):
    out=[a[i]-b[i] for i in range(5)]
    if min(out)<-1e-12:raise ValueError('Specified collection cut exceeds available component mass')
    return [max(0,x) for x in out]

def condition(feed,target_M=.05,base_M=1.):
    """Scenario NH3(aq) addition to analytical free-HCl target; no pH prediction."""
    acid=feed[2]/(MW[2]/1000)
    if feed[1]>1e-12 or feed[3]>1e-12:raise ValueError('Conditioner boundary assumes HCl carryover only')
    def at(v):
        base=solution(v,{'Ammonia':base_M}) if v>0 else [0.]*5
        out=mix_mass(feed,base)
        n=min(acid,base[3]/(MW[3]/1000))
        out[2]-=n*MW[2]/1000;out[3]-=n*MW[3]/1000
        out[4]+=n*(MW[2]+MW[3])/1000
        return base,out,(acid-n)/(liquid_volume(out)/1000)
    lo=0.;hi=acid/base_M*1000
    if at(0)[2]<target_M-1e-10:raise ValueError('Initial acid below specified target')
    for k in range(65):
        mid=(lo+hi)/2
        if at(mid)[2]>target_M:lo=mid
        else:hi=mid
    base,out,c=at((lo+hi)/2)
    return base,out,c

def ergun(masses,flow_mL_min,bed_mL,diameter_mm,particle_um=100.,voidage=.4):
    """Clean spherical fixed-bed screening; no resin compression/dispersion claim."""
    rho=liquid_density(masses);mu=liquid_viscosity(masses)
    area=pi*(diameter_mm*.001)**2/4
    length=bed_mL*1e-6/area;u=flow_mL_min*1e-6/60/area;dp=particle_um*1e-6;e=voidage
    dP=length*(150*mu*(1-e)**2*u/(e**3*dp**2)+1.75*rho*(1-e)*u*u/(e**3*dp))
    return dP

def ph_recovery(ph):
    if ph<1.5 or ph>2.5:raise ValueError('Measured pH interpolation restricted to 1.5--2.5')
    if ph<=2:return .95+.06*(ph-1.5)
    return .98-.64*(ph-2)
