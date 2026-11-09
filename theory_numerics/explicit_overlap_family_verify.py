"""Independent, deterministic whole-grid two-kink overlap witness verification.
Model: accepted charge measure rho0 * (r/(omega+1))**2 * x/4,
G_M=0, G_Q=1-t L min(x,X), G_1=1+t L[min(x,X)-2 min(x,h)],
X=2(sqrt(1+(9 sigma)^2)-1), h=X/32.
If W_C(0)>0>W_C(1), IVT proves a continuous-class exact overlap.
"""
from scipy.integrate import quad
from scipy.optimize import brentq
import json,math,csv
from pathlib import Path
SIGMAS=(.02,.05,.08,.1,.15,.2,.3,.4,.5)
LS=(.25,.5,1,2,3,4,6,8,12,16,24)
OUT=Path(__file__).parent

def x_of_r(r):
    return 2*math.expm1(0.5*math.log1p(r*r))
def weight(r,s):
    om=math.hypot(1,r)
    x=x_of_r(r)
    return math.exp(-0.5*(r/s)**2)/(2*om)*(r/(om+1))**2*(x/4)

def profiles(x,s,L,t):
    X=x_of_r(9*s);h=X/32
    return (1+t*L*(min(x,X)-2*min(x,h)),1-t*L*min(x,X))

def brackets(s):
    X=x_of_r(9*s);h=X/32
    r_h=math.sqrt(h*(1+h/4));r_X=9*s
    return (0,r_h,r_X,14*s)

def integ(f,s):
    bb=brackets(s)
    return sum(quad(f,a,b,epsabs=5e-13,epsrel=5e-12,limit=150)[0] for a,b in zip(bb[:-1],bb[1:]))

def WC(s,L,t):
    den=integ(lambda r:weight(r,s),s)
    num=integ(lambda r:weight(r,s)*((lambda a,b,x:a*b+x*b*b/4)(*profiles(x_of_r(r),s,L,t),x_of_r(r))),s)
    return num/den

def main():
    # prior exact status grid from manuscript; classification labels include only 2 unresolved
    red=[]
    certified=[]
    unresolved={(0.08,24),(0.10,16)}
    # exact formulas for pointwise envelope to discover candidate labels independently
    def h(x,L):
        a0,a1=1-L*x,1+L*x
        lo,hi=1-L*x,1+L*x
        if x<1e-13:
            return lo
        v=math.inf
        for a in (a0,a1):
            b=max(lo,min(hi,-2*a/x))
            for bb in (lo,hi,b):
                v=min(v,a*bb+x*bb*bb/4)
        return v
    def BC(s,L):
        den=integ(lambda r:weight(r,s),s)
        return integ(lambda r:weight(r,s)*h(x_of_r(r),L),s)/den
    for s in SIGMAS:
        for L in LS:
            b=BC(s,L)
            if b>0:certified.append((s,L,b))
            elif (s,L) not in unresolved:
                w0=WC(s,L,0);w1=WC(s,L,1)
                root=brentq(lambda t:WC(s,L,t),0,1,xtol=4e-14) if w0>0>w1 else float('nan')
                red.append({'sigma_k_over_m':s,'L':L,'B_C':b,'W_C_t0':w0,'W_C_t1':w1,'t_root':root,'max_abs_slope':root*L,'X':x_of_r(9*s),'h':x_of_r(9*s)/32})
    assert len(certified)==64,(len(certified),'expected 64')
    assert len(red)==33,(len(red),'expected 33')
    assert all(r['W_C_t0']>0>r['W_C_t1'] and 0<r['t_root']<1 for r in red)
    minimal=min(red,key=lambda q:-q['W_C_t1'])
    data={'grid':{'certified':len(certified),'constructive_overlap':len(red),'unresolved':2},
          'witness':{'G_M':'0','G_1':'1 + t L(min(x,X)-2 min(x,h))','G_Q':'1 - t L min(x,X)','X':'2(sqrt(1+(9 sigma)^2)-1)','h':'X/32'},
          'least_negative_target':minimal,'overlaps':red}
    (OUT/'explicit_overlap_family_results.json').write_text(json.dumps(data,indent=2))
    with open(OUT/'explicit_overlap_family_results.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(red[0]));w.writeheader();w.writerows(red)
    print('certified',len(certified),'overlap',len(red),'unresolved',2)
    print('LEAST_NEGATIVE',minimal)
    print('MIN_NEGATIVE',min(r['W_C_t1'] for r in red))

if __name__=='__main__':main()
