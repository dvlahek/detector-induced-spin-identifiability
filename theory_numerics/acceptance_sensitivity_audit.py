"""Sensitivity of the sufficient charge certificate and two-kink overlaps
for A_C^(alpha)=(2m/(E_k+m))^alpha. Four fixed architecture choices,
NOT a certificate uniform in alpha. Reproducible with Python+SciPy.
"""
from math import exp, hypot, sqrt
from pathlib import Path
from scipy.integrate import quad
import json, csv

SIGMAS=(.02,.05,.08,.1,.15,.2,.3,.4,.5)
LS=(.25,.5,1,2,3,4,6,8,12,16,24)
ALPHAS=(0,1,2,3)


def x_of_r(r): return 2*(hypot(1,r)-1)

def weight(r,s,alpha):
    omega=hypot(1,r)
    return exp(-0.5*(r/s)**2)/(2*omega)*r*r*x_of_r(r)*(2/(omega+1))**alpha

def envelope(x,L):
    if x<1.e-16: return 1.0
    lo,hi=1-L*x,1+L*x
    return min(a*b+x*b*b/4 for a in (lo,hi)
               for b in (lo,hi,max(lo,min(hi,-2*a/x))))

def witness(x,L,s):
    X=x_of_r(9*s); h=X/32
    a=1+L*(min(x,X)-2*min(x,h))
    b=1-L*min(x,X)
    return a*b+x*b*b/4

def avg(fn,s,L,alpha):
    X=x_of_r(9*s); h=X/32
    r_h=sqrt(h*(1+h/4))
    intervals=(0,r_h,9*s,14*s)
    numerator=sum(quad(lambda r:weight(r,s,alpha)*fn(x_of_r(r),L,s), a,b,
                       epsabs=2e-12,epsrel=2e-12)[0]
                  for a,b in zip(intervals[:-1],intervals[1:]))
    denominator=quad(lambda r:weight(r,s,alpha),0,14*s,
                     epsabs=2e-12,epsrel=2e-12)[0]
    return numerator/denominator

def run(outdir=Path(__file__).parent):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    rows=[]; summary=[]; original_overlaps=None
    for alpha in ALPHAS:
        positive=overlap=unresolved=0; overlaps=set()
        for s in SIGMAS:
            for L in LS:
                bc=avg(lambda x,L,s: envelope(x,L),s,L,alpha)
                wc=avg(witness,s,L,alpha)
                if bc>1e-9: status='certified';positive+=1
                elif wc<-1e-9: status='overlap';overlap+=1;overlaps.add((s,L))
                else: status='unresolved';unresolved+=1
                rows.append({'alpha':alpha,'sigma_k_over_m':s,'L':L,
                             'BC':f'{bc:.13g}','W_C_t1':f'{wc:.13g}','status':status})
        assert positive+overlap+unresolved==99
        if original_overlaps is None: original_overlaps=overlaps
        assert overlaps==original_overlaps
        ref=next(r for r in rows if r['alpha']==alpha and r['sigma_k_over_m']==.5 and r['L']==1)
        summary.append({'alpha':alpha,'certified':positive,'overlap':overlap,
                        'unresolved':unresolved,'B_C_0p5_L1':float(ref['BC'])})
    assert [(r['certified'],r['overlap'],r['unresolved']) for r in summary]==[(63,33,3),(63,33,3),(64,33,2),(64,33,2)]
    with open(outdir/'acceptance_sensitivity_99point.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    (outdir/'acceptance_sensitivity_summary.json').write_text(json.dumps({'summary':summary,'method':'fixed-alpha positive acceptance, knot-split quadrature, full-Lipschitz sufficient envelope, two-kink constructive witness'},indent=2))
    print('PASS acceptance sensitivity:',summary)

if __name__=='__main__':run()
