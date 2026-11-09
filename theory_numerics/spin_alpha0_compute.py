#!/usr/bin/env python3
"""Independent α=0 recalculation from the displayed formulas (no imported repository scripts).
Recomputes normalized measures, full Lipschitz envelopes, half-strip distance,
explicit upper profile, all covariance shapes, Sachs coordinate audit, delta,
and deterministic two-kink signed witnesses. Finite-basis QCQP/SOS handled separately.
"""
import math,json,sys,time
from pathlib import Path
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import quad
OUT=Path(__file__).resolve().parent / 'reproduced_alpha0';OUT.mkdir(exist_ok=True)
SIGMAS=(.02,.05,.08,.10,.15,.20,.30,.40,.50)
LS=(.25,.5,1.,2.,3.,4.,6.,8.,12.,16.,24.)
DELTAS=(0,.05,.1,.25,.5,.75,.9)
REF=(.10,2.)

def measures(s, n=6001, tail=12, alpha=0):
    k=np.linspace(-tail*s,tail*s,n)
    omega=np.hypot(1.,k);x=2.*(omega-1.)
    rho=np.exp(-.5*(k/s)**2)/(2*omega)
    wm=rho*k*k;wc=rho*k*k*x*(2/(omega+1))**alpha
    return k,x,wm/np.trapezoid(wm,k),wc/np.trapezoid(wc,k)

def box_extrema(x,g,L,c=.25,Delta=0):
    alo,ahi=1-L*x,1+L*x
    blo,bhi=1-g-Delta-L*x,1-g+Delta+L*x
    cand=[]
    for a in (alo,ahi):
        for b in (blo,bhi):cand.append(a*b+c*x*b*b)
        star=np.divide(-a,2*c*x,out=np.zeros_like(x),where=x>1.e-15)
        b=np.clip(star,blo,bhi)
        cand.append(a*b+c*x*b*b)
    return np.minimum.reduce(cand),np.maximum.reduce(cand)

def rect(M,g,L,coord='G1',Delta=0):
    k,x,mu_m,mu_c=M
    mm=float(np.trapezoid(mu_m*np.maximum(abs(g)-L*x,0)**2,k))
    lo,hi=box_extrema(x,g,L,c=(.25 if coord=='G1' else 1./12),Delta=Delta)
    return (mm,float(np.trapezoid(mu_c*lo,k)),float(np.trapezoid(mu_c*hi,k)))

def distance_strip(precision,mm,lo,hi):
    # project onto M>=mm and lo <= C<=hi; reliable 1D candidate and boundaries
    Q=precision
    cases=[]
    # M=mm plus optimal C
    c=float(np.clip(-Q[0,1]*mm/Q[1,1],lo,hi));cases.append((mm,c))
    # C=endpoint plus optimal M
    for c in (lo,hi):cases.append((max(mm,-Q[0,1]*c/Q[0,0]),c))
    if lo<=0<=hi:cases.append((mm,0.))
    return min(math.sqrt(max(0,np.array([m,c])@Q@np.array([m,c]))) for m,c in cases)

def d_class(s=.1,L=2.,ratio=1.,corr=0.,coord='G1',Delta=0,n=6001):
    M=measures(s,n=n)
    cov=np.array([[1.,corr*ratio],[corr*ratio,ratio**2]])
    q=np.linalg.inv(cov)
    def fun(g):return distance_strip(q,*rect(M,g,L,coord,Delta))
    gs=np.linspace(-4,4,161);vals=np.asarray([fun(float(g)) for g in gs]);idx=int(vals.argmin())
    sol=minimize_scalar(fun,bounds=(gs[max(0,idx-1)],gs[min(len(gs)-1,idx+1)]),method='bounded',options={'xatol':1e-10})
    return float(sol.fun),float(sol.x)

def upper(s=.1,L=2.,n=20001):
    k,x,mu_m,mu_c=measures(s,n=n)
    def f(g,coords=False):
        gm=np.sign(g)*np.maximum(abs(g)-L*x,0)
        g1=1-L*x;gq=1-g-L*x
        wm=float(np.trapezoid(mu_m*gm*gm,k))
        wc=float(np.trapezoid(mu_c*(gq*g1+x*gq*gq/4),k))
        return (math.hypot(wm,wc),wm,wc) if coords else math.hypot(wm,wc)
    sol=minimize_scalar(f,bounds=(-4,4),method='bounded',options={'xatol':1e-12})
    d,wm,wc=f(float(sol.x),True)
    return {'distance':d,'g':float(sol.x),'W_M':wm,'W_C':wc}

def envelope(s,L,coord='G1',Delta=0,n=6001):
    M=measures(s,n=n); k,x,_,mu_c=M
    lo,_=box_extrema(x,0.,L,c=.25 if coord=='G1' else 1/12,Delta=Delta)
    return float(np.trapezoid(mu_c*lo,k))

def componentwise_certificate(s,L,n=6001):
    k,x,mu_m,mu_c=measures(s,n=n)
    xm=float(np.trapezoid(mu_m*x,k)); x2m=float(np.trapezoid(mu_m*x*x,k))
    x1=float(np.trapezoid(mu_c*x,k));x2=float(np.trapezoid(mu_c*x*x,k));x3=float(np.trapezoid(mu_c*x**3,k))
    # For all admitted g, either WM zero excluded by magnetic box or WC zero excluded by charge box
    # Exact interval search for worst g via fine global sampling plus smooth subdivision
    def fails(g):
        errm=2*abs(g)*L*xm+L*L*x2m
        errc=(L+abs(1-g)*L)*x1 +L*L*x2+.25*(1-g)**2*x1+.5*abs(1-g)*L*x2+.25*L*L*x3
        return (g*g<=errm) and abs(1-g)<=errc
    gs=np.linspace(-4,4,8001)
    return not any(fails(float(g)) for g in gs)

def x_r(r):return 2*math.expm1(.5*math.log1p(r*r))
def w_r(r,s,alpha=0):
    om=math.hypot(1,r)
    return math.exp(-.5*(r/s)**2)/(2*om)*r*r*x_r(r)*(2/(om+1))**alpha

def signed_witness(s,L,alpha=0):
    X=x_r(9*s);h=X/32
    r_h=math.sqrt(h*(1+h/4));breaks=(0.,r_h,9*s,14*s)
    def weighted(t):
        def integrand(r):
            x=x_r(r);g1=1+t*L*(min(x,X)-2*min(x,h));gq=1-t*L*min(x,X)
            return w_r(r,s,alpha)*(g1*gq+x*gq*gq/4)
        num=sum(quad(integrand,a,b,epsabs=5e-13,epsrel=5e-12,limit=150)[0] for a,b in zip(breaks[:-1],breaks[1:]))
        den=sum(quad(lambda r:w_r(r,s,alpha),a,b,epsabs=5e-13,epsrel=5e-12)[0] for a,b in zip(breaks[:-1],breaks[1:]))
        return num/den
    wc0,wc1=weighted(0.),weighted(1.)
    t=brentq(weighted,0,1,xtol=4e-14) if wc0>0>wc1 else None
    return {'sigma_k_over_m':s,'L':L,'W_C_t0':wc0,'W_C_t1':wc1,'t_root':t,'max_abs_slope':(t*L if t is not None else None)}

def main():
    start=time.time()
    logs={}
    logging=lambda message: print(message,flush=True)
    logging('Calculating alpha0 99-cell envelope and componentwise...')
    fixture=Path(__file__).resolve().parent / 'acceptance_sensitivity_99point.csv'
    old=pd.read_csv(fixture) if fixture.exists() else None
    results=[];counts={}
    for s in SIGMAS:
        for L in LS:
            b=envelope(s,L)
            cb=componentwise_certificate(s,L)
            b_sachs=envelope(s,L,'GC')
            if old is not None:
                oldrow=old[(old.alpha==0)&np.isclose(old.sigma_k_over_m,s)&np.isclose(old.L,L)].iloc[0]
                if abs(b-oldrow.BC)>1.e-5:raise RuntimeError(('BC mismatch',s,L,b,oldrow.BC))
            if b>1e-9:
                status='certified'
            else:
                witness=signed_witness(s,L,alpha=0)
                status='overlap' if witness['W_C_t0']>0>witness['W_C_t1'] else 'unresolved'
            results.append({'sigma_k_over_m':s,'L':L,'B_C':b,'componentwise_certified':cb,'Sachs_B_C':b_sachs,'status':status})
    df=pd.DataFrame(results)
    df.to_csv(OUT/'classification_99_alpha0.csv',index=False)
    logging('classification '+str(df.status.value_counts().to_dict())+' componentwise='+str(df.componentwise_certified.sum())+' Sachs='+str((df.Sachs_B_C>0).sum()))
    logging('Reference whole class distances and bracket...')
    lower,g=d_class(n=20001);up=upper();ls,gs=d_class(coord='GC',n=20001)
    logs['reference']={'D_lower':lower,'g_lower':g,'upper':up,'Sachs_lower':ls,'g_sachs':gs}
    logging(str(logs['reference']))
    logging('Covariance table...')
    rows=[]
    for ratio in (.5,1,2,4):
        for corr in (0,.5,-.5):
            d,g=d_class(ratio=ratio,corr=corr,n=6001)
            rows.append({'r_sigma':ratio,'rho':corr,'lower_bound':d,'epsilon_3sigma':d/3,'epsilon_5sigma':d/5,'g':g})
    pd.DataFrame(rows).to_csv(OUT/'covariance_alpha0.csv',index=False)
    logging('Delta robustness 99x7...')
    delta_rows=[]
    for dlt in DELTAS:
        cert=0;over=33
        for row in results:
            if envelope(row['sigma_k_over_m'],row['L'],Delta=dlt,n=4001)>1e-9:cert+=1
        bd=envelope(.1,2.,Delta=dlt,n=20001);dr,gr=d_class(Delta=dlt,n=6001)
        delta_rows.append({'Delta':dlt,'certified':cert,'overlap':over,'unresolved':99-cert-over,'BC_ref':bd,'D_ref':dr,'g_D_ref':gr})
    pd.DataFrame(delta_rows).to_csv(OUT/'delta_alpha0.csv',index=False)
    logging(str(delta_rows))
    logging('Signed witnesses 33...')
    wrows=[]
    for row in results:
        if row['status']=='overlap':
            v=signed_witness(row['sigma_k_over_m'],row['L']);assert v['W_C_t0']>0>v['W_C_t1'];wrows.append(v)
    pd.DataFrame(wrows).to_csv(OUT/'explicit_overlaps_alpha0.csv',index=False)
    logging('Minimum margin witness: '+str(max(wrows,key=lambda x:x['W_C_t1'])))
    logging('Distance map 99...')
    distances=[]
    for row in results:
        d,g=d_class(row['sigma_k_over_m'],row['L'],n=3001)
        distances.append(d)
    df['distance_lower_bound']=distances
    df.to_csv(OUT/'classification_99_alpha0.csv',index=False)
    logging('99 distance results finished in %.1f s'%(time.time()-start))
    logs['counts']=df.status.value_counts().to_dict()
    logs['componentwise_count_estimate']=int(df.componentwise_certified.sum())
    logs['sachs_count']=int((df.Sachs_B_C>0).sum())
    logs['delta']=delta_rows
    logs['witness_critical']=max(wrows,key=lambda x:x['W_C_t1'])
    logs['runtime_s']=round(time.time()-start,2)
    (OUT/'numerical_summary_alpha0.json').write_text(json.dumps(logs,indent=2))
    return df

if __name__=='__main__':main()
