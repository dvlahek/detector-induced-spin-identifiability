#!/usr/bin/env python3
"""Seeded high-resolution adversarial search for the primary alpha=0 charge response.
Numerical search is not a whole-class separation proof. Every saved negative profile
is checked by direct independently integrated spectral response.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.integrate import quad

CASES=[(0.15,8.0,'known_overlap'),(0.08,24.0,'unresolved'),(0.10,16.0,'unresolved'),(0.50,1.0,'unresolved')]

def matrices(sigma,nseg=48,grid=2601,tail=9.0):
    k=np.linspace(-tail*sigma,tail*sigma,grid)
    e=np.hypot(1.0,k)
    x=2.0*(e-1.0)
    rho=np.exp(-0.5*(k/sigma)**2)/(2*e)*k**2*x
    rho/=np.trapezoid(rho,k)
    xmax=float(x.max())
    knots=np.linspace(0,xmax,nseg+1)
    delta=knots[1]-knots[0]
    phi=np.clip(x[:,None]-knots[:-1],0,delta)
    weighted=rho[:,None]*phi
    au=np.trapezoid(weighted,k,axis=0)
    av=np.trapezoid(weighted*(1+x[:,None]/2),k,axis=0)
    H=np.empty((nseg,nseg));Q=np.empty((nseg,nseg))
    for i in range(nseg):
        H[i,:]=np.trapezoid(weighted[:,i,None]*phi,k,axis=0)
        Q[i,:]=np.trapezoid(weighted[:,i,None]*phi*x[:,None]/4,k,axis=0)
    Q=(Q+Q.T)/2
    c=float(np.trapezoid(rho*(1+x/4),k))
    return c,au,av,H,Q,knots

def obj_q(q,u,v):
    c,a,b,H,Q,_=q
    return float(c+a@u+b@v+u@H@v+v@Q@v)

def optimize_case(sigma,L,nseg=48,starts=200,seed=20261109,max_cycles=60):
    q=matrices(sigma,nseg)
    c,a,b,H,Q,knots=q
    rng=np.random.default_rng(seed)
    candidates=[np.zeros(nseg),np.full(nseg,L),np.full(nseg,-L)]
    candidates.extend(rng.uniform(-L,L,size=(starts,nseg)))
    best=(np.inf,None,None,np.inf,None)
    st=time.time()
    for j,v0 in enumerate(candidates):
        v=v0.copy(); prev_u=None
        for it in range(max_cycles):
            u=np.where(a+H@v>=0,-L,L)
            linear=b+H.T@u
            def f(z):return c+a@u+linear@z+z@Q@z
            def jac(z):return linear+2*Q@z
            result=minimize(f,v,method='L-BFGS-B',jac=jac,bounds=[(-L,L)]*nseg,
                options={'ftol':1e-12,'gtol':1e-9,'maxiter':400,'maxls':40})
            vn=np.asarray(result.x)
            if prev_u is not None and np.array_equal(u,prev_u) and np.max(np.abs(vn-v))<2e-7:
                v=vn;break
            prev_u=u.copy();v=vn
        val=obj_q(q,u,v)
        grad=b+H.T@u+2*Q@v
        kkt=float(np.max(abs(v-np.clip(v-grad,-L,L))))
        if val<best[0]:best=(val,u.copy(),v.copy(),kkt,j)
    value,u,v,kkt,index=best
    # Re-evaluate best piecewise linear profiles on a finer independent grid.
    def response(n,tail=12):
        kk=np.linspace(-tail*sigma,tail*sigma,n)
        ee=np.hypot(1,kk); xx=2*(ee-1)
        wt=np.exp(-.5*(kk/sigma)**2)/(2*ee)*kk**2*xx
        ph=np.clip(xx[:,None]-knots[:-1],0,knots[1]-knots[0]); g1=1+ph@u;gq=1+ph@v
        return float(np.trapezoid(wt*(g1*gq+xx*gq*gq/4),kk)/np.trapezoid(wt,kk))
    checks=[response(3001),response(9001),response(20001)]
    payload={'sigma_k_over_m':sigma,'L':L,'n_segments':nseg,'random_starts':starts,'seed':seed,
             'best_charge_grid':value,'W_C_n3001':checks[0],'W_C_n9001':checks[1],'W_C_n20001':checks[2],
             'projected_v_KKT_residual':kkt,'best_start_index':index,'max_abs_u':float(max(abs(u))),
             'max_abs_v':float(max(abs(v))),'all_slope_bounds_ok':bool(max(abs(u).max(),abs(v).max())<=L+1e-9),
             'xmax_9sigma':float(knots[-1]),'u':u.tolist(),'v':v.tolist(),'elapsed_seconds':round(time.time()-st,2)}
    return payload

def run(out,starts):
    out.mkdir(parents=True,exist_ok=True)
    rows=[];allcases=[]
    for i,(s,L,label) in enumerate(CASES):
        v=optimize_case(s,L,starts=starts,seed=20261109+97*i)
        v['case']=label
        print(i+1,len(CASES),s,L,'value=',v['W_C_n20001'],'sec=',v['elapsed_seconds'],flush=True)
        rows.append({k:v[k] for k in ('case','sigma_k_over_m','L','n_segments','random_starts','seed','best_charge_grid','W_C_n3001','W_C_n9001','W_C_n20001','projected_v_KKT_residual','max_abs_u','max_abs_v','all_slope_bounds_ok')})
        allcases.append(v)
    pd.DataFrame(rows).to_csv(out/'alpha0_48segments.csv',index=False)
    (out/'alpha0_48segments_profiles.json').write_text(json.dumps({'cases':allcases},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=Path('theory_numerics/reproduced_alpha0'));p.add_argument('--starts',type=int,default=200)
    args=p.parse_args();run(args.out,args.starts)
