#!/usr/bin/env python3
"""Robustness audit for g_M + g_Q = 1 + delta.

The audit treats delta as an additional nuisance parameter with |delta|<=Delta.
It preserves the detector measures and anchored Lipschitz class used by the
manuscript. A uniform whole-class certificate follows from the enlarged
pointwise interval

    G_Q(x) in [1-Delta-Lx, 1+Delta+Lx]

once magnetic positivity has reduced a possible protected-origin overlap to
G_M=0 and g_M=0.

The script also checks that the compact-prior magnetic margins are unchanged,
computes the static protected distance, and evaluates the unit-covariance lower
enclosure at the manuscript reference point.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize_scalar

SIGMAS=(0.02,0.05,0.08,0.10,0.15,0.20,0.30,0.40,0.50)
SLOPES=(0.25,0.5,1.0,2.0,3.0,4.0,6.0,8.0,12.0,16.0,24.0)
DELTAS=(0.0,0.05,0.10,0.25,0.50,0.75,0.90)

def detector_measures(sigma: float, grid_n: int=20001, tail: float=12.0):
    k=np.linspace(-tail*sigma,tail*sigma,grid_n)
    e=np.sqrt(1.0+k*k)
    x=2.0*(e-1.0)
    rho0=np.exp(-0.5*(k/sigma)**2)/(2.0*e)
    rho_m=rho0*k*k
    rho_c=rho0*(k/(e+1.0))**2*(x/4.0)
    return k,x,rho_m/np.trapezoid(rho_m,k),rho_c/np.trapezoid(rho_c,k)

def average(k: np.ndarray, density: np.ndarray, values: np.ndarray) -> float:
    return float(np.trapezoid(density*values,k))

def robust_charge_bounds(x: np.ndarray, g: float, slope: float,
                         Delta: float) -> tuple[np.ndarray,np.ndarray]:
    a_lo,a_hi=1.0-slope*x,1.0+slope*x
    b_lo=1.0-g-Delta-slope*x
    b_hi=1.0-g+Delta+slope*x
    candidates=[]
    for a in (a_lo,a_hi):
        candidates.append(a*b_lo+0.25*x*b_lo*b_lo)
        candidates.append(a*b_hi+0.25*x*b_hi*b_hi)
        stationary=np.zeros_like(x)
        nonzero=x>1.0e-15
        np.divide(-2.0*a,x,out=stationary,where=nonzero)
        stationary=np.clip(stationary,b_lo,b_hi)
        candidates.append(a*stationary+0.25*x*stationary*stationary)
    candidates.append(np.where((b_lo<=0.0)&(b_hi>=0.0),0.0,np.nan))
    stack=np.asarray(candidates)
    return np.nanmin(stack,axis=0),np.nanmax(stack,axis=0)

def robust_origin_envelope(sigma: float, slope: float, Delta: float,
                           grid_n: int=20001) -> float:
    k,x,_,nu_c=detector_measures(sigma,grid_n=grid_n)
    lower,_=robust_charge_bounds(x,g=0.0,slope=slope,Delta=Delta)
    return average(k,nu_c,lower)

def magnetic_margin(sigma: float, slope: float, g0: float=1.0,
                    grid_n: int=20001) -> float:
    k,x,nu_m,_=detector_measures(sigma,grid_n=grid_n)
    return average(k,nu_m,np.maximum(g0-slope*x,0.0)**2)

def static_distance(Delta: float) -> tuple[float,float,float]:
    def objective(g: float) -> float:
        delta=float(np.clip(g-1.0,-Delta,Delta))
        charge=1.0+delta-g
        return math.sqrt(g**4+charge**2)
    res=minimize_scalar(objective,bounds=(-4.0,4.0),method='bounded',
                        options={'xatol':1.0e-13})
    g=float(res.x)
    d=float(np.clip(g-1.0,-Delta,Delta))
    return float(res.fun),g,d

def half_strip_distance(magnetic: float, c_lo: float, c_hi: float) -> float:
    charge=0.0 if c_lo<=0.0<=c_hi else min((c_lo,c_hi),key=abs)
    return math.hypot(magnetic,charge)

def reference_distance(Delta: float, sigma: float=0.10, slope: float=2.0,
                       grid_n: int=20001) -> tuple[float,float]:
    k,x,nu_m,nu_c=detector_measures(sigma,grid_n=grid_n)
    def objective(g: float) -> float:
        magnetic=average(k,nu_m,np.maximum(abs(g)-slope*x,0.0)**2)
        lo,hi=robust_charge_bounds(x,g,slope,Delta)
        c_lo=average(k,nu_c,lo)
        c_hi=average(k,nu_c,hi)
        return half_strip_distance(magnetic,c_lo,c_hi)
    anchors=np.linspace(-4.0,4.0,401)
    vals=np.asarray([objective(float(g)) for g in anchors])
    i=int(np.argmin(vals))
    lo=anchors[max(0,i-2)]
    hi=anchors[min(len(anchors)-1,i+2)]
    res=minimize_scalar(objective,bounds=(lo,hi),method='bounded',
                        options={'xatol':1.0e-11})
    return float(res.fun),float(res.x)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--status',type=Path,required=True)
    ap.add_argument('--magnetic-reference',type=Path,required=True)
    ap.add_argument('--outdir',type=Path,required=True)
    args=ap.parse_args()
    args.outdir.mkdir(parents=True,exist_ok=True)

    old=pd.read_csv(args.status)
    overlap={(float(r.sigma_k_over_m),float(r.L))
             for r in old.itertuples(index=False)
             if r.status_publication=='explicit_overlap'}

    grid=[]
    for Delta in DELTAS:
        for sigma in SIGMAS:
            for slope in SLOPES:
                bound=robust_origin_envelope(sigma,slope,Delta)
                if bound>1.0e-10:
                    status='guaranteed'
                elif (sigma,slope) in overlap:
                    status='explicit_overlap'
                else:
                    status='unresolved'
                grid.append({'Delta':Delta,'sigma_k_over_m':sigma,'L':slope,
                             'robust_charge_lower_bound':bound,'status':status})
    frame=pd.DataFrame(grid)
    frame.to_csv(args.outdir/'delta_robustness_grid.csv',index=False)

    summaries=[]
    for Delta,sub in frame.groupby('Delta',sort=True):
        counts=sub.status.value_counts().to_dict()
        ref=sub[(sub.sigma_k_over_m==0.10)&(sub.L==2.0)].iloc[0]
        sd,sg,sdelta=static_distance(float(Delta))
        rd,rg=reference_distance(float(Delta))
        summaries.append({
            'Delta':float(Delta),
            'guaranteed':int(counts.get('guaranteed',0)),
            'explicit_overlap':int(counts.get('explicit_overlap',0)),
            'unresolved':int(counts.get('unresolved',0)),
            'reference_charge_lower_bound':float(ref.robust_charge_lower_bound),
            'reference_distance_lower_bound':rd,
            'reference_distance_anchor_gM':rg,
            'static_min_distance':sd,
            'static_min_gM':sg,
            'static_min_delta':sdelta,
        })
    summary=pd.DataFrame(summaries)
    summary.to_csv(args.outdir/'delta_robustness_summary.csv',index=False)

    threshold=float(brentq(lambda d: robust_origin_envelope(0.10,2.0,d,grid_n=10001),
                            0.75,0.95,xtol=1.0e-10))

    ref=pd.read_csv(args.magnetic_reference)
    computed=[magnetic_margin(float(r.sigma_k_over_m),float(r.L))
              for r in ref.itertuples(index=False)]
    value_col='magnetic_only_lower_bound'
    err=np.max(np.abs(np.asarray(computed)-ref[value_col].to_numpy(dtype=float)))

    payload={
        'relation':'g_M + g_Q = 1 + delta',
        'standard_TGC_mapping':'delta = 2 lambda_gamma',
        'Delta_values':list(DELTAS),
        'reference_point':{'sigma_k_over_m':0.10,'L':2.0},
        'reference_certificate_zero_at_Delta':threshold,
        'magnetic_prior_K':'[1,3]',
        'magnetic_margin_max_abs_regression_error':float(err),
        'delta_zero_status_counts':summary.iloc[0][['guaranteed','explicit_overlap','unresolved']].to_dict(),
    }
    (args.outdir/'delta_robustness_summary.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(summary.to_string(index=False))
    print(json.dumps(payload,indent=2))

if __name__=='__main__':
    main()
