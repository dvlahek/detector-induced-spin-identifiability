"""Re-run published estimator-space protocol with the unfiltered longitudinal charge measure.
Uses identical random seed, cloud size, basis, and MC conventions as v1.8.0.
The explicit distance upper profile is recomputed from new reference measure.
"""
from pathlib import Path
import argparse
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
OUT=Path(__file__).resolve().parent / 'reproduced_alpha0'; OUT.mkdir(parents=True,exist_ok=True)
SEED=20260816;SIGMA=.1;L=2.;N=8;GRID_N=2601;KMAX_SIGMA=9.;INFERENCE=30000;TRUTH=10000;MC=8000
EPS=(.5,.3,.2,.15)

def spectral(alpha=0):
    k=np.linspace(-KMAX_SIGMA*SIGMA,KMAX_SIGMA*SIGMA,GRID_N)
    omega=np.hypot(1,k);x=2*(omega-1)
    rho=np.exp(-.5*(k/SIGMA)**2)/(2*omega)
    m=rho*k*k;c=rho*k*k*x*(2/(omega+1))**alpha
    return k,x,m/np.trapezoid(m,k),c/np.trapezoid(c,k)

def moments(spec):
    k,x,m,c=spec
    xmax=x.max();knots=np.linspace(0,xmax,N+1);d=knots[1]-knots[0]
    phi=np.column_stack([np.clip(x-knots[i],0,d) for i in range(N)])
    mm=np.trapezoid(m[:,None]*phi,k,axis=0)
    hmm=np.trapezoid(m[:,None,None]*phi[:,:,None]*phi[:,None,:],k,axis=0)
    cm=np.trapezoid(c[:,None]*phi,k,axis=0)
    cx=float(np.trapezoid(c*x,k))
    cxm=np.trapezoid(c[:,None]*x[:,None]*phi,k,axis=0)
    hc=np.trapezoid(c[:,None,None]*phi[:,:,None]*phi[:,None,:],k,axis=0)
    hxc=np.trapezoid(c[:,None,None]*x[:,None,None]*phi[:,:,None]*phi[:,None,:],k,axis=0)
    return mm,hmm,cm,cx,cxm,hc,hxc

def cloud(n,rng,mom):
    # Same RNG calling pattern as original; vectorized only during bilinear evaluation
    pars=np.empty((n,1+3*N))
    for i in range(n):
        pars[i,0]=rng.uniform(-4,4)
        pars[i,1:]=rng.uniform(-L,L,size=3*N)
    g=pars[:,0];u=pars[:,1:1+N];m=pars[:,1+N:1+2*N];v=pars[:,1+2*N:]
    mm,hmm,cm,cx,cxm,hc,hxc=mom
    wm=g*g+2*g*(m@mm)+np.einsum('ni,ij,nj->n',m,hmm,m)
    b=1-g
    wc=b+.25*cx*b*b+b*(u@cm)+v@cm+.5*b*(v@cxm)+np.einsum('ni,ij,nj->n',u,hc,v)+.25*np.einsum('ni,ij,nj->n',v,hxc,v)
    return np.column_stack([wm,wc])

def auc(null,alt):
    z=np.r_[null,alt];ranks=pd.Series(z).rank(method='average').to_numpy();n0=len(null);n1=len(alt)
    return (ranks[n0:].sum()-n1*(n1+1)/2)/(n0*n1)

def run(alpha=0):
    spec=spectral(alpha);mom=moments(spec)
    rng=np.random.default_rng(SEED+3333)
    inf=cloud(INFERENCE,rng,mom);tree=cKDTree(inf)
    truth=cloud(TRUTH,rng,mom)
    near_rng=np.random.default_rng(SEED+4444)
    boundary=pd.read_json(OUT/'numerical_summary_alpha0.json',typ='series')['reference']['upper'] if alpha==0 else None
    if boundary is None:
        boundary_f=np.array([0.2677438808009649,0.30289269495740834])
    else: boundary_f=np.array([boundary['W_M'],boundary['W_C']])
    print('near boundary',boundary_f,flush=True)
    rows=[]
    for epsilon in EPS:
        null=np.empty(MC);alt=np.empty(MC)
        for t in range(MC):
            null_obs=rng.multivariate_normal(np.zeros(2),epsilon**2*np.eye(2));null_dist=tree.query(null_obs,k=1)[0]
            null[t]=null_obs@null_obs-null_dist**2
            true=truth[rng.integers(0,len(truth))]
            obs=true+rng.multivariate_normal(np.zeros(2),epsilon**2*np.eye(2));alt_dist=tree.query(obs,k=1)[0]
            alt[t]=obs@obs-alt_dist**2
        cut=float(np.quantile(null,.95))
        near_obs=boundary_f+near_rng.multivariate_normal(np.zeros(2),epsilon**2*np.eye(2),size=MC)
        origin2=np.einsum('ij,ij->i',near_obs,near_obs)
        model_dist=tree.query(near_obs,k=1)[0];score=origin2-model_dist**2
        rows.append(dict(epsilon=epsilon,cloud_averaged_spin1_power=float(np.mean(alt>cut)),near_boundary_profile_power=float(np.mean(score>cut)),cloud_averaged_roc_auc=float(auc(null,alt))))
        print(rows[-1],flush=True)
    return pd.DataFrame(rows)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--alpha',type=int,default=0);a=ap.parse_args()
    frame=run(alpha=a.alpha);frame.to_csv(OUT/f'estimator_power_alpha{a.alpha}.csv',index=False)
