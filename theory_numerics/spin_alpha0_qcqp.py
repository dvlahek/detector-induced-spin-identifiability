import itertools,time,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.optimize import minimize
OUT=Path(__file__).resolve().parent / 'reproduced_alpha0';OUT.mkdir(exist_ok=True)
TAB=[(.05,24),(.08,12),(.08,16),(.08,24),(.10,6),(.10,8),(.10,12),(.10,16),(.15,3),(.15,4),(.15,6),(.15,8),(.20,2),(.20,3),(.20,4),(.30,1),(.30,2),(.40,.5),(.40,1),(.50,.5),(.50,1)]

def qcqp(s,N=8,alpha=0,grid=2601):
    k=np.linspace(-9*s,9*s,grid);e=np.hypot(1,k);x=2*(e-1);r0=np.exp(-.5*(k/s)**2)/(2*e)
    mu=r0*k*k*x*(2/(e+1))**alpha
    mu/=np.trapezoid(mu,k)
    xmax=float(x.max());knots=np.linspace(0,xmax,N+1);delta=knots[1]-knots[0]
    phi=np.column_stack([np.clip(x-knots[i],0,delta) for i in range(N)])
    c0=float(np.trapezoid(mu*(1+x/4),k))
    au=np.trapezoid(mu[:,None]*phi,k,axis=0)
    av=np.trapezoid(mu[:,None]*phi*(1+x[:,None]/2),k,axis=0)
    H=np.trapezoid(mu[:,None,None]*phi[:,:,None]*phi[:,None,:],k,axis=0)
    Q=np.trapezoid(mu[:,None,None]*(x/4)[:,None,None]*phi[:,:,None]*phi[:,None,:],k,axis=0)
    return c0,au,av,H,.5*(Q+Q.T)

def solve(q,u,L,vstart=None):
    c,au,av,H,Q=q
    linear=av+H.T@u
    def f(v):return float(c+au@u+linear@v+v@Q@v)
    def jac(v):return linear+2*Q@v
    res=minimize(f,np.zeros_like(u) if vstart is None else vstart,jac=jac,method='L-BFGS-B',bounds=[(-L,L)]*len(u),options={'ftol':1e-15,'gtol':1e-12,'maxiter':3000,'maxls':80})
    v=res.x; grad=jac(v)
    kkt=np.max(np.abs(v-np.clip(v-grad,-L,L)))
    return f(v),v,float(kkt)

def global_min(q,L):
    N=len(q[1]);best=(float('inf'),None,None,None)
    for signs in itertools.product((-1.,1.),repeat=N):
        u=L*np.array(signs)
        value,v,kkt=solve(q,u,L)
        if value<best[0]:best=(value,u,v,kkt)
    return {'minimum':best[0],'max_kkt_best':best[3],'u':best[1].tolist(),'v':best[2].tolist()}

def sos(q,L):
    c,au,av,H,Q=q;n=len(au)
    lin=np.r_[au,av]
    R=np.block([[np.zeros((n,n)),.5*H],[.5*H.T,Q]])
    R=(R+R.T)/2
    shift=max(1.e-6,-np.linalg.eigvalsh(R)[0]+1.e-4)
    lam0=np.full(2*n,shift)
    def obj(lam):
        A=R+np.diag(lam)
        mineig=np.linalg.eigvalsh(A)[0]
        if mineig<=1.e-13:return 1.e12+1.e12*(1.e-13-mineig),np.full(2*n,-1.e6)
        y=np.linalg.solve(A,lin)
        value=L*L*np.sum(lam)+.25*lin@y-c
        grad=L*L-.25*y*y
        return float(value),np.asarray(grad)
    try:
        res=minimize(lambda xx:obj(xx)[0],lam0,jac=lambda xx:obj(xx)[1],method='L-BFGS-B',bounds=[(0,None)]*(2*n),options={'ftol':1.e-14,'gtol':1.e-10,'maxiter':10000,'maxls':100})
        lam=np.maximum(res.x,0)
        A=R+np.diag(lam);y=np.linalg.solve(A,lin)
        gamma=float(c-L*L*np.sum(lam)-.25*lin@y)-1.e-10
        S=np.empty((2*n+1,2*n+1));S[0,0]=c-gamma-L*L*np.sum(lam);S[0,1:]=.5*lin;S[1:,0]=.5*lin;S[1:,1:]=A
        return {'gamma':gamma,'minimum_S_eigenvalue':float(np.linalg.eigvalsh(S)[0]),'success':bool(res.success),'lam':lam.tolist()}
    except Exception as e:return {'error':str(e)}

def main():
    rows=[];t0=time.time()
    for i,(s,L) in enumerate(TAB):
        q=qcqp(s)
        v=global_min(q,L)
        sosr=sos(q,L)
        rows.append({'sigma_k_over_m':s,'L':L,'q8_min':v['minimum'],'kkt':v['max_kkt_best'],'SOS_gamma':sosr.get('gamma',float('nan')),'SOS_min_eigenvalue':sosr.get('minimum_S_eigenvalue',float('nan')),'SOS_optimizer_success':sosr.get('success',False)})
        print(i+1,'/',len(TAB),s,L,'q=',round(v['minimum'],8),'sos=',round(sosr.get('gamma',np.nan),8),'in %.1fs'%(time.time()-t0),flush=True)
    df=pd.DataFrame(rows)
    df.to_csv(OUT/'qcqp_sos_alpha0.csv',index=False)
    print('qcqp all positive',bool((df.q8_min>0).all()),'sos positive',int(((df.SOS_gamma>0)&(df.SOS_min_eigenvalue>=-1e-7)).sum()),'runtime',time.time()-t0)
if __name__=='__main__':main()
