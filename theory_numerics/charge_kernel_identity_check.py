"""Self-test of the chosen positive charge acceptance and its exact normalized density.
The raw longitudinal-current rank-2 coefficient is (e^2*k^2*x/3)
(G1*GQ+x/4*GQ^2). After acceptance A_C=[2m/(E+m)]^2
its nonnegative k-weight is exactly 16 times the pre-normalization
weight printed in the paper (the constant 16 cancels after normalization).
"""
import math,json
from pathlib import Path
rs=[0.001,0.03,0.1,0.3,1.0,3.0]
sigma=.15
items=[]
for r in rs:
    om=math.hypot(1,r)
    x=2*(om-1)
    rho=math.exp(-r*r/(2*sigma*sigma))/(2*om)
    A=(2/(om+1))**2
    paper=rho*(r/(om+1))**2*(x/4)
    accepted=rho*(r*r*x)*A
    ratio=accepted/paper
    assert math.isclose(ratio,16,rel_tol=1e-12),ratio
    assert 0<A<=1
    G1=1.2;GQ=-0.3
    J_long=r*(G1+x*GQ/2)
    J_perp=r*G1
    difference=J_long**2-J_perp**2
    predicted=r*r*x*(G1*GQ+x*GQ*GQ/4)
    assert math.isclose(difference,predicted,abs_tol=1e-13,rel_tol=1e-10)
    items.append({'r':r,'acceptance':A,'accepted_over_paper':ratio,'effect_anisotropy':difference,'predicted':predicted})
result={'PASS':True,'m':1,'sigma_k_over_m':sigma,
        'statement':'A_C=(2/(omega+1))^2, 0<A_C<=1, weighted raw charge measure =16 times paper numerator',
        'samples':items}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2))
print('PASS positive acceptance and density identity',len(rs),'points')
if __name__=='__main__':pass
