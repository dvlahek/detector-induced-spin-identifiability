"""Independent knot-split quadrature of saved 48-segment alpha=0 profiles."""
from pathlib import Path
import json,math,csv
import numpy as np
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parent
cases=json.loads((ROOT/'results'/'alpha0'/'alpha0_48segments_profiles.json').read_text())['cases']
rows=[]
for a in cases:
 sig=float(a['sigma_k_over_m']);L=float(a['L']);N=int(a['n_segments'])
 u=np.array(a['u']);v=np.array(a['v']);k9=9*sig
 xmax=2*(math.sqrt(1+k9*k9)-1)
 knot=np.linspace(0,xmax,N+1)
 split=[0.0]+[math.sqrt((1+z/2)**2-1) for z in knot[1:-1] if math.sqrt((1+z/2)**2-1)<12*sig]+[12*sig]
 def value(k,mode):
  e=math.sqrt(1+k*k);x=2*(e-1)
  phi=np.clip(x-knot[:-1],0,knot[1]-knot[0])
  g1=1+float(phi@u);gq=1+float(phi@v)
  weight=math.exp(-.5*(k/sig)**2)/(2*e)*k*k*x
  if mode=='numerator':return weight*(g1*gq+x*gq*gq/4)
  return weight
 def integ(mode):
  acc=0
  for k0,k1 in zip(split[:-1],split[1:]):
   acc+=quad(lambda k:value(k,mode),k0,k1,epsabs=1e-13,epsrel=2e-10,limit=100)[0]
  return acc
 wc=integ('numerator')/integ('denominator')
 row={'case':a['case'],'sigma_k_over_m':sig,'L':L,'W_C_quad':wc,'W_C_grid':a['W_C_n20001'],'difference':wc-a['W_C_n20001'],'sign_ok':(wc<0)==(a['W_C_n20001']<0),'slopes_within_L':a['all_slope_bounds_ok']}
 rows.append(row)
 print(row)
assert all(row['sign_ok'] and row['slopes_within_L'] for row in rows)
output=ROOT/'reproduced_alpha0'/'alpha0_48segments_quadrature.csv'
output.parent.mkdir(parents=True,exist_ok=True)
with output.open('w',newline='') as out:
 writer=csv.DictWriter(out,rows[0].keys());writer.writeheader();writer.writerows(rows)

# Independently validate the cited saved profile integrals.
reference=ROOT/'results'/'alpha0'/'alpha0_48segments_quadrature.csv'
with reference.open(newline='') as handle:
 frozen=list(csv.DictReader(handle))
assert len(rows)==len(frozen)==4
for row,record in zip(rows,frozen):
 assert row['case']==record['case']
 assert abs(row['W_C_quad']-float(record['W_C_quad']))<1e-7
 assert row['sign_ok'] and row['slopes_within_L']
print('PASS alpha0 48-segment frozen-profile quadrature')
