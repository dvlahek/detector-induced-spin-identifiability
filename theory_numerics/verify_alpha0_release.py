#!/usr/bin/env python3
"""Check the frozen primary alpha=0 data against independently regenerated outputs.
Run after spin_alpha0_compute.py, spin_alpha0_qcqp.py, spin_alpha0_power.py.
The audit checks numerical results, not validated interval proofs.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

BASE=Path(__file__).resolve().parent
FROZEN=BASE/'results'/'alpha0'
ACTUAL=BASE/'reproduced_alpha0'
ATOL={'classification_99_alpha0.csv':1e-7,'delta_alpha0.csv':2e-7,
      'covariance_alpha0.csv':2e-7,'explicit_overlaps_alpha0.csv':3e-7,
      'qcqp_sos_alpha0.csv':1e-5,'estimator_power_alpha0.csv':0.02}

def verify_pair(filename, atol):
    expected=pd.read_csv(FROZEN/filename)
    produced=pd.read_csv(ACTUAL/filename)
    assert expected.shape==produced.shape,(filename,expected.shape,produced.shape)
    assert list(expected.columns)==list(produced.columns),filename
    for column in expected.columns:
        if pd.api.types.is_numeric_dtype(expected[column]):
            if filename=='qcqp_sos_alpha0.csv' and column=='SOS_min_eigenvalue':
                # Semidefinite margins may be near machine precision and solver dependent.
                assert np.all(np.isfinite(produced[column])),(filename,column)
            else:
                assert np.allclose(expected[column],produced[column],rtol=0.01,atol=atol,equal_nan=True),(filename,column)
        else:
            assert expected[column].astype(str).tolist()==produced[column].astype(str).tolist(),(filename,column)

for n,tol in ATOL.items(): verify_pair(n,tol)
cls=pd.read_csv(ACTUAL/'classification_99_alpha0.csv')
assert cls.status.value_counts().to_dict()=={'certified':63,'overlap':33,'unresolved':3}
assert int(cls.componentwise_certified.sum())==46
assert int((cls.Sachs_B_C>0).sum())==62
D=pd.read_csv(ACTUAL/'delta_alpha0.csv').sort_values('Delta')
assert D.certified.astype(int).tolist()==[63,61,61,59,55,45,31]
Q=pd.read_csv(ACTUAL/'qcqp_sos_alpha0.csv')
assert len(Q)==21 and Q.q8_min.gt(0).all()
assert int(((Q.SOS_gamma>0)&(Q.SOS_min_eigenvalue>=-1e-7)).sum())>=12
w=pd.read_csv(ACTUAL/'explicit_overlaps_alpha0.csv')
assert len(w)==33 and w.W_C_t0.gt(0).all() and w.W_C_t1.lt(0).all()
sos=json.loads((FROZEN/'sos_witness_alpha0.json').read_text())
assert len(sos['records'])==21 and sos['positive_feasible']>=12
print('PASS: alpha=0 frozen numerical tables and signed witnesses reproduced.')
