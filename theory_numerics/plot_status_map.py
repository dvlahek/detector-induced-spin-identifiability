#!/usr/bin/env python3
"""Recreate the publication finite-bandwidth status map from the archived table."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

ROOT = Path(__file__).resolve().parent
inp = ROOT / "results" / "S1_slope_bandwidth_certificates_monotone.csv"
out = ROOT / "figures" / "S1_status_map_reproduced.png"

df = pd.read_csv(inp)
sigmas = [0.02,0.05,0.08,0.10,0.15,0.20,0.30,0.40,0.50]
Ls = [0.25,0.5,1,2,3,4,6,8,12,16,24]
code = {"explicit_overlap":0, "unresolved":1, "guaranteed":2}
mat = np.empty((len(Ls), len(sigmas)))
for iy,L in enumerate(Ls):
    for ix,s in enumerate(sigmas):
        row = df[(df.sigma_k_over_m==s) & (df.L==L)].iloc[0]
        mat[iy,ix] = code[row.status_monotone]

fig, ax = plt.subplots(figsize=(8.2,5.8))
cmap = ListedColormap([plt.cm.viridis(0.0), plt.cm.viridis(0.5), plt.cm.viridis(1.0)])
norm = BoundaryNorm([-0.5,0.5,1.5,2.5], cmap.N)
im = ax.imshow(mat, origin='lower', aspect='auto', cmap=cmap, norm=norm)
ax.set_xticks(range(len(sigmas)), [str(x) for x in sigmas])
ax.set_yticks(range(len(Ls)), [str(x) for x in Ls])
ax.set_xlabel(r"kernel width $\sigma_k/m$")
ax.set_ylabel(r"slope bound $L$")
ax.set_title("Finite-bandwidth certificate status")
cbar = fig.colorbar(im, ax=ax, ticks=[0,1,2])
cbar.ax.set_yticklabels(["explicit overlap", "unresolved", "guaranteed"])
fig.tight_layout()
fig.savefig(out, dpi=220, bbox_inches='tight')
print(out)
