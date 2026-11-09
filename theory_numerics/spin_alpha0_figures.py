import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm,LogNorm
from pathlib import Path
p=Path(__file__).resolve().parent / 'reproduced_alpha0';p.mkdir(parents=True,exist_ok=True);df=pd.read_csv(p/'classification_99_alpha0.csv')
sigmas=sorted(df.sigma_k_over_m.unique());Ls=sorted(df.L.unique())
cats=np.empty((len(Ls),len(sigmas)),int);dists=np.empty_like(cats,dtype=float);sachs=np.empty_like(cats,dtype=bool)
for _,r in df.iterrows():
 i=Ls.index(r.L);j=sigmas.index(r.sigma_k_over_m)
 cats[i,j]={'certified':0,'overlap':2,'unresolved':1}[r.status]
 dists[i,j]=r.distance_lower_bound
 sachs[i,j]=r.Sachs_B_C>0
assert np.sum(cats==0)==63 and np.sum(cats==2)==33 and np.sum(cats==1)==3
fig,ax=plt.subplots(figsize=(8.2,4.2),dpi=180)
cm=ListedColormap(['#4ea87c','#d6d9dc','#d96969']);norm=BoundaryNorm([-.5,.5,1.5,2.5],3)
im=ax.imshow(cats,origin='lower',aspect='auto',cmap=cm,norm=norm)
ax.set_xticks(range(len(sigmas)),[f'{s:g}' for s in sigmas]);ax.set_yticks(range(len(Ls)),[f'{l:g}' for l in Ls]);ax.set_xlabel(r'Gaussian width $\sigma_k/m$');ax.set_ylabel(r'Lipschitz bound $L$');ax.set_title('Unfiltered charge response: full-class classification')
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor='#4ea87c',label='whole-class certified'),Patch(facecolor='#d96969',label='explicit overlap'),Patch(facecolor='#d6d9dc',label='unresolved')],loc='upper left',fontsize=8,framealpha=.95)
fig.tight_layout();fig.savefig(p/'figure2a.png',dpi=180);plt.close(fig)
# Numerical distance lower enclosure is zero for cells outside certified set.
m=np.ma.masked_less_equal(dists,1e-10)
fig,ax=plt.subplots(figsize=(8.2,4.2),dpi=180)
im=ax.imshow(m,origin='lower',aspect='auto',cmap='viridis',norm=LogNorm(vmin=1e-9,vmax=max(1.,float(dists.max()))))
ax.set_facecolor('#d6d9dc');ax.set_xticks(range(len(sigmas)),[f'{s:g}' for s in sigmas]);ax.set_yticks(range(len(Ls)),[f'{l:g}' for l in Ls]);ax.set_xlabel(r'Gaussian width $\sigma_k/m$');ax.set_ylabel(r'Lipschitz bound $L$');ax.set_title(r'Whole-class lower distance $D_L(I)$ ($\alpha=0$)')
c=fig.colorbar(im,ax=ax,shrink=.86);c.set_label(r'$D_L(I)$')
fig.tight_layout();fig.savefig(p/'figure2b.png',dpi=180);fig.savefig(p/'figureS1.png',dpi=180);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(9,4.2),dpi=180,sharey=True)
for ax,mask,title in zip(axes,(cats==0,sachs),[r'Lipschitz in $(G_1,G_M,G_Q)$',r'Lipschitz in $(G_C,G_M,G_Q)$']):
 ax.imshow(mask.astype(int),origin='lower',aspect='auto',cmap=ListedColormap(['#d6d9dc','#4ea87c']),vmin=0,vmax=1)
 ax.set_xticks(range(len(sigmas)),[f'{s:g}' for s in sigmas]);ax.set_yticks(range(len(Ls)),[f'{l:g}' for l in Ls]);ax.set_xlabel(r'$\sigma_k/m$');ax.set_title(title,fontsize=10)
axes[0].set_ylabel(r'Lipschitz bound $L$')
fig.legend(handles=[Patch(facecolor='#4ea87c',label='origin excluded'),Patch(facecolor='#d6d9dc',label='sufficient test inconclusive')],loc='lower center',bbox_to_anchor=(.5,-.04),ncol=2,fontsize=9)
fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(p/'figureS2.png',dpi=180,bbox_inches='tight');plt.close(fig)
print('Images generated; certified',np.sum(cats==0),'Sachs',sachs.sum())
