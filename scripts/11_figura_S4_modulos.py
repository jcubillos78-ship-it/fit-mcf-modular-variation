"""Figura S4: % de proteinas con cada modulo Pfam por genero (S1, N=300)."""
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
SOFT=plt.colormaps['YlOrRd']
plt.rcParams['svg.fonttype']='none'; plt.rcParams['font.family']=['Arial','Helvetica','Liberation Sans','DejaVu Sans']
x=pd.read_excel('Supplementary_Material.xlsx',sheet_name='S1_Accesiones')
mods=['Mcf1-N','Activador','HopA1','C58','C80','Poro']
M=pd.DataFrame({m:x['Módulos Pfam detectados'].str.contains(m).groupby(x['Género']).mean()*100 for m in mods}).loc[['Pseudomonas','Photorhabdus','Xenorhabdus']]
fig,ax=plt.subplots(figsize=(5.2,2.6)); im=ax.imshow(M.values,cmap=SOFT,vmin=0,vmax=100,aspect='auto')
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        pass
ax.set_xticks(range(6)); ax.set_xticklabels(mods,fontsize=8); ax.set_yticks(range(3)); ax.set_yticklabels(M.index,fontsize=8,style='italic')
ax.set_xticks(np.arange(-.5,6),minor=True); ax.set_yticks(np.arange(-.5,3),minor=True); ax.grid(which='minor',color='white',lw=1.5); ax.tick_params(which='both',length=0)
for s in ax.spines.values(): s.set_visible(False)
cb=fig.colorbar(im,ax=ax,fraction=0.04,pad=0.03); cb.set_label('% de proteínas (n = 100 por género)',fontsize=8); cb.ax.tick_params(labelsize=7)
fig.savefig('Figures_Vector/Figure_S4.svg',format='svg',bbox_inches='tight'); fig.savefig('/tmp/s4.png',dpi=130,bbox_inches='tight'); print(M.round(0))
