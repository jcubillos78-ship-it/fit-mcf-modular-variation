"""Figure4.svg: A) triada C58 (render ChimeraX); B) mapa de calor por categorias (Tabla 4)."""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from PIL import Image
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
SOFT=plt.colormaps['YlOrRd']
plt.rcParams['svg.fonttype']='none'; plt.rcParams['font.family']=['Arial','Helvetica','Liberation Sans','DejaVu Sans']
CAT={'req':('Fit requerido','#0072B2'),'con':('Fit contribuye','#56B4E9'),'nor':('Fit no requerido','#E69F00'),
     'alt':('Mortalidad alta (sin mutante)','#009E73'),'sd':('Sin datos','#F0F0F0')}
cols=['G. mellonella\n(inyección)','S. littoralis\n(oral)','D. melanogaster\n(oral)','T. solanivora\n(vía n. i.)']
rows=[('P. protegens/P. chlororaphis\n(17 cepas Fit⁺)',['req','sd','sd','sd']),
      ('P. protegens CHA0,\nP. chlororaphis PCL1391',['sd','con','sd','sd']),
      ('P. protegens Pf-5',['sd','sd','nor','sd']),
      ('P. chlororaphis PCL1606',['con','sd','sd','sd']),
      ('P. protegens CHA0 y 59C',['sd','sd','sd','alt']),
      ('Xenorhabdus spp. (C80)',['sd','sd','sd','sd'])]
W,H=7.6,4.0
fig=plt.figure(figsize=(W,H))
axA=fig.add_axes([0.0,0.10,0.35,0.80]); axA.axis('off')
axA.imshow(np.asarray(Image.open('render/Figure4A_triada_3D.png')),interpolation='none')
fig.text(0.005,0.975,'A',fontsize=12,fontweight='bold',va='top')
fig.text(0.38,0.975,'B',fontsize=12,fontweight='bold',va='top')
fig.text(0.175,0.07,'Tríada catalítica de C58 (distancias en Å)',fontsize=8,ha='center',va='top')
axB=fig.add_axes([0.60,0.30,0.395,0.47]); n=len(rows); m=len(cols)
TXT={'req':'Requerido','con':'Contribuye','nor':'No requerido','alt':'>75 %','sd':''}
# reordenar columnas: Lepidoptera (G.m., S.l., T.s.) | Diptera (D.m.)
order=[0,1,3,2]; cols=[cols[k] for k in order]; rows=[(l,[v[k] for k in order]) for l,v in rows]
cols=['G. mellonella\nInyección','S. littoralis\nOral','T. solanivora\nNo indicada','D. melanogaster\nOral']
# quitar fila sin ningun dato (Xenorhabdus) y declararla en la leyenda de la figura
rows=[r for r in rows if any(v!='sd' for v in r[1])]; n=len(rows)
import matplotlib.cm as cm
from matplotlib.colors import Normalize
cmap=SOFT; nm=Normalize(0,2); SC={'nor':0,'con':1,'req':2}
for i,(lab,vals) in enumerate(rows):
    for j,v in enumerate(vals):
        y=n-1-i
        if v=='sd':
            axB.add_patch(Rectangle((j,y),0.94,0.9,facecolor='#F2F2F2',edgecolor='white',lw=0.8)); continue
        if v=='alt':
            axB.add_patch(Rectangle((j,y),0.94,0.9,facecolor='white',edgecolor='#555555',lw=0.8,hatch='////')); continue
        c=cmap(nm(SC[v])); axB.add_patch(Rectangle((j,y),0.94,0.9,facecolor=c,edgecolor='white',lw=0.8))
    axB.text(-0.08,n-1-i+0.45,lab,fontsize=8,ha='right',va='center')
cax=fig.add_axes([0.60,0.085,0.28,0.02]); import numpy as np
cax.imshow(np.linspace(0,2,256)[None,:],aspect='auto',cmap=cmap,extent=[0,2,0,1]); cax.set_yticks([]); cax.set_xticks([0,1,2]); cax.set_xticklabels(['0 No requerido','1 Contribuye','2 Requerido'],fontsize=7); cax.tick_params(length=2,pad=1)
for s_ in cax.spines.values(): s_.set_visible(False)
fig.text(0.60,0.125,'Contribución de Fit a la virulencia',fontsize=7.5)
axB.set_xlim(0,m); axB.set_ylim(0,n)
axB.set_xticks([j+0.47 for j in range(m)]); axB.set_xticklabels(cols,fontsize=7.5,ha='center')
axB.xaxis.tick_top(); axB.set_yticks([]); axB.tick_params(length=0,pad=3)
for s_ in axB.spines.values(): s_.set_visible(False)
# cabeceras de orden
axB.annotate('',xy=(0,n+0.95),xytext=(2.94,n+0.95),xycoords='data',arrowprops=dict(arrowstyle='-',lw=0.8),annotation_clip=False)
axB.annotate('',xy=(3,n+0.95),xytext=(3.94,n+0.95),xycoords='data',arrowprops=dict(arrowstyle='-',lw=0.8),annotation_clip=False)
axB.text(1.47,n+1.05,'Lepidoptera',ha='center',fontsize=8,fontweight='bold',clip_on=False)
axB.text(3.47,n+1.05,'Diptera',ha='center',fontsize=8,fontweight='bold',clip_on=False)
fig.text(0.60,0.255,'Gris: sin ensayo publicado.\nRayado: mortalidad alta sin mutante (no puntuable).\nXenorhabdus (C80): sin ensayos, no se muestra.',fontsize=6.5,color='#555555',va='top')
fig.savefig('Figures_Vector/Figure4.svg',format='svg'); fig.savefig('/tmp/fig4_prev.png',dpi=130); print('ok')
