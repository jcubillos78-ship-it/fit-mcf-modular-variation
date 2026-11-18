"""Panel de sintenia del locus fit/mcf (datos exportados de clinker: sintenia_3cepas.html) redibujado en vectorial."""
import json, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Polygon
plt.rcParams['svg.fonttype']='none'; plt.rcParams['font.family']=['Arial','Helvetica','Liberation Sans','DejaVu Sans']
d=json.load(open('data/clinker_sintenia_3cepas.json'))
name={'P_protegens_Pf5_cut':'P. protegens Pf-5','P_protegens_CHA0_cut':'P. protegens CHA0','P_chlororaphis_3084_cut':'P. chlororaphis 30-84'}
g2grp={u:gr['label'] for gr in d['groups'] for u in gr['genes']}
# color: toxina (Group 8), resto de grupos conservados en gris azulado, sin homologo en blanco
COL0={'Group 8':'#D55E00','Group 6':'#0072B2','Group 7':'#56B4E9','Group 9':'#009E73','Group 10':'#CC79A7','Group 11':'#E69F00','Group 12':'#F0E442'}
GC=['#7B3FB5','#A63DB5','#E8438A','#F0505A','#F08A3A','#E8A92E','#B5D334','#86E84A','#4BE87A','#2FD7A0','#1CB8C8','#3A8EE8','#5A5CD0']
COL={f'Group {i}':c for i,c in enumerate(GC)}
DEF='#BBBBBB'
LAB={'Group 0':'transportador','Group 1':'PelD','Group 2':'HEAT','Group 3':'PelF','Group 4':'PelG','Group 5':'proteína periplásmica','Group 6':'T1SS permeasa/ATPasa','Group 7':'T1SS HlyD','Group 8':'mcf/fit (TcdA/TcdB)','Group 9':'TolC','Group 10':'ATPasa / PAS','Group 11':'LysR','Group 12':'regulador de respuesta'}
fig,ax=plt.subplots(figsize=(7.2,3.3))
Y=[2.0,1.0,0.0]; H=0.22; pos={}
for y,c in zip(Y,d['clusters']):
    L=c['loci'][0]['genes']; flip=False
    ax.plot([0,40],[y,y],color='#888888',lw=0.6,zorder=1)
    ax.text(-0.8,y,name[c['name']],ha='right',va='center',fontsize=8,style='italic')
    for g in L:
        s,e=g['start']/1000,g['end']/1000
        st=g['strand']
        if flip: s,e=40-e,40-s; st=-st
        col=COL.get(g2grp.get(g['uid']),DEF if g['uid'] in g2grp else 'white')
        hl=min(0.9,(e-s)*0.5)
        if st==1: pts=[(s,y-H),(e-hl,y-H),(e,y),(e-hl,y+H),(s,y+H)]
        else: pts=[(e,y-H),(s+hl,y-H),(s,y),(s+hl,y+H),(e,y+H)]
        ax.add_patch(Polygon(pts,closed=True,fc=col,ec='#333333',lw=0.5,zorder=3)); pos[g['uid']]=(s,e,y,st)
for l in d['links']:
    a,b=pos.get(l['query']['uid']),pos.get(l['target']['uid'])
    if not a or not b or abs(a[2]-b[2])<0.1: continue
    if a[2]<b[2]: a,b=b,a  # a arriba
    ya,yb=a[2]-H,b[2]+H
    gv=str(round(1-l['identity'],3)); ax.add_patch(Polygon([(a[0],ya),(a[1],ya),(b[1],yb),(b[0],yb)],closed=True,fc=plt.colormaps['Greys'](l['identity']),ec='none',alpha=0.9,zorder=2))
ax.set_xlim(-14,41); ax.set_ylim(-0.5,2.5); ax.axis('off')
ax.plot([30,40],[-0.45,-0.45],color='black',lw=1); ax.text(35,-0.52,'10 kb',ha='center',va='top',fontsize=8)
import numpy as np
cb=ax.inset_axes([0.42,-0.005,0.25,0.03]); cb.imshow(np.linspace(0,1,256)[None,:],aspect='auto',cmap='Greys',extent=[0,100,0,1]); cb.set_yticks([]); cb.set_xticks([0,50,100]); cb.tick_params(labelsize=8,length=2); [x.set_visible(False) for x in cb.spines.values()]; cb.set_title('Identidad de aminoácidos (%)',fontsize=8,pad=2,loc='left')
import matplotlib.patches as mp
h=[mp.Patch(fc=COL[k],ec='#333',lw=0.5,label=f'{k.split()[1]}: {LAB[k]}') for k in LAB]+[mp.Patch(fc=DEF,ec='#333',lw=0.5,label='gen sin grupo asignado'),mp.Patch(fc='white',ec='#333',lw=0.5,label='sin homólogo')]
ax.legend(handles=h,loc='upper center',bbox_to_anchor=(0.42,-0.11),ncol=3,fontsize=8,frameon=False,handlelength=1.2)
fig.savefig('render/sintenia_B.svg',format='svg',bbox_inches='tight'); fig.savefig('/tmp/sint.png',dpi=130,bbox_inches='tight'); print('ok')
