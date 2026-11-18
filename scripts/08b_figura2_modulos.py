"""Figure2.svg: tres arboles ML por modulo (A extremo N, B poro, C proteasa C58) en estilo iTOL: sin texto en las puntas,
franjas de genero y arquitectura, soporte UFBoot 10 000."""
import re
from io import StringIO
from Bio import Phylo
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
plt.rcParams['svg.fonttype']='none'; plt.rcParams['font.family']=['Arial','Helvetica','Liberation Sans','DejaVu Sans']
def cons(path):
    t=open(path).read()
    if path.endswith('.iqtree'):
        i=t.index('Consensus tree in newick format:'); t=t[i:].split('\n',1)[1].strip().split('\n')[0]
    return Phylo.read(StringIO(t.strip()),'newick')
GC={'Ps':'#009E73','Ph':'#E69F00','Xe':'#CC79A7'}; GN={'Ps':'Pseudomonas','Ph':'Photorhabdus','Xe':'Xenorhabdus'}
AC={'N-A-P-T':'#0072B2','H-A-P-T':'#56B4E9','N-C-T':'#D55E00','A-P-T':'#999999'}
AL={'N-A-P-T':'Mcf1 (N-A-P-T)','H-A-P-T':'Mcf2 (H-A-P-T)','N-C-T':'Mcf3 (N-C-T)','A-P-T':'A-P-T'}
def draw(ax,path,cap,title,letter):
    t=cons(path); ps=[x for x in t.get_terminals() if x.name.startswith('Ps_')]; t.root_with_outgroup(t.common_ancestor(ps))
    def lay(c,x):
        c._x=x
        for k in c.clades:
            bl=k.branch_length or 0; k._clip=bl>cap; lay(k,x+min(bl,cap))
    lay(t.root,0); tips=t.get_terminals()
    for i,tp in enumerate(tips): tp._y=i
    def ys(c):
        if c.is_terminal(): return c._y
        v=[ys(k) for k in c.clades]; c._y=sum(v)/len(v); return c._y
    ys(t.root); mx=max(tp._x for tp in tips)
    def rec(c):
        for k in c.clades:
            ax.plot([c._x,k._x],[k._y,k._y],color='black',lw=0.5,solid_capstyle='butt')
            if k._clip:
                xm=(c._x+k._x)/2
                for d in (-0.012*mx,0.012*mx): ax.plot([xm+d-0.006*mx,xm+d+0.006*mx],[k._y-0.9,k._y+0.9],color='black',lw=0.5)
            if k.clades and k.confidence is not None:
                v=float(k.confidence)
                if v>=70: ax.scatter([k._x],[k._y],s=5 if v>=95 else 3.5,c='black' if v>=95 else '#aaaaaa',zorder=3,linewidths=0)
            rec(k)
        if c.clades:
            y=[k._y for k in c.clades]; ax.plot([c._x]*2,[min(y),max(y)],color='black',lw=0.5)
    rec(t.root)
    w=mx*0.07; x1=mx*1.03; x2=x1+w*1.15
    for tp in tips:
        sp,arch,acc=tp.name.split('__'); g=sp.split('_')[0]
        ax.add_patch(plt.Rectangle((x1,tp._y-0.5),w,1.0,fc=GC[g],ec='none',lw=0))
        ax.add_patch(plt.Rectangle((x2,tp._y-0.5),w,1.0,fc=AC[arch],ec='none',lw=0))
    ax.set_xlim(-0.02*mx,x2+w*1.1); ax.set_ylim(len(tips)+6,-1); ax.axis('off')
    ax.plot([0,0.1],[len(tips)+2.5]*2,color='black',lw=0.8); ax.text(0.05,len(tips)+3.5,'0,1',fontsize=8,ha='center',va='top')
    ax.text(-0.02*mx,-4,title,fontsize=8.5,fontweight='bold',va='bottom',ha='left',clip_on=False)
fig=plt.figure(figsize=(7.2,6.4))
specs=[('data/extremoN_B10k.contree',0.5,'Extremo N (n = 65)','A'),('data/poro_B10k.contree',0.4,'Poro (n = 111)','B'),('data/proteasaC58_B10k.iqtree',0.4,'Proteasa C58 (n = 100)','C')]
for i,(p,cap,ti,L) in enumerate(specs):
    ax=fig.add_axes([0.01+i*0.33,0.20,0.31,0.72]); draw(ax,p,cap,ti,L); fig.text(0.01+i*0.33,0.995,L,fontsize=12,fontweight='bold',va='top')
h=[Patch(fc=GC[k],label=GN[k]) for k in GC]+[Patch(fc=AC[k],label=AL[k]) for k in AC]
fig.legend(handles=h[:3],loc='lower left',bbox_to_anchor=(0.01,0.0),ncol=1,fontsize=8,frameon=False,title='Franja interna: género',title_fontsize=8,alignment='left',handlelength=1.1,columnspacing=1.0)
fig.legend(handles=h[3:],loc='lower left',bbox_to_anchor=(0.30,0.0),ncol=1,fontsize=8,frameon=False,title='Franja externa: arquitectura',title_fontsize=8,alignment='left',handlelength=1.1,columnspacing=1.0)
fig.legend(handles=[Line2D([],[],marker='o',ls='',color='black',ms=3,label='UFBoot ≥ 95'),Line2D([],[],marker='o',ls='',color='#aaaaaa',ms=3,label='UFBoot 70-94')],loc='lower right',bbox_to_anchor=(0.995,0.0),fontsize=8,frameon=False,handletextpad=0.2)
fig.savefig('Figures_Vector/Figure2.svg',format='svg'); fig.savefig('/tmp/f2.png',dpi=130); print('ok')
