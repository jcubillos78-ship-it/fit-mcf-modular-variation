"""Figuras suplementarias S1-S3: arboles de consenso (UFBoot 10 000) por modulo, vectorizados."""
import re, sys
from io import StringIO
from Bio import Phylo
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype']='none'; plt.rcParams['font.family']=['Arial','Helvetica','Liberation Sans','DejaVu Sans']
def cons(path):
    t=open(path).read()
    if path.endswith('.iqtree'):
        i=t.index('Consensus tree in newick format:'); t=t[i:].split('\n',1)[1].strip().split('\n')[0]
    return Phylo.read(StringIO(t.strip()),'newick')
GEN={'Ps':'Ps.','Ph':'Ph.','Xe':'Xe.'}; COL={'Ps':'#009E73','Ph':'#E69F00','Xe':'#CC79A7'}
ARC={'N-A-P-T':'#0072B2','H-A-P-T':'#56B4E9','N-C-T':'#CC79A7','A-P-T':'#999999'}
def disp(n):
    sp,arch,acc=n.split('__'); g=sp.split('_')[0]; rest=sp.split('_',1)[1].replace('_',' ')
    parts=rest.split(' ',1); rest=parts[1] if parts[0] in ('Ps','Ph','Xe') and len(parts)>1 else rest
    return g,f'{GEN[g]} {rest}',arch,acc.replace('v1','').replace('WP_','WP_')
def fig(path,out,cap,title,ht):
    t=cons(path); L=t.get_terminals()
    ps=[x for x in L if x.name.startswith('Ps_')]; t.root_with_outgroup(t.common_ancestor(ps))
    def lay(c,x):
        c._x=x
        for k in c.clades:
            bl=k.branch_length or 0; k._clip=bl>cap; lay(k,x+min(bl,cap))
    lay(t.root,0)
    tips=t.get_terminals()
    for i,tp in enumerate(tips): tp._y=i
    def ys(c):
        if c.is_terminal(): return c._y
        v=[ys(k) for k in c.clades]; c._y=sum(v)/len(v); return c._y
    ys(t.root)
    mx=max(tp._x for tp in tips)
    f,ax=plt.subplots(figsize=(7.2,ht)); ax.axis('off'); ax.set_xlim(-0.02,mx*1.0+0.7*mx+0.3); ax.set_ylim(len(tips)+2,-2)
    def rec(c):
        for k in c.clades:
            ax.plot([c._x,k._x],[k._y,k._y],color='black',lw=0.6,solid_capstyle='butt')
            if k._clip:
                xm=(c._x+k._x)/2
                for d in (-0.01*mx,0.01*mx): ax.plot([xm+d-0.005*mx,xm+d+0.005*mx],[k._y-0.4,k._y+0.4],color='black',lw=0.6)
            if k.clades and k.confidence is not None:
                v=float(k.confidence)
                if v>=70: ax.scatter([k._x],[k._y],s=7 if v>=95 else 5,c='black' if v>=95 else '#aaaaaa',zorder=3,linewidths=0)
            rec(k)
        if c.clades:
            y=[k._y for k in c.clades]; ax.plot([c._x]*2,[min(y),max(y)],color='black',lw=0.6)
    rec(t.root)
    f.canvas.draw(); rend=f.canvas.get_renderer(); inv=ax.transData.inverted()
    fs=4.6 if len(tips)>80 else 5.6
    for tp in tips:
        g,sp,arch,acc=disp(tp.name)
        x0=tp._x+0.012*mx
        ax.scatter([x0],[tp._y],s=5,marker='s',c=ARC[arch],linewidths=0,zorder=3)
        a=ax.text(x0+0.015*mx,tp._y,sp,style='italic',fontsize=fs,va='center',ha='left',color=COL[g])
        bb=a.get_window_extent(rend); w=inv.transform((bb.x1,0))[0]-inv.transform((bb.x0,0))[0]
        ax.text(x0+0.015*mx+w+0.01*mx,tp._y,acc,fontsize=fs,va='center',ha='left',color='#333333')
    L0=0.1; ax.plot([0,L0],[len(tips)+0.8]*2,color='black',lw=0.8); ax.text(L0/2,len(tips)+1.1,'0,1',fontsize=5.5,ha='center',va='top')
    ax.text(0,-1.3,title,fontsize=7.5,fontweight='bold',va='center')
    # leyenda
    for i,(k,c) in enumerate(ARC.items()):
        ax.scatter([mx*1.25+i*0.0],[len(tips)-6+i*1.6],s=8,marker='s',c=c,linewidths=0); ax.text(mx*1.25+0.015*mx,len(tips)-6+i*1.6,k,fontsize=5.5,va='center')
    ax.scatter([mx*1.25],[len(tips)-6+6.4+0],s=7,c='black',linewidths=0); ax.text(mx*1.25+0.015*mx,len(tips)-6+6.4,'UFBoot ≥ 95',fontsize=5.5,va='center')
    ax.scatter([mx*1.25],[len(tips)-6+8.0],s=5,c='#aaaaaa',linewidths=0); ax.text(mx*1.25+0.015*mx,len(tips)-6+8.0,'UFBoot 70-94',fontsize=5.5,va='center')
    f.savefig(out,format='svg',bbox_inches='tight'); plt.close(f); print(out)
fig('data/extremoN_B10k.contree','Figures_Vector/Figura_S1_arbol_extremo_N.svg',0.5,'Extremo N de tipo Mcf1 (65 secuencias)',7.5)
fig('data/poro_B10k.contree','Figures_Vector/Figura_S2_arbol_poro.svg',0.4,'Módulo de poro (111 secuencias)',10.5)
fig('data/proteasaC58_B10k.iqtree','Figures_Vector/Figura_S3_arbol_proteasa_C58.svg',0.4,'Proteasa C58 (100 secuencias)',9.5)
