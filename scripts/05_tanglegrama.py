"""Figura 3: tanglegrama entre el arbol genomico podado y el arbol del modulo de poro (seis cepas)."""
import re, sys
from io import StringIO
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype']='none'
plt.rcParams['font.family']=['Arial','Helvetica','Liberation Sans','DejaVu Sans']
exec(open('scripts/04_congruencia.py').read().split('# opciones')[0])   # reutiliza gmap, prune, splits, load_genome
NEM=sys.argv[1] if len(sys.argv)>1 else 'N-C-T'
nem={'N-C-T':'Xe_Xe_nematophila__N-C-T__WP_493910284v1','H-A-P-T':'Xe_Xe_nematophila__H-A-P-T__WP_493907850v1'}[NEM]
def load_poro():
    t=Phylo.read(P,'newick')
    keep={'Ps_Ps_protegens__N-A-P-T__WP_232642130v1':'Pf-5/CHA0','Ps_Ps_chlororaphis__N-A-P-T__WP_097304090v1':'30-84',
          'Ph_Ph_laumondii__N-A-P-T__WP_058589554v1':'TT01','Xe_Xe_bovienii__N-C-T__WP_038188846v1':'SS-2004',nem:'ATCC 19061'}
    t=prune(t,set(keep))
    for x in t.get_terminals(): x.name=keep[x.name]
    tip=[x for x in t.get_terminals() if x.name=='Pf-5/CHA0'][0]
    tip.name=None; tip.clades=[Clade(branch_length=0.0,name='CHA0'),Clade(branch_length=0.0,name='Pf-5')]
    return t
def reroot(t):
    ps=[x for x in t.get_terminals() if x.name in ('CHA0','Pf-5','30-84')]
    t.root_with_outgroup(t.common_ancestor(ps)); 
    # raiz binaria: si la raiz tiene 3 hijos, se resuelve con el grupo externo
    return t
ORDER=['CHA0','Pf-5','30-84','TT01','ATCC 19061','SS-2004']
def rotate(c):
    if c.is_terminal(): return ORDER.index(c.name)
    for k in c.clades: rotate(k)
    c.clades.sort(key=lambda k:min(ORDER.index(x.name) for x in k.get_terminals()))
    return 0
def layout(t,cap):
    pos={}; 
    tips=t.get_terminals()
    def walk(c,x):
        c._x=x
        for k in c.clades:
            bl=k.branch_length or 0.0
            k._clip=bl>cap
            walk(k,x+min(bl,cap))
    walk(t.root,0.0)
    for tp in tips: tp._y=ORDER.index(tp.name)
    def ys(c):
        if c.is_terminal(): return c._y
        v=[ys(k) for k in c.clades]; c._y=sum(v)/len(v); return c._y
    ys(t.root); return t
def draw(ax,t,xs,scale,flip=False):
    global t_cur; t_cur=t
    sgn=-1 if flip else 1
    def X(c): return xs+sgn*c._x*scale
    def rec(c):
        for k in c.clades:
            ax.plot([X(c),X(k)],[k._y,k._y],color='black',lw=1.0,solid_capstyle='butt')
            if k.clades:
                S=frozenset(x.name for x in k.get_terminals()); key=min(S,frozenset(ORDER)-S,key=lambda q:(len(q),sorted(q)))
                v=SUP[id(t_cur)].get(key)
                if v is not None:
                    xm2=X(k)-sgn*0.004 if k._clip else (X(c)+X(k))/2
                    ax.text(xm2,k._y-0.14,v,fontsize=8,ha='center',va='bottom')
            if k._clip:
                xm=(X(c)+X(k))/2
                for d in (-0.012,0.012): ax.plot([xm+d-0.006*sgn*-1,xm+d+0.006*sgn*-1],[k._y-0.16,k._y+0.16],color='black',lw=1.0)
            rec(k)
        if c.clades:
            ys=[k._y for k in c.clades]
            ax.plot([X(c),X(c)],[min(ys),max(ys)],color='black',lw=1.0)
    rec(t.root)
def support_labels(ax,t,xs,scale,flip,nodes):
    pass
# soporte por biparticion (de los arboles completos; ver 04_congruencia.py)
K=lambda *a:frozenset(a)
SUP_G={K('CHA0','Pf-5'):'93',K('CHA0','Pf-5','30-84'):'100',K('ATCC 19061','SS-2004'):'100'}
SUP_P={K('CHA0','Pf-5','30-84'):'100',K('ATCC 19061','SS-2004'):'100'}
g=reroot(load_genome()); p=reroot(load_poro())
rotate(g.root); rotate(p.root)
SUP={id(g):SUP_G,id(p):SUP_P}
CAP=0.4
layout(g,CAP); layout(p,CAP)
# cruces
order_l=[x.name for x in g.get_terminals()]; order_r=[x.name for x in p.get_terminals()]
cross=sum(1 for i in range(6) for j in range(i+1,6) if (g.find_any(name=ORDER[i])._y-g.find_any(name=ORDER[j])._y)*(p.find_any(name=ORDER[i])._y-p.find_any(name=ORDER[j])._y)<0)
gs=splits(load_genome()); ps=splits(p)
rf=len(gs^ps)
print('RF',rf,'cruces',cross)
lab={'CHA0':('P. protegens',' CHA0'),'Pf-5':('P. protegens',' Pf-5'),'30-84':('P. chlororaphis',' 30-84'),
     'TT01':('P. laumondii',' TT01'),'ATCC 19061':('X. nematophila',' ATCC 19061'),'SS-2004':('X. bovienii',' SS-2004')}
col={'CHA0':'#009E73','Pf-5':'#009E73','30-84':'#009E73','TT01':'#E69F00','ATCC 19061':'#CC79A7','SS-2004':'#CC79A7'}
fig,ax=plt.subplots(figsize=(7.2,3.5)); ax.set_xlim(0,1); ax.set_ylim(5.7,-0.9); ax.axis('off')
xl,xr=0.02,0.98
sgl=0.27/max(tp._x for tp in g.get_terminals()); sgr=0.27/max(tp._x for tp in p.get_terminals())
draw(ax,g,xl,sgl); draw(ax,p,xr,sgr,flip=True)
fig.canvas.draw(); rend=fig.canvas.get_renderer(); inv=ax.transData.inverted()
def wid(t):
    bb=t.get_window_extent(rend); return inv.transform((bb.x1,0))[0]-inv.transform((bb.x0,0))[0]
for name in ORDER:
    yl=g.find_any(name=name)._y; yr=p.find_any(name=name)._y
    xtl=xl+g.find_any(name=name)._x*sgl; xtr=xr-p.find_any(name=name)._x*sgr
    sp,st=lab[name]
    t1=ax.text(0,yl,sp,style='italic',fontsize=8,va='center',ha='left'); w1=wid(t1)
    t2=ax.text(0,yl,st.strip(),fontsize=8,va='center',ha='left'); w2=wid(t2)
    tot=w1+0.008+w2; x0=0.5-tot/2
    t1.set_x(x0); t2.set_x(x0+w1+0.008)
    ax.plot([xtl+0.006,x0-0.008],[yl,yl],color=col[name],lw=1.4)
    ax.plot([x0+tot+0.008,xtr-0.006],[yr,yr],color=col[name],lw=1.4)
ax.text(xl,-0.75,'Filogenia genómica (100 genes de copia única)',fontsize=8.5,ha='left',va='center',fontweight='bold')
ax.text(xr,-0.75,'Módulo de poro (JTT+I+G4)',fontsize=8.5,ha='right',va='center',fontweight='bold')
for xs,sc,fl in ((xl,sgl,False),(xr,sgr,True)):
    L=0.1*sc; x0=xs if not fl else xs-L
    ax.plot([x0,x0+L],[5.3,5.3],color='black',lw=1.0); ax.text(x0+L/2,5.45,'0,1',fontsize=8,ha='center',va='top')
ax.text(0.5,5.45,f'Robinson-Foulds = {rf}; cruces = {cross}',fontsize=8,ha='center',va='top')
out=f'Figures_Vector/Figura_3_tanglegrama_{NEM}.svg'
fig.savefig(out,format='svg',bbox_inches='tight'); print(out)
