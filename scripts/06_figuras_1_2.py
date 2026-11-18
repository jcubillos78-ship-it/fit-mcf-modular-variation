"""Figura 1 (arbol genomico) y Figura 2 (arquitecturas modulares). SVG con texto editable, Arial/Helvetica."""
import re, glob, collections, statistics
from io import StringIO
from Bio import Phylo
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams['svg.fonttype']='none'
plt.rcParams['font.family']=['Arial','Helvetica','Liberation Sans','DejaVu Sans']
# ---------- Figura 1 ----------
txt=open('data/genomico_nombres.nwk').read()
txt=re.sub(r'"([^"]*)"',lambda m:"'"+m.group(1)+"'",txt)
t=Phylo.read(StringIO(txt),'newick')
for x in t.get_terminals(): x.name=x.name.strip("'")
ps=[x for x in t.get_terminals() if x.name.startswith('Pseudomonas') and 'aeruginosa' in x.name]
out=t.common_ancestor([x for x in t.get_terminals() if x.name.startswith('Pseudomonas')])
t.root_with_outgroup(out)
CAP=0.3
def lay(c,x):
    c._x=x
    for k in c.clades:
        bl=k.branch_length or 0; k._clip=bl>CAP; lay(k,x+min(bl,CAP))
lay(t.root,0)
tips=t.get_terminals()
for i,tp in enumerate(tips): tp._y=i
def ys(c):
    if c.is_terminal(): return c._y
    v=[ys(k) for k in c.clades]; c._y=sum(v)/len(v); return c._y
ys(t.root)
fig,ax=plt.subplots(figsize=(7.2,3.6)); ax.axis('off'); ax.set_ylim(len(tips)+0.4,-1); ax.set_xlim(-0.01,1.15)
sc=1.0
def rec(c):
    for k in c.clades:
        ax.plot([c._x*sc,k._x*sc],[k._y,k._y],color='black',lw=1.0,solid_capstyle='butt')
        if k.clades and k.confidence is not None:
            ax.text((c._x+k._x)/2*sc if (k._x-c._x)>0.03 else k._x*sc-0.008,k._y-0.15,f'{k.confidence:.0f}',fontsize=8,ha='center' if (k._x-c._x)>0.03 else 'right',va='bottom')
        if k._clip:
            xm=(c._x+k._x)/2*sc
            for d in (-0.012,0.012): ax.plot([xm+d-0.006,xm+d+0.006],[k._y-0.2,k._y+0.2],color='black',lw=1.0)
        rec(k)
    if c.clades:
        y=[k._y for k in c.clades]; ax.plot([c._x*sc]*2,[min(y),max(y)],color='black',lw=1.0)
rec(t.root)
fig.canvas.draw(); rend=fig.canvas.get_renderer()
colg={'Pseudomonas':'#009E73','Photorhabdus':'#E69F00','Xenorhabdus':'#CC79A7'}
sp=re.compile(r'^(\w+ [\w\.]+)(.*)$')
DISP={'Xenorhabdus nematophila ATCC 19061':('X. nematophila','ATCC 19061'),'Xenorhabdus bovienii SS-2004':('X. bovienii','SS-2004'),
'Photorhabdus luminescens subsp. laumondii TTO1':('P. laumondii','TT01'),'Pseudomonas aeruginosa PAO1':('P. aeruginosa','PAO1'),
'Pseudomonas putida KT2440':('P. putida','KT2440'),'Pseudomonas fluorescens Pf0-1':('P. fluorescens','Pf0-1'),
'Pseudomonas protegens Cab57':('P. protegens','Cab57'),'Pseudomonas protegens Pf-5':('P. protegens','Pf-5'),
'Pseudomonas protegens CHA0 anot1':('P. protegens','CHA0'),
'Pseudomonas chlororaphis subsp. chlororaphis strain DSM 50083T':('P. chlororaphis subsp. chlororaphis','DSM 50083T'),
'Pseudomonas chlororaphis subsp. aureofaciens 30-84':('P. chlororaphis subsp. aureofaciens','30-84')}
inv=ax.transData.inverted()
for tp in tips:
    n=tp.name; gen=n.split()[0]; spn,strn=DISP[n]; focal=(strn=='CHA0')
    tx=tp._x*sc+0.012
    a=ax.text(tx,tp._y,spn,style='italic',fontsize=8,va='center',ha='left',color=colg[gen],fontweight='bold' if focal else 'normal')
    bb=a.get_window_extent(rend)
    w=inv.transform((bb.x1,0))[0]-inv.transform((bb.x0,0))[0]
    ax.text(tx+w+0.012,tp._y,strn,fontsize=8,va='center',ha='left',fontweight='bold' if focal else 'normal')
L=0.1; ax.plot([0,L],[len(tips)-0.2]*2,color='black',lw=1); ax.text(L/2,len(tips)-0.0,'0,1',fontsize=8,ha='center',va='top')
fig.savefig('Figures_Vector/Figura_1_arbol_genomico.svg',format='svg',bbox_inches='tight')
plt.close(fig)
# ---------- Figura 2 ----------
rows=collections.defaultdict(list); length={}
for f in glob.glob('data/*iprscan6*.tsv'):
    for l in open(f):
        c=l.rstrip('\n').split('\t'); length[c[0]]=int(c[2])
        if c[3]=='Pfam': rows[c[0]].append((c[4],int(c[6]),int(c[7])))
genus={}
for g in ('Pseudomonas','Photorhabdus','Xenorhabdus'):
    for l in open(f'data/fasta_{g}.fasta'):
        if l.startswith('>'): genus[l[1:].split()[0]]=g
code={'PF30645':'N','PF17914':'H','PF30646':'A','PF30660':'P','PF11713':'C','PF12920':'T'}
def cls(S):
    if 'C' in S and 'N' in S: return 'N-C-T'
    if 'H' in S and 'A' in S and 'P' in S: return 'H-A-P-T'
    if 'N' in S and 'A' in S and 'P' in S: return 'N-A-P-T'
    if 'A' in S and 'P' in S and 'N' not in S and 'H' not in S: return 'A-P-T'
arch={a:cls({code[h[0]] for h in v if h[0] in code}) for a,v in rows.items()}
spec=[('N-A-P-T','Mcf1 / Fit'),('H-A-P-T','Mcf2'),('N-C-T','Mcf3'),('A-P-T','Sin módulo N-terminal')]
cols={'N':'#0072B2','H':'#56B4E9','A':'#E69F00','P':'#D55E00','C':'#CC79A7','T':'#009E73'}
names={'N':'Mcf1-N (PF30645)','H':'HopA1 (PF17914)','A':'Unión al activador (PF30646)','P':'Proteasa C58 (PF30660)','C':'Peptidasa C80 (PF11713)','T':'Poro TcdA/TcdB (PF12920)'}
fig,ax=plt.subplots(figsize=(7.2,3.8)); ax.axis('off')
ymax=3500; xs=0.30; W=0.45
for i,(a,lab) in enumerate(spec):
    prots=[p for p,v in arch.items() if v==a]
    med=statistics.median(length[p] for p in prots)
    pure=[p for p in prots if 'C' not in {code[h[0]] for h in rows[p] if h[0] in code}] or prots
    rep=min(pure,key=lambda p:abs(length[p]-med)) if a!='N-C-T' else min(prots,key=lambda p:abs(length[p]-med))
    y=i*1.0
    ax.text(0.0,y,a,fontsize=8.5,fontweight='bold',va='center',ha='left'); ax.text(0.0,y+0.27,lab,fontsize=8,va='center',ha='left')
    L=length[rep]; ax.add_patch(Rectangle((xs,y-0.03),L/ymax*W,0.06,color='#b0b0b0',lw=0))
    dom=collections.defaultdict(list)
    for (acc,s,e) in rows[rep]:
        if acc in code: dom[code[acc]].append((s,e))
    for k,v in dom.items():
        s=min(x[0] for x in v); e=max(x[1] for x in v)
        ax.add_patch(Rectangle((xs+s/ymax*W,y-0.17),(e-s)/ymax*W,0.34,color=cols[k],lw=0))
    n=collections.Counter(genus[p] for p in prots)
    ax.text(xs+W+0.03,y,f"Ps {n.get('Pseudomonas',0)} · Ph {n.get('Photorhabdus',0)} · Xe {n.get('Xenorhabdus',0)}",fontsize=8,va='center')
    ax.text(xs+L/ymax*W+0.008,y-0.28,f'{int(med):,}'.replace(',',' ')+' aa',fontsize=8,va='center',color='#444444')
# eje
ax.plot([xs,xs+W],[3.75,3.75],color='black',lw=0.8)
for v in range(0,3501,500):
    ax.plot([xs+v/ymax*W]*2,[3.75,3.82],color='black',lw=0.8); ax.text(xs+v/ymax*W,3.9,f'{v}',fontsize=8,ha='center',va='top')
ax.text(xs+W/2,4.3,'Posición (aminoácidos)',fontsize=8,ha='center',va='top')
for i,(k,v) in enumerate(names.items()):
    ax.add_patch(Rectangle((0.0+ (i%3)*0.34,4.75+(i//3)*0.28),0.02,0.16,color=cols[k],lw=0)); ax.text(0.03+(i%3)*0.34,4.83+(i//3)*0.28,v,fontsize=8,va='center')
ax.set_xlim(-0.02,1.0); ax.set_ylim(5.5,-0.6)
fig.savefig('Figures_Vector/Figura_2_arquitecturas.svg',format='svg',bbox_inches='tight')
print('ok')
