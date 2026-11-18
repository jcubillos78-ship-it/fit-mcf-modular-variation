"""Congruencia genoma-toxina: poda a seis cepas, Robinson-Foulds sin raíz, cruces del tanglegrama."""
import itertools, json, re
from io import StringIO
from Bio import Phylo
G='data/genomico_nombres.nwk'; P='data/poro_B10k.contree'
gmap={'Pseudomonas protegens CHA0 anot1':'CHA0','Pseudomonas protegens Pf-5':'Pf-5',
 'Pseudomonas chlororaphis subsp. aureofaciens 30-84':'30-84',
 'Pseudomonas protegens Cab57':None,'Pseudomonas chlororaphis subsp. chlororaphis strain DSM 50083T':None,
 'Pseudomonas fluorescens Pf0-1':None,'Pseudomonas putida KT2440':None,'Pseudomonas aeruginosa PAO1':None,
 'Photorhabdus luminescens subsp. laumondii TTO1':'TT01','Xenorhabdus nematophila ATCC 19061':'ATCC 19061',
 'Xenorhabdus bovienii SS-2004':'SS-2004'}
def prune(t,keep):
    for tip in [x for x in t.get_terminals() if x.name not in keep]: t.prune(tip)
    return t
def splits(t):
    L=frozenset(x.name for x in t.get_terminals()); out=set()
    for c in t.find_clades(terminal=False):
        S=frozenset(x.name for x in c.get_terminals())
        if 2<=len(S)<=len(L)-2:
            if 'CHA0' in S: S=L-S if False else S
            out.add(min(S,L-S,key=lambda s:(len(s),sorted(s))))
    return out
def load_genome():
    txt=open(G).read()
    txt=re.sub(r'"([^"]*)"',lambda m:"'"+m.group(1)+"'",txt)
    t=Phylo.read(StringIO(txt),'newick')
    for x in t.get_terminals(): x.name=gmap[x.name.strip("'")] or x.name
    return prune(t,set(v for v in gmap.values() if v))
# opciones de copia de la toxina en poro: Pf-5 y CHA0 comparten proteína (P. protegens)
P_PROT='Ps_Ps_protegens__N-A-P-T__WP_232642130v1'; P_CHL='Ps_Ps_chlororaphis__N-A-P-T__WP_097304090v1'
P_TT='Ph_Ph_laumondii__N-A-P-T__WP_058589554v1'; P_SB='Xe_Xe_bovienii__N-C-T__WP_038188846v1'
opts={'nematophila N-C-T (WP_493910284)':'Xe_Xe_nematophila__N-C-T__WP_493910284v1',
      'nematophila H-A-P-T (WP_493907850)':'Xe_Xe_nematophila__H-A-P-T__WP_493907850v1'}
res={}
g=load_genome(); gs=splits(g)
print('genomico podado:',[sorted(s) for s in gs])
for k,nem in opts.items():
    t=Phylo.read(P,'newick')
    t=prune(t,{P_PROT,P_CHL,P_TT,P_SB,nem})
    ren={P_PROT:'Pf-5/CHA0',P_CHL:'30-84',P_TT:'TT01',P_SB:'SS-2004',nem:'ATCC 19061'}
    for x in t.get_terminals(): x.name=ren[x.name]
    # duplica Pf-5 y CHA0 (proteina identica) como hermanas con rama 0
    tip=[x for x in t.get_terminals() if x.name=='Pf-5/CHA0'][0]
    tip.name=None
    from Bio.Phylo.BaseTree import Clade
    tip.clades=[Clade(branch_length=0.0,name='CHA0'),Clade(branch_length=0.0,name='Pf-5')]
    ps=splits(t)
    rf=len(gs^ps)
    print(k,'poro podado:',[sorted(s) for s in ps],'RF=',rf)
    res[k]=rf
