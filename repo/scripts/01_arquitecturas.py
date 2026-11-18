"""Clasifica las 300 proteinas por arquitectura modular a partir de las salidas de InterProScan 6 (TSV)."""
import glob, collections
CODE={'PF30645':'N','PF17914':'H','PF30646':'A','PF30660':'P','PF11713':'C','PF12920':'T'}
def cls(S):
    if 'C' in S and 'N' in S: return 'N-C-T'
    if 'H' in S and 'A' in S and 'P' in S: return 'H-A-P-T'
    if 'N' in S and 'A' in S and 'P' in S: return 'N-A-P-T'
    if 'A' in S and 'P' in S: return 'A-P-T'
    return 'No asignada'
mods=collections.defaultdict(set)
for f in sorted(set(glob.glob('data/*iprscan6*.tsv'))):
    for l in open(f):
        c=l.rstrip('\n').split('\t')
        if c[3]=='Pfam' and c[4] in CODE: mods[c[0]].add(CODE[c[4]])
genus={}
for g in ('Pseudomonas','Photorhabdus','Xenorhabdus'):
    for l in open(f'data/fasta_{g}.fasta'):
        if l.startswith('>'): genus[l[1:].split()[0]]=g
for (a,g),n in sorted(collections.Counter((cls(mods[p]),genus[p]) for p in genus).items()): print(a,g,n,sep='\t')
