"""Supplementary_Material.xlsx: S1 accesiones, S2 dominios InterProScan, S3 soporte 1000 vs 10000, S4 arbol genomico."""
import glob, re, collections, statistics
from io import StringIO
from Bio import Phylo
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter
code={'PF30645':'Mcf1-N','PF17914':'HopA1','PF30646':'Activador','PF30660':'C58','PF11713':'C80','PF12920':'Poro'}
gen={}; desc={}
for g in ('Pseudomonas','Photorhabdus','Xenorhabdus'):
    for l in open(f'data/fasta_{g}.fasta'):
        if l.startswith('>'):
            a=l[1:].split()[0]; gen[a]=g; desc[a]=l[1:].strip()
rows=collections.defaultdict(list); length={}; hits=[]
for f in sorted(set(glob.glob('data/*iprscan6*.tsv'))):
    for l in open(f):
        c=l.rstrip('\n').split('\t'); length[c[0]]=int(c[2])
        if c[3] in ('Pfam','CDD','PROSITE profiles') and float(c[8]) if c[8] not in ('-',) else True:
            hits.append((c[0],c[3],c[4],c[5],int(c[6]),int(c[7]),c[8]))
        if c[3]=='Pfam': rows[c[0]].append(c[4])
def cls(S):
    if 'C80' in S and 'Mcf1-N' in S: return 'N-C-T (Mcf3)'
    if 'HopA1' in S and 'Activador' in S and 'C58' in S: return 'H-A-P-T (Mcf2)'
    if 'Mcf1-N' in S and 'Activador' in S and 'C58' in S: return 'N-A-P-T (Mcf1)'
    if 'Activador' in S and 'C58' in S: return 'A-P-T'
    return 'No asignada'
wb=Workbook(); ws=wb.active; ws.title='S1_Accesiones'
ws.append(['Accesión','Género','Descripción (cabecera RefSeq)','Longitud (aa)','Módulos Pfam detectados','Arquitectura'])
for a in sorted(gen,key=lambda x:(gen[x],x)):
    S={code[p] for p in rows[a] if p in code}
    ws.append([a,gen[a],desc[a].split(' ',1)[1] if ' ' in desc[a] else '',length[a],', '.join(sorted(S)),cls(S)])
ws2=wb.create_sheet('S2_Dominios_InterProScan')
ws2.append(['Accesión','Base de datos','Firma','Descripción','Inicio','Fin','Valor E'])
for h in sorted(hits): ws2.append(list(h))
# S3
def splits(path):
    t=open(path).read(); i=t.index('Tree in newick format:'); s=t[i:].split('\n',2)[2].strip().split('\n')[0]
    tr=Phylo.read(StringIO(s),'newick'); allL=frozenset(x.name for x in tr.get_terminals()); d={}
    for c in tr.find_clades(terminal=False):
        if c is tr.root: continue
        S=frozenset(x.name for x in c.get_terminals())
        if len(S)<2 or len(allL-S)<2: continue
        ref=sorted(allL)[0]
        if ref in S: S=allL-S
        v=str(c.name if c.name else c.confidence); d[S]=float(v.split('/')[-1])
    return d
ws3=wb.create_sheet('S3_Soporte_1000_vs_10000')
ws3.append(['Módulo','Taxones en la rama','Taxón de referencia del clado','UFBoot 1 000','UFBoot 10 000','Diferencia (10 000 − 1 000)','Cambia de categoría (≥95 / 70-94 / <70)'])
cat=lambda v:'≥95' if v>=95 else ('70-94' if v>=70 else '<70')
for mod,a,b in (('Extremo N','data/extremoN_1000.iqtree','data/extremoN_B10k.iqtree'),('Proteasa C58','data/proteasaC58_1000.iqtree','data/proteasaC58_B10k.iqtree')):
    A=splits(a); B=splits(b)
    for S in sorted(set(A)&set(B),key=len):
        ws3.append([mod,len(S),sorted(S)[0],A[S],B[S],B[S]-A[S],'Sí' if cat(A[S])!=cat(B[S]) else 'No'])
    print(mod,len(A),len(B),'comunes',len(set(A)&set(B)))
ws4=wb.create_sheet('S4_Arbol_genomico')
ws4.append(['Parámetro','Valor'])
st=open('data/ARBOL_GENOMICO.analysisStats').read().split('\n'); pg=[l for l in st if l.startswith('PGFams')][0].split('\t')[1].split(',')
for k,v in (('Genomas',11),('Familias de copia única recuperadas',785),('Familias retenidas',100),('Aminoácidos alineados',48086),('Nucleótidos alineados',144258),
            ('Programa','RAxML 8.2.12 (GTRCAT particionado; -f a -x 12345 -N 100)'),('Duración del trabajo (s)',12821)): ws4.append([k,v])
ws4.append([]); ws4.append(['PGFam retenida (n.º)','Identificador'])
for i,p in enumerate(pg,1): ws4.append([i,p])
for w in wb.worksheets:
    for c in w[1]: c.font=Font(bold=True,name='Arial',size=10); c.fill=PatternFill('solid',fgColor='D9D9D9')
    w.freeze_panes='A2'
    for col in range(1,w.max_column+1): w.column_dimensions[get_column_letter(col)].width=min(60,max(12,max(len(str(c.value or '')) for c in w[get_column_letter(col)][:60])+2))
    for row in w.iter_rows(min_row=2):
        for c in row: c.font=Font(name='Arial',size=10)
wb.save('Supplementary_Material.xlsx'); print([ (w.title,w.max_row) for w in wb.worksheets],len(pg))
