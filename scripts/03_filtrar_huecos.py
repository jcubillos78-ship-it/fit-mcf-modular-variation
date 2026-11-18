"""Elimina columnas con mas del 50 % de huecos de un alineamiento FASTA. Uso: python 03_filtrar_huecos.py in.fasta out.fasta"""
import sys
from Bio import AlignIO
aln=AlignIO.read(sys.argv[1],'fasta'); n=len(aln); keep=[i for i in range(aln.get_alignment_length()) if sum(r.seq[i]=='-' for r in aln)/n<=0.5]
out=aln[:,0:0]
for i in keep: out=out+aln[:,i:i+1]
AlignIO.write(out,sys.argv[2],'fasta'); print(len(keep),'columnas retenidas')
