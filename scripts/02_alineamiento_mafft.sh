#!/usr/bin/env bash
# Alineamiento de cada conjunto con MAFFT 7.505, modo L-INS-i. Entrada: FASTA por modulo (poro, extremoN, proteasaC58).
set -euo pipefail
for m in poro extremoN proteasaC58; do
  mafft --localpair --maxiterate 1000 --thread -1 "data/${m}.fasta" > "data/alineamiento_${m}.fasta"
done
