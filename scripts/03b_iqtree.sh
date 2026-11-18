#!/usr/bin/env bash
# Filogenias por modulo con el modelo fijado por BIC (JTT+I+G4). Seleccion previa: iqtree -s ALN -m MFP -mset LG,WAG,JTT
set -euo pipefail
for m in poro extremoN proteasaC58; do
  iqtree3 -s "data/alineamiento_${m}.fasta" -m JTT+I+G4 -B 10000 -alrt 1000 -bnni --seed 1 -T AUTO --prefix "${m}_B10k"
done
