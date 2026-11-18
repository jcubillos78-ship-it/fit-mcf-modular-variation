# Variación modular de la toxina Fit/Mcf y especificidad de hospedero en *Pseudomonas*, *Photorhabdus* y *Xenorhabdus*

Código, datos de entrada y figuras vectoriales del manuscrito enviado a *Current Microbiology*.
Autor: Juan David Cubillos Martínez (Universidad Antonio Nariño, Bogotá). Docente: Javier Vanegas Guerrero.

## Estructura
- `scripts/` scripts numerados (ejecutar desde la raíz: `python -I scripts/NN_nombre.py`)
- `data/` alineamientos y árboles IQ-TREE, tablas InterProScan, árbol genómico, datos de clinker, metadatos de los 11 genomas (BV-BRC)
- `Figures_Vector/` Figure1.svg a Figure4.svg y Figure_S1 a S4
- `render/` render 3D de la tríada (ChimeraX) y paneles intermedios
- `Supplementary_Material.xlsx` tablas S1 a S4
- `environment.yml` entorno conda

## Scripts
| Script | Qué hace |
|---|---|
| 01_arquitecturas.py | Clasifica 300 proteínas por arquitectura |
| 02_alineamiento_mafft.sh, 03_filtrar_huecos.py, 03b_iqtree.sh | Alineamiento MAFFT L-INS-i, filtro de huecos > 50 % e IQ-TREE 3 (10 000 réplicas) |
| 04_congruencia.py | Poda a seis cepas y distancia de Robinson-Foulds |
| 05_tanglegrama.py | Panel A de la Figura 3 |
| 06_figuras_1_2.py | Paneles de la Figura 1 |
| 07c_triada_chimerax.cxc | Render 3D de la tríada catalítica (ChimeraX) |
| 07d_figura4_compuesta.py | Figura 4 (panel A: render; panel B: mapa de calor) |
| 08b_figura2_modulos.py | Figura 2 (tres árboles por módulo) |
| 08_arboles_modulos.py | Figuras S1 a S3 (árboles con accesiones) |
| 09_suplementario.py | Supplementary_Material.xlsx |
| 11_figura_S4_modulos.py | Figura S4 |
| 12_sintenia_clinker.py | Panel B de la Figura 3 (datos de clinker) |
| 13_componer_figuras.py | Ensambla Figure1.svg y Figure3.svg |

Pasos realizados en servidores web, sin script: BLASTP (NCBI), InterProScan 6 (EMBL-EBI), CodonTree/RAxML (BV-BRC) y modelo ESMFold.

## Licencia
Código bajo MIT; datos y figuras bajo CC BY 4.0.
