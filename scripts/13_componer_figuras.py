"""Compone Figure1.svg (A arbol genomico, B mapa de arquitecturas) y Figure3.svg (A tanglegrama, B sintenia) con letras de panel."""
import re
def info(p):
    t=open(p).read(); m=re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"',t); inner=t[t.index('<svg'):]
    body=re.sub(r'^<svg[^>]*>','',inner,1).rsplit('</svg>',1)[0]
    return float(m.group(1)),float(m.group(2)),body
def comp(panels,out,gap=14):
    W=max(info(p)[0] for p in panels)+22; y=4; parts=[]
    for L,p in zip('ABCD',panels):
        w,h,b=info(p)
        parts.append(f'<text x="2" y="{y+11:.1f}" font-family="Arial, Helvetica, sans-serif" font-size="12" font-weight="bold">{L}</text>')
        parts.append(f'<svg x="{18+(W-22-w)/2:.1f}" y="{y:.1f}" width="{w}" height="{h}" viewBox="0 0 {w} {h}">{b}</svg>'); y+=h+gap
    H=y
    open(out,'w').write(f'<?xml version="1.0" encoding="utf-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{W:.1f}pt" height="{H:.1f}pt" viewBox="0 0 {W:.1f} {H:.1f}" version="1.1">{"".join(parts)}</svg>')
comp(['render/paneles/Figura_1_arbol_genomico.svg','render/paneles/Figura_2_arquitecturas.svg'],'Figures_Vector/Figure1.svg')
comp(['render/paneles/Figura_3_tanglegrama_N-C-T.svg','render/sintenia_B.svg'],'Figures_Vector/Figure3.svg')
