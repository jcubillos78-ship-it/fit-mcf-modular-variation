"""Render 3D (PyMOL, ray-tracing, 300 DPI) de la triada catalitica C58 del modelo ESMFold de FitD (fragmento 1275-1600).
Numeracion local = global - 1274: Cys1400=126, His1489=215, Asp1508=234. Uso: pymol -cq scripts/07b_render_triada_pymol.py"""
from pymol import cmd
cmd.load('data/FitD_1275-1600.pdb','fitd')
cmd.bg_color('white'); cmd.set('ray_opaque_background',1); cmd.set('antialias',2)
cmd.set('ray_trace_mode',1); cmd.set('ray_trace_gain',0.0); cmd.set('ray_shadow',0)
cmd.hide('everything'); cmd.show('cartoon','fitd'); cmd.set('cartoon_transparency',0.35); cmd.set('cartoon_fancy_helices',1); cmd.set('ray_shadows',0); cmd.set('depth_cue',0)
# pLDDT (escala 0-1 en la columna B) -> paleta AlphaFold apta para daltonismo parcial
cmd.color('0x0053D6','fitd and b>=0.9'); cmd.color('0x65CBF3','fitd and b<0.9 and b>=0.7')
cmd.color('0xFFDB13','fitd and b<0.7 and b>=0.5'); cmd.color('0xFF7D45','fitd and b<0.5')
tri='fitd and resi 126+215+234 and not (name N+C+O)'
cmd.show('sticks',tri); cmd.set('stick_radius',0.18)
cmd.color('0x2B2B2B',tri+' and elem C'); cmd.color('0x2166AC',tri+' and elem N'); cmd.color('0xD6301F',tri+' and elem O'); cmd.color('0xC9A400',tri+' and elem S')
cmd.set('stick_radius',0.22)
cmd.distance('d1','fitd and resi 126 and name SG','fitd and resi 215 and name ND1')
cmd.distance('d2','fitd and resi 215 and name NE2','fitd and resi 234 and name OD1')
cmd.set('dash_color','black'); cmd.set('dash_gap',0.25); cmd.set('dash_radius',0.07)
cmd.set('label_font_id',7); cmd.set('label_size',-1.6); cmd.set('label_color','black'); cmd.set('label_distance_digits',2)
cmd.pseudoatom('l1',pos=[0,0,0]); cmd.hide('everything','l1')
for nm,res,at,txt in (('L1',126,'CB','Cys1400'),('L2',215,'CB','His1489'),('L3',234,'CB','Asp1508')):
    cmd.label(f'fitd and resi {res} and name {at}',f'"{txt}"')
cmd.set('label_position',[1.5,1.2,2.5])
cmd.set('label_color','black','d1'); cmd.set('label_color','black','d2')
cmd.hide('cartoon','fitd'); cmd.select('vecino','byres (fitd within 7 of ('+tri+'))'); cmd.show('cartoon','vecino and not resi 126+215+234')
cmd.set('cartoon_transparency',0.65,'vecino'); cmd.set('cartoon_tube_radius',0.12); cmd.set('cartoon_loop_radius',0.12)
cmd.set('label_outline_color','white'); cmd.set('label_shadow_mode',2)
cmd.orient(tri); cmd.turn('x',-35); cmd.zoom(tri,7); cmd.clip('slab',25)
cmd.png('render/Figure4A_triada_3D.png',width=2400,height=1800,dpi=300,ray=1)
print('d1',cmd.get_distance('fitd and resi 126 and name SG','fitd and resi 215 and name ND1'),'d2',cmd.get_distance('fitd and resi 215 and name NE2','fitd and resi 234 and name OD1'))
