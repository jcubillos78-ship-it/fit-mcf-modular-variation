"""Figura 4: triada catalitica de C58 en el modelo ESMFold (FitD 1275-1600), proyeccion 2D vectorial."""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype']='none'; plt.rcParams['font.family']=['Arial','Helvetica','Liberation Sans','DejaVu Sans']
want={126:('Cys',1400),215:('His',1489),234:('Asp',1508)}
atoms={}; plddt={}; bf={126:[],215:[],234:[]}
for l in open('data/FitD_1275-1600.pdb'):
    if l.startswith('ATOM'):
        r=int(l[22:26]); 
        if r in want:
            n=l[12:16].strip(); atoms[(r,n)]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]); bf[r].append(float(l[60:66]))
plddt={r:float(np.mean(v)) for r,v in bf.items()}  # pLDDT del residuo = media de sus atomos (0-1)
sc=100
side={126:['CA','CB','SG'],215:['CA','CB','CG','ND1','CE1','NE2','CD2'],234:['CA','CB','CG','OD1','OD2']}
bonds={126:[('CA','CB'),('CB','SG')],215:[('CA','CB'),('CB','CG'),('CG','ND1'),('ND1','CE1'),('CE1','NE2'),('NE2','CD2'),('CD2','CG')],234:[('CA','CB'),('CB','CG'),('CG','OD1'),('CG','OD2')]}
pts=np.array([atoms[(r,n)] for r in side for n in side[r]]); c=pts.mean(0)
u,s,vt=np.linalg.svd(pts-c); e1,e2=vt[0],vt[1]
P=lambda r,n:((atoms[(r,n)]-c)@e1,(atoms[(r,n)]-c)@e2)
d=lambda a,b:np.linalg.norm(atoms[a]-atoms[b])
def pc(v):
    v=v*(100 if sc==100 else 1)
    return '#0053D6' if v>=90 else '#65CBF3' if v>=70 else '#FFDB13' if v>=50 else '#FF7D45'
el={'C':'#555555','N':'#2166ac','O':'#d6301f','S':'#c9a400'}
fig,ax=plt.subplots(figsize=(5.2,4.2)); ax.set_aspect('equal'); ax.axis('off')
for r in side:
    col=pc(plddt[r])
    for a,b in bonds[r]:
        (x1,y1),(x2,y2)=P(r,a),P(r,b); ax.plot([x1,x2],[y1,y2],color=col,lw=4,solid_capstyle='round',zorder=1)
    for n in side[r]:
        x,y=P(r,n); e=n[0]; ax.scatter([x],[y],s=170 if e!='C' else 90,c=el[e],edgecolors='white',linewidths=0.8,zorder=3)
        if e!='C': ax.text(x,y,e,color='white',fontsize=6.5,ha='center',va='center',zorder=4,fontweight='bold')
def hb(a,b,lab,dx=0,dy=0):
    (x1,y1),(x2,y2)=P(*a),P(*b); ax.plot([x1,x2],[y1,y2],color='black',lw=1.2,ls=(0,(4,3)),zorder=2)
    ax.text((x1+x2)/2+dx,(y1+y2)/2+dy,lab,fontsize=8.5,ha='center',va='center',bbox=dict(fc='white',ec='none',pad=1.2),zorder=5)
dA=d((126,'SG'),(215,'ND1')); dB=d((215,'NE2'),(234,'OD1'))
hb((126,'SG'),(215,'ND1'),f'{dA:.2f} Å'.replace('.',','),dy=0.0)
hb((215,'NE2'),(234,'OD1'),f'{dB:.2f} Å'.replace('.',','),dy=0.0)
for r,(aa,g) in want.items():
    xs=[P(r,n)[0] for n in side[r]]; ys=[P(r,n)[1] for n in side[r]]
    ax.text(np.mean(xs),max(ys)+0.55 if False else np.mean(ys)-0.0,'',fontsize=1)
ann={126:'Cys1400',215:'His1489',234:'Asp1508'}
for r,t in ann.items():
    x,y=P(r,{126:'SG',215:'CE1',234:'CG'}[r]); ax.annotate(f'{t}\npLDDT {plddt[r]*(100 if sc==100 else 1):.1f}'.replace('.',','),(x,y),xytext=((x+1.7,y+0.1) if r==126 else (x+(0.9 if r==234 else -0.9),y+1.5)),fontsize=8.5,ha='center',arrowprops=dict(arrowstyle='-',color='#888888',lw=0.7))
# leyenda pLDDT
for i,(lb,cl) in enumerate((('≥ 90','#0053D6'),('70-90','#65CBF3'),('50-70','#FFDB13'))):
    ax.add_patch(plt.Rectangle((0.02+i*0.2,0.02),0.03,0.03,transform=ax.transAxes,color=cl)); ax.text(0.06+i*0.2,0.035,lb,transform=ax.transAxes,fontsize=7,va='center')
ax.text(0.02,0.085,'Color de enlaces: pLDDT',transform=ax.transAxes,fontsize=7,va='center')
fig.savefig('Figures_Vector/Figura_4_triada_C58.svg',format='svg',bbox_inches='tight')
print('SG-ND1',dA,'NE2-OD1',dB,{r:plddt[r] for r in plddt})
