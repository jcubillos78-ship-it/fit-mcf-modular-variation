"""Cover_Letter.pdf (1 pagina)."""
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
st=ParagraphStyle('b',fontName='Helvetica',fontSize=10.3,leading=14,spaceAfter=7)
doc=SimpleDocTemplate('Cover_Letter.pdf',pagesize=letter,leftMargin=1*inch,rightMargin=1*inch,topMargin=0.8*inch,bottomMargin=0.8*inch,title='Carta de presentación',author='Juan David Cubillos Martínez')
P=lambda t:Paragraph(t,st)
body=[
P('Bogotá D.C., 18 de noviembre de 2026'),
P('Señores Editores<br/><i>Current Microbiology</i>'),
P('<b>Asunto: envío del manuscrito «Variación modular de la toxina Fit/Mcf en <i>Pseudomonas</i>, <i>Photorhabdus</i> y <i>Xenorhabdus</i>: arquitecturas de dominios, historias filogenéticas discordantes y especificidad de hospedero»</b>'),
P('Estimados editores:'),
P('Les presento para su consideración un artículo de investigación original sobre la toxina insecticida Fit/Mcf, presente en pseudomonas benéficas de plantas y en los simbiontes entomopatógenos <i>Photorhabdus</i> y <i>Xenorhabdus</i>. Estas cepas difieren mucho en los insectos que logran matar, y el trabajo pregunta si esa diferencia se acompaña de variación en la composición modular de la propia toxina.'),
P('Con 300 proteínas homólogas anotadas con InterProScan, el trabajo encuentra cuatro arquitecturas de dominios: el módulo de poro está presente en prácticamente todas, mientras que los módulos efectores del extremo amino y la proteasa se intercambian. Las filogenias independientes de tres módulos (poro, extremo amino y proteasa C58, con 10 000 réplicas de bootstrap) muestran historias discordantes dentro de una misma molécula: el poro de la variante con peptidasa C80 forma un grupo con soporte 100 que no se reproduce en el extremo amino. La comparación con una filogenia genómica de once genomas es congruente para seis cepas solo si se elige la copia de la toxina de arquitectura equivalente (distancia de Robinson-Foulds 0 frente a 2), y el trabajo declara ese límite. Un modelo estructural del dominio C58 de FitD de <i>P. protegens</i> CHA0 sitúa la tríada Cys1400, His1489 y Asp1508 a 2,82 y 2,86 Å, en disposición compatible con un sitio activo funcional. Una matriz que cruza cepa, arquitectura, insecto y vía de exposición resume la evidencia publicada y muestra que la especificidad de hospedero aún no puede asignarse a una arquitectura concreta.'),
P('Creemos que el trabajo interesa a los lectores de la revista por tres razones: documenta con filogenias explícitas el mosaico modular propuesto para esta toxina, ofrece una anotación reproducible de sus variantes y delimita con franqueza qué conclusiones permiten los datos actuales. Las secuencias proceden de bases públicas y todos los alineamientos, árboles, scripts y el archivo de entorno están en un repositorio con identificador permanente.'),
P('El manuscrito es original, no se ha publicado ni está en evaluación en otra revista, y el autor declara no tener conflictos de intereses. El uso de asistentes de inteligencia artificial generativa se declara en el manuscrito conforme a la política de la revista.'),
P('Agradezco su atención y quedo atento a sus comentarios.'),
Spacer(1,6),
P('Atentamente,<br/><br/><b>Juan David Cubillos Martínez</b><br/>Programa de Biotecnología, Universidad Antonio Nariño, Bogotá D.C., Colombia<br/>jcubillos@uan.edu.co'),
]
doc.build(body)
