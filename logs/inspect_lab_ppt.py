from pathlib import Path
import zipfile
import xml.etree.ElementTree as E

p = Path('C:/Users/정호영/Downloads/GWAS_First_Lab_Presentation_CRC_2.pptx')
ns = {'a':'http://schemas.openxmlformats.org/drawingml/2006/main','p':'http://schemas.openxmlformats.org/presentationml/2006/main'}
with zipfile.ZipFile(p) as z:
    print('SIZE', E.fromstring(z.read('ppt/presentation.xml')).find('p:sldSz',ns).attrib)
    names = sorted([n for n in z.namelist() if n.startswith('ppt/slides/slide') and n.endswith('.xml')], key=lambda n:int(n.rsplit('slide',1)[1].split('.')[0]))
    for n in names:
        root = E.fromstring(z.read(n))
        print(n, ' | '.join(t.text or '' for t in root.findall('.//a:t',ns)))
    for i in (9,10,11):
        name = f'ppt/slides/slide{i}.xml'
        if name not in z.namelist(): continue
        print('\nLAYOUT', i)
        root = E.fromstring(z.read(name))
        for s in root.findall('.//p:spTree/*',ns):
            print(E.tostring(s, encoding='unicode')[:5500])
    print('FONT', E.fromstring(z.read('ppt/theme/theme1.xml')).find('.//a:fontScheme',ns).attrib)
