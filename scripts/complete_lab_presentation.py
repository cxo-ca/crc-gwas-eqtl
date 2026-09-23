"""Append an editable tenth slide while preserving original slides byte-for-byte."""
import copy
import hashlib
import json
from pathlib import Path
import zipfile
import xml.etree.ElementTree as E

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path('C:/Users/정호영/Downloads/GWAS_First_Lab_Presentation_CRC_2.pptx')
OUT = ROOT / 'reports/presentation'
DEST = OUT / 'GWAS_First_Lab_Presentation_CRC_completed.pptx'
NS = {'p':'http://schemas.openxmlformats.org/presentationml/2006/main',
      'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
      'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
REL = 'http://schemas.openxmlformats.org/package/2006/relationships'
CT = 'http://schemas.openxmlformats.org/package/2006/content-types'
for prefix, uri in NS.items(): E.register_namespace(prefix, uri)


def tag(n):
    prefix, name = n.split(':')
    return '{' + NS[prefix] + '}' + name


def sub(parent, element_name, **attrs):
    return E.SubElement(parent, tag(element_name), {k:str(v) for k,v in attrs.items()})


def xml(root):
    return E.tostring(root, encoding='utf-8', xml_declaration=True)


def emu(pt): return str(round(pt * 12700))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    result = json.loads((ROOT / 'results/final/final_summary.json').read_text(encoding='utf-8'))
    assert result['candidate_genes'] == 102 and result['analyzed_groups'] == 144
    assert len(result['prior_stable_genes']) == 9
    root = E.Element(tag('p:sld'))
    cs = sub(root, 'p:cSld', name='실제 분석 진행 결과')
    bg = sub(sub(cs, 'p:bg'), 'p:bgPr')
    sub(sub(bg, 'a:solidFill'), 'a:srgbClr', val='FFFFFF'); sub(bg, 'a:effectLst')
    tree = sub(cs, 'p:spTree')
    nv = sub(tree, 'p:nvGrpSpPr'); sub(nv, 'p:cNvPr', id=1, name=''); sub(nv, 'p:cNvGrpSpPr'); sub(nv, 'p:nvPr')
    xf = sub(sub(tree, 'p:grpSpPr'), 'a:xfrm')
    for n, values in [('off', {'x':0,'y':0}),('ext', {'cx':0,'cy':0}),('chOff', {'x':0,'y':0}),('chExt', {'cx':0,'cy':0})]:
        sub(xf, 'a:' + n, **values)
    shapes = []

    def shape(x,y,w,h,text='',size=18,bold=False,color='172D43',fill=None,name=''):
        ident = len(shapes) + 2
        s = sub(tree, 'p:sp')
        nv = sub(s, 'p:nvSpPr'); sub(nv, 'p:cNvPr', id=ident, name=name or f'result-{ident}')
        sub(nv, 'p:cNvSpPr', txBox='1' if text else '0'); sub(nv, 'p:nvPr')
        sp = sub(s, 'p:spPr'); xf = sub(sp, 'a:xfrm')
        sub(xf, 'a:off', x=emu(x), y=emu(y)); sub(xf, 'a:ext', cx=emu(w), cy=emu(h))
        sub(sub(sp, 'a:prstGeom', prst='rect'), 'a:avLst')
        if fill: sub(sub(sp, 'a:solidFill'), 'a:srgbClr', val=fill)
        else: sub(sp, 'a:noFill')
        sub(sub(sp, 'a:ln', w=0), 'a:noFill')
        if text:
            body = sub(s, 'p:txBody')
            bp = sub(body, 'a:bodyPr', wrap='none', anchor='ctr', lIns=0,tIns=0,rIns=0,bIns=0)
            sub(bp, 'a:noAutofit'); sub(body, 'a:lstStyle')
            for line in text.split('\n'):
                p = sub(body, 'a:p'); pp = sub(p, 'a:pPr', algn='l')
                sub(pp,'a:buNone')
                r = sub(p,'a:r'); rp = sub(r,'a:rPr', lang='ko-KR', sz=round(size*100), b=int(bold))
                sub(sub(rp,'a:solidFill'),'a:srgbClr',val=color)
                for n in ('latin','ea','cs'): sub(rp,'a:'+n,typeface='Noto Sans KR')
                sub(r,'a:t').text=line
                sub(p,'a:endParaRPr', lang='ko-KR', sz=round(size*100))
        shapes.append(dict(id=ident,x=x,y=y,w=w,h=h,text=text,size=size,bold=bold,color=color,fill=fill))

    shape(48,29,864,55,'실제 분석 진행 결과: 후보 재정리까지 완료',31.5,True)
    shape(49,87,864,24,'공개 자료 확보 → 기준 정렬 → 실제 coloc → 누락 점검·최종 보고서',17,color='556675')
    cards = [
        ('01  데이터 확보','100개 유전자','36개 구간 · 대장 2조직','GWAS·eQTL 양쪽 자료 확보'),
        ('02  기준 정렬','111,807 SNP','GRCh38·대립유전자 대조','중복·결측·모호한 변이 QC'),
        ('03  실제 coloc','144개 조합','3개 prior 설정 · 총 432회','전체 사후확률 독립 검증'),
        ('04  후보 재정리','9개 유전자','세 prior 모두 PP.H4 ≥ 0.8','조건부 지지 · 11개 조합'),
    ]
    for i,(label,metric,line1,line2) in enumerate(cards):
        x=48+i*219
        shape(x,130,207,132,fill='EAF4FB' if i==3 else 'F2F6FA')
        shape(x,130,207,3,fill='1675B8')
        shape(x+13,142,182,24,label,16,True,color='1675B8')
        shape(x+13,172,182,35,metric,24,True)
        shape(x+13,213,184,20,line1,13.5)
        shape(x+13,236,184,19,line2,12.9,color='556675')
    shape(49,282,830,26,'조건부 지지 후보 9개',20,True,color='1675B8')
    shape(49,316,850,25,'ACTR1B · LRP1 · COLCA2 · LMOD1 · ATF1',19)
    shape(49,344,850,25,'COLCA1 · LINC02044 · MZF1-AS1 · RP11-351I24.1',19)
    shape(48,389,864,104,fill='F2F6FA')
    shape(61,397,840,22,'누락 점검과 해석 범위',16.5,True,color='1675B8')
    shape(61,425,840,19,'최강 GWAS 변이 누락 137개 조합: eQTL 원본 부재 131 · cis 범위 밖 6',15)
    shape(61,447,840,19,'점검 범위 내 처리 누락 0건 · 모호한 lead 4개 제외 · X염색체 3개 조합 보류',14.5)
    shape(61,469,840,19,'원인 유전자 확정 아님: LD·복수 신호·독립 검증 필요',15,True)
    shape(49,509,818,16,'GCST90129505 · GTEx v8  |  분석·누락 점검·보고서 완료  ·  2026.09.23',10,color='556675')
    shape(876,508.5,36,18,'10',11.25,color='556675')
    sub(sub(root, 'p:clrMapOvr'), 'a:masterClrMapping')
    notes_text = [
        '9페이지의 초기 계획 이후 수행한 실제 분석 결과입니다. 분석 범위는 36개 상염색체 거리 구간으로 확대됐습니다.',
        'GCST90129505의 대장암 GWAS와 GTEx v8 Colon_Sigmoid/Colon_Transverse를 사용했습니다. 102개 발견 후보 중 100개 유전자·144개 유전자–조직 조합에 양쪽 자료를 확보했습니다.',
        '확보 자료는 GWAS 298,636행, eQTL 876,402행입니다. GRCh38 위치·변이 ID·대립유전자를 확인하고 중복·결측·모호한 SNP를 처리했습니다. 고유 공통 SNP 111,807개, 유전자·조직별 반복을 포함한 공통 통계 630,083행입니다.',
        'coloc 5.2.3의 단일 원인 변이 ABF 방법을 사용했습니다. 요청 p1=p2=10^-4, p12=10^-6, 10^-5, 10^-4로 432회 실행했습니다. 29개 조합에서 SNP 수에 따른 prior 자동 조정이 있었으며 실제 적용값을 기록했습니다. beta/SE² 및 공식 정규화 발현행렬에서 직접 측정한 sdY를 사용했습니다.',
        '기본 prior에서 18개 유전자·23개 조합이 PP.H4 ≥ 0.8을 충족했고, 시험한 세 prior 모두에서 유지한 것은 9개 유전자·11개 조합입니다. ACTR1B와 ATF1은 두 조직에서 기준을 충족했습니다. 조직 간 공유 표본 때문에 독립 재현은 아닙니다.',
        '전체 432개 posterior를 별도 Python 계산으로 검증했습니다. 최대 절대 차이는 약 5×10^-16입니다. 이는 계산 재현성 확인이며 생물학적 검증을 뜻하지 않습니다.',
        '누락 점검은 SHA-256 기록 181건과 원본 색인 좌표 조회 96건을 확인했습니다. 최강 GWAS 변이가 없는 137개 조합 중 131개는 사용한 eQTL 원본에 위치가 없었고 6개는 유전자 TSS ±1 Mb 밖이었습니다. 최강 eQTL 변이 누락 4개 조합도 해당 GWAS 원본에 위치가 없었습니다. 확인 범위에서 처리 누락 0건, 복구 행 0개입니다.',
        'A/T·C/G lead rs2581817, rs209142, rs9271363, rs7398375는 의도적으로 제외했습니다. SHROOM2 두 조직과 GPR143 Transverse는 X염색체 GWAS 부재로 보류했습니다. 점검한 137개 GWAS peak의 EAF는 모두 NA여서 빈도로 문제를 해결하지 못했습니다.',
        '102개 후보의 최종 분류: 세 prior 지지·조건부 9개, 기본 prior 지지·민감 9개, 큰 prior에서만 지지 7개, 기준 미충족 75개, 자료 부재 2개입니다. 모든 후보와 147개 대상 조합을 최종 보고서에 보존했습니다.',
        '현재 결과는 공통 SNP 집합, 거리로 정한 구간, 단일 원인 변이 가정에 조건부입니다. LD/재조합 경계·복수 신호·원본 극단 beta의 품질·독립 reference FASTA·indel 복구·독립 및 기능 검증은 남아 있습니다. 인과적 매개나 원인 유전자 확정으로 해석하지 않습니다.',
        '분석 출처: reports/final_report.html; results/final/candidates.tsv; results/coloc/20260923_094957_real/.',
        '원본: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90129001-GCST90130000/GCST90129505/harmonised/',
        'eQTL: https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/',
    ]
    with zipfile.ZipFile(SOURCE) as zin:
        changed = {}
        presentation = E.fromstring(zin.read('ppt/presentation.xml'))
        ids = presentation.find('p:sldIdLst',NS)
        assert len(ids)==9
        rels = E.fromstring(zin.read('ppt/_rels/presentation.xml.rels'))
        rid='rId'+str(max(int(r.get('Id')[3:]) for r in rels)+1)
        E.SubElement(rels,'{'+REL+'}Relationship',Id=rid,Type=NS['r']+'/slide',Target='slides/slide10.xml')
        E.SubElement(ids,tag('p:sldId'),{'id':str(max(int(s.get('id')) for s in ids)+1),tag('r:id'):rid})
        changed['ppt/presentation.xml']=xml(presentation)
        changed['ppt/_rels/presentation.xml.rels']=xml(rels)
        changed['ppt/slides/slide10.xml']=xml(root)
        slide_rels=E.fromstring(zin.read('ppt/slides/_rels/slide9.xml.rels'))
        for r in slide_rels:
            if r.get('Type').endswith('/notesSlide'): r.set('Target','../notesSlides/notesSlide10.xml')
        changed['ppt/slides/_rels/slide10.xml.rels']=xml(slide_rels)
        notes=E.fromstring(zin.read('ppt/notesSlides/notesSlide9.xml'))
        for s in notes.findall('.//p:sp',NS):
            ph=s.find('.//p:ph',NS)
            if ph is not None and ph.get('type')=='body':
                tx=s.find('p:txBody',NS)
                for p in list(tx.findall('a:p',NS)): tx.remove(p)
                for line in notes_text:
                    p=sub(tx,'a:p'); r=sub(p,'a:r'); sub(r,'a:t').text=line
        changed['ppt/notesSlides/notesSlide10.xml']=xml(notes)
        nr=E.fromstring(zin.read('ppt/notesSlides/_rels/notesSlide9.xml.rels'))
        for r in nr:
            if r.get('Type').endswith('/slide'):r.set('Target','../slides/slide10.xml')
        changed['ppt/notesSlides/_rels/notesSlide10.xml.rels']=xml(nr)
        types=E.fromstring(zin.read('[Content_Types].xml'))
        for part,kind in [('slides/slide10','slide'),('notesSlides/notesSlide10','notesSlide')]:
            E.SubElement(types,'{'+CT+'}Override',PartName='/ppt/'+part+'.xml',ContentType='application/vnd.openxmlformats-officedocument.presentationml.'+kind+'+xml')
        changed['[Content_Types].xml']=xml(types)
        app=E.fromstring(zin.read('docProps/app.xml'))
        ans={'ep':'http://schemas.openxmlformats.org/officeDocument/2006/extended-properties','vt':'http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes'}
        app.find('ep:Slides',ans).text='10'; app.find('ep:Notes',ans).text='10'
        hp=app.findall('.//ep:HeadingPairs/vt:vector/vt:variant/vt:i4',ans); hp[-1].text='10'
        titles=app.find('ep:TitlesOfParts/vt:vector',ans)
        E.SubElement(titles,'{'+ans['vt']+'}lpstr').text='실제 분석 진행 결과'
        titles.set('size',str(len(titles)))
        changed['docProps/app.xml']=xml(app)
        with zipfile.ZipFile(DEST,'w',zipfile.ZIP_DEFLATED) as zout:
            for name in zin.namelist():zout.writestr(name,changed.pop(name,zin.read(name)))
            for name,data in changed.items():zout.writestr(name,data)
        with zipfile.ZipFile(DEST) as zout:
            assert zout.testzip() is None
            for i in range(1,10):
                for name in (f'ppt/slides/slide{i}.xml',f'ppt/notesSlides/notesSlide{i}.xml'):
                    assert zout.read(name)==zin.read(name)
            for name in zout.namelist():
                if name.endswith(('.xml','.rels')): E.fromstring(zout.read(name))
    for s in shapes:
        assert s['x']>=0 and s['y']>=0 and s['x']+s['w']<=960 and s['y']+s['h']<=540
    (OUT/'slide10_layout.json').write_text(json.dumps(shapes,ensure_ascii=False,indent=2),encoding='utf-8')
    (OUT/'slide10_speaker_notes.txt').write_text('\n\n'.join(notes_text),encoding='utf-8')
    print(DEST)
    print('PASS: 10 slides, original slides/notes 1-9 unchanged, ZIP/XML valid, new shapes within slide bounds')


if __name__=='__main__':main()
