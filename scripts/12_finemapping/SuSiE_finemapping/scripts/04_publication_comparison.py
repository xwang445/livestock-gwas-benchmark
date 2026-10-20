#!/usr/bin/env python3
"""Matched-locus publication comparison; only reads already evaluated results."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/susie_publication_mpl')
from pathlib import Path
import hashlib, json, zipfile
from xml.sax.saxutils import escape
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
W=Path(__file__).resolve().parents[1]
OUT=W/'publication'; OUT.mkdir(exist_ok=True)
source=W/'results/comparison_common_converged_loci.csv'
D=pd.read_csv(source)
D['fine_mapping']=D.analysis.map({'Standard SuSiE-RSS':'SuSiE','LOCO-adjusted SuSiE':'SuSiE-adj'})
assert D.fine_mapping.notna().all() and len(D)==40
assert not D.duplicated(['fine_mapping','method','phenotype']).any()
METHODS=['BOLT-LMM','REGENIE-LOCO']; FM=['SuSiE','SuSiE-adj']
COL={'SuSiE':'#0072B2','SuSiE-adj':'#D55E00'}
METRICS=['coverage','power','median_CS_size']
for method in METHODS:
    a=D[(D.method==method)&(D.fine_mapping=='SuSiE')].set_index('phenotype').sort_index()
    b=D[(D.method==method)&(D.fine_mapping=='SuSiE-adj')].set_index('phenotype').sort_index()
    assert len(a)==len(b)==10
    assert a.index.equals(b.index)
    assert (a.n_loci==b.n_loci).all() and (a.n_causal_analyzed==b.n_causal_analyzed).all()
D.to_csv(OUT/'source_data.csv',index=False)
rows=[]; numeric=[]
for method in METHODS:
    for fm in FM:
        q=D[(D.method==method)&(D.fine_mapping==fm)]
        rows.append({'GWAS input':method,'Fine-mapping':fm,'Locus fits':f'{q.n_loci.sum():,}',
            '95% CSs':f'{q.n_CS.sum():,}',
            'Coverage (%)':f'{q.coverage.mean():.2f} ± {q.coverage.std(ddof=1):.2f}',
            'Power (%)':f'{q.power.mean():.2f} ± {q.power.std(ddof=1):.2f}',
            'CS size (SNPs)':f'{q.median_CS_size.mean():.2f} ± {q.median_CS_size.std(ddof=1):.2f}'})
        r=dict(GWAS_input=method,fine_mapping=fm,n_phenotypes=len(q),n_locus_fits=int(q.n_loci.sum()),n_credible_sets=int(q.n_CS.sum()))
        for metric in METRICS:
            r[metric+'_mean']=float(q[metric].mean());r[metric+'_sd']=float(q[metric].std(ddof=1))
        numeric.append(r)
T=pd.DataFrame(rows); T.to_csv(OUT/'Table_1.csv',index=False)
pd.DataFrame(numeric).to_csv(OUT/'Table_1_numeric.csv',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':9,
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'xtick.major.width':.7,
    'ytick.major.width':.7,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
fig,axs=plt.subplots(1,3,figsize=(7.2,3.05))
for ax,key,title,ylabel,letter in zip(axs,METRICS,['Coverage','Power','Resolution'],
    ['Credible sets with a causal SNP (%)','Analyzed causal SNPs recovered (%)','Median SNPs per credible set'],'ABC'):
    for i,method in enumerate(METHODS):
        for j,fm in enumerate(FM):
            q=D[(D.method==method)&(D.fine_mapping==fm)].copy()
            q['number']=q.phenotype.str.replace('pheno','',regex=False).astype(int)
            y=q.sort_values('number')[key].to_numpy(); x=i+(j-.5)*.30
            ax.scatter(x+np.linspace(-.055,.055,len(y)),y,s=12,facecolors='white',edgecolors=COL[fm],lw=.65,alpha=.8,zorder=3)
            ax.errorbar(x,y.mean(),yerr=y.std(ddof=1),fmt='_',markersize=11,markeredgewidth=1.7,
                color=COL[fm],elinewidth=1.2,capsize=3,capthick=1.0,zorder=4)
    ax.set_xticks([0,1]); ax.set_xticklabels(['BOLT-LMM','REGENIE\n(LOCO)'])
    ax.set_xlim(-.5,1.5);ax.set_title(title,pad=10);ax.set_ylabel(ylabel,labelpad=4)
    ax.text(-.20,1.06,letter,transform=ax.transAxes,fontweight='bold',fontsize=11)
    ax.yaxis.grid(True,color='#E8E8E8',lw=.5);ax.set_axisbelow(True)
    if key=='coverage':
        ax.set_ylim(0,103);ax.set_yticks([0,20,40,60,80,100]);ax.axhline(95,color='#777777',ls='--',lw=.7)
        ax.text(1.42,96.5,'95% nominal',ha='right',fontsize=6.7,color='#666666')
    elif key=='power': ax.set_ylim(0,25);ax.set_yticks([0,5,10,15,20,25])
    else: ax.set_ylim(0,8);ax.set_yticks([0,2,4,6,8])
handles=[Line2D([0],[0],color=COL[fm],marker='_',lw=1.2,markersize=10,label=fm) for fm in FM]
fig.legend(handles=handles,ncol=2,loc='lower center',bbox_to_anchor=(.5,.065),frameon=False,handlelength=1.4,columnspacing=2)
fig.text(.5,.025,'Points: 10 phenotypes  |  Bars: mean ± SD  |  Matched loci',ha='center',fontsize=7,color='#555555')
fig.subplots_adjust(left=.085,right=.985,bottom=.28,top=.88,wspace=.51)
for ext in ['pdf','svg','png','tiff']:
    kw={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
    fig.savefig(OUT/f'Figure_1.{ext}',dpi=600,**kw)
plt.close(fig)
figure_caption='''Figure 1. Fine-mapping performance of SuSiE and LOCO SuSiE-adj using BOLT-LMM and REGENIE-LOCO association statistics. (A) Empirical coverage of nominal 95% credible sets (CSs), defined as the percentage of reported CSs containing at least one simulated causal SNP; the dashed line marks 95%. (B) Conditional fine-mapping power, defined as the percentage of unique causal SNPs analyzed within the matched converged loci that occur in at least one reported CS. SNPs are deduplicated within each GWAS method and phenotype. (C) Median number of SNPs per reported CS, calculated separately for each phenotype. Open points represent the 10 simulated phenotypes; central horizontal bars and whiskers represent the mean and sample standard deviation across phenotypes. All 6,732 BOLT-LMM and 6,569 REGENIE-LOCO locus fits across the 10 phenotypes converged in both analyses. All locus types, including null-region loci, are included; locus-fit counts count a genomic region again when analyzed for a different phenotype. The two fine-mapping pipelines share loci, SNPs, harmonized z-scores, L=10, 95% CSs, and minimum absolute within-CS correlation of 0.5. SuSiE uses ordinary in-sample LD, N=20,000, and fixed residual variance. SuSiE-adj refers to the LOCO-adjusted analysis, with method-specific adjusted LD, locus-specific effective sample size, and estimated residual variance. This is a pipeline comparison, not an isolated test of LD adjustment. Power is conditional on the analyzed causal SNPs and is not genome-wide discovery power. Causal status is defined by the simulated set of 10,000 SNPs and used only for evaluation. The same simulated causal set and genotype panel are shared across phenotypes.\n'''
table_caption='''Table 1. Summary of SuSiE and LOCO SuSiE-adj fine-mapping performance. Locus fits and 95% credible sets (CSs) are totals across the 10 phenotypes. Coverage and power are percentages reported as mean ± sample SD across phenotypes (n=10), with each phenotype weighted equally. CS size is the mean ± SD of the phenotype-specific median numbers of SNPs per CS, not the median of all pooled CSs. Coverage is the fraction of CSs containing at least one true causal SNP. Power is the fraction of analyzed unique causal SNPs recovered in any CS; its denominator is limited to matched converged loci. All reported fits converged in both analyses. SuSiE uses ordinary LD, N=20,000, and fixed residual variance; SuSiE-adj denotes LOCO-adjusted LD, effective sample size, and estimated residual variance. The methods share GWAS inputs, loci, and other fine-mapping settings. See Figure 1 for phenotype-level values.\n'''
(OUT/'Figure_1_caption.txt').write_text(figure_caption)
(OUT/'Table_1_caption.txt').write_text(table_caption)
# Print-quality table preview; CSV, LaTeX, and Word versions are editable.
fig,ax=plt.subplots(figsize=(10,3.1));ax.axis('off')
ax.text(0,.98,'Table 1. SuSiE and LOCO SuSiE-adj fine-mapping performance',fontsize=11,fontweight='bold',va='top')
labels=['GWAS input','Fine-mapping','Locus fits','95% CSs','Coverage (%)','Power (%)','CS size (SNPs)']
tab=ax.table(cellText=T.values,colLabels=labels,cellLoc='center',colLoc='center',bbox=[0,.37,1,.48],colWidths=[.16,.13,.10,.10,.17,.17,.17])
tab.auto_set_font_size(False);tab.set_fontsize(8.5)
for (r,c),cell in tab.get_celld().items():
    cell.set_edgecolor('#BBBBBB');cell.set_linewidth(.4)
    if r==0:cell.set_facecolor('#EDEFF1');cell.set_text_props(fontweight='bold')
    elif r in [3,4]:cell.set_facecolor('#F7F7F7')
ax.text(0,.27,'Counts are totals across 10 phenotypes; performance metrics are mean ± SD (n = 10).\n'
    'Coverage: CSs containing a causal SNP. Power: analyzed causal SNPs recovered in a CS.\n'
    'CS size: phenotype-specific median SNP count. All fits converged in both analyses.\n'
    'SuSiE-adj uses LOCO-adjusted LD and estimated residual variance; SuSiE uses ordinary LD and fixed residual variance.',
    fontsize=8,va='top',linespacing=1.6)
fig.subplots_adjust(left=.025,right=.975,top=.96,bottom=.04)
for ext in ['pdf','png']:fig.savefig(OUT/f'Table_1.{ext}',dpi=300)
plt.close(fig)
# LaTeX uses booktabs and a minipage; no external file dependencies.
tex=['% Requires \\usepackage{booktabs}',r'\begin{table}[htbp]',r'\centering',r'\small',
     r'\caption{Summary of SuSiE and LOCO SuSiE-adj fine-mapping performance.}',r'\label{tab:susie-comparison}',
     r'\begin{tabular}{llrrrrr}',r'\toprule',r'GWAS input & Fine-mapping & Locus fits & 95\% CSs & Coverage (\%) & Power (\%) & CS size (SNPs) \\',r'\midrule']
for i,row in T.iterrows():
    tex.append(' & '.join(str(v).replace('±',r'$\pm$') for v in row.values)+r' \\')
    if i==1:tex.append(r'\addlinespace')
tex += [r'\bottomrule',r'\end{tabular}',r'\vspace{0.5em}',r'\begin{minipage}{0.98\linewidth}',r'\footnotesize',
    table_caption.strip().replace('Table 1. Summary of SuSiE and LOCO SuSiE-adj fine-mapping performance. ','').replace('%',r'\%').replace('±',r'$\pm$'),
    r'\end{minipage}',r'\end{table}']
(OUT/'Table_1.tex').write_text('\n'.join(tex)+'\n')
# Minimal valid Office Open XML manuscript table; standard library only.
def para(txt,bold=False,size=20):
    return '<w:p><w:r><w:rPr>'+('<w:b/>' if bold else '')+f'<w:sz w:val="{size}"/></w:rPr><w:t xml:space="preserve">'+escape(str(txt))+'</w:t></w:r></w:p>'
widths=[2200,1800,1300,1300,2100,2100,2100]
xmlrows=[]
for ri,row in enumerate([labels]+T.values.tolist()):
    cells=[]
    for ci,value in enumerate(row):
        cells.append(f'<w:tc><w:tcPr><w:tcW w:w="{widths[ci]}" w:type="dxa"/>'+('<w:shd w:fill="EDEFF1"/>' if ri==0 else '')+'</w:tcPr>'+para(value,ri==0,18)+'</w:tc>')
    xmlrows.append('<w:tr>'+('<w:trPr><w:tblHeader/></w:trPr>' if ri==0 else '')+''.join(cells)+'</w:tr>')
borders=''.join(f'<w:{name} w:val="single" w:sz="4" w:color="BBBBBB"/>' for name in ['top','left','bottom','right','insideH','insideV'])
table='<w:tbl><w:tblPr><w:tblW w:w="0" w:type="auto"/><w:tblBorders>'+borders+'</w:tblBorders></w:tblPr><w:tblGrid>'+''.join(f'<w:gridCol w:w="{x}"/>' for x in widths)+'</w:tblGrid>'+''.join(xmlrows)+'</w:tbl>'
doc='<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body>'+para('Table 1. SuSiE and LOCO SuSiE-adj fine-mapping performance',True,24)+table+para(table_caption.strip(),False,18)+para('Figure 1 caption',True,24)+para(figure_caption.strip(),False,18)+'<w:sectPr><w:pgSz w:w="15840" w:h="12240" w:orient="landscape"/><w:pgMar w:top="720" w:right="720" w:bottom="720" w:left="720"/></w:sectPr></w:body></w:document>'
with zipfile.ZipFile(OUT/'Table_1_and_captions.docx','w',zipfile.ZIP_DEFLATED) as z:
    z.writestr('[Content_Types].xml','<?xml version="1.0" encoding="UTF-8"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/></Types>')
    z.writestr('_rels/.rels','<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>')
    z.writestr('word/document.xml',doc)
meta={'source':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'figure_width_inches':7.2,'figure_height_inches':3.05,'raster_dpi':600,'n_phenotypes':10,'variation':'sample SD across phenotypes; not standard error or confidence interval','analysis':'LOCO SuSiE-adj versus standard SuSiE-RSS','refitting':False}
(OUT/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
print(T.to_string(index=False))
print('Saved publication files:',OUT)
