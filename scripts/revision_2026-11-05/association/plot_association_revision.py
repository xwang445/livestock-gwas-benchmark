#!/usr/bin/env python3
"""Draw review figures and editable Supplementary Table S1 from canonical TSVs."""
import argparse, copy, os, zipfile
from pathlib import Path
from xml.sax.saxutils import escape
os.environ.setdefault('MPLCONFIGDIR','/tmp/gwas_revision_mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ORDER=['BOLT-LMM','REGENIE_LOCO','REGENIE_WO_LOCO','fastGWA','FarmCPU_NYC','FarmCPU_Kansas','BLINK_NYC','BLINK_Kansas','SLEMM','GCTA_MLMA']
LABEL=['BOLT-LMM','REGENIE (LOCO)','REGENIE (no LOCO)','fastGWA','FarmCPU (NYC)','FarmCPU (Kansas)','BLINK (NYC)','BLINK (Kansas)','SLEMM','GCTA MLMA']
SHORT=['BOLT','REG-L','REG-W','fast','FC-N','FC-K','BL-N','BL-K','SLEMM','GCTA']
COLORS=dict(zip(ORDER,['#0072B2','#D55E00','#E69F00','#009E73','#CC79A7','#7B3294','#56B4E9','#999933','#222222','#A6761D']))

def docx_table(rows,path,caption,notes):
    """Minimal OOXML; table cells and caption remain native editable Word text."""
    def p(s,bold=False):return '<w:p><w:r>'+('<w:rPr><w:b/></w:rPr>' if bold else '')+'<w:t xml:space="preserve">'+escape(str(s))+'</w:t></w:r></w:p>'
    body=p('Supplementary Table S1. Regional distribution of significant loci',True)+p(caption)
    body+='<w:tbl><w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="0" w:type="auto"/></w:tblPr>'
    for row in rows:
        body+='<w:tr>'+''.join('<w:tc><w:tcPr><w:tcW w:w="2000" w:type="dxa"/></w:tcPr>'+p(c)+'</w:tc>' for c in row)+'</w:tr>'
    body+='</w:tbl>'+''.join(p(x) for x in notes)
    body+='<w:sectPr><w:pgSz w:w="16838" w:h="11906" w:orient="landscape"/><w:pgMar w:top="720" w:bottom="720" w:left="720" w:right="720"/></w:sectPr>'
    document='<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body>'+body+'</w:body></w:document>'
    types='<?xml version="1.0"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/></Types>'
    rels='<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>'
    with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml',types);z.writestr('_rels/.rels',rels);z.writestr('word/document.xml',document)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--tables',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--write-table-docx',action='store_true');a=ap.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    def read(x):return pd.read_csv(a.tables/x,sep='\t',keep_default_na=False,na_values=['NA'])
    by=read('association_by_phenotype.tsv');by=by[by.configuration!='SLEMM_previous']
    summary=read('association_summary.tsv').set_index('configuration').reindex(ORDER)
    incomplete=[m for m in ORDER if summary.loc[m,'complete_ten_phenotypes']!=True]
    available=set(zip(by.configuration,by.phenotype))
    gcta_missing=[str(ph) for ph in range(1,11) if ('GCTA_MLMA',ph) not in available]
    gcta_coverage_note=('GCTA phenotypes '+', '.join(gcta_missing)+' are pending.'
                        if gcta_missing else 'All ten configurations have complete data for ten phenotypes.')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
    def save(fig,name):
        for ext in ['pdf','svg','png']:fig.savefig(a.output/(name+'.'+ext),dpi=300,facecolor='white',bbox_inches='tight')
        plt.close(fig)
    def note(fig):
        pending=[]
        for m in incomplete:
            missing=[str(ph) for ph in range(1,11) if (m,ph) not in available]
            phlabel='2–10' if missing==[str(ph) for ph in range(2,11)] else ', '.join(missing)
            pending.append(LABEL[ORDER.index(m)]+' phenotypes '+phlabel+' remain pending')
        fig.text(.5,.008,'Review version: '+('; '.join(pending)+'.' if pending else 'all configurations complete.'),ha='center',fontsize=8,color='#666666')
    c=by.pivot(index='phenotype',columns='configuration',values='loci').reindex(index=range(1,11),columns=ORDER)
    s=by.pivot(index='phenotype',columns='configuration',values='significant_snps').reindex(index=range(1,11),columns=ORDER)
    fig,ax=plt.subplots(figsize=(12.7,8));fig.subplots_adjust(left=.09,right=.91,bottom=.23,top=.94)
    cmap=copy.copy(plt.get_cmap('plasma'));cmap.set_bad('#dddddd');norm=LogNorm(vmin=1,vmax=np.nanmax(c.values))
    im=ax.imshow(np.ma.masked_invalid(np.maximum(c.values,1)),aspect='auto',cmap=cmap,norm=norm)
    for i in range(10):
        for j in range(10):
            v=c.iloc[i,j];txt='pending' if pd.isna(v) else '%d\n(%d)'%(v,s.iloc[i,j])
            ax.text(j,i,txt,ha='center',va='center',fontsize=9,color='#444444' if pd.isna(v) or v>700 else 'white')
    ax.set_xticks(range(10));ax.set_xticklabels(LABEL,rotation=42,ha='right');ax.set_yticks(range(10));ax.set_yticklabels(range(1,11));ax.set_ylabel('Phenotype')
    ax.set_title('Significant loci and SNPs',loc='left',pad=12)
    ax.set_xticks(np.arange(-.5,10,1),minor=True);ax.set_yticks(np.arange(-.5,10,1),minor=True);ax.grid(which='minor',color='white',lw=.6);ax.tick_params(which='minor',length=0)
    fig.colorbar(im,ax=ax,fraction=.03,pad=.03,label='500 kb loci');note(fig);save(fig,'Figure_3_revision_review')
    ov=read('pairwise_locus_overlap.tsv');fig,axes=plt.subplots(5,2,figsize=(12.5,17));fig.subplots_adjust(left=.09,right=.96,top=.94,bottom=.13,wspace=.18,hspace=.32)
    cmap=copy.copy(plt.get_cmap('Blues'));cmap.set_bad('#dddddd')
    for ph,ax in enumerate(axes.flat,1):
        z=ov[ov.phenotype==ph].pivot(index='configuration_a',columns='configuration_b',values='jaccard').reindex(index=ORDER,columns=ORDER)
        im=ax.imshow(np.ma.masked_invalid(z.values),cmap=cmap,vmin=0,vmax=1,aspect='auto')
        for i in range(10):
            for j in range(10):
                v=z.iloc[i,j];pending=(ORDER[i],ph) not in available or (ORDER[j],ph) not in available
                ax.text(j,i,'—' if pending else ('NA' if pd.isna(v) else '%.2f'%v),ha='center',va='center',fontsize=6.7,color='white' if v>.6 else '#333333')
        ax.set_title('%s  Phenotype %d'%('abcdefghij'[ph-1],ph),loc='left',fontsize=10)
        ax.set_xticks(range(10));ax.set_xticklabels(SHORT if ph>=9 else ['']*10,rotation=45,ha='right',fontsize=8)
        ax.set_yticks(range(10));ax.set_yticklabels(SHORT if ph%2 else ['']*10,fontsize=8);ax.tick_params(length=0)
    fig.suptitle('Pairwise overlap of significant loci',fontsize=13,y=.98)
    fig.text(.5,.068,'REG-L/W: REGENIE with/without LOCO; FC: FarmCPU; N/K: NYC/Kansas. NA: either set empty; —: analysis pending.',ha='center',fontsize=8)
    cax=fig.add_axes([.39,.036,.24,.009]);fig.colorbar(im,cax=cax,orientation='horizontal',label='Jaccard index');note(fig);save(fig,'Supplementary_phenotype_overlap_review')
    # A page-sized review layout. The historical ten-panel layout is retained
    # separately; this mean matrix does not assert manuscript comment mapping.
    ovs=read('pairwise_locus_overlap_summary.tsv')
    z=ovs.pivot(index='configuration_a',columns='configuration_b',values='mean_jaccard').reindex(index=ORDER,columns=ORDER)
    nn=ovs.pivot(index='configuration_a',columns='configuration_b',values='n_defined_phenotypes').reindex(index=ORDER,columns=ORDER)
    complete_pair=ovs.pivot(index='configuration_a',columns='configuration_b',values='complete_ten_phenotypes').reindex(index=ORDER,columns=ORDER)
    fig,ax=plt.subplots(figsize=(7.2,6.9));fig.subplots_adjust(left=.16,right=.92,top=.91,bottom=.25)
    im=ax.imshow(np.ma.masked_invalid(z.values),cmap=cmap,vmin=0,vmax=1,aspect='auto')
    for i in range(10):
        for j in range(10):
            v=z.iloc[i,j];nvalid=nn.iloc[i,j]
            txt=('pending\nn=%d'%nvalid) if not complete_pair.iloc[i,j] else ('NA\nn=0' if pd.isna(v) else '%.3f\nn=%d'%(v,nvalid))
            ax.text(j,i,txt,ha='center',va='center',fontsize=8,color='white' if v>.6 else '#333333',linespacing=1.15)
    ax.set_xticks(range(10));ax.set_xticklabels(SHORT,rotation=45,ha='right',fontsize=8)
    ax.set_yticks(range(10));ax.set_yticklabels(SHORT,fontsize=8);ax.tick_params(length=0)
    ax.set_xticks(np.arange(-.5,10,1),minor=True);ax.set_yticks(np.arange(-.5,10,1),minor=True);ax.grid(which='minor',color='white',lw=.6);ax.tick_params(which='minor',length=0)
    ax.set_title('Mean overlap of significant loci',loc='left',fontsize=12,pad=12)
    cax=fig.add_axes([.26,.14,.46,.012]);fig.colorbar(im,cax=cax,orientation='horizontal',label='Mean phenotype-specific Jaccard index')
    fig.text(.53,.055,'n = phenotypes with nonempty locus sets for both configurations.\nREG-L/W: REGENIE with/without LOCO; FC: FarmCPU; N/K: NYC/Kansas.',ha='center',va='top',fontsize=7,linespacing=1.4)
    fig.text(.53,.008,'Review layout; GCTA aggregates withheld until ten phenotypes are complete.' if incomplete else 'Review layout; all ten configurations complete.',ha='center',fontsize=7,color='#666666')
    save(fig,'Figure_5_revision_review')
    spatial=read('spatial_by_phenotype.tsv');complete=[m for m in ORDER if summary.loc[m,'complete_ten_phenotypes']==True]
    fig,axes=plt.subplots(2,2,figsize=(14,10));fig.subplots_adjust(left=.075,right=.98,top=.93,bottom=.17,hspace=.6,wspace=.22)
    for label,ax in zip('abcd',axes.flat):ax.text(-.12,1.08,label,transform=ax.transAxes,fontsize=15,fontweight='bold')
    x=np.arange(10)
    for ax,col,factor,title,ylabel in [(axes[0,0],'central_null_spillover_rate',100,'Central-null windows detected','Percent of 1,480 central-null windows'),(axes[0,1],'lambda_gc_central_null_historical',1,'Central-null genomic inflation',r'$\lambda_{GC}$')]:
        vals=summary[col+'_mean'].values*factor;err=summary[col+'_sd'].values*factor
        ax.bar(x,vals,yerr=err,color=[COLORS[m] for m in ORDER],capsize=2)
        ax.set_xticks(x);ax.set_xticklabels(LABEL,rotation=45,ha='right',fontsize=8);ax.set_ylabel(ylabel);ax.set_title(title,fontsize=11)
        ax.set_xlim(-.6,len(ORDER)-.4)
        for m in incomplete:ax.text(ORDER.index(m),0,'pending',rotation=90,ha='center',va='bottom',fontsize=8,color='#666666')
        if factor==1:ax.axhline(1,color='#777777',ls='--',lw=.8)
    for ax,ms,title in [(axes[1,0],complete,'All complete configurations'),(axes[1,1],[m for m in complete if m not in ['BOLT-LMM','REGENIE_LOCO']],'Configurations excluding BOLT-LMM and REGENIE LOCO')]:
        for m in ms:
            z=spatial[spatial.configuration==m].groupby('distance_index').observed_expected_equal_nonempty_bins.mean()
            ax.plot((z.index-.5)*.05,z.values,color=COLORS[m],lw=1.5,label=LABEL[ORDER.index(m)])
        ax.set_title(title,fontsize=10);ax.set_xlabel('Distance from nearest effect boundary (Mb)');ax.set_ylabel('Observed / expected significant SNPs');ax.set_xlim(0,2);ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.grid(alpha=.2)
    handles,labels=axes[1,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.53,.035),ncol=5,frameon=False,fontsize=8)
    note(fig);save(fig,'Figure_6_revision_review')
    rows=[['Configuration','Loci, mean ± SD','Effect loci (%)','Null loci (%)','Central-null loci (%)','Defined percentage phenotypes','Status']]
    for m,label in zip(ORDER,LABEL):
        r=summary.loc[m]
        if not r.complete_ten_phenotypes:rows.append([label,'pending','pending','pending','pending','pending','%d/10 available'%r.n_phenotypes]);continue
        fmt=lambda k,f=1:'%.2f ± %.2f'%(r[k+'_mean']*f,r[k+'_sd']*f)
        rows.append([label,fmt('loci'),fmt('effect_locus_fraction',100),fmt('null_locus_fraction',100),fmt('central_null_locus_fraction',100),str(int(r.effect_locus_fraction_n_defined)),'10/10 complete'])
    caption=('A locus is a non-overlapping, chromosome-specific 500 kb genomic window containing at least one SNP with P ≤ 5 × 10⁻⁸. Windows are anchored at zero, [start,end), with terminal fragments retained. Counts are mean ± sample SD across all ten phenotypes. Percentages are the equally weighted mean ± sample SD of phenotype-specific percentages: detected loci of the stated region divided by all detected loci. Phenotypes with zero detected loci have undefined percentages and are excluded from percentage summaries; their number is reported. Central-null is a subset of null, not a third disjoint category. These percentages are not the central-null window detection rate in Figure 6a, which has a fixed denominator of 1,480 windows per phenotype. Pooled percentages are separate columns in the source TSV. GCTA aggregates remain pending until all ten phenotypes are complete.')
    notes=['Review status: source statistics are verified; manuscript and response-letter quantitative reconciliation is recorded separately. Final ten-phenotype GCTA integration remains pending.',
        'Source: association_by_phenotype.tsv and association_summary.tsv. No effect-allele orientation or generating causal-effect size is inferred from the reported association coefficient.']
    if a.write_table_docx:docx_table(rows,a.output/'Supplementary_Table_S1_review.docx',caption,notes)
    (a.output/'captions_and_definitions.txt').write_text('Figure 3. Significant loci and SNPs\nCells show detected 500 kb loci and significant SNPs in parentheses; P ≤ 5 × 10⁻⁸. Zero-anchored half-open windows retain chromosome-end fragments. GCTA phenotypes 2–10 are pending.\n\nFigure 5. Pairwise overlap of significant loci\nJaccard indices use the inclusive-threshold locus sets of Figure 3. As in the historical plotting convention, comparisons are undefined (NA) if either set is empty; pending fits are also NA. A set-theoretic alternative with zero for one empty set is supplied separately in the canonical table.\n\nFigure 6. Regional spillover and genomic inflation\n(a) Detected central-null 500 kb windows divided by 1,480 eligible windows, mean ± sample SD across ten phenotypes. (b) Historical central-null genomic inflation: median finite one-df association chi-square divided by 0.4549364, then mean ± SD across phenotypes; the P-inverse-survival estimate is supplied separately. (c,d) Significant/(tested × 5 × 10⁻⁸) is computed within nonempty 50 kb null windows, averaged across windows at mirrored distance indices, then equally across phenotypes. Panel d excludes BOLT-LMM and REGENIE LOCO. Complete configurations only are aggregated; GCTA is pending. Closed versus half-open endpoint conventions have no effect because no tested marker lies at a relevant boundary.\n\nTable S1\n'+caption+'\n')
    power_file=a.tables/'power_and_null_snp_rate_summary.tsv'
    if power_file.exists():
        power=read('power_and_null_snp_rate_summary.tsv')
        fig,axes=plt.subplots(1,2,figsize=(12,5));fig.subplots_adjust(left=.07,right=.98,top=.9,bottom=.27,wspace=.27)
        for ax,region,title,ylabel in [(axes[0],'causal_positions','Causal-SNP detection','Fraction of 10,000 causal SNPs detected'),(axes[1],'central_null','Central-null SNP associations','Fraction of tested central-null SNPs significant')]:
            for m in ORDER+['SLEMM_previous']:
                z=power[(power.configuration==m)&(power.region==region)&(power.n_phenotypes==10)].sort_values('alpha')
                if z.empty:continue
                ax.plot(z.alpha,z.mean_rate,color=COLORS.get(m,'#222222'),ls='--' if m=='SLEMM_previous' else '-',lw=1.4,label='SLEMM previous' if m=='SLEMM_previous' else LABEL[ORDER.index(m)])
            ax.set_xscale('log');ax.set_title(title);ax.set_xlabel('Association P threshold');ax.set_ylabel(ylabel);ax.axvline(5e-8,color='#777777',ls=':',lw=.8);ax.grid(alpha=.2)
        axes[0].text(-.12,1.08,'a',transform=axes[0].transAxes,fontweight='bold',fontsize=15);axes[1].text(-.12,1.08,'b',transform=axes[1].transAxes,fontweight='bold',fontsize=15)
        handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.53,.04),ncol=5,frameon=False,fontsize=8);note(fig);save(fig,'Supplementary_power_revision_review')
    if not incomplete:
        cp=a.output/'captions_and_definitions.txt';txt=cp.read_text().replace('GCTA phenotypes 2–10 are pending.','All ten configurations have complete data for ten phenotypes.').replace('Complete configurations only are aggregated; GCTA is pending.','All ten configurations are aggregated over ten phenotypes.').replace('GCTA aggregates remain pending until all ten phenotypes are complete.','GCTA aggregates include all ten phenotypes.');cp.write_text(txt)
    cp=a.output/'captions_and_definitions.txt';txt=cp.read_text().replace('pending fits are also NA.','in the phenotype-specific supplement, pending fits are marked by a dash.').replace('Closed versus half-open endpoint conventions have no effect because no tested marker lies at a relevant boundary.','No tested marker lies at a 500 kb or central-null boundary. Spatial 50 kb bins preserve the historical closed endpoints, duplicating three shared-boundary markers per phenotype.');cp.write_text(txt)
    cp=a.output/'captions_and_definitions.txt';txt=cp.read_text().replace('Figure 5. Pairwise overlap of significant loci\n','Figure 5. Mean overlap of significant loci\nThe compact matrix shows the equal-phenotype mean Jaccard index, conditional on both locus sets being nonempty. Each cell states the number of defined phenotypes (n), excluding phenotypes with an empty locus set for either configuration. Pending cells also show n, not the number of completed fits. Averages are withheld until all ten phenotype outputs are available for both configurations. All ten phenotype-specific matrices are retained as Supplementary_phenotype_overlap_review.\n').replace('GCTA phenotypes 2–10 are pending.',gcta_coverage_note);cp.write_text(txt)

if __name__=='__main__':main()
