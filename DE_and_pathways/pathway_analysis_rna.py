#!/usr/bin/env python3
#%%
from gprofiler import GProfiler
import pandas as pd
from gtfparse import read_gtf
import seaborn as sns
import matplotlib.pyplot as plt
from venn import venn
#%%
dg_alb=pd.read_csv('analysis/noiseq/diff_gene_saur_50cpm.csv').rename({'Unnamed: 0':'gid'},axis=1)
dg_aur=pd.read_csv('analysis/noiseq/diff_gene_calb_50cpm.csv').rename({'Unnamed: 0':'gid'},axis=1)
dg_lps=pd.read_csv('analysis/noiseq/diff_gene_lps_50cpm.csv').rename({'Unnamed: 0':'gid'},axis=1)
dg_polyic=pd.read_csv('analysis/noiseq/diff_gene_polyic_50cpm.csv').rename({'Unnamed: 0':'gid'},axis=1)
#%%
for d in [dg_alb,dg_aur,dg_lps,dg_polyic]:
    d['gene']=d.gid.apply(lambda x: x.split('.',1)[0])

#%% split up and downregulated subsets and combine afterwards
pathways=[]
for stimulus in [dg_lps,dg_aur,dg_alb,dg_polyic]:
    dg_up_stim,dg_down_stim=stimulus[stimulus.M>0],stimulus[stimulus.M<0]
    stim_genep_up=GProfiler(return_dataframe=True).profile(organism='hsapiens',query=list(dg_up_stim['gene'].unique()))
    stim_genep_down=GProfiler(return_dataframe=True).profile(organism='hsapiens',query=list(dg_down_stim['gene'].unique()))
    stim_genep_up=stim_genep_up[(stim_genep_up['term_size']>100)&(stim_genep_up['term_size']<500)]
    stim_genep_down=stim_genep_down[(stim_genep_down['term_size']>100)&(stim_genep_down['term_size']<500)]
    stim_genep_up['enrichment']=(stim_genep_up['intersection_size']/stim_genep_up['query_size'])/(stim_genep_up['term_size']/stim_genep_up['effective_domain_size'])
    stim_genep_down['enrichment']=-(stim_genep_down['intersection_size']/stim_genep_down['query_size'])/(stim_genep_down['term_size']/stim_genep_down['effective_domain_size'])
    pathways.append(pd.concat([stim_genep_up,stim_genep_down]))
lps_genep,saur_genep,calb_genep,polyic_genep=pathways

#%% venn diagram
venn({'c.alb':set(calb_genep['native'].unique()),'s.aur':set(saur_genep['native'].unique()),'polyic':set(polyic_genep['native'].unique()),'lps':set(lps_genep['native'].unique())})#,cmap=[])
#%%
#polyic specific figure
polyic_only_g=set(polyic_genep['native'].unique()).difference(set(calb_genep['native'].unique())).difference(set(saur_genep['native'].unique())).difference(set(lps_genep['native'].unique()))
pset=polyic_genep[polyic_genep['native'].isin(polyic_only_g)][['enrichment','name']]
main_fig_set=['release of sequestered calcium ion into cytosol','T cell receptor signaling pathway','regulation of sequestering of calcium ion','viral gene expression']
pset[pset.name.isin(main_fig_set)].set_index('name').sort_values(by='enrichment',key=abs).plot.barh(figsize=(8,2),fontsize=22,legend=False,xlim=(0,8),color='#35b779',alpha=0.8)
plt.xlabel('Enrichment', fontsize=22)
plt.ylabel('Pathway', fontsize=22)
plt.tight_layout()
plt.savefig('polyic_rna_pathways.pdf')
