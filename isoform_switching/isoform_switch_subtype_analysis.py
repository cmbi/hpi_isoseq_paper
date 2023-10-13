#!/usr/bin/env python3
#%%
import pandas as pd
from gprofiler import GProfiler
import seaborn as sns
sns.set_palette("colorblind")
sns.set_context("paper")
sns.set_style("white")
import matplotlib.pyplot as plt
from gtfparse import read_gtf
import numpy as np
plt.rcParams.update({'font.size': 24})
#%% read files from isoformswitchanalyzeR
switches=pd.read_csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/switch_consequence.csv')
switches=switches.dropna(subset=['isoformsDifferent'])
switches_diff=switches[switches['isoformsDifferent']]
switches_diff.switchConsequence.value_counts().plot.bar(figsize=(8,8))
#%% isoform switches can have multiple "consequences"
switches_diff[['condition_2','isoformUpregulated','isoformDownregulated']].groupby(['condition_2','isoformUpregulated','isoformDownregulated']).agg(len).reset_index().rename({0:'num_cons'},axis=1).sort_values('num_cons',ascending=False).drop_duplicates(subset=['isoformUpregulated','isoformDownregulated'])['num_cons'].value_counts().plot.bar()
plt.xlabel('max # consequences in an IS')
plt.ylabel('# isoform switches')

#%% number of switches per gene
ispg=switches_diff[['gene_id','condition_2','isoformUpregulated','isoformDownregulated']].drop_duplicates().groupby(['condition_2','gene_id']).agg(list)
ispg['num_switches']=ispg.isoformUpregulated.str.len()
ispg.num_switches.value_counts().plot.bar(xlabel='# isoform switches per gene')

# %% abundance analyses - % isoform usage in CIS
isofeat=pd.read_csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/isoform_features.csv')
isofeat['abs_dIF']=isofeat.dIF.abs()
isofeat[(isofeat.switchConsequencesGene==True)&(isofeat['isoform_id'].isin(switches['isoformUpregulated'].unique()))&(isofeat.abs_dIF>0.1)]['abs_dIF'].plot.kde(xlim=(0,1),label='Consequential',figsize=(8,8))
isofeat[(isofeat.switchConsequencesGene==True)&(~isofeat['isoform_id'].isin(switches['isoformUpregulated'].unique()))&(isofeat.abs_dIF>0.1)]['abs_dIF'].plot.kde(xlim=(0,1),label='Non consequential')
plt.xlabel('dIF')
plt.legend()

#%% how often is a isoform switch involving a novel transcript - Figure 4D
def typeswitch(up,down):
    if '.n' in up and '.n' in down:
        return('both novel')
    elif '.n' not in up and '.n' in down:
        return('novel down')
    elif '.n' in up and '.n' not in down:
        return('novel up')
    elif '.n' not in up and '.n' not in down:
        return('both known')
df=switches_diff[['condition_2','isoformUpregulated','isoformDownregulated']].drop_duplicates()
df['type_switch']=df.apply(lambda x: typeswitch(x['isoformUpregulated'],x['isoformDownregulated']),axis=1)
df.type_switch.value_counts().plot.bar(fontsize=22)
plt.ylabel('# isoform switches',fontsize=22)
plt.tight_layout()
# plt.savefig('/mnt/xomics/renees/data/host_pathogen_PID/figures/novel_is.pdf')

#%% what effects of switch involving a novel transcript - Figure 4D
plt.style.use('default')
normed_effects=pd.crosstab(switches_diff.merge(df)['type_switch'],switches_diff.merge(df)['switchConsequence'],normalize='index').T
normed_effects['max']=normed_effects.max(axis=1)
normed_effects[normed_effects['max']>0.06].drop(columns=['max'])[['novel down','both known', 'novel up', 'both novel']].plot.bar(fontsize=22,figsize=(15,15))
plt.legend(fontsize=22,bbox_to_anchor=(1.0, 1.0))
plt.xlabel('Switch Consequence',fontsize=22)
plt.tight_layout()
plt.savefig('/mnt/xomics/renees/data/host_pathogen_PID/figures/novel_cis.pdf')

#%%
domain=pd.read_csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/domain_analysis.csv')[['isoform_id','hmm_name']].drop_duplicates()
domainswitch_gain=switches_diff[switches_diff['switchConsequence']=='Domain gain']
domainswitch_gain=domainswitch_gain.merge(domainswitch_gain.merge(domain,left_on='isoformUpregulated',right_on='isoform_id')[['isoform_id','hmm_name']].groupby('isoform_id').agg(lambda x: set(x)),left_on='isoformUpregulated',right_index=True,how='left',suffixes=('','_up')).merge(domainswitch_gain.merge(domain[['isoform_id','hmm_name']],left_on='isoformDownregulated',right_on='isoform_id')[['isoform_id','hmm_name']].groupby('isoform_id').agg(lambda x: set(x)),left_on='isoformDownregulated',right_index=True,how='left',suffixes=('','_down'))
domainswitch_gain['hmm_name_down']=domainswitch_gain['hmm_name_down'].fillna('').apply(lambda x: set(x))
domainswitch_gain.apply(lambda x: list(x['hmm_name'].difference(x['hmm_name_down'])),axis=1).explode().value_counts().head(20).plot.bar(figsize=(15,10),fontsize=22)
pd.DataFrame(domainswitch_gain.apply(lambda x: list(x['hmm_name'].difference(x['hmm_name_down'])),axis=1).explode().value_counts().head(20)).merge(pd.read_csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/Pfam-A.hmm.csv'),left_index=True,right_on='ID').to_csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/domainsgained.csv',index=False) # look at the names of the gained domains in an excel
domainswitch_gain['gene']=domainswitch_gain['gene_id'].apply(lambda x: x.split('.',1)[0])
dgain_genep=GProfiler(return_dataframe=True).profile(organism='hsapiens',query=list(domainswitch_gain['gene'].unique()))
dgain_genep=dgain_genep[(dgain_genep['term_size']>100)&(dgain_genep['term_size']<500)]
dgain_genep['enrichment']=(dgain_genep['intersection_size']/dgain_genep['query_size'])/(dgain_genep['term_size']/dgain_genep['effective_domain_size'])
dgain_genep[['enrichment','name']].set_index('name').sort_values(by='enrichment',ascending=False).head(20).plot.barh(figsize=(5,8),fontsize=22)

#%% Pathway analyses on isoform switch sub-types
domainswitch_loss=switches_diff[switches_diff['switchConsequence']=='Domain loss']
domainswitch_loss=domainswitch_loss.merge(domainswitch_loss.merge(domain,left_on='isoformDownregulated',right_on='isoform_id')[['isoform_id','hmm_name']].groupby('isoform_id').agg(lambda x: set(x)),left_on='isoformDownregulated',right_index=True,how='left',suffixes=('','_down')).merge(domainswitch_loss.merge(domain[['isoform_id','hmm_name']],left_on='isoformUpregulated',right_on='isoform_id')[['isoform_id','hmm_name']].groupby('isoform_id').agg(lambda x: set(x)),left_on='isoformUpregulated',right_index=True,how='left',suffixes=('','_up'))
domainswitch_loss['hmm_name_up']=domainswitch_loss['hmm_name_up'].fillna('').apply(lambda x: set(x))
domainswitch_loss.apply(lambda x: list(x['hmm_name'].difference(x['hmm_name_up'])),axis=1).explode().value_counts().head(20).plot.bar(figsize=(8,10),fontsize=22)
domainswitch_loss['gene']=domainswitch_loss['gene_id'].apply(lambda x: x.split('.',1)[0])
dloss_genep=GProfiler(return_dataframe=True).profile(organism='hsapiens',query=list(domainswitch_loss['gene'].unique()))
#dloss_genep=dloss_genep[(dloss_genep['term_size']>100)&(dloss_genep['term_size']<500)]
dloss_genep['enrichment']=(dloss_genep['intersection_size']/dloss_genep['query_size'])/(dloss_genep['term_size']/dloss_genep['effective_domain_size'])
dloss_genep[['enrichment','name']].set_index('name').sort_values(by='enrichment',ascending=False).head(20).plot.barh(figsize=(8,8),fontsize=22)

nmdinsensitive=switches_diff[switches_diff.switchConsequence=='NMD insensitive']
nmdinsensitive['gene']=nmdinsensitive['gene_id'].apply(lambda x: x.split('.',1)[0])
nmd_genep=GProfiler(return_dataframe=True).profile(organism='hsapiens',query=list(nmdinsensitive['gene'].unique()))
nmd_genep=nmd_genep[(nmd_genep['term_size']>100)&(nmd_genep['term_size']<500)]
nmd_genep['enrichment']=(nmd_genep['intersection_size']/nmd_genep['query_size'])/(nmd_genep['term_size']/nmd_genep['effective_domain_size'])
nmd_genep[['enrichment','name']].set_index('name').sort_values(by='enrichment',ascending=False).head(20).plot.barh(figsize=(10,10),fontsize=22)

orflen=switches_diff[switches_diff.switchConsequence=='ORF is longer']
orflen['gene']=orflen['gene_id'].apply(lambda x: x.split('.',1)[0])
orflen_genep=GProfiler(return_dataframe=True).profile(organism='hsapiens',query=list(orflen['gene'].unique()))
orflen_genep=orflen_genep[(orflen_genep['term_size']>100)&(orflen_genep['term_size']<500)]
orflen_genep['enrichment']=(orflen_genep['intersection_size']/orflen_genep['query_size'])/(orflen_genep['term_size']/orflen_genep['effective_domain_size'])
orflen_genep[['enrichment','name']].set_index('name').sort_values(by='enrichment',ascending=False).head(20).plot.barh(figsize=(10,10),fontsize=22)

#%%

ir_diff=switches_diff[switches_diff['featureCompared']=='intron_retention']
ir_diff['gene']=ir_diff['gene_id'].apply(lambda x: x.split('.',1)[0])
ir_genep=GProfiler(return_dataframe=True).profile(organism='hsapiens',query=list(ir_diff['gene'].unique()))
ir_genep=ir_genep[(ir_genep['term_size']>100)&(ir_genep['term_size']<500)]
ir_genep['enrichment']=(ir_genep['intersection_size']/ir_genep['query_size'])/(ir_genep['term_size']/ir_genep['effective_domain_size'])
ir_genep[['enrichment','name']].set_index('name').sort_values(by='enrichment',ascending=False).head(20).plot.barh(figsize=(10,10),fontsize=22)
#%% split up loss and gain from above
ir_loss=switches_diff[switches_diff['switchConsequence']=='Intron retention loss']
ir_loss['gene']=ir_loss['gene_id'].apply(lambda x: x.split('.',1)[0])
ir_genep=GProfiler(return_dataframe=True).profile(organism='hsapiens',query=list(ir_loss['gene'].unique()))
ir_genep=ir_genep[(ir_genep['term_size']>100)&(ir_genep['term_size']<500)]
ir_genep['enrichment']=(ir_genep['intersection_size']/ir_genep['query_size'])/(ir_genep['term_size']/ir_genep['effective_domain_size'])
ir_genep[['enrichment','name']].set_index('name').sort_values(by='enrichment',ascending=False).head(20).plot.barh(figsize=(10,10),fontsize=22)
len(ir_loss.gene.unique())
ir_loss.shape[0]


#%%
ir_gain=switches_diff[switches_diff['switchConsequence']=='Intron retention gain']
ir_gain['gene']=ir_gain['gene_id'].apply(lambda x: x.split('.',1)[0])
ir_genep=GProfiler(return_dataframe=True).profile(organism='hsapiens',query=list(ir_gain['gene'].unique()))
ir_genep=ir_genep[(ir_genep['term_size']>100)&(ir_genep['term_size']<500)]
ir_genep['enrichment']=(ir_genep['intersection_size']/ir_genep['query_size'])/(ir_genep['term_size']/ir_genep['effective_domain_size'])
ir_genep[['enrichment','name']].set_index('name').sort_values(by='enrichment',ascending=False).head(20).plot.barh(figsize=(10,10),fontsize=22)

#%%
#what about those IS that didn't have DE on gene level
switches=switches.dropna(subset=['switchConsequence'])
s=pd.concat([switches[(switches['condition_2']=='calb')&(~switches['gene_id'].isin(dg_alb.gid.unique()))][['gene_id']].drop_duplicates(),switches[(switches['condition_2']=='saur')&(~switches['gene_id'].isin(dg_aur.gid.unique()))][['gene_id']].drop_duplicates(),switches[(switches['condition_2']=='lps')&(~switches['gene_id'].isin(dg_lps.gid.unique()))][['gene_id']].drop_duplicates(),switches[(switches['condition_2']=='polyic')&(~switches['gene_id'].isin(dg_polyic.gid.unique()))][['gene_id']].drop_duplicates()]).drop_duplicates()
s[s.gene_id.isin(switches.gene_id.unique())].shape[0]

# %% stat test for above
import numpy as np
from scipy.stats import chi2_contingency
chi2_contingency(pd.DataFrame(switches_diff_de.switchConsequence.value_counts()).rename({'switchConsequence':'DE'},axis=1).merge(pd.DataFrame(switches_diff_node.switchConsequence.value_counts()).rename({'switchConsequence':'no DE'},axis=1),how='outer',left_index=True,right_index=True).fillna(0))

# %%
#how often was a switch more specific to one condition
paired=switches_diff[['isoformUpregulated','isoformDownregulated','condition_2']].drop_duplicates().groupby(['isoformUpregulated','isoformDownregulated']).agg(list)
paired['num_conds']=paired.condition_2.str.len()
paired.num_conds.plot.hist() # num conditions per isoform switch
paired.condition_2.explode().value_counts().plot.bar() # num isoform switches per condition

# %% on a gene level
paired_gene=switches_diff[['gene_id','condition_2']].drop_duplicates().groupby('gene_id').agg(list)
paired_gene['num_conds']=paired_gene.condition_2.str.len()
paired_gene.num_conds.plot.hist() # num conditions per isoform switch
paired_gene.condition_2.explode().value_counts().plot.bar() # num isoform switches per condition

swan=pd.concat([pd.read_csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/diff_exp/isoquant_based/genes_diff_isoform_10dpi_calb.csv'),pd.read_csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/diff_exp/isoquant_based/genes_diff_isoform_10dpi_saur.csv'),pd.read_csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/diff_exp/isoquant_based/genes_diff_isoform_10dpi_lps.csv'),pd.read_csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/diff_exp/isoquant_based/genes_diff_isoform_10dpi_polyic.csv')])
pd.DataFrame(pd.DataFrame(swan.gid.value_counts()).gid.value_counts()).rename({'gid':'All IS'},axis=1).merge(pd.DataFrame(paired.num_conds.value_counts()).rename({'num_conds':'CIS'},axis=1),left_index=True,right_index=True).sort_index().plot.bar(xlabel='Number of conditions with same IS',ylabel='Occurences')


# %% abundance analyses - % mean gene abundance in CIS
isofeat[(isofeat.switchConsequencesGene==True)&(isofeat['isoform_id'].isin(switches['isoformUpregulated'].unique()))&(isofeat.abs_dIF>0.1)]['gene_overall_mean'].plot.kde(xlim=(0,500),label='Consequential',figsize=(8,8))
isofeat[(isofeat.switchConsequencesGene==True)&(~isofeat['isoform_id'].isin(switches['isoformUpregulated'].unique()))&(isofeat.abs_dIF>0.1)]['gene_overall_mean'].plot.kde(xlim=(0,500),label='Non consequential')
plt.xlabel('mean gene abundance across conditions')
plt.legend()
#%% abundance analyses - mean gene abundance in CIS

# %% overall gene expression for consequential 
sd_extra=switches_diff.merge(isofeat[['gene_ref','gene_overall_mean','gene_log2_fold_change']].drop_duplicates()).merge(isofeat[['isoform_id','condition_1','condition_2','iso_log2_fold_change','IF1','IF2']].rename({'isoform_id':'isoformUpregulated'},axis=1).drop_duplicates())



# %% connect to immunity
immune=pd.read_table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/go_immune.tsv',names=['uniprot','gene_name','type'])
df=read_gtf('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoquant/five_conditions/combined/00_m64102e_220102_030451.flnc.transcript_models.gtf')[['gene_id']].drop_duplicates()
df['gene']=df['gene_id'].apply(lambda x: x.split('.',1)[0])
df=df.merge(pd.read_table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/biomart_mapping.tsv')[['Gene stable ID','Gene name']].drop_duplicates(),left_on='gene',right_on='Gene stable ID').dropna()
immune_switches=switches_diff[switches_diff.gene_name.isin(immune.gene_name.unique())]
nonimmune_switches=switches_diff[~switches_diff.gene_name.isin(immune.gene_name.unique())]
#immune_switches[['condition_2','gene_name']].drop_duplicates().to_csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/immune_switches.tsv',sep='\t',index=False)
# for n in immune_switches.gene_name.dropna().unique():
#     path='/mnt/xomics/renees/data/host_pathogen_PID/figures/immune_CIS/'+n
#     os.mkdir(path)
# %% immune v non immune domains
domain=pd.read_csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/domain_analysis.csv')[['isoform_id','hmm_name']].drop_duplicates()
domainswitch_gain=immune_switches[immune_switches['switchConsequence']=='Domain gain']
domainswitch_gain=domainswitch_gain.merge(domainswitch_gain.merge(domain,left_on='isoformUpregulated',right_on='isoform_id')[['isoform_id','hmm_name']].groupby('isoform_id').agg(lambda x: set(x)),left_on='isoformUpregulated',right_index=True,how='left',suffixes=('','_up')).merge(domainswitch_gain.merge(domain[['isoform_id','hmm_name']],left_on='isoformDownregulated',right_on='isoform_id')[['isoform_id','hmm_name']].groupby('isoform_id').agg(lambda x: set(x)),left_on='isoformDownregulated',right_index=True,how='left',suffixes=('','_down'))
domainswitch_gain['hmm_name_down']=domainswitch_gain['hmm_name_down'].fillna('').apply(lambda x: set(x))
domainswitch_gain_ni=nonimmune_switches[nonimmune_switches['switchConsequence']=='Domain gain']
domainswitch_gain_ni=domainswitch_gain_ni.merge(domainswitch_gain_ni.merge(domain,left_on='isoformUpregulated',right_on='isoform_id')[['isoform_id','hmm_name']].groupby('isoform_id').agg(lambda x: set(x)),left_on='isoformUpregulated',right_index=True,how='left',suffixes=('','_up')).merge(domainswitch_gain_ni.merge(domain[['isoform_id','hmm_name']],left_on='isoformDownregulated',right_on='isoform_id')[['isoform_id','hmm_name']].groupby('isoform_id').agg(lambda x: set(x)),left_on='isoformDownregulated',right_index=True,how='left',suffixes=('','_down'))
domainswitch_gain_ni['hmm_name_down']=domainswitch_gain_ni['hmm_name_down'].fillna('').apply(lambda x: set(x))
domains_comp=pd.DataFrame(domainswitch_gain.apply(lambda x: list(x['hmm_name'].difference(x['hmm_name_down'])),axis=1).explode().value_counts(normalize=True)).rename({0:'immune_related'},axis=1).merge(pd.DataFrame(domainswitch_gain_ni.apply(lambda x: list(x['hmm_name'].difference(x['hmm_name_down'])),axis=1).explode().value_counts(normalize=True)).rename({0:'non_immune_related'},axis=1),how='outer',left_index=True,right_index=True).fillna(0)
domains_comp['max']=domains_comp.max(axis=1)
domains_comp.sort_values(by='max',ascending=False).drop(columns=['max']).head(40).plot.bar(figsize=(17,8))
