#!/usr/bin/env python3
#%%
from gprofiler import GProfiler
import pandas as pd
from gtfparse import read_gtf
import seaborn as sns
import matplotlib.pyplot as plt
from venn import venn
from Bio import SeqIO
from collections import defaultdict

#%% 
bfc=pd.read_table('uniprot_results/Task2CalibrationTask/BayesianFoldChangeAnalysis.tsv')
calb=bfc[['Protein Group', 'Gene', 'Organism', 'Control Condition','Treatment Condition', 'Null Hypothesis Width','Protein Log2 Fold-Change', 'Uncertainty in Protein Log2 Fold-Change','Standard Deviation of Peptide Log2 Fold-Changes','Protein Intensity in Control Condition','Protein Intensity in Treatment Condition', 'Number of Peptides','Number of Control Condition Measurements','Number of Treatment Condition Measurements', 'Control Measurements','Treatment Measurements', 'Bayes Factor', 'Posterior Error Probability','False Discovery Rate']]
saur=bfc[['Protein Group.3', 'Gene.3', 'Organism.3', 'Control Condition.3','Treatment Condition.3', 'Null Hypothesis Width.3','Protein Log2 Fold-Change.3', 'Uncertainty in Protein Log2 Fold-Change.3','Standard Deviation of Peptide Log2 Fold-Changes.3','Protein Intensity in Control Condition.3','Protein Intensity in Treatment Condition.3', 'Number of Peptides.3','Number of Control Condition Measurements.3','Number of Treatment Condition Measurements.3', 'Control Measurements.3','Treatment Measurements.3', 'Bayes Factor.3', 'Posterior Error Probability.3','False Discovery Rate.3']]
saur.columns=saur.columns.str.strip('.3')
polyic=bfc[['Protein Group.2', 'Gene.2', 'Organism.2', 'Control Condition.2','Treatment Condition.2', 'Null Hypothesis Width.2','Protein Log2 Fold-Change.2', 'Uncertainty in Protein Log2 Fold-Change.2','Standard Deviation of Peptide Log2 Fold-Changes.2','Protein Intensity in Control Condition.2','Protein Intensity in Treatment Condition.2', 'Number of Peptides.2','Number of Control Condition Measurements.2','Number of Treatment Condition Measurements.2', 'Control Measurements.2','Treatment Measurements.2', 'Bayes Factor.2', 'Posterior Error Probability.2','False Discovery Rate.2']]
polyic.columns=polyic.columns.str.strip('.2')
lps=bfc[['Protein Group.1', 'Gene.1', 'Organism.1', 'Control Condition.1','Treatment Condition.1', 'Null Hypothesis Width.1','Protein Log2 Fold-Change.1', 'Uncertainty in Protein Log2 Fold-Change.1','Standard Deviation of Peptide Log2 Fold-Changes.1','Protein Intensity in Control Condition.1','Protein Intensity in Treatment Condition.1', 'Number of Peptides.1','Number of Control Condition Measurements.1','Number of Treatment Condition Measurements.1', 'Control Measurements.1','Treatment Measurements.1', 'Bayes Factor.1', 'Posterior Error Probability.1','False Discovery Rate.1']]
lps.columns=lps.columns.str.strip('.1')
bfc=pd.concat([calb,saur,polyic,lps])
bfc=bfc[bfc.Organism=='Homo sapiens']
bfc=bfc[bfc['False Discovery Rate']<0.05]
bfc=bfc.drop_duplicates(subset=['Control Measurements','Treatment Measurements','Treatment Condition'])
uniprot_to_ensg=pd.read_table('uniprot/uniprot_to_ensg.tsv',names=['Protein Group','Gene stable ID']).dropna()
uniprot_to_ensg['Gene stable ID']=uniprot_to_ensg['Gene stable ID'].apply(lambda x: x.split('.',1)[0])
bfc=bfc.merge(uniprot_to_ensg).merge(pd.read_table('host_pathogen_PID/pbmc_isoseq/analysis/biomart_mapping.tsv')[['Gene stable ID','Gene name']].drop_duplicates())
bfc['Gene name']=bfc['Gene name'].fillna('None')
bfc=bfc[~bfc['Gene name'].str.contains('KRT')]
#%% method to reduce error from shared peptide quantification
bfc=bfc.drop_duplicates(subset=['Control Measurements','Treatment Measurements','Treatment Condition'])
#%% remove uniprot contaminants
contam=pd.read_csv('host_pathogen_PID/pbmc_secretome/generated_data/databases/contaminantsidentifiers.tsv',delim_whitespace=True,usecols=[0],names=['prot'])
contam['prot']=contam['prot'].str[1:]
bfc=bfc[~bfc['Protein Group'].isin(contam.prot.unique())]

#%% go enrichment
pathways=[]
for stimulus in [bfc[bfc['Treatment Condition']=='LPS'],bfc[bfc['Treatment Condition']=='SAUR'],bfc[bfc['Treatment Condition']=='CALB'],bfc[bfc['Treatment Condition']=='POLYIC']]:
    dg_up_stim,dg_down_stim=stimulus[stimulus['Protein Log2 Fold-Change']>0],stimulus[stimulus['Protein Log2 Fold-Change']<0]
    stim_genep_up=GProfiler(return_dataframe=True).profile(organism='hsapiens',query=list(dg_up_stim['Gene stable ID'].unique()))
    stim_genep_down=GProfiler(return_dataframe=True).profile(organism='hsapiens',query=list(dg_down_stim['Gene stable ID'].unique()))
    stim_genep_up=stim_genep_up[(stim_genep_up['term_size']>100)&(stim_genep_up['term_size']<500)&(stim_genep_up['source'].str.contains('GO'))]
    stim_genep_down=stim_genep_down[(stim_genep_down['term_size']>100)&(stim_genep_down['term_size']<500)&(stim_genep_down['source'].str.contains('GO'))]
    stim_genep_up['enrichment']=(stim_genep_up['intersection_size']/stim_genep_up['query_size'])/(stim_genep_up['term_size']/stim_genep_up['effective_domain_size'])
    stim_genep_down['enrichment']=-(stim_genep_down['intersection_size']/stim_genep_down['query_size'])/(stim_genep_down['term_size']/stim_genep_down['effective_domain_size'])
    pathways.append(pd.concat([stim_genep_up,stim_genep_down]))
lps_genep,saur_genep,calb_genep,polyic_genep=pathways

#%% pathway analysis for protein cluster 4
clust=pd.read_csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_secretome/generated_data/metamorpheus/analysis_files/ids_cluster4.txt',names=['gid'])
pathways=GProfiler(return_dataframe=True).profile(organism='hsapiens',query=list(clust.gid.unique()))
pathways=pathways[(pathways['term_size']>100)&(pathways['term_size']<500)&(pathways['source']=='GO:BP')]
pathways['enrichment']=(pathways['intersection_size']/pathways['query_size'])/(pathways['term_size']/pathways['effective_domain_size'])
pathways[['enrichment','name']].set_index('name').sort_values(by='enrichment').tail(20).plot.barh(figsize=(10,10),fontsize=22,legend=False,alpha=0.8)
plt.xlabel('Enrichment', fontsize=22)
plt.ylabel('Pathway', fontsize=22)
plt.tight_layout()
plt.savefig('/mnt/xomics/renees/data/host_pathogen_PID/figures/clus4_gobp.pdf')

# %% make pathway analysis plots
calb_only_g=set(calb_genep['native'].unique()).difference(set(polyic_genep['native'].unique())).difference(set(saur_genep['native'].unique())).difference(set(lps_genep['native'].unique()))
plot=calb_genep[calb_genep['native'].isin(calb_only_g)][['enrichment','name']].set_index('name').sort_values(by='enrichment',key=abs).tail(20).sort_values(by='enrichment').plot.barh(figsize=(10,10),fontsize=22,legend=False,color='purple',alpha=0.8)
plot.axvline(x=0,color='black')

lps_only_g=set(lps_genep['native'].unique()).difference(set(polyic_genep['native'].unique())).difference(set(saur_genep['native'].unique())).difference(set(calb_genep['native'].unique()))
plot=lps_genep[lps_genep['native'].isin(lps_only_g)][['enrichment','name']].set_index('name').sort_values(by='enrichment',key=abs).tail(20).sort_values(by='enrichment').plot.barh(figsize=(10,10),fontsize=22,legend=False,color='#fde725',alpha=0.8)
plot.axvline(x=0,color='black')
plt.xlabel('Enrichment', fontsize=22)
plt.ylabel('Pathway', fontsize=22)

polyic_only_g=set(polyic_genep['native'].unique()).difference(set(lps_genep['native'].unique())).difference(set(saur_genep['native'].unique())).difference(set(calb_genep['native'].unique()))
plot=polyic_genep[polyic_genep['native'].isin(polyic_only_g)][['enrichment','name']].set_index('name').sort_values(by='enrichment',key=abs).tail(20).sort_values(by='enrichment').plot.barh(figsize=(10,10),fontsize=22,legend=False,color='#35b779',alpha=0.8)
plot.axvline(x=0,color='black')

saur_only_g=set(saur_genep['native'].unique()).difference(set(lps_genep['native'].unique())).difference(set(polyic_genep['native'].unique())).difference(set(calb_genep['native'].unique()))
plot=saur_genep[saur_genep['native'].isin(saur_only_g)][['enrichment','name']].set_index('name').sort_values(by='enrichment',key=abs).tail(20).sort_values(by='enrichment').plot.barh(figsize=(10,10),fontsize=22,legend=False,alpha=0.8)
plot.axvline(x=0,color='black')

albpoly_g=set(calb_genep['native'].unique()).intersection(set(polyic_genep['native'].unique())).difference(set(lps_genep['native'].unique())).difference(set(saur_genep['native'].unique()))
interall_g=set(calb_genep['native'].unique()).intersection(set(saur_genep['native'].unique())).intersection(set(polyic_genep['native'].unique())).intersection(set(lps_genep['native'].unique()))

apg=genep[genep['native'].isin(albpoly_g)].sort_values(by='max').drop(columns=['native','max']).reset_index()
apg=apg[['name','c.alb','s.aur','polyic','lps']]
plot=apg.set_index('name').plot.barh(figsize=(10,10),fontsize=22,legend=False,color=['#440154','#31688e','#35b779','#fde725'],alpha=0.8)
plot.axvline(x=0,color='black')
plt.xlabel('Enrichment', fontsize=22)
plt.ylabel('Pathway', fontsize=22)
