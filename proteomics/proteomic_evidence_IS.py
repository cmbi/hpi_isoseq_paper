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
from pyteomics.parser import cleave
from scipy import stats

#%% 
# read in isoform switches
switches=pd.read_csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/switch_consequence.csv')
switches=switches.dropna(subset=['isoformsDifferent'])
switches_diff=switches[switches['isoformsDifferent']]
pep_switches_up=switches_diff[['condition_2','switchConsequence','isoformUpregulated','isoformDownregulated']].groupby(['isoformUpregulated','isoformDownregulated','condition_2']).agg(lambda x: set(x)).reset_index().merge(all_prot.explode('tid'),left_on='isoformUpregulated',right_on='tid')
pep_switches_down=switches_diff[['condition_2','switchConsequence','isoformUpregulated','isoformDownregulated']].groupby(['isoformUpregulated','isoformDownregulated','condition_2']).agg(lambda x: set(x)).reset_index().merge(all_prot.explode('tid'),left_on='isoformDownregulated',right_on='tid')

#%% check evidence IR loss in proteome
bfc=pd.read_table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_secretome/generated_data/metamorpheus/uniprot_results/Task2CalibrationTask/BayesianFoldChangeAnalysis.tsv')
calb=bfc[['Protein Group', 'Gene', 'Organism', 'Control Condition','Treatment Condition', 'Null Hypothesis Width','Protein Log2 Fold-Change', 'Uncertainty in Protein Log2 Fold-Change','Standard Deviation of Peptide Log2 Fold-Changes','Protein Intensity in Control Condition','Protein Intensity in Treatment Condition', 'Number of Peptides','Number of Control Condition Measurements','Number of Treatment Condition Measurements', 'Control Measurements','Treatment Measurements', 'Bayes Factor', 'Posterior Error Probability','False Discovery Rate']]
saur=bfc[['Protein Group.3', 'Gene.3', 'Organism.3', 'Control Condition.3','Treatment Condition.3', 'Null Hypothesis Width.3','Protein Log2 Fold-Change.3', 'Uncertainty in Protein Log2 Fold-Change.3','Standard Deviation of Peptide Log2 Fold-Changes.3','Protein Intensity in Control Condition.3','Protein Intensity in Treatment Condition.3', 'Number of Peptides.3','Number of Control Condition Measurements.3','Number of Treatment Condition Measurements.3', 'Control Measurements.3','Treatment Measurements.3', 'Bayes Factor.3', 'Posterior Error Probability.3','False Discovery Rate.3']]
saur.columns=saur.columns.str.strip('.3')
polyic=bfc[['Protein Group.2', 'Gene.2', 'Organism.2', 'Control Condition.2','Treatment Condition.2', 'Null Hypothesis Width.2','Protein Log2 Fold-Change.2', 'Uncertainty in Protein Log2 Fold-Change.2','Standard Deviation of Peptide Log2 Fold-Changes.2','Protein Intensity in Control Condition.2','Protein Intensity in Treatment Condition.2', 'Number of Peptides.2','Number of Control Condition Measurements.2','Number of Treatment Condition Measurements.2', 'Control Measurements.2','Treatment Measurements.2', 'Bayes Factor.2', 'Posterior Error Probability.2','False Discovery Rate.2']]
polyic.columns=polyic.columns.str.strip('.2')
lps=bfc[['Protein Group.1', 'Gene.1', 'Organism.1', 'Control Condition.1','Treatment Condition.1', 'Null Hypothesis Width.1','Protein Log2 Fold-Change.1', 'Uncertainty in Protein Log2 Fold-Change.1','Standard Deviation of Peptide Log2 Fold-Changes.1','Protein Intensity in Control Condition.1','Protein Intensity in Treatment Condition.1', 'Number of Peptides.1','Number of Control Condition Measurements.1','Number of Treatment Condition Measurements.1', 'Control Measurements.1','Treatment Measurements.1', 'Bayes Factor.1', 'Posterior Error Probability.1','False Discovery Rate.1']]
lps.columns=lps.columns.str.strip('.1')
bfc=pd.concat([calb,saur,polyic,lps])
bfc=bfc[bfc.Organism=='Homo sapiens']
bfc=bfc[bfc['False Discovery Rate']<0.2]
bfc=bfc.drop_duplicates(subset=['Control Measurements','Treatment Measurements','Treatment Condition'])
uniprot_to_ensg=pd.read_table('/mnt/xomics/renees/data/uniprot/uniprot_to_ensg.tsv',names=['Protein Group','Gene stable ID']).dropna()
uniprot_to_ensg['Gene stable ID']=uniprot_to_ensg['Gene stable ID'].apply(lambda x: x.split('.',1)[0])
bfc=bfc.merge(uniprot_to_ensg)
uniprot_to_enst=pd.read_table('/mnt/xomics/renees/data/uniprot/HUMAN_9606_idmapping_enst.dat',usecols=[0,2],names=['Protein Group','Transcript stable ID']).dropna()
bfc=bfc.merge(uniprot_to_enst)
irloss=switches[(switches['switchConsequence']=='Intron retention loss')&(switches['isoformUpregulated'].str.contains('ENST'))]
irloss['Transcript stable ID']=irloss['isoformUpregulated'].apply(lambda x: x.split('.',1)[0])
irloss=irloss.merge(bfc)

#%% 
# check Bayesian DE analysis for support
bfc=pd.read_table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_secretome/generated_data/metamorpheus/LRP_results/Task2CalibrationTask/BayesianFoldChangeAnalysis.tsv')
calb=bfc[['Protein Group', 'Gene', 'Organism', 'Control Condition','Treatment Condition', 'Null Hypothesis Width','Protein Log2 Fold-Change', 'Uncertainty in Protein Log2 Fold-Change','Standard Deviation of Peptide Log2 Fold-Changes','Protein Intensity in Control Condition','Protein Intensity in Treatment Condition', 'Number of Peptides','Number of Control Condition Measurements','Number of Treatment Condition Measurements', 'Control Measurements','Treatment Measurements', 'Bayes Factor', 'Posterior Error Probability','False Discovery Rate']]
calb['tid']=calb['Protein Group'].apply(lambda y: name_conversion_dict[y] if y in name_conversion_dict else y)
saur=bfc[['Protein Group.3', 'Gene.3', 'Organism.3', 'Control Condition.3','Treatment Condition.3', 'Null Hypothesis Width.3','Protein Log2 Fold-Change.3', 'Uncertainty in Protein Log2 Fold-Change.3','Standard Deviation of Peptide Log2 Fold-Changes.3','Protein Intensity in Control Condition.3','Protein Intensity in Treatment Condition.3', 'Number of Peptides.3','Number of Control Condition Measurements.3','Number of Treatment Condition Measurements.3', 'Control Measurements.3','Treatment Measurements.3', 'Bayes Factor.3', 'Posterior Error Probability.3','False Discovery Rate.3']]
saur.columns=saur.columns.str.strip('.3')
saur['tid']=saur['Protein Group'].apply(lambda y: name_conversion_dict[y] if y in name_conversion_dict else y)
polyic=bfc[['Protein Group.2', 'Gene.2', 'Organism.2', 'Control Condition.2','Treatment Condition.2', 'Null Hypothesis Width.2','Protein Log2 Fold-Change.2', 'Uncertainty in Protein Log2 Fold-Change.2','Standard Deviation of Peptide Log2 Fold-Changes.2','Protein Intensity in Control Condition.2','Protein Intensity in Treatment Condition.2', 'Number of Peptides.2','Number of Control Condition Measurements.2','Number of Treatment Condition Measurements.2', 'Control Measurements.2','Treatment Measurements.2', 'Bayes Factor.2', 'Posterior Error Probability.2','False Discovery Rate.2']]
polyic.columns=polyic.columns.str.strip('.2')
polyic['tid']=polyic['Protein Group'].apply(lambda y: name_conversion_dict[y] if y in name_conversion_dict else y)
lps=bfc[['Protein Group.1', 'Gene.1', 'Organism.1', 'Control Condition.1','Treatment Condition.1', 'Null Hypothesis Width.1','Protein Log2 Fold-Change.1', 'Uncertainty in Protein Log2 Fold-Change.1','Standard Deviation of Peptide Log2 Fold-Changes.1','Protein Intensity in Control Condition.1','Protein Intensity in Treatment Condition.1', 'Number of Peptides.1','Number of Control Condition Measurements.1','Number of Treatment Condition Measurements.1', 'Control Measurements.1','Treatment Measurements.1', 'Bayes Factor.1', 'Posterior Error Probability.1','False Discovery Rate.1']]
lps.columns=lps.columns.str.strip('.1')
lps['tid']=lps['Protein Group'].apply(lambda y: name_conversion_dict[y] if y in name_conversion_dict else y)
#bfc=pd.concat([calb,saur,polyic,lps])

calb_switches=switches_diff[switches_diff['condition_2']=='calb'][['switchConsequence','isoformUpregulated','isoformDownregulated']].groupby(['isoformUpregulated','isoformDownregulated']).agg(lambda x: set(x)).reset_index()
calb_switches['isoformUpregulated']=calb_switches['isoformUpregulated'].str.upper()
calb_switches['isoformDownregulated']=calb_switches['isoformDownregulated'].str.upper()
calb_switches.merge(calb[calb['False Discovery Rate']<0.05],left_on='isoformUpregulated',right_on='tid')
calb_switches.merge(calb[calb['False Discovery Rate']<0.05],left_on='isoformDownregulated',right_on='tid')

saur_switches=switches_diff[switches_diff['condition_2']=='saur'][['switchConsequence','isoformUpregulated','isoformDownregulated']].groupby(['isoformUpregulated','isoformDownregulated']).agg(lambda x: set(x)).reset_index()
saur_switches['isoformUpregulated']=saur_switches['isoformUpregulated'].str.upper()
saur_switches['isoformDownregulated']=saur_switches['isoformDownregulated'].str.upper()
saur_switches.merge(saur[saur['False Discovery Rate']<0.05],left_on='isoformUpregulated',right_on='tid')
saur_switches.merge(saur[saur['False Discovery Rate']<0.05],left_on='isoformDownregulated',right_on='tid')

lps_switches=switches_diff[switches_diff['condition_2']=='lps'][['switchConsequence','isoformUpregulated','isoformDownregulated']].groupby(['isoformUpregulated','isoformDownregulated']).agg(lambda x: set(x)).reset_index()
lps_switches['isoformUpregulated']=lps_switches['isoformUpregulated'].str.upper()
lps_switches['isoformDownregulated']=lps_switches['isoformDownregulated'].str.upper()
lps_switches.merge(lps[lps['False Discovery Rate']<0.05],left_on='isoformUpregulated',right_on='tid')
lps_switches.merge(lps[lps['False Discovery Rate']<0.05],left_on='isoformDownregulated',right_on='tid')


polyic_switches=switches_diff[switches_diff['condition_2']=='polyic'][['switchConsequence','isoformUpregulated','isoformDownregulated']].groupby(['isoformUpregulated','isoformDownregulated']).agg(lambda x: set(x)).reset_index()
polyic_switches['isoformUpregulated']=polyic_switches['isoformUpregulated'].str.upper()
polyic_switches['isoformDownregulated']=polyic_switches['isoformDownregulated'].str.upper()
polyic_switches.merge(polyic[polyic['False Discovery Rate']<0.05],left_on='isoformUpregulated',right_on='tid')
polyic_switches.merge(polyic[polyic['False Discovery Rate']<0.05],left_on='isoformDownregulated',right_on='tid')
