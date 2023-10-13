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
import numpy as np
from pyteomics.parser import cleave
#%%
#functions
def give_ids_from_fasta(ff):
    pathogens=[]
    with open(ff) as patho:
        for precord in SeqIO.parse(patho,"fasta"):
            pathogens.append(precord.id)
    return(pathogens)

def digest(protein_sequence):
    '''digests the protein sequences with trypsin
    '''
    seq_cut = cleave(protein_sequence, rule='[KR]', min_length=6, missed_cleavages=2)
    plist=set() #to prevent duplicates
    for peptide in seq_cut:
        if peptide in plist or len(peptide) < 6 or len(peptide) > 40:
            continue
        plist.add(peptide)
    return(plist)

def refine_source(sources):
    if len(set(sources))==1:
        return(sources[0])
    return('Multi-mapping')

def novel_only(proteins):
    listp=proteins.split('|') if '|' in proteins else [proteins]
    for p in listp:
        if 'TRANSCRIPT' not in p:
            return(False)
    return(True)

def find_secr(tid):
    for t in tid:
        if t not in secreted_all:
            return(False)
    return(True)
def find_secr_novel(tid):
    for t in tid:
        if t not in secreted_novel:
            return(False)
    return(True)

def rename_ids(index):
    i=index.split('|') if '|' in index else [index]
    new=[]
    for y in i:
        new.append(name_conversion_dict[y] if y in name_conversion_dict else y)
    return(new)

#%%
#important definitions
filename_key={'Intensity_B09905_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-1-calib':'rpmi1-1','Intensity_B09907_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-2-calib':'lps1','Intensity_B09909_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-3-calib':'saur1','Intensity_B09911_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-4-calib':'polyic1','Intensity_B09913_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-5-calib':'calb1','Intensity_B09915_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-6-calib':'rpmi1-2','Intensity_B09917_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-7-calib':'rpmi2-1','Intensity_B09919_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-8-calib':'lps2','Intensity_B09921_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-9-calib':'saur2','Intensity_B09923_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-10-calib':'polyic2','Intensity_B09925_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-11-calib':'calb2','Intensity_B09925_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-12-calib':'rpmi2-2','Intensity_B10105_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-13-calib':'rpmi3-1','Intensity_B10107_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-14-calib':'lps3','Intensity_B10109_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-15-calib':'saur3','Intensity_B10111_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-16-calib':'polyic3','Intensity_B10113_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-17-calib':'calb3','Intensity_B10115_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-18-calib':'rpmi3-2','Intensity_B10117_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-19-calib':'rpmi4-1','Intensity_B10119_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-20-calib':'lps4','Intensity_B10121_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-21-calib':'saur4','Intensity_B10123_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-22-calib':'polyic4','Intensity_B10125_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-23-calib':'calb4','Intensity_B10127_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-24-calib':'rpmi4-2','Intensity_B10129_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-25-calib':'lps5','Intensity_B10131_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-26-calib':'polyic5','Intensity_B10133_Bp_WNE17_trap9_PRC-5493_LennartMartens_shotgun_5095_1-calib':'rpmi5-1','Intensity_B10135_Bp_WNE17_trap9_PRC-5493_LennartMartens_shotgun_5095_2-calib':'calb5','Intensity_B10139_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun_5095_3-calib':'saur5','Intensity_B10141_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun_5095_4-calib':'rpmi5-2'}
filename_key_short={'B09905_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-1-calib':'rpmi1-1','B09907_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-2-calib':'lps1','B09909_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-3-calib':'saur1','B09911_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-4-calib':'polyic1','B09913_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-5-calib':'calb1','B09915_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-6-calib':'rpmi1-2','B09917_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-7-calib':'rpmi2-1','B09919_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-8-calib':'lps2','B09921_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-9-calib':'saur2','B09923_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-10-calib':'polyic2','B09925_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-11-calib':'calb2','B09925_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-12-calib':'rpmi2-2','B10105_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-13-calib':'rpmi3-1','B10107_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-14-calib':'lps3','B10109_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-15-calib':'saur3','B10111_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-16-calib':'polyic3','B10113_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-17-calib':'calb3','B10115_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-18-calib':'rpmi3-2','B10117_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-19-calib':'rpmi4-1','B10119_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-20-calib':'lps4','B10121_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-21-calib':'saur4','B10123_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-22-calib':'polyic4','B10125_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-23-calib':'calb4','B10127_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-24-calib':'rpmi4-2','B10129_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-25-calib':'lps5','B10131_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-26-calib':'polyic5','B10133_Bp_WNE17_trap9_PRC-5493_LennartMartens_shotgun_5095_1-calib':'rpmi5-1','B10135_Bp_WNE17_trap9_PRC-5493_LennartMartens_shotgun_5095_2-calib':'calb5','B10139_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun_5095_3-calib':'saur5','B10141_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun_5095_4-calib':'rpmi5-2'}
prot=pd.read_table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_secretome/generated_data/metamorpheus/LRP_results/Task2CalibrationTask/QuantifiedProteins.tsv').rename(filename_key,axis=1)
prot=prot[prot.Organism.isna()]
subcloc=pd.read_table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_secretome/generated_data/subcellular_location.tsv')
biomart=pd.read_table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/mart_export.txt')

name_conversion_dict = {}
with open('/mnt/xomics/renees/data/ref_grch38/gencode/gencode.v39.pc_translations.fa') as fp:
  for record in SeqIO.parse(fp,"fasta"):
    ident = record.id.split('|')
    name_conversion_dict[ident[5]]=ident[1]


#%% pie charts ORF vs peptide
# first read in the data

data = defaultdict(list)
with open('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_secretome/generated_data/databases/LRP/hybrid_db_plus_pathogens.fasta') as fp:
  for record in SeqIO.parse(fp,"fasta"):
    ident = record.id
    sequence = str(record.seq)
    data['id'].append(ident)
    data['sequence'].append(sequence)

db = pd.DataFrame.from_dict(data)
#%% add sources
pathogens=give_ids_from_fasta('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_secretome/generated_data/databases/pathogens_calb_saur.fasta')
db['source']=db['id'].apply(lambda x: 'Pathogen' if x in pathogens else x.split('|',1)[0])
db['source']=db['source'].str.replace('gc','Gencode')
db['source']=db['source'].str.replace('pb','Novel')

#make ORF pie chart
db.source.value_counts().plot.pie(colors=['#0173b2','#d55e00','#cc78bc','#ca9161'])

#%% recalculate FDR and check for pathogen proteins
psm=pd.read_table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_secretome/generated_data/metamorpheus/LRP_results/Task5SearchTask/AllPSMs.psmtsv')#,usecols=['File Name','Scan Number', 'Scan Retention Time','Base Sequence','Ambiguity Level','PSM Count (unambiguous, <0.01 q-value)','Protein Accession','Protein Name', 'Gene Name', 'Organism Name', 'Mods','Contaminant','Decoy','Cumulative Target', 'Cumulative Decoy','Cumulative Target Notch', 'Cumulative Decoy Notch','QValue Notch', 'PEP', 'PEP_QValue'])
psm=psm[psm.Contaminant=='N']
psm['File Name']=psm['File Name'].replace({'B09905_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-1-calib':'rpmi1-1','B09907_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-2-calib':'lps1','B09909_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-3-calib':'saur1','B09911_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-4-calib':'polyic1','B09913_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-5-calib':'calb1','B09915_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-6-calib':'rpmi1-2','B09917_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-7-calib':'rpmi2-1','B09919_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-8-calib':'lps2','B09921_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-9-calib':'saur2','B09923_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-10-calib':'polyic2','B09925_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-11-calib':'calb2','B09925_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun-12-calib':'rpmi2-2','B10105_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-13-calib':'rpmi3-1','B10107_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-14-calib':'lps3','B10109_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-15-calib':'saur3','B10111_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-16-calib':'polyic3','B10113_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-17-calib':'calb3','B10115_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-18-calib':'rpmi3-2','B10117_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-19-calib':'rpmi4-1','B10119_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-20-calib':'lps4','B10121_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-21-calib':'saur4','B10123_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-22-calib':'polyic4','B10125_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-23-calib':'calb4','B10127_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-24-calib':'rpmi4-2','B10129_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-25-calib':'lps5','B10131_Bp_WNE19_trap10_PRC-5493_LennartMartens_shotgun-26-calib':'polyic5','B10133_Bp_WNE17_trap9_PRC-5493_LennartMartens_shotgun_5095_1-calib':'rpmi5-1','B10135_Bp_WNE17_trap9_PRC-5493_LennartMartens_shotgun_5095_2-calib':'calb5','B10139_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun_5095_3-calib':'saur5','B10141_Bp_WNE17_trap10_PRC-5493_LennartMartens_shotgun_5095_4-calib':'rpmi5-2'})
psm_pat=psm[~psm['Organism Name'].isna()]
psm_pat=psm_pat[~psm_pat['Organism Name'].str.contains('\|')]
calb_psms=psm_pat[psm_pat['Organism Name'].str.contains('Candida')]
calb_psms['decoy']=calb_psms.Decoy.apply(lambda x: True if x=='Y' else False)
df=calculate_qvalues(calb_psms,decoy_col='decoy',score_col='Score')
df=df[df['q_value']<0.01]#[['File Name','PSM Count (unambiguous, <0.01 q-value)']].groupby('File Name').sum().reset_index()
df['condition']=df['File Name'].str.replace('-', '').str.replace('\d+','')
df.groupby('condition').size()

saur_psms=psm_pat[psm_pat['Organism Name'].str.contains('aureus')]
saur_psms['decoy']=saur_psms.Decoy.apply(lambda x: True if x=='Y' else False)
df=calculate_qvalues(saur_psms,decoy_col='decoy',score_col='Score')
df=df[df['q_value']<0.01]
df['condition']=df['File Name'].str.replace('-', '').str.replace('\d+','')
df.groupby('condition').size()

# %% QC: check for secreted proteins
prot['Gene name']=prot['Gene Name'].str.strip('primary:')
peps=pd.read_table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_secretome/generated_data/metamorpheus/LRP_results/Task5SearchTask/AllQuantifiedPeptides.tsv')#.rename(filename_key,axis=1).iloc[:,:35]
peps['Gene name']=peps['Gene Names'].str.strip('primary:')
# how many ID'ed peps are secreted
peps.merge(subcloc,on='Gene name')['Extracellular location'].value_counts() #3730 of the 16798 peptides are predicted to be secreted (protein atlas)
peps.merge(subcloc,on='Gene name').shape[0]
# how many proteins are secreted
prot.merge(subcloc,on='Gene name')['Extracellular location'].value_counts()
prot.merge(subcloc,on='Gene name').shape[0] #404 of 5483 are predicted secreted