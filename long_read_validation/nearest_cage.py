#!/usr/bin/env python3
# %%
import pandas as pd
from gtfparse import read_gtf

def is_within_peak(chrom,tss,strand,peaks):
    chr_df=peaks[peaks['chrom']==chrom]
    chr_df['within_cage']=chr_df.apply(lambda x: 'yes' if x['start'] <= tss <= x['end'] else 'no', axis=1)
    return(True if 'yes' in chr_df.within_cage.unique() else False)

def nearest_cage(tss,chrom,dict_sites):
    distances=[abs(site-tss) for site in dict_sites[chrom]]
    return(min(distances))
# %%
cage_peaks=pd.read_table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/cage_fantom5/monocytes.hg38.clusters.bed',usecols=[0,1,2],names=['chrom','start','end'])
transcripts=read_gtf('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoquant/five_conditions/combined/00_m64102e_220102_030451.flnc.transcript_models.gtf')[['seqname','feature','start','end','strand','transcript_id']]
transcripts=transcripts[(transcripts.feature=='transcript')&(transcripts.transcript_id.str.contains('transcript'))]
gencode_tss=read_gtf('/mnt/xomics/renees/data/ref_grch38/gencode/gencode.v39.annotation.gtf')[['seqname','feature','start','end','strand']]
gencode=gencode[gencode.feature=='transcript']
# %%
##first check for transcripts that have TSS *inside* a cage peak
pos=transcripts[transcripts.strand=='+'].merge(gencode[gencode.strand=='+'][['seqname','start']],how='left',indicator=True)
pos=pos[pos['_merge']!='both']
pos['inside_cage']=pos.apply(lambda x: is_within_peak(x['seqname'],x['start'],x['strand'],cage_peaks),axis=1)
neg=transcripts[transcripts.strand=='-'].merge(gencode[gencode.strand=='-'][['seqname','end']],how='left',indicator=True)
neg=neg[neg['_merge']!='both']
neg['inside_cage']=neg.apply(lambda x: is_within_peak(x['seqname'],x['end'],x['strand'],cage_peaks),axis=1)
# %%
##for those that remain, check nearby cage peaks
cage_peaks['combined']=cage_peaks[['start','end']].values.tolist()
cage=cage_peaks[['chrom','combined']].explode('combined').groupby('chrom').agg(list).to_dict()['combined']
pos_outside=pos[~pos.inside_cage]
pos_outside['min_dist']=pos_outside.apply(lambda x: nearest_cage(x['start'],x['seqname'],cage),axis=1)
neg_outside=neg[~neg.inside_cage]
neg_outside['min_dist']=neg_outside.apply(lambda x: nearest_cage(x['end'],x['seqname'],cage),axis=1)
