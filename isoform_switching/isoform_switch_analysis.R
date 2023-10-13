if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
    BiocManager::install()
}
library(IsoformSwitchAnalyzeR)
library(tidyverse)

##### isoquant #####
isoquant_expression <- readr::read_tsv(file = '/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoquant/five_conditions/combined/00_m64102e_220102_030451.flnc.transcript_model_grouped_counts_tmm.tsv')
names(isoquant_expression)[names(isoquant_expression) == "#feature_id"] <- "isoform_id"


myDesign <- data.frame(
  sampleID = c('calb','lps','polyic','rpmi','saur'),
  condition = c('calb','lps','polyic','rpmi','saur')
)

comparisons <- data.frame(
  condition_1 = c('rpmi','rpmi','rpmi','rpmi'),
  condition_2 = c('calb','lps','polyic','saur')
)

isoquantSwitchList <- importRdata(
  isoformCountMatrix   = isoquant_expression,
  designMatrix         = myDesign,
  comparisonsToMake    = comparisons,
  addIFmatrix          = FALSE,
  isoformExonAnnoation = '/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoquant/five_conditions/combined/00_m64102e_220102_030451.flnc.transcript_models.gtf',
  isoformNtFasta       = '/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoquant/five_conditions/combined/00_m64102e_220102_030451.flnc.transcript_models.fa'
)


swan_polyic <- read.table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/diff_exp/isoquant_based/genes_diff_isoform_polyic.csv',sep=',',header=TRUE,colClasses=(c("character","NULL","NULL","numeric")))
swan_polyic$condition_1 = 'rpmi'
swan_polyic$condition_2 = 'polyic'
names(swan_polyic)[names(swan_polyic) == "gid"] <- "gene_id"
swan_calb <- read.table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/diff_exp/isoquant_based/genes_diff_isoform_calb.csv',sep=',',header=TRUE,colClasses=(c("character","NULL","NULL","numeric")))
swan_calb$condition_1 = 'rpmi'
swan_calb$condition_2 = 'calb'
names(swan_calb)[names(swan_calb) == "gid"] <- "gene_id"
swan_saur <- read.table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/diff_exp/isoquant_based/genes_diff_isoform_saur.csv',sep=',',header=TRUE,colClasses=(c("character","NULL","NULL","numeric")))
swan_saur$condition_1 = 'rpmi'
swan_saur$condition_2 = 'saur'
names(swan_saur)[names(swan_saur) == "gid"] <- "gene_id"
swan_lps <- read.table('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/diff_exp/isoquant_based/genes_diff_isoform_lps.csv',sep=',',header=TRUE,colClasses=(c("character","NULL","NULL","numeric")))
swan_lps$condition_1 = 'rpmi'
swan_lps$condition_2 = 'lps'
names(swan_lps)[names(swan_lps) == "gid"] <- "gene_id"

swan=rbind(swan_calb,swan_saur,swan_polyic,swan_lps)

##Filter
isoquantSwitchListFilteredStrict <- preFilter(
    switchAnalyzeRlist = isoquantSwitchList,
    geneExpressionCutoff = 10,
    isoformExpressionCutoff = 0,
    removeSingleIsoformGenes = TRUE
)
#isoquantSwitchListFilteredStrict$gene_switch_q_value <- 1

## import isoform switching data in from swan and gene name mapping
gnames=read.csv('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoquant/five_conditions/combined/gene_mapping.csv')
#names(swan)[names(swan) == "gid"] <- "gene_id"
isoquantSwitchListFilteredStrict$isoformFeatures = merge(isoquantSwitchListFilteredStrict$isoformFeatures,swan)
isoquantSwitchListFilteredStrict$isoformFeatures$gene_switch_q_value = isoquantSwitchListFilteredStrict$isoformFeatures$adj_p_val
isoquantSwitchListFilteredStrict$isoformFeatures=isoquantSwitchListFilteredStrict$isoformFeatures[1:(length(isoquantSwitchListFilteredStrict$isoformFeatures)-1)]

isoquantSwitchListFilteredStrict$isoformFeatures = merge(isoquantSwitchListFilteredStrict$isoformFeatures,gnames,on='gene_id')
isoquantSwitchListFilteredStrict$isoformFeatures$gene_name = isoquantSwitchListFilteredStrict$isoformFeatures$Gene.name
isoquantSwitchListFilteredStrict$isoformFeatures=isoquantSwitchListFilteredStrict$isoformFeatures[1:(length(isoquantSwitchListFilteredStrict$isoformFeatures)-1)]


### If analysing (some) novel isoforms (else use CDS from ORF as explained in importRdata() )
isoquantSwitchListFilteredStrict <- addORFfromGTF( isoquantSwitchListFilteredStrict, '/mnt/xomics/renees/data/ref_grch38/gencode/gencode.v39.annotation.gtf' )
isoquantSwitchListFilteredStrict <- analyzeNovelIsoformORF( isoquantSwitchListFilteredStrict, analysisAllIsoformsWithoutORF = TRUE )

### Extract Sequences
isoquantSwitchListFilteredStrict <- extractSequence(
    isoquantSwitchListFilteredStrict, 
    pathToOutput = '/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/',
    onlySwitchingGenes = FALSE,
    alsoSplitFastaFile = TRUE
)

isoquantSwitchListAnalyzed <- analyzePFAM( #Added domain information to 15118 (74.1%) transcripts
    switchAnalyzeRlist   = isoquantSwitchListFilteredStrict,
    pathToPFAMresultFile = '/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/isoformSwitchAnalyzeR_hmmscanwebserver.txt',
    showProgress=FALSE
)

isoquantSwitchListAnalyzed <- analyzeCPAT( # Added coding potential to 20401 (100%) transcripts
    switchAnalyzeRlist   = isoquantSwitchListAnalyzed,
    pathToCPATresultFile = '/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/isoformSwitchAnalyzeR_cpat2',
    codingCutoff         = 0.725, # the coding potential cutoff we suggested for human
    removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
)

isoquantSwitchListAnalyzed <- analyzeSignalP( #run with signalp5.Added signal peptide information to 1360 (6.67%) transcripts
    switchAnalyzeRlist       = isoquantSwitchListAnalyzed,
    pathToSignalPresultFile  = '/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/signalp5_summary.signalp5'
)

isoquantSwitchListAnalyzed <- analyzeAlternativeSplicing(
    switchAnalyzeRlist = isoquantSwitchListAnalyzed,
    onlySwitchingGenes=FALSE,
    quiet=TRUE,
    alpha=1
)
#table( isoquantSwitchListAnalyzed$AlternativeSplicingAnalysis$IR )

isoquantSwitchListAnalyzed <- analyzeSwitchConsequences(
    isoquantSwitchListAnalyzed,
    consequencesToAnalyze = c('intron_retention', 'coding_potential','ORF_seq_similarity','NMD_status','domains_identified','signal_peptide_identified'),
    #alpha=1
)

## post analysis and plotting

#if wish to look at a particular subset
isoquantSwitchListAnalyzedSubset <- subsetSwitchAnalyzeRlist(
  isoquantSwitchListAnalyzed, 
  isoquantSwitchListAnalyzed$isoformFeatures$condition_2 == 'calb'
)
isoquantSwitchListAnalyzedSubset

#list top most significant switches
extractTopSwitches(
  isoquantSwitchListAnalyzedSubset, 
  filterForConsequences = TRUE, 
  n = 5, 
  sortByQvals = TRUE
)

#plot gene- and condition-specific switching
switchPlot(isoquantSwitchListAnalyzed, gene='NFKB1',condition1 = 'rpmi',condition2 = 'calb')
switchPlot(isoquantSwitchListAnalyzedSubset, gene='CFLAR')

pdf(file='/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/CIS/NFKB1.pdf',width = 8, height = 6)
switchPlot(isoquantSwitchListAnalyzed, gene='NFKB1',condition1 = 'rpmi',condition2 = 'calb')
dev.off()

pdf(file='/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/CIS/CASP1.pdf',width = 8, height = 6)
switchPlot(isoquantSwitchListAnalyzed, gene='CASP1',condition1 = 'rpmi',condition2 = 'polyic')
dev.off()


#genome wide plots
pdf(file='/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/splicingsummary.pdf',width = 8, height = 6)
extractSplicingSummary(isoquantSwitchListAnalyzed)
dev.off()

pdf(file='/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/splicingenrichment.pdf',width = 10, height = 6)
extractSplicingEnrichment(
  isoquantSwitchListAnalyzed,
  splicingToAnalyze='all'
)
dev.off()

pdf(file='/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/consequenceenrichment.pdf',width = 10, height = 8)
extractConsequenceEnrichment(
  isoquantSwitchListAnalyzed,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = TRUE,
  returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)
dev.off()

write.csv(extractConsequenceEnrichment(
  isoquantSwitchListAnalyzed,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = TRUE,
  returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
),'/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/cons_enrichment.csv', row.names=FALSE,quote=FALSE)

#my own plots - export to python

write.csv(isoquantSwitchListAnalyzed$isoformFeatures, '/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/isoform_features.csv', row.names=FALSE,quote=FALSE)
write.csv(isoquantSwitchListAnalyzed$domainAnalysis, '/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/domain_analysis.csv', row.names=FALSE,quote=FALSE)
write.csv(isoquantSwitchListAnalyzed$signalPeptideAnalysis, '/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/signalpep_analysis.csv', row.names=FALSE,quote=FALSE)
write.csv(isoquantSwitchListAnalyzed$switchConsequence, '/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/switch_consequence.csv', row.names=FALSE,quote=FALSE)
write.csv(isoquantSwitchListAnalyzed$AlternativeSplicingAnalysis, '/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/AS_analysis.csv', row.names=FALSE,quote=FALSE)

#save all immune-related switch plots to files

immune=readr::read_tsv(file = '/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/immune_switches.tsv')

for (i in unique(isoquantSwitchListAnalyzed$isoformFeatures$gene_name)){
  dir.create(paste('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/CIS/',i,'/',sep=''))
}

for (i in 1:nrow(isoquantSwitchListAnalyzed$isoformFeatures)) {
  print(paste(isoquantSwitchListAnalyzed$isoformFeatures$gene_name[i],'&',isoquantSwitchListAnalyzed$isoformFeatures$condition_2[i],sep=' '))
  png(file=paste('/mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/isoformswitchR/old/CIS/',isoquantSwitchListAnalyzed$isoformFeatures$gene_name[i],'/',isoquantSwitchListAnalyzed$isoformFeatures$gene_name[i],'_',isoquantSwitchListAnalyzed$isoformFeatures$condition_2[i],'.png',sep=''),width = 1058, height = 777,res=130)
  switchPlot(isoquantSwitchListAnalyzed, gene=isoquantSwitchListAnalyzed$isoformFeatures$gene_name[i], condition1 = 'rpmi', condition2 = isoquantSwitchListAnalyzed$isoformFeatures$condition_2[i])
  dev.off()
}

