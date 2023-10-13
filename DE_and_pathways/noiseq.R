library(NOISeq)
library(tidyverse)

gene = read.table("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/isoquant/five_conditions/combined/LR_gene_counts_50cpm.tsv",row.names = 1)
colnames(gene)=c("calb","lps","polyic","rpmi","saur")

transcript = read.table("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/isoquant/five_conditions/combined/LR_transcript_counts.tsv",row.names = 1)
colnames(transcript)=c("calb","lps","polyic","rpmi","saur")

#for the grouped set
transcript = read.table("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/3_prime_validation/grouped_rawcounts_filtered.tsv", header = TRUE)
transcript=column_to_rownames(transcript,var="rep")


#biomart= read.table("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/biomart_mapping.tsv",sep='\t',header=T,fill = T)
#genebiotypes=unique(biomart[c('Gene.stable.ID.version','Gene.type','Gene...GC.content')])
#rownames(genebiotypes) <- genebiotypes$'Gene.stable.ID.version'
#gccontent=as.matrix(genebiotypes['Gene...GC.content'])
#genebiotypes=as.matrix(genebiotypes['Gene.type'])
#rownames(biomart)=biomart$Transcript.stable.ID.version
#transcriptbiotypes=as.matrix(biomart['Transcript.type'])
#transcriptlength=as.matrix(biomart$Transcript.length..including.UTRs.and.CDS.)


myfactors = data.frame(Tissue=c("calb","lps","polyic","rpmi","saur"),
                       TissueRun = c("calb_1","lps_1","polyic_1","rpmi_1","saur_1"),
                       Run = c(rep("R1", 5)))
myfactors

#for gene
mydata <- readData(data=gene,factors=myfactors)
mydata
calb <- noiseq(mydata, factor = "Tissue", conditions=c("calb","rpmi"),replicates = "no",norm = 'tmm')
saur <- noiseq(mydata, factor = "Tissue", conditions=c("saur","rpmi"),replicates = "no",norm = 'tmm')
lps <- noiseq(mydata, factor = "Tissue", conditions=c("lps","rpmi"),replicates = "no",norm = 'tmm')
polyic <- noiseq(mydata, factor = "Tissue", conditions=c("polyic","rpmi"),replicates = "no",norm = 'tmm')

lps_genes=degenes(lps,q=0.95,M=NULL)
calb_genes=degenes(calb,q=0.95,M=NULL)
saur_genes=degenes(saur,q=0.95,M=NULL)
polyic_genes=degenes(polyic,q=0.95,M=NULL)
write.csv(lps_genes,"~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/noiseq/diff_gene_lps_50cpm.csv", quote = FALSE)
write.csv(calb_genes,"~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/noiseq/diff_gene_calb_50cpm.csv", quote = FALSE)
write.csv(saur_genes,"~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/noiseq/diff_gene_saur_50cpm.csv", quote = FALSE)
write.csv(polyic_genes,"~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/noiseq/diff_gene_polyic_50cpm.csv", quote = FALSE)



#for transcript
mydata <- readData(data=transcript,factors=myfactors)
mydata
calb <- noiseq(mydata, factor = "Tissue", conditions=c("calb","rpmi"),replicates = "no",norm = 'tmm')
saur <- noiseq(mydata, factor = "Tissue", conditions=c("saur","rpmi"),replicates = "no",norm = 'tmm')
lps <- noiseq(mydata, factor = "Tissue", conditions=c("lps","rpmi"),replicates = "no",norm = 'tmm')
polyic <- noiseq(mydata, factor = "Tissue", conditions=c("polyic","rpmi"),replicates = "no",norm = 'tmm')

lps_transcripts=degenes(lps,q=0.95,M=NULL)
calb_transcripts=degenes(calb,q=0.95,M=NULL)
saur_transcripts=degenes(saur,q=0.95,M=NULL)
polyic_transcripts=degenes(polyic,q=0.95,M=NULL)
write.csv(lps_transcripts,"~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/noiseq/diff_transcript_lps.csv", quote = FALSE)
write.csv(calb_transcripts,"~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/noiseq/diff_transcript_calb.csv", quote = FALSE)
write.csv(saur_transcripts,"~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/noiseq/diff_transcript_saur.csv", quote = FALSE)
write.csv(polyic_transcripts,"~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/noiseq/diff_transcript_polyic.csv", quote = FALSE)

