## RNA analysis

library(tidyverse)

transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(t_df)
}
#transcript_expression=read_tsv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/talon/min5_talon_abundance_filtered.tsv",col_select = c("annot_gene_id","annot_gene_name","annot_transcript_id")) %>% inner_join(read_csv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/TMN_5cond_Qfilter_PI.csv"))
transcript_expression=read_tsv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/isoquant/five_conditions/combined/LR_transcript_counts_withinfo.tsv")

gene_expression=read_csv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/TMN_genesum_5cond_Qfilter.csv")
top_calb=read_csv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/SR_top_genes/top_log2_CAlb_24hr.csv") #can use padj later if desired. starting with log2fc
top_saur=read_csv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/SR_top_genes/top_log2_SAur_24hr.csv")
top_lps=read_csv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/SR_top_genes/top_log2_LPS_24hr.csv")
top_polyic=read_csv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/SR_top_genes/top_log2_PolIC_24hr.csv")
#lr=read_tsv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/talon/min5_talon_abundance_filtered.tsv",col_select = c("annot_gene_id","annot_gene_name","annot_transcript_id")) # import gene symbols for later
lr=read_tsv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/isoquant/LR_transcript_counts_withinfo.tsv",col_select = c("gene_id","Gene.name","transcript_id")) # import gene symbols for later
lr_transcript= lr[!duplicated(lr$annot_transcript_id), ] 
lr=read_tsv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/talon/min5_talon_abundance_filtered.tsv",col_select = c("annot_gene_id","annot_gene_name")) # import gene symbols for later
lr= lr[!duplicated(lr$annot_gene_id), ] 

#pid only analysis - not informative
pid_panel=read_csv("~/Desktop/narrativum/unsolved_cases/PID_genelist.csv",col_names = c("annot_gene_name","blank"),col_select = c("annot_gene_name"), trim_ws = TRUE) %>% left_join(lr)#import pid genes for comparison
pid_expr=pid_panel %>% inner_join(gene_expression) %>% select(-annot_gene_name) %>% column_to_rownames(var="annot_gene_id")
#pid_expr=na_if(pid_expr,0) #cannot do PCA with infinite values
#pid_expr = pid_expr %>% drop_na() # get rid of columns where genes were not found expressed in a condition
pid_expr=transpose_df(log10(pid_expr + 0.5))
pid_expr=column_to_rownames(pid_expr,var="rowname")
x.pca <- prcomp(pid_expr,scale=TRUE,center=TRUE)

#subset analysis - gene level
gene_expression$Transcript=sub("\\.[0-9]+","",gene_expression$annot_gene_id)
calb_expr=gene_expression %>% inner_join(top_calb[, "Transcript"]) %>% column_to_rownames(var="annot_gene_id")%>% select(-Transcript)
calb_expr= log10(calb_expr+0.5)
calb_expr=calb_expr[,2:5]-calb_expr[,1] 
calb_expr= calb_expr %>% rownames_to_column("annot_gene_id") %>% inner_join(lr) %>% select(-annot_gene_id) %>% column_to_rownames(var="annot_gene_name")#avoid log10 errors by adding 0.5 to whole df
heatmap(as.matrix(calb_expr),scale="column")
saur_expr=gene_expression %>% inner_join(top_saur[, "Transcript"]) %>% column_to_rownames(var="annot_gene_id")%>% select(-Transcript)
saur_expr= log10(saur_expr+0.5) # normalize log 10
saur_expr=saur_expr[,2:5]-saur_expr[,1] # subtract RPMI from rest
saur_expr= saur_expr %>% rownames_to_column("annot_gene_id") %>% inner_join(lr) %>% select(-annot_gene_id) %>% column_to_rownames(var="annot_gene_name")#avoid log10 errors by adding 0.5 to whole df
heatmap(as.matrix(saur_expr),scale="column")
lps_expr=gene_expression %>% inner_join(top_lps[, "Transcript"]) %>% column_to_rownames(var="annot_gene_id")%>% select(-Transcript)
lps_expr= log10(lps_expr+0.5)
lps_expr=lps_expr[,2:5]-lps_expr[,1] 
lps_expr= lps_expr %>% rownames_to_column("annot_gene_id") %>% inner_join(lr) %>% select(-annot_gene_id) %>% column_to_rownames(var="annot_gene_name")#avoid log10 errors by adding 0.5 to whole df
heatmap(as.matrix(lps_expr),scale="column")
polyic_expr=gene_expression %>% inner_join(top_polyic[, "Transcript"]) %>% column_to_rownames(var="annot_gene_id")%>% select(-Transcript)
polyic_expr= log10(polyic_expr+0.5)
polyic_expr=polyic_expr[,2:5]-polyic_expr[,1] 
polyic_expr= polyic_expr %>% rownames_to_column("annot_gene_id") %>% inner_join(lr) %>% select(-annot_gene_id) %>% column_to_rownames(var="annot_gene_name")#avoid log10 errors by adding 0.5 to whole df
heatmap(as.matrix(polyic_expr),scale="column")


#all gene analysis
rna=gene_expression %>% column_to_rownames(var="annot_gene_id") %>% select(-Transcript)
rna=log10(rna+0.5) #cannot do PCA with infinite values
rna = rna[,2:5]-rna[,1] # subtract RPMI values
rna=transpose_df(rna)
rna=column_to_rownames(rna,var="rowname")
x.pca <- prcomp(rna,scale=TRUE,center=TRUE)

#separate category where RPMI=0
#TODO

#all transcript analysis
isoform=transcript_expression %>% select(c("transcript_id","rpmi","saur","lps","calb","polyic")) %>% column_to_rownames(var="transcript_id") %>% rename("s.aur" = "saur","c.alb"="calb")
isoform[isoform==0]=NA
isoform=isoform[complete.cases(isoform),]
isoform=transpose_df(log10(isoform))
isoform=column_to_rownames(isoform,var="rowname")
x.pca <- prcomp(isoform,scale=TRUE,center=TRUE)

transcript_expr=transcript_expression %>% column_to_rownames(var="annot_transcript_id")
transcript_expr_abs=log10(transcript_expr[,3:7]+0.5)
transcript_expr_abs=transcript_expr_abs[,2:5]-transcript_expr_abs[,1]
transcript_expr_abs=transcript_expr_abs %>% bind_cols(transcript_expr[,c(2,8:12)]) %>% filter(rpmi_pi==0) #solely isoforms that are not expressed in control!!
transcript_expr_abs$max = apply(transcript_expr_abs[,1:4], 1, FUN = max)
transcript_expr_abs_hm=head(transcript_expr_abs[order(-transcript_expr_abs$max),],1000) %>% select(-max) 
heatmap(as.matrix(transcript_expr_abs_hm[,1:4]),scale="column",labRow=transcript_expr_abs$annot_gene_name)
heatmap(as.matrix(transcript_expr_abs_hm[,6:10]),scale="column",labRow=transcript_expr_abs$annot_gene_name)


calb_iso=transcript_expression %>% filter(calb_pi-rpmi_pi>30) %>% arrange(desc(calb)) %>% head(100)
saur_iso=transcript_expression %>% filter(saur_pi-rpmi_pi>30) %>% arrange(desc(saur))
lps_iso=transcript_expression %>% filter(lps_pi-rpmi_pi>30) %>% arrange(desc(lps))
polyic_iso=transcript_expression %>% filter(polyic_pi-rpmi_pi>30) %>% arrange(desc(polyic))
isoform=head(calb_iso,100) %>% bind_rows(head(saur_iso,100)) %>% bind_rows(head(lps_iso,100)) %>% bind_rows(head(polyic_iso,100)) %>% distinct() %>% column_to_rownames(var="annot_transcript_id")

isoform=transcript_expression %>% select(c("annot_transcript_id","rpmi","saur","lps","calb","polyic")) %>% column_to_rownames(var="annot_transcript_id")
isoform=log10(isoform+0.5)
isoform=isoform[,2:5]-isoform[,1] 
isoform$max = apply(isoform, 1, FUN = max)
isoform_hm=head(isoform[order(-isoform$max),],50) %>% select(-max)
heatmap(as.matrix(isoform_hm),scale="column")

#PCA plots and info
var_explained <- x.pca$sdev^2/sum(x.pca$sdev^2)
screeplot(x.pca) #first PC responsible for 50% of the variance
group.colors <- c(rpmi="#949494",lps = "#cc78bc", c.alb= "#ca9161", s.aur ="#d55e00", polyic = "#ece133")
pca_df=x.pca$x %>% as.data.frame 
pca_df$stimulus <- rownames(pca_df)
# %>%
#  rownames_to_column("stimulus") %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(aes(color=stimulus),size=10) + theme_bw(base_size=24,base_family = 'sans') + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + theme(legend.position="top") + scale_fill_manual(values=group.colors)
ggplot()+geom_point(data=pca_df, aes(x=PC1, y=PC2, color=stimulus), size=10) +scale_color_manual(values = group.colors)+ theme_bw(base_size=24,base_family = 'sans') + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + theme(legend.position="top")
# observe PCA
PC1 <- x.pca$rotation[,1]
PC1_scores <- abs(PC1)
PC1_scores_ordered <- sort(PC1_scores, decreasing = TRUE)
#PC1_scores_ordered <- rownames_to_column(as.data.frame(PC1_scores_ordered), var = "annot_gene_id") %>% as_tibble() %>% left_join(lr %>% select(c("annot_gene_id","annot_gene_name")) %>% distinct()) #put a gene symbol to the important ones
isoform_PC1 <- rownames_to_column(as.data.frame(PC1_scores_ordered), var = "annot_transcript_id") %>% as_tibble() %>% left_join(lr) #put a gene symbol to the important ones
venn.diagram(
  x = list(isoform_PC1_head$annot_gene_id, PC1_scores_ordered_head$annot_gene_id),
  filename = '~/Desktop/overlap.png',
  category.names = c("Transcript level" , "Gene level"),
  output=TRUE
  )

#heatmaps
pid=transpose_df(pid_expr) %>% inner_join(lr,by=c("rowname"="annot_gene_id")) %>% select(-rowname) %>% column_to_rownames(var="annot_gene_name")
#pid$row_std = apply(pid, 1, sd)
#pid_hm=head(pid[order(-pid$row_std),],30) %>% select(-row_std)
pid$max = apply(pid, 1, FUN = max)
pid$diff = pid$max - pid$rpmi # look for the ones with the biggest difference max-rpmi
pid_hm=head(pid[order(-pid$diff),],50) %>% select(-c(max,diff))
heatmap(as.matrix(pid_hm),scale="column")

all_rna=transpose_df(rna) %>% inner_join(lr,by=c("rowname"="annot_gene_id")) %>% select(-rowname) %>% distinct(annot_gene_name, .keep_all = TRUE) %>% column_to_rownames(var="annot_gene_name")
all_rna$max = apply(all_rna, 1, FUN = max)
all_rna_hm=head(all_rna[order(-all_rna$max),],1000) %>% select(-max)
heatmap(as.matrix(all_rna_hm),scale="column")

# subset diff expressed genes 
de_calb=read_csv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/diff_exp/genes_diff_isoform_calb_Qfilter.csv") #%>% column_to_rownames(var="annot_transcript_id")
de_calb=log10(de_calb)
de_calb_hm=head(de_calb[order(-de_calb$row_std),],30) %>% select(-row_std)
heatmap(as.matrix(de_calb_hm),scale="column")
de_saur=read_csv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/diff_exp/genes_diff_isoform_saur_Qfilter.csv") #%>% column_to_rownames(var="annot_transcript_id")
de_saur$row_std = apply(de_saur, 1, sd)
de_saur_hm=head(de_saur[order(-de_saur$row_std),],30) %>% select(-row_std)
heatmap(as.matrix(de_saur_hm),scale="column")
de_polyic=read_csv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/diff_exp/genes_diff_isoform_polyic_Qfilter.csv") #%>% column_to_rownames(var="annot_transcript_id")
de_polyic$row_std = apply(de_polyic, 1, sd)
de_polyic_hm=head(de_polyic[order(-de_polyic$row_std),],30) %>% select(-row_std)
heatmap(as.matrix(de_polyic_hm),scale="column")
de_lps=read_csv("~/Desktop/narrativum/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/diff_exp/genes_diff_isoform_lps_Qfilter.csv") #%>% column_to_rownames(var="annot_transcript_id")
de_lps$row_std = apply(de_lps, 1, sd)
de_lps_hm=head(de_lps[order(-de_lps$row_std),],30) %>% select(-row_std)
heatmap(as.matrix(de_lps_hm),scale="column")

intersection=de_calb %>% inner_join(de_saur[,'annot_transcript_id']) %>% inner_join(de_polyic[,'annot_transcript_id']) %>% inner_join(de_lps[,'annot_transcript_id']) %>% column_to_rownames(var="annot_transcript_id")
heatmap(as.matrix(intersection[,11:15]),scale="column",labRow=intersection$annot_gene_name)

interesting_group1=de_calb %>% inner_join(de_saur[,'annot_transcript_id']) %>% inner_join(de_polyic[,'annot_transcript_id']) %>% anti_join(de_lps[,'annot_transcript_id']) %>% column_to_rownames(var="annot_transcript_id")
heatmap(as.matrix(interesting_group1[,11:15]),scale="column",labRow=interesting_group1$annot_gene_name)

group1=de_polyic %>% anti_join(de_lps[,'annot_transcript_id']) %>% anti_join(de_calb[,'annot_transcript_id']) %>% anti_join(de_saur[,'annot_transcript_id']) #%>% column_to_rownames(var="annot_transcript_id")
#group1_hm=head(group1[order(-group1$polyic),],50)
group2=de_calb %>% anti_join(de_lps[,'annot_transcript_id']) %>% anti_join(de_polyic[,'annot_transcript_id']) %>% anti_join(de_saur[,'annot_transcript_id']) #%>% column_to_rownames(var="annot_transcript_id")
#group2_hm=head(group1[order(-group1$calb),],50)
group3=de_saur %>% anti_join(de_lps[,'annot_transcript_id']) %>% anti_join(de_calb[,'annot_transcript_id']) %>% anti_join(de_polyic[,'annot_transcript_id']) #%>% column_to_rownames(var="annot_transcript_id")
#group3_hm=head(group1[order(-group1$saur),],50)
group4=de_lps %>% anti_join(de_polyic[,'annot_transcript_id']) %>% anti_join(de_calb[,'annot_transcript_id']) %>% anti_join(de_saur[,'annot_transcript_id']) #%>% column_to_rownames(var="annot_transcript_id")
#group4_hm=head(group1[order(-group1$lps),],50)
combi=group1 %>% bind_rows(group2) %>% bind_rows(group3) %>% bind_rows(group4) %>% distinct() %>% column_to_rownames(var="annot_transcript_id")
combi$max = apply(combi[,3:7], 1, FUN = max)
combi$diff = combi$max - combi$rpmi # look for the ones with the biggest difference max-min
combi_hm=head(combi[order(-combi$diff),],60) %>% select(-c(max,diff))
combi_hm=head(combi[order(combi$rpmi_pi),],60) %>% select(-c(max,diff))
heatmap(as.matrix(combi_hm[,11:15]),scale="column",labRow=combi_hm$annot_gene_name)

heatmap(as.matrix(group4[,12:16]),scale="column",labRow=combi_hm$annot_gene_name)

intersection$novel <- "known"
intersection$novel[grepl("novelT",intersection$annot_transcript_name,fixed=TRUE)] <- "novel"
combi$novel <- "known"
combi$novel[grepl("novelT",combi$annot_transcript_name,fixed=TRUE)] <- "novel"
