library(tidyverse)
library(factoextra)
library(ggrepel)
library(ggpubr)
library(ComplexHeatmap)
library(cluster)
# Read in the proteomics data from TSV file
p<-read.table("~/Desktop/narrativum/host_pathogen_PID/pbmc_secretome/generated_data/metamorpheus/uniprot_results/Task2CalibrationTask/QuantifiedProteins.tsv",header=T,sep='\t')#,row.names=1)
proteins_tb= as_tibble(p)
proteins_tb=proteins_tb[proteins_tb$Organism=="",]
proteins_tb=proteins_tb %>% 
  rename(
    rpmi_1=Intensity_B09905_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun.1.calib,
    lps_1=Intensity_B09907_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun.2.calib,
    saur_1=Intensity_B09909_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun.3.calib,
    polyic_1=Intensity_B09911_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun.4.calib,
    calb_1=Intensity_B09913_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun.5.calib,
    #rpmi_1.2=Intensity_B09915_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun.6.calib,
    rpmi_2=Intensity_B09917_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun.7.calib,
    lps_2=Intensity_B09919_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun.8.calib,
    saur_2=Intensity_B09921_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun.9.calib,
    polyic_2=Intensity_B09923_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun.10.calib,
    calb_2=Intensity_B09925_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun.11.calib,
    #rpmi_2.2=Intensity_B09925_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun.12.calib,
    rpmi_3=Intensity_B10105_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.13.calib,
    lps_3=Intensity_B10107_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.14.calib,
    saur_3=Intensity_B10109_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.15.calib,
    polyic_3=Intensity_B10111_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.16.calib,
    calb_3=Intensity_B10113_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.17.calib,
    #rpmi_3.2=Intensity_B10115_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.18.calib,
    rpmi_4.1=Intensity_B10117_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.19.calib,
    lps_4=Intensity_B10119_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.20.calib,
    saur_4=Intensity_B10121_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.21.calib,
    polyic_4=Intensity_B10123_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.22.calib,
    calb_4=Intensity_B10125_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.23.calib,
    #rpmi_4.2=Intensity_B10127_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.24,
    lps_5=Intensity_B10129_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.25.calib,
    polyic_5=Intensity_B10131_Bp_WNE19_trap10_PRC.5493_LennartMartens_shotgun.26.calib,
    rpmi_5.1=Intensity_B10133_Bp_WNE17_trap9_PRC.5493_LennartMartens_shotgun_5095_1.calib,
    calb_5=Intensity_B10135_Bp_WNE17_trap9_PRC.5493_LennartMartens_shotgun_5095_2.calib,
    saur_5=Intensity_B10139_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun_5095_3.calib
    #rpmi_5.2=Intensity_B10141_Bp_WNE17_trap10_PRC.5493_LennartMartens_shotgun_5095_4.calib
    )
proteins_tb=select(proteins_tb, -c(Gene.Name,Organism,X)) #%>% column_to_rownames('Protein.Groups')
#proteins_tb=na_if(proteins_tb,"")
proteins_tb[proteins_tb==0] <- NA
proteins_tb = proteins_tb %>% drop_na() # get rid of columns where protein groups were not found in a condition
proteins_tb=column_to_rownames(proteins_tb,var="Protein.Groups") #set protein groups as index

transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(t_df)
}

proteins_tb=transpose_df(proteins_tb)
proteins_tb=column_to_rownames(proteins_tb,var="rowname")

# define x.pca to receive the object of the princomp function applied to x.mat.
x.pca <- prcomp(proteins_tb,scale=TRUE,center=TRUE)
# print out the x.pca dataset.
x.pca

var_explained <- x.pca$sdev^2/sum(x.pca$sdev^2)


# what are the names of the objects in the x.pca principal component object?
names(x.pca)
attributes(x.pca)
x.pca$sdev     #std of PC
x.pca$sdev^2   #eigenvalues
x.pca$rotation #PC loadings 
screeplot(x.pca)

x.pca$x   #PC scores, values of each obs for each PC

#  compute a summary of the x.pca object.
summary(x.pca)

# plot the PCs
x.pca$x %>% 
  as.data.frame %>%
  rownames_to_column("stimulus_individual") %>%
  separate(stimulus_individual, c("stimulus","individual"),extra="drop") %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(aes(color=individual),size=5) + theme_bw(base_size=20) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + theme(legend.position="top")


#  compute the covariance matrix for x.mat.
x.cov <- cov(proteins_tb)

#  run the eigen function on the covariance of proteins_tb.
eigen(x.cov)

#  compare this to what you got from princomp.
x.pca$sdev
x.pca$rotation
#the sqaure of the singular values equal th eigenvalues
x.pca$sdev^2

#  compute the correlations between the principal components scores and the
#  original data. 
#  here we create a new matrix, xxx by binding the columns of proteins_tb to the 
#  columns of x.pca$scores - the pca scores of proteins_tb.
xxx <- cbind(proteins_tb,x.pca$x)
#  compute the correlation between the original variables and pca scores.
cor(xxx)


## heatmaps proteins - figure 7
df=read.table("~/Desktop/narrativum/host_pathogen_PID/pbmc_secretome/generated_data/metamorpheus/analysis_files/uniprot_bfc_wsha_filt_sqrtnorm_znorm_indivnorm.tsv",header=T,sep='\t')
df=column_to_rownames(df,var="Protein.Group")
conditions=colnames(df)
conditions=str_remove(conditions,"\\d")
cond=HeatmapAnnotation(condition= conditions,col = list(condition = c(rpmi="#949494","calb" = "#440154", "saur" = "#31688e", "lps" = "#fde725","polyic" = "#35b779")))
Heatmap(df, name = "mat", row_km = 4,column_km = 3,show_row_names = FALSE,bottom_annotation = cond)

fviz_nbclust(df, kmeans, method = "wss")
set.seed(82)
km=kmeans(df,4,nstart = 100)
names(km$cluster[km$cluster == 4]) %>% clipr::write_clip()
pdf(file = "~/Desktop/heatmap_prots.pdf", width= 7, height = 7)
Heatmap(df, name = "mat", row_split = km$cluster,column_km = 3,show_row_names = FALSE,bottom_annotation = cond)
dev.off()


##volcano plots - not in manuscript
bayesianfc <- read.table("~/Desktop/narrativum/host_pathogen_PID/pbmc_secretome/generated_data/metamorpheus/uniprot_results/Task2CalibrationTask/BayesianFoldChangeAnalysisProcessed.tsv",header=T,sep='\t')
bayesianfc <- as_tibble(bayesianfc)
calb <- bayesianfc %>% select(-contains(".1")) %>% select(-contains(".2")) %>% select(-contains(".3")) %>% select(-contains("X"))
calb=calb[(calb$Organism=="Homo sapiens"),]
saur <- bayesianfc %>% select(contains(".3")) %>% select(-contains("X"))
colnames(saur) <- gsub('.{2}$', '', colnames(saur))
saur=saur[(saur$Organism=="Homo sapiens"),]
polyic <- bayesianfc %>% select(contains(".2")) %>% select(-contains("X"))
colnames(polyic) <- gsub('.{2}$', '', colnames(polyic))
polyic=polyic[(polyic$Organism=="Homo sapiens"),]
lps <- bayesianfc %>% select(contains(".1")) %>% select(-contains("X"))
colnames(lps) <- gsub('.{2}$', '', colnames(lps))
lps=lps[(lps$Organism=="Homo sapiens"),]

calb$genename=str_replace_all(calb$Gene,"primary:","")
saur$genename=str_replace_all(saur$Gene,"primary:","")
lps$genename=str_replace_all(lps$Gene,"primary:","")
polyic$genename=str_replace_all(polyic$Gene,"primary:","")
calb$genename[duplicated(calb$genename)] <- ""
saur$genename[duplicated(saur$genename)] <- ""
lps$genename[duplicated(lps$genename)] <- ""
polyic$genename[duplicated(polyic$genename)] <- ""

sig=rbind(calb[calb$False.Discovery.Rate<0.05,],saur[saur$False.Discovery.Rate<0.05,],lps[lps$False.Discovery.Rate<0.05,],polyic[polyic$False.Discovery.Rate<0.05,])


calb_plot <- calb %>% 
  mutate(calb = ifelse(Protein.Log2.Fold.Change > 1 & Posterior.Error.Probability < 0.05, 
                       "up",
                       ifelse(Protein.Log2.Fold.Change < -1 & Posterior.Error.Probability < 0.05, 
                              "down",
                              "ns"))) %>% 
  ggplot(aes(x = Protein.Log2.Fold.Change, y = -log10(False.Discovery.Rate))) +
  geom_point(aes(color = factor(calb)), size = 0.5) +
  scale_color_manual(
    values = c(
      ns = "black", 
      up = "blue", 
      down = "red"
    ), 
    labels = c(
      ns = "NS", 
      up = paste("PEP < 0.05 & Log2FoldChange > 1"), 
      down = paste("PEP < 0.05 & Log2FoldChange < -1"))
  ) + 
  theme_minimal() + 
  facet_grid(~Treatment.Condition) +
  geom_text_repel(aes(label = ifelse(calb != "ns",Gene.name, "")),box.padding = 0.7, max.overlaps = Inf) + 
  labs(x="log2FC",y="-log10 Pval", color = "DEG") +
  scale_y_continuous(limits = c(0, 3))+
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  theme(legend.position = "right", 
        aspect.ratio = 0.75,
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_rect(fill="white"),
        panel.grid.major.x = element_line(size = 0.25, color = "grey"),
        panel.grid.major.y = element_line(size = 0.25, color = "grey"),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color ="black", size = 14),
        axis.title.y = element_text(color = "black", size = 14),
        plot.title = element_text(color = "black", size = 16),
        plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm"),
        plot.tag = element_text(color ="black", size= 20, face="bold")
  )

saur_plot <- saur %>% 
  mutate(saur = ifelse(Protein.Log2.Fold.Change > 1 & Posterior.Error.Probability < 0.05, 
                       "up",
                       ifelse(Protein.Log2.Fold.Change < -1 & Posterior.Error.Probability < 0.05, 
                              "down",
                              "ns"))) %>% 
  ggplot(aes(x = Protein.Log2.Fold.Change, y = -log10(False.Discovery.Rate))) +
  geom_point(aes(color = factor(saur)), size = 0.5) +
  scale_color_manual(
    values = c(
      ns = "black", 
      up = "blue", 
      down = "red"
    ), 
    labels = c(
      ns = "NS", 
      up = paste("PEP < 0.05 & Log2FoldChange > 1"), 
      down = paste("PEP < 0.05 & Log2FoldChange < -1"))
  ) + 
  theme_minimal() + 
  facet_grid(~Treatment.Condition) +
  geom_text_repel(aes(label = ifelse(saur != "ns",Gene.name, "")),box.padding = 0.7, max.overlaps = Inf) + 
  labs(x="log2FC",y="-log10 Pval", color = "DEG") +
  scale_y_continuous(limits = c(0, 3))+
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  theme(legend.position = "right", 
        aspect.ratio = 0.75,
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_rect(fill="white"),
        panel.grid.major.x = element_line(size = 0.25, color = "grey"),
        panel.grid.major.y = element_line(size = 0.25, color = "grey"),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color ="black", size = 14),
        axis.title.y = element_text(color = "black", size = 14),
        plot.title = element_text(color = "black", size = 16),
        plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm"),
        plot.tag = element_text(color ="black", size= 20, face="bold")
  )

lps_plot <- lps %>% 
  mutate(lps = ifelse(Protein.Log2.Fold.Change > 1 & Posterior.Error.Probability < 0.05, 
                      "up",
                      ifelse(Protein.Log2.Fold.Change < -1 & Posterior.Error.Probability < 0.05, 
                             "down",
                             "ns"))) %>% 
  ggplot(aes(x = Protein.Log2.Fold.Change, y = -log10(False.Discovery.Rate))) +
  geom_point(aes(color = factor(lps)), size = 0.5) +
  scale_color_manual(
    values = c(
      ns = "black", 
      up = "blue", 
      down = "red"
    ), 
    labels = c(
      ns = "NS", 
      up = paste("PEP < 0.05 & Log2FoldChange > 1"), 
      down = paste("PEP < 0.05 & Log2FoldChange < -1"))
  ) + 
  theme_minimal() + 
  facet_grid(~Treatment.Condition) +
  geom_text_repel(aes(label = ifelse(lps != "ns",Gene.name, "")),box.padding = 0.7, max.overlaps = Inf) + 
  labs(x="log2FC",y="-log10 Pval", color = "DEG") +
  scale_y_continuous(limits = c(0, 3))+
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  theme(legend.position = "right", 
        aspect.ratio = 0.75,
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_rect(fill="white"),
        panel.grid.major.x = element_line(size = 0.25, color = "grey"),
        panel.grid.major.y = element_line(size = 0.25, color = "grey"),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color ="black", size = 14),
        axis.title.y = element_text(color = "black", size = 14),
        plot.title = element_text(color = "black", size = 16),
        plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm"),
        plot.tag = element_text(color ="black", size= 20, face="bold")
  )

polyic_plot <- polyic %>% 
  mutate(polyic = ifelse(Protein.Log2.Fold.Change > 1 & Posterior.Error.Probability < 0.05, 
                         "up",
                         ifelse(Protein.Log2.Fold.Change < -1 & Posterior.Error.Probability < 0.05, 
                                "down",
                                "ns"))) %>% 
  ggplot(aes(x = Protein.Log2.Fold.Change, y = -log10(False.Discovery.Rate))) +
  geom_point(aes(color = factor(polyic)), size = 0.5) +
  scale_color_manual(
    values = c(
      ns = "black", 
      up = "blue", 
      down = "red"
    ), 
    labels = c(
      ns = "NS", 
      up = paste("PEP < 0.05 & Log2FoldChange > 1"), 
      down = paste("PEP < 0.05 & Log2FoldChange < -1"))
  ) + 
  theme_minimal() + 
  facet_grid(~Treatment.Condition) +
  geom_text_repel(aes(label = ifelse(polyic != "ns",Gene.name, "")),box.padding = 0.7, max.overlaps = Inf) + 
  labs(x="log2FC",y="-log10 Pval", color = "DEG") +
  scale_y_continuous(limits = c(0, 3))+
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  theme(legend.position = "right", 
        aspect.ratio = 0.75,
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_rect(fill="white"),
        panel.grid.major.x = element_line(size = 0.25, color = "grey"),
        panel.grid.major.y = element_line(size = 0.25, color = "grey"),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color ="black", size = 14),
        axis.title.y = element_text(color = "black", size = 14),
        plot.title = element_text(color = "black", size = 16),
        plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm"),
        plot.tag = element_text(color ="black", size= 20, face="bold")
  )

ggarrange(saur_plot, calb_plot,lps_plot,polyic_plot,nrow=2, ncol=2,common.legend = TRUE, legend = "bottom")
