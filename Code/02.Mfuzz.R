no_function()
masstools::setwd_project()
rm(list = ls())
library(stringr)
library(tidyverse)
source("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/0.code-tools.R")


setwd("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/5.all_omics/1.data_combination_1")

#### 1.data combination ####

process_omics_data <- function(omics_type, tissue_types, file_path) {
  omics_expression_data <- NULL
  
  for (t in tissue_types) {
    # Load the data for the current tissue
    file_path_prefix <- paste0("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/",file_path,"/2.linear_regression/1.linear_changing/object_cross_section_")
    load(paste0(file_path_prefix, t, ".rdata"))
    inf <- object_cross_section[["inf"]]
    
    if (file_path=="3.Single-cell_deconvolution"){t=str_split_fixed(t, "_", 2)[,2]
    }else if (file_path=="4.Single-cell_transcriptomics"){t=str_split_fixed(t, "_", 2)[,2]
    }else {t=str_split_fixed(t, "_", 2)[,1]}
    
    # Normalize expression data
    expression_data <- apply(object_cross_section[["expression"]], 2, function(x) (x - min(x)) / (max(x) - min(x)) * 10)  # Normalize range
    if (file_path=="4.Single-cell_transcriptomics") {
      colnames(expression_data) <- paste("JH", "diet1", inf$day, str_split_fixed(inf$sample_name, "_", 3)[, 3], sep = "_")
    } else {colnames(expression_data) <- paste(inf$breed, inf$diet, inf$day, str_split_fixed(inf$old_name, "_", 3)[, 3], sep = "_")}
    rownames(expression_data) <- paste(omics_type, t, rownames(expression_data), sep = ":")
    
    # Group by day and compute mean expression data for each group
    day_ages <- gsub(".*_(\\d+)_.*", "\\1", colnames(expression_data))
    expression_data_list <- lapply(unique(day_ages), function(dayage) {
      cols_for_dayage <- colnames(expression_data)[day_ages == dayage]
      expression_data[, cols_for_dayage, drop=FALSE]
    })
    names(expression_data_list) <- unique(day_ages)
    
    # Compute mean for each day group
    if (file_path=="4.Single-cell_transcriptomics"){
      mean_expression_data_list <- lapply(expression_data_list, function(x) apply(x, 1, function(row_data) mean(row_data, trim=0.1, na.rm=T)))
    }else {mean_expression_data_list <- lapply(expression_data_list, function(x) apply(x, 1, mean))}
    mean_data_matrix <- do.call(cbind, mean_expression_data_list) %>% as.data.frame()
    
    # Ensure all required day columns are present (0, 60, 90, 180, 240)
    missing_columns <- setdiff(c("0", "60", "90", "180", "240"), colnames(mean_data_matrix))
    if (length(missing_columns) > 0) {
      for (col in missing_columns) {
        mean_data_matrix[[col]] <- NA
        print(t)
      }
    }
    
    # Reorder columns
    mean_data_matrix <- mean_data_matrix[, c("0", "60", "90", "180", "240")]
    print(mean_data_matrix %>% filter(grepl("ENSSSCG00000007507", rownames(.))))
    
    # Combine the data across tissues
    omics_expression_data <- rbind(omics_expression_data, mean_data_matrix)
  }
  
  return(omics_expression_data)
}

#### Transcriptomics Data
tissues <- c("duodenum", "jejunum", "ileum", "cecum", "colon")
transcriptomics_expression_data <- process_omics_data("RNA_Muscle", tissues, "1.Transcriptomics")

#### Untargeted Metabolomics Data
tissues <- c("duodenum", "jejunum", "ileum", "cecum", "colon")
metabolomics_expression_data <- process_omics_data("Meta_Muscle", tissues, "2.Untargeted_metabolomics")
metabolomics_expression_data <- metabolomics_expression_data %>% filter(if_all(everything(), ~ !is.na(.)))

#### scDeconv Data
tissues <- c("cecum", "colon")
scDeconv.Colonocytes_expression_data <- process_omics_data("scDeconv_Colonocytes", paste0("Colonocytes_",tissues), "3.Single-cell_deconvolution")

tissues <- c("duodenum", "jejunum", "ileum")
scDeconv.Enterocytes_expression_data <- process_omics_data("scDeconv_Enterocytes", paste0("Enterocytes_",tissues), "3.Single-cell_deconvolution")

tissues <- c("duodenum", "jejunum", "ileum", "cecum", "colon")
scDeconv.Tlymphocytes_expression_data <- process_omics_data("scDeconv_Tlymphocytes", paste0("Tlymphocytes_",tissues), "3.Single-cell_deconvolution")
scDeconv.Bcells_expression_data <- process_omics_data("scDeconv_Bcells", paste0("Bcells_",tissues), "3.Single-cell_deconvolution")
scDeconv.Plasma_expression_data <- process_omics_data("scDeconv_Plasma", paste0("Plasma_",tissues), "3.Single-cell_deconvolution")
scDeconv.Epithelial_expression_data <- process_omics_data("scDeconv_Epithelial", paste0("Epithelial_",tissues), "3.Single-cell_deconvolution")
scDeconv.Mesenchymal_expression_data <- process_omics_data("scDeconv_Mesenchymal", paste0("Mesenchymal_",tissues), "3.Single-cell_deconvolution")
scDeconv.Myeloid_expression_data <- process_omics_data("scDeconv_Myeloid", paste0("Myeloid_",tissues), "3.Single-cell_deconvolution")
scDeconv.Neuronal_expression_data <- process_omics_data("scDeconv_Neuronal", paste0("Neuronal_",tissues), "3.Single-cell_deconvolution")
scDeconv.Endothelial_expression_data <- process_omics_data("scDeconv_Endothelial", paste0("Endothelial_",tissues), "3.Single-cell_deconvolution")
scDeconv.Goblet_expression_data <- process_omics_data("scDeconv_Goblet", paste0("Goblet_",tissues), "3.Single-cell_deconvolution")
scDeconv.BEST4enterocytes_expression_data <- process_omics_data("scDeconv_BEST4enterocytes", paste0("BEST4enterocytes_",tissues), "3.Single-cell_deconvolution")
scDeconv.TA_expression_data <- process_omics_data("scDeconv_TA", paste0("TA_",tissues), "3.Single-cell_deconvolution")
scDeconv.Stem_expression_data <- process_omics_data("scDeconv_Stem", paste0("Stem_",tissues), "3.Single-cell_deconvolution")
scDeconv.Progenitor_expression_data <- process_omics_data("scDeconv_Progenitor", paste0("Progenitor_",tissues), "3.Single-cell_deconvolution")
scDeconv.EECs_expression_data <- process_omics_data("scDeconv_EECs", paste0("EECs_",tissues), "3.Single-cell_deconvolution")
scDeconv.Tuft_expression_data <- process_omics_data("scDeconv_Tuft", paste0("Tuft_",tissues), "3.Single-cell_deconvolution")

#### scRNA Data
tissues <- c("cecum", "colon")
scRNA.Colonocytes_expression_data <- process_omics_data("scRNA_Colonocytes", paste0("Colonocytes_",tissues), "4.Single-cell_transcriptomics")

tissues <- c("duodenum", "jejunum", "ileum")
scRNA.Enterocytes_expression_data <- process_omics_data("scRNA_Enterocytes", paste0("Enterocytes_",tissues), "4.Single-cell_transcriptomics")
scRNA.Myeloid_expression_data <- process_omics_data("scRNA_Myeloid", paste0("Myeloid_",tissues), "4.Single-cell_transcriptomics")

tissues <- c("duodenum", "jejunum", "ileum", "cecum", "colon")
scRNA.Tlymphocytes_expression_data <- process_omics_data("scRNA_Tlymphocytes", paste0("Tlymphocytes_",tissues), "4.Single-cell_transcriptomics")
scRNA.Epithelial_expression_data <- process_omics_data("scRNA_Epithelial", paste0("Epithelial_",tissues), "4.Single-cell_transcriptomics")

tissues <- c("duodenum", "colon")
scRNA.Goblet_expression_data <- process_omics_data("scRNA_Goblet", paste0("Goblet_",tissues), "4.Single-cell_transcriptomics")

tissues <- c("jejunum", "colon")
scRNA.BEST4enterocytes_expression_data <- process_omics_data("scRNA_BEST4enterocytes", paste0("BEST4enterocytes_",tissues), "4.Single-cell_transcriptomics")

tissues <- c("ileum", "cecum")
scRNA.Bcells_expression_data <- process_omics_data("scRNA_Bcells", paste0("Bcells_",tissues), "4.Single-cell_transcriptomics")

tissues <- c("ileum")
scRNA.Plasma_expression_data <- process_omics_data("scRNA_Plasma", paste0("Plasma_",tissues), "4.Single-cell_transcriptomics")

#### combine all omics together ####
temp_data_mean_list <- list(
  Transcriptomics          = transcriptomics_expression_data,
  unMetabolomics           = metabolomics_expression_data,
  scDeconv.Tlymphocytes    = scDeconv.Tlymphocytes_expression_data,
  scDeconv.Bcells          = scDeconv.Bcells_expression_data,
  scDeconv.Plasma          = scDeconv.Plasma_expression_data,
  scDeconv.Epithelial      = scDeconv.Epithelial_expression_data,
  scDeconv.Mesenchymal     = scDeconv.Mesenchymal_expression_data,
  scDeconv.Myeloid         = scDeconv.Myeloid_expression_data,
  scDeconv.Neuronal        = scDeconv.Neuronal_expression_data,
  scDeconv.Endothelial     = scDeconv.Endothelial_expression_data,
  scDeconv.Colonocytes     = scDeconv.Colonocytes_expression_data,
  scDeconv.Enterocytes     = scDeconv.Enterocytes_expression_data,
  scDeconv.Goblet          = scDeconv.Goblet_expression_data,
  scDeconv.BEST4enterocytes= scDeconv.BEST4enterocytes_expression_data,
  scDeconv.TA              = scDeconv.TA_expression_data,
  scDeconv.Stem            = scDeconv.Stem_expression_data,
  scDeconv.Progenitor      = scDeconv.Progenitor_expression_data,
  scDeconv.EECs            = scDeconv.EECs_expression_data,
  scDeconv.Tuft            = scDeconv.Tuft_expression_data,
  scRNA.Tlymphocytes       = scRNA.Tlymphocytes_expression_data,
  scRNA.Bcells             = scRNA.Bcells_expression_data,
  scRNA.Plasma             = scRNA.Plasma_expression_data,
  scRNA.Epithelial         = scRNA.Epithelial_expression_data,
  scRNA.Myeloid            = scRNA.Myeloid_expression_data,
  scRNA.Colonocytes        = scRNA.Colonocytes_expression_data,
  scRNA.Enterocytes        = scRNA.Enterocytes_expression_data,
  scRNA.Goblet             = scRNA.Goblet_expression_data,
  scRNA.BEST4enterocytes   = scRNA.BEST4enterocytes_expression_data
)
save(temp_data_mean_list,file="temp_data_mean_list.rdata")

na_summary <- sapply(temp_data_mean_list, function(df) {sum(rowSums(is.na(df)) > 0)})
na_summary

temp_data_mean <- do.call(rbind, temp_data_mean_list) %>% as.data.frame() 
save(temp_data_mean,file="temp_data_mean.rdata")


#### 2.k-means clustering ####
setwd("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/5.all_omics/2.k-means_clustering_1")

library(tibble)
library(dplyr)
library(tidyr)
library(data.table)
library(purrr)
library(stringr)
library(Mfuzz)
library(corrplot)
library(ComplexHeatmap)
library(ggplot2)
library(ggsci)
library(cowplot)

source("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/0.code-tools.R")

#### load data
load("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/5.all_omics/1.data_combination/temp_data_mean.rdata")
temp_data_mean <- as.matrix(temp_data_mean)

#### building object
data <- new('ExpressionSet',exprs=temp_data_mean)
data.s <- standardise(data)
save(data.s,file="data.s.rdata")
#load("data.s.rdata")
m1 <- mestimate(data.s)

#### clustering
# distance_k_number
Min_centroid_distance <- Dmin(data.s, m=m1, crange=seq(2, 30, 2), repeats=3, visu=TRUE)
Min_centroid_distance <- data.frame(distance=Min_centroid_distance, k=seq(2, 30, 2))
write.csv(Min_centroid_distance,"Min_centroid_distance.csv",row.names=F)

cluster_number <- 12 #k
c <- mfuzz(data.s, c=cluster_number, m=m1) 
save(c,file="c_mfuzz.rdata")

#### only molecules with memberships above 0.4 were retained within each cluster for further analysis.
#load("c_mfuzz.rdata")
membership_cutoff <- 0.4
center <- get_mfuzz_center(data=data.s, c=c, membership_cutoff=membership_cutoff)
rownames(center) <- paste("Cluster", rownames(center), sep=' ')
center_cor <- cor(t(center),method="spearman") %>% as.data.frame()
write.csv(center_cor,file="center_cor.csv")
corrplot(center_cor, type="full", diag=TRUE, order="hclust", hclust.method="ward.D", number.cex=.7, addCoef.col="black",
         col=colorRampPalette(colors=rev(RColorBrewer::brewer.pal(n=11, name="Spectral")))(n=100))

#### plot for each cluster
cluster_info <- data.frame(variable_id=names(c$cluster),c$membership,cluster=c$cluster,stringsAsFactors=FALSE) %>%
  as.data.frame() %>% {rownames(.) <- NULL; .} %>% 
  arrange(cluster)
table(cluster_info$cluster)
write.csv(cluster_info,file="cluster_info.csv")

temp_data <- exprs(data.s) %>% as.data.frame()