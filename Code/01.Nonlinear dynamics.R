#### 0.BayesPrism ####

setwd("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/3.Single-cell_deconvolution/0.BayesPrism")

library(SingleCellExperiment)
suppressWarnings(library(BayesPrism))
library(InstaPrism)
library(biomaRt)
library(arrow)
library(dplyr)
library(stringr)
library(tidyr)

#### 谱系的BayesPrism反卷积 ####

# 获取单细胞数据（需要EMSEMBLE ID）
load("sc.dat.filtered.rdata") 
sc.stat <- plot.scRNA.outlier(input=sc.dat, cell.type.labels=cell.type.labels, species="hs", return.raw=TRUE)

# 获取转录组数据（需要EMSEMBLE ID）
load("bk.dat.rdata")
bk.stat <- plot.bulk.outlier(bulk.input=bk.dat, sc.input=sc.dat, cell.type.labels=cell.type.labels, species="hs", return.raw=TRUE)

# 绘制批量和单细胞数据对比图
plot.bulk.vs.sc(sc.input=sc.dat.filtered, bulk.input=bk.dat, pdf.prefix="bk.vs.sc")

# 选择相关性最高的组别
sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered, gene.type="protein_coding") #protein_coding=10320

# 构造一个Prism对象
myPrism <- new.prism(reference=sc.dat.filtered, mixture=bk.dat, input.type="count.matrix", 
                     cell.type.labels=cell.type.labels, cell.state.labels=cell.type.labels, 
                     key=NULL, outlier.cut=0.01, outlier.fraction=0.1)
save(myPrism, file="myPrism.rdata")

# (可选)使用instaprism进行加速
InstaPrism.res.initial <- InstaPrism(input_type='prism', prismObj=myPrism, n.core=16)
save(InstaPrism.res.initial, file="InstaPrism.res.initial.rdata")

load("InstaPrism.res.initial.rdata")
theta <- t(InstaPrism.res.initial@Post.ini.ct@theta)
write.csv(theta, file="InstaPrism_fraction_CellLineage.csv",row.names=T) # 获取细胞比例,每行加和为1

ct=unique(cell.type.labels)
for (i in ct){
  Z.ct <- reconstruct_Z_ct_initial(InstaPrism_obj=InstaPrism.res.initial, cell.type.of.interest=i) # 获取细胞的表达量
  write.csv(Z.ct, file=paste0('InstaPrism_exp_decon_CellLineage_',i,'.csv'), row.names=T)
}

# (可选)原生bayesprism慢速
bp.res <- run.prism(prism=myPrism, n.cores=50)
save(bp.res, file="bp.res.rdata")

theta <- get.fraction(bp=bp.res, which.theta="final", state.or.type="type") 
write.csv(theta, file="BayesPrism_fraction_CellLineage.csv",row.names=T) # 获取细胞比例,每行加和为1

ct=unique(cell.type.labels)
for (i in ct){
  Z.ct <- get.exp(bp=bp.res, state.or.type="type", cell.name=i)  # 获取细胞的表达量
  write.csv(Z.ct, file=paste0('BayesPrism_exp_decon_CellLineage_',i,'.csv'), row.names=T)
}


#### 反卷积准确性（与真实数据进行比较） ####
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)

# 计算单细胞数据的基因平均表达
genes <- read.csv("gene.csv")

compute_group_means <- function(sc.dat, cell.type.labels) {
  sc.dat <- as.data.frame(sc.dat) %>% 
    select(genes$gene) %>% 
    mutate(celltype=cell.type.labels, 
           group=paste0(str_split_fixed(rownames(.),"_",3)[,2],"_",str_split_fixed(rownames(.),"_",3)[,3])) %>%
    gather(key="gene", value="expression", -celltype, -group)
  result <- sc.dat %>% group_by(celltype, group, gene) %>%
    summarize(mean_expression=mean(expression, na.rm=TRUE))
  result_list <- list()
  for (i in unique(sc.dat$celltype)){
    result_list[[i]] <- result %>% filter(celltype==i) %>% pivot_wider(names_from=group, values_from=mean_expression)
  }
  return(result_list)
}

load("sc.dat.filtered.rdata")

load("CellLineage.labels.rdata")
CellLineage.labels <- cell.type.labels
load("CellType.labels.rdata")
CellType.labels <- cell.type.labels

group_means_result <- compute_group_means(sc.dat.filtered, CellLineage.labels)
save(group_means_result, file="scRNA_CellLineage.means.rdata")
group_means_result <- compute_group_means(sc.dat.filtered, CellType.labels)
save(group_means_result, file="scRNA_CellType.means.rdata")

load("sc.dat.filtered_SI.rdata")
load("CellType_SI.labels.rdata")
CellType.labels <- cell.type.labels
group_means_result <- compute_group_means(sc.dat.filtered, CellType.labels)
save(group_means_result, file="scRNA_CellType_SI.means.rdata")

load("sc.dat.filtered_LI.rdata")
load("CellType_LI.labels.rdata")
CellType.labels <- cell.type.labels
group_means_result <- compute_group_means(sc.dat.filtered, CellType.labels)
save(group_means_result, file="scRNA_CellType_LI.means.rdata")


# 相关性分析
sample <- c("JH_0_DU_2","JH_0_JE_2","JH_0_CE_2","JH_0_CO_2",
            "JH_60_DU_1","JH_60_JE_1","JH_60_IL_1","JH_60_CE_1",
            "JH_90_DU_3","JH_90_JE_3","JH_90_IL_3","JH_90_CE_3","JH_90_CO_3",
            "JH_180_DU_4","JH_180_JE_4","JH_180_IL_4","JH_180_CE_4","JH_180_CO_4",
            "JH_240_DU_2","JH_240_JE_2","JH_240_IL_2","JH_240_CE_2","JH_240_CO_2")
column_name <- c("DU_0","JE_0","CE_0","CO_0",
                 "DU_60","JE_60","IL_60","CE_60",
                 "DU_90","JE_90","IL_90","CE_90","CO_90",
                 "DU_180","JE_180","IL_180","CE_180","CO_180",
                 "DU_240","JE_240","IL_240","CE_240","CO_240")

df <- NULL
load("scRNA_CellLineage.means.rdata")
CellLineage.labels <- group_means_result
celltypes_1 <- c("T_ILC_NKcells","Bcells","Plasma","Epithelial","Mesenchymal","Myeloid","Neuronal","Endothelial")
for (c in celltypes_1){
  filename_1 <- paste0('InstaPrism_exp_decon_CellLineage_',c,'.csv')
  sc.ct <- CellLineage.labels[[c]] %>% as.data.frame(.) %>% column_to_rownames(var=names(.)[2]) %>% select(any_of(column_name))
  Z.ct <- read.csv(filename_1,row.names=1) %>% select(any_of(sample))
  Z.ct_column_name <- paste0(str_split_fixed(colnames(Z.ct),"_",4)[,3],"_",str_split_fixed(colnames(Z.ct),"_",4)[,2])
  Z.ct <- Z.ct %>% rename_with(~ Z.ct_column_name) %>% select(any_of(colnames(sc.ct)))
  Z.ct <- Z.ct[order(rownames(Z.ct)), ]
  
  cor_matrix <- cor(Z.ct, sc.ct, method="spearman")
  diag_data <- diag(cor_matrix)
  diag_df <- data.frame(celltype=c, sample=names(diag_data), correlation=diag_data, row.names=NULL)
  df <- rbind(df,diag_df)
}

load("scRNA_CellType.means.rdata")
CellType.labels <- group_means_result
celltypes_2 <- c("Colonocytes","Enterocytes","Goblet","BEST4enterocytes","TA","Stem","Progenitor","EECs","Tuft")
for (c in celltypes_2){
  filename_2 <- paste0('InstaPrism_exp_decon_CellType_',c,'.csv')
  Z.ct <- read.csv(filename_2,row.names=1) %>% select(any_of(sample))
  Z.ct_column_name <- paste0(str_split_fixed(colnames(Z.ct),"_",4)[,3],"_",str_split_fixed(colnames(Z.ct),"_",4)[,2])
  Z.ct <- Z.ct %>% rename_with(~ Z.ct_column_name) %>% select(any_of(colnames(sc.ct)))
  Z.ct <- Z.ct[order(rownames(Z.ct)), ]
  
  cor_matrix <- cor(Z.ct, sc.ct, method="spearman")
  diag_data <- diag(cor_matrix)
  diag_df <- data.frame(celltype=c, sample=names(diag_data), correlation=diag_data, row.names=NULL)
  df <- rbind(df,diag_df)
}

load("scRNA_CellType_SI.means.rdata")
CellType.labels <- group_means_result
c="Enterocytes"
filename_3 <- paste0('InstaPrism_exp_decon_CellType.SI_',c,'.csv')
Z.ct <- read.csv(filename_3,row.names=1) %>% select(any_of(sample)) %>% filter(rownames(.) %in% rownames(sc.ct))
Z.ct_column_name <- paste0(str_split_fixed(colnames(Z.ct),"_",4)[,3],"_",str_split_fixed(colnames(Z.ct),"_",4)[,2])
Z.ct <- Z.ct %>% rename_with(~ Z.ct_column_name) %>% select(any_of(colnames(sc.ct)))
Z.ct <- Z.ct[order(rownames(Z.ct)), ]
sc.ct.new <- sc.ct %>% select(any_of(colnames(Z.ct))) %>% filter(rownames(.) %in% rownames(Z.ct))
cor_matrix <- cor(Z.ct, sc.ct.new, method="spearman")
diag_data <- diag(cor_matrix)
diag_df <- data.frame(celltype=paste0(c,".SI"), sample=names(diag_data), correlation=diag_data, row.names=NULL)
df <- rbind(df,diag_df)

load("scRNA_CellType_LI.means.rdata")
CellType.labels <- group_means_result
c="Colonocytes"
filename_3 <- paste0('InstaPrism_exp_decon_CellType.LI_',c,'.csv')
Z.ct <- read.csv(filename_3,row.names=1) %>% select(any_of(sample)) %>% filter(rownames(.) %in% rownames(sc.ct))
Z.ct_column_name <- paste0(str_split_fixed(colnames(Z.ct),"_",4)[,3],"_",str_split_fixed(colnames(Z.ct),"_",4)[,2])
Z.ct <- Z.ct %>% rename_with(~ Z.ct_column_name) %>% select(any_of(colnames(sc.ct)))
Z.ct <- Z.ct[order(rownames(Z.ct)), ]
sc.ct.new <- sc.ct %>% select(any_of(colnames(Z.ct))) %>% filter(rownames(.) %in% rownames(Z.ct))
cor_matrix <- cor(Z.ct, sc.ct.new, method="spearman")
diag_data <- diag(cor_matrix)
diag_df <- data.frame(celltype=paste0(c,".LI"), sample=names(diag_data), correlation=diag_data, row.names=NULL)
df <- rbind(df,diag_df)

unique(df$celltype)
df$celltype <- factor(df$celltype,
                      levels=c("Bcells","Endothelial","Epithelial","Mesenchymal","Myeloid","Neuronal","Plasma","T_ILC_NKcells",
                               "BEST4enterocytes","Colonocytes","Colonocytes.LI","EECs","Enterocytes","Enterocytes.SI","Goblet","Progenitor","Stem","TA","Tuft"),
                      labels=c("B cells","Endothelial","Epithelial","Mesenchymal","Myeloid","Neuronal","Plasma","T lymphocytes",
                               "BEST4 enterocytes","Colonocytes","Colonocytes.LI","EECs","Enterocytes","Enterocytes.SI","Goblet","Progenitor","Stem","TA","Tuft"))
write.csv(df,"Spearman_corr_scRNA&bulk.csv")



#### 1.data_preparation ####

library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(impute)        # 用于KNN填补缺失值
library(limma)         # 用于批次效应校正(removeBatchEffect)
library(sva)           # 用于批次效应校正(ComBat)

#### 转录组数据 ####
setwd("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/1.Transcriptomics/1.data_preparation")

# 加载函数
plot_pca_batch_correction <- function(data, batch, combat_data, limma_data) {
  pca_data <- prcomp(t(data))
  pca_data_df <- as.data.frame(pca_data$x) %>% mutate(batch=batch)
  pca_data_plot <- ggplot(pca_data_df, aes(x=PC1,y=PC2,color=batch)) + geom_point() + theme_minimal() + labs(title="PCA plot before batch effect correction")
  
  pca_combat <- prcomp(t(combat_data))
  pca_combat_df <- as.data.frame(pca_combat$x) %>% mutate(batch = batch)
  pca_combat_plot <- ggplot(pca_combat_df, aes(x=PC1,y=PC2,color=batch)) + geom_point() + theme_minimal() + labs(title="PCA plot after ComBat batch effect correction")
  
  pca_limma <- prcomp(t(limma_data))
  pca_limma_df <- as.data.frame(pca_limma$x) %>% mutate(batch = batch)
  pca_limma_plot <- ggplot(pca_limma_df, aes(x=PC1,y=PC2,color=batch)) + geom_point() + theme_minimal() + labs(title="PCA plot after Limma batch effect correction")
  
  combined_plot <- pca_data_plot + pca_combat_plot + pca_limma_plot + plot_layout(ncol=3)
  return(combined_plot)
}

# 生成数据集
rna_data <- read.csv("/disk213/xieqq/JINHUA138/Transcriptome_analysis.2.RNAseq/gene_count_matrix.csv",row.names=1)
inf <- read.table("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/1.Transcriptomics/sample_name_rna.txt",header=T)
colnames(rna_data) <- sapply(colnames(rna_data), function(x) {
  new_name <- inf$sample_name[inf$old_name == x]
  if (length(new_name) > 0) {return(new_name)} else {return(x)}
})
rna_data_list <- list()
rna_data_list[["Duodenum"]] <- rna_data %>% select(contains("_DU_"))
rna_data_list[["Jejunum"]] <- rna_data %>% select(contains("_JE_"))
rna_data_list[["Ileum"]] <- rna_data %>% select(contains("_IL_"))
rna_data_list[["Cecum"]] <- rna_data %>% select(contains("_CE_"))
rna_data_list[["Colon"]] <- rna_data %>% select(contains("_CO_"))
save(rna_data_list, file="rna_data_list.rdata")

# 处理数据集
for (i in c("Duodenum","Jejunum","Ileum","Cecum","Colon")){
  rna_data_tissue <- rna_data_list[[i]]
  rna_data_dataset <- list()

  # 过滤掉低表达基因
  rna_data_filtered <- rna_data_tissue[rowMeans(rna_data_tissue) > 1, ]
  rna_data_dataset[["expression"]] <- rna_data_filtered
  
  # 混杂因素和批次效应
  inf <- data.frame(sample_name=colnames(rna_data_filtered)) %>% 
    left_join(.,read.table("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/1.Transcriptomics/sample_name_rna.txt",header=T),by="sample_name")
  rna_data_dataset[["inf"]] <- inf
  
  # 提取注释信息
  rna_annotation <- data.frame(ensenbleID=str_split_fixed(rownames(rna_data_filtered), "\\|", 2)[, 1], genename=str_split_fixed(rownames(rna_data_filtered), "\\|", 2)[, 2])
  rna_data_dataset[["annotation"]] <- rna_annotation
  
  # 评估批次效应并采用多种方法校正（最后选择ComBat）
  batch <- factor(inf$batch)
  combat_data <- sva::ComBat(dat=rna_data_filtered, batch=batch)         # 使用ComBat函数
  limma_data <- limma::removeBatchEffect(rna_data_filtered, batch=batch) # 使用removeBatchEffect函数
  combined_plot <- plot_pca_batch_correction(rna_data_filtered, batch, combat_data, limma_data)
  ggsave(paste0("pca_batch_effects_",tolower(i),".pdf"), combined_plot, width=16, height=5)
  
  # 输出数据
  rna_data_dataset[["combat_data"]] <- combat_data
  rna_data_dataset[["limma_data"]] <- limma_data
  save(rna_data_dataset, file=paste0("rna_data_dataset_",tolower(i),".rdata"))
}

for (i in c("Duodenum","Jejunum","Ileum","Cecum","Colon")){
  load(paste0("rna_data_dataset_",tolower(i),".rdata"))
  rna_data_filtered <- rna_data_dataset[["expression"]]
  inf <- data.frame(sample_name=colnames(rna_data_filtered)) %>% 
    left_join(.,read.table("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/1.Transcriptomics/sample_name_rna.txt",header=T),by="sample_name")
  rna_data_dataset[["inf"]] <- inf
  save(rna_data_dataset, file=paste0("rna_data_dataset_",tolower(i),".rdata"))
}


#### 2.Linear regression changing ####

#### linear changing ####
library(limma)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(plyr)
library(tibble)
library(lme4)
source("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/0.code-tools.R")

setwd("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/1.Transcriptomics/2.linear_regression/1.linear_changing")

for (tissue in c("duodenum","jejunum","ileum","cecum","colon")){
  # load data
  load(paste0("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/1.Transcriptomics/1.data_preparation/rna_data_dataset_",tissue,".rdata"))
  object_cross_section <- rna_data_dataset
  expression_data <- object_cross_section[["combat_data"]] %>% 
    apply(1, function(x) {(x-mean(x))/sd(x)}) %>% t() %>% as.data.frame()  #z-score标准化
  sample_info <- object_cross_section[["inf"]]
  
  # adjust sex, breed and diet ethnicity
  expression_data <- lm_adjust(expression_data=expression_data, sample_info=sample_info, threads=3)
  temp_object <- object_cross_section
  temp_object[["expression"]] <- expression_data
  object_cross_section <- temp_object
  save(object_cross_section, file=paste0("object_cross_section_",tissue,".rdata"))
  
  # Spearman correlation
  cor_data <- seq_len(nrow(expression_data)) %>%
    purrr::map(function(i) {
      value <- as.numeric(expression_data[i, , drop=TRUE])
      cor_result <- cor.test(value, temp_object[["inf"]]$day, method="spearman")
      data.frame(variable_id=rownames(expression_data)[i], cor_p=cor_result$p.value, spearman_cor=cor_result$estimate)
    }) %>%
    dplyr::bind_rows() %>%
    as.data.frame() %>% 
    mutate(cor_p_adjust=p.adjust(cor_p, method="fdr")) %>% 
    {rownames(.) <- NULL; .}
  write.csv(cor_data, file=paste0("correlation_with_age_",tissue,".csv"))
  
  # permutation to get the p value
  permutation_cor_list <- list()
  for (idx in 1:100) {
    #cat(idx, " ")
    permutation_cor_data <- seq_len(nrow(expression_data)) %>%
      purrr::map(function(i) {
        value <- as.numeric(expression_data[i, , drop=TRUE])
        cor_result <- cor.test(value, sample(temp_object[["inf"]]$day, replace=FALSE), method="spearman") # 随机打乱年龄标签
        data.frame(variable_id=rownames(expression_data)[i], cor_p=cor_result$p.value, spearman_cor=cor_result$estimate)
      }) %>%
      dplyr::bind_rows() %>% 
      as.data.frame()
    permutation_cor_list[[as.character(idx)]] <- permutation_cor_data
  }
  save(permutation_cor_list, file=paste0("permutation_cor_with_age_",tissue,".rdata"))
  
  permutaton_cor_data_all <- permutation_cor_list[[1]] %>% select(-cor_p)
  for (idx in 2:100) {
    permutation_cor_data <- permutation_cor_list[[idx]]
    permutaton_cor_data_all <- cbind(permutaton_cor_data_all, permutation_cor_data[, 3, drop=FALSE])
  }
  rownames(permutaton_cor_data_all) <- NULL
  permutaton_cor_data_all <- permutaton_cor_data_all %>% as.data.frame() %>% column_to_rownames(var="variable_id")
  
  permutated_p_value <- 1:nrow(cor_data) %>%
    purrr::map(function(idx) {
      #cat(idx, " ")
      original_cor <- cor_data$spearman_cor[idx]
      permutation_cor <- sample(as.numeric(permutaton_cor_data_all[idx,]), 10000, replace=TRUE)
      if (original_cor > 0) {sum(permutation_cor > original_cor) / 10000} else{sum(permutation_cor < original_cor) / 10000}
    }) %>%
    unlist() %>%
    as.data.frame() %>% 
    mutate(permutated_p_adjust=p.adjust(permutated_p_value, method="fdr")) %>% 
    {rownames(.) <- NULL; .}
  perm_df <- permutated_p_value
  cor_data <- cor_data %>%
    mutate(permutated_p_value=perm_df$permutated_p_value, 
           permutated_p_adjust=perm_df$permutated_p_adjust)
  write.csv(cor_data, file=paste0("permutation_cor_with_age_",tissue,".csv"))
}

##linear mixed model
for (tissue in c("duodenum","jejunum","ileum","cecum","colon")){
  load(paste0("object_cross_section_",tissue,".rdata"))
  expression_data <- object_cross_section[["expression"]]
  sample_info <- object_cross_section[["inf"]]
  linear_mix_model_data <- seq_len(nrow(expression_data)) %>% 
    purrr::map(function(i) {
      #cat(i, " ")
      value <- as.numeric(expression_data[i, , drop=TRUE])
      temp_data <- data.frame(sample_info, value) %>% mutate(day=as.numeric(as.character(day)))
      lm_result <- glm(formula=value ~ day, data=temp_data)
      lm_result <- lm_result %>% broom::tidy()
      data.frame(lm_p=lm_result$p.value[2], coefficient=lm_result$estimate[2])
    }) %>%
    dplyr::bind_rows() %>%
    as.data.frame()
  save(linear_mix_model_data, file=paste0("linear_mix_model_with_age_",tissue,".rdata"))
  temp_object <- data.frame(variable_id=rownames(expression_data)) %>%
    dplyr::mutate(lm_p=linear_mix_model_data$lm_p, 
                  coefficient=linear_mix_model_data$coefficient,
                  lm_p_adjust=p.adjust(lm_p, method="fdr"))
  write.csv(temp_object, file=paste0("linear_mix_model_with_age_",tissue,".csv"))
}

#### dysregulated molecules during aging ####
library(dplyr)
library(plyr)
library(pbapply)

setwd("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/1.Transcriptomics/2.linear_regression/2.dysregulated_molecules")

for (tissue in c("duodenum","jejunum","ileum","cecum","colon")){
  # load data
  load(paste0("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/1.Transcriptomics/2.linear_regression/1.linear_changing/object_cross_section_",tissue,".rdata"))
  subject_data <- object_cross_section[["expression"]]
  subject_data2 <- subject_data %>% t() %>% as.data.frame() %>%
    mutate(age_range=object_cross_section[["inf"]]$day) %>% 
    plyr::dlply(.variables=.(age_range)) %>%
    lapply(., function(x) {x <- x %>% select(-age_range)})
  
  # two.sided t.test
  fc_p_value <- pbapply::pblapply(subject_data2[-1], function(x) {
      p_value <- lapply(1:ncol(x), function(idx) {t.test(x[, idx], subject_data2[[1]][, idx], paired=FALSE)$p.value}) %>%
        unlist() %>%
        p.adjust(method="fdr")
      data.frame(p_value, variable_id=rownames(subject_data),stringsAsFactors=FALSE)
    })
  save(fc_p_value, file=paste0("fc_p_value_",tissue,".rdata"))  #find all the peaks in different time points
  
  marker_each_point <- lapply(fc_p_value, function(x) {
      idx1 <- which(x$p_value < 0.05)
      gene1 <- try(data.frame(x[idx1,], stringsAsFactors=FALSE), silent=TRUE)
      if (class(gene1) == "try-error") {gene1 <- NULL}
      gene1
    })
  save(marker_each_point, file=paste0("marker_each_point_",tissue,".rdata")) # find markers for each time points
  
  # permutation to get the p value
  original_difference <- 2:length(subject_data2) %>%
    purrr::map(function(idx) {
      x <- subject_data2[[idx]]
      difference <- 1:ncol(x) %>% purrr::map(function(i) {abs(mean(subject_data2[[1]][, i]) -mean(x[, i]))}) %>% unlist()}) %>%
    do.call(cbind, .) %>%
    as.data.frame()
  rownames(original_difference) <- colnames(subject_data2[[1]])
  colnames(original_difference) <- names(subject_data2)[-1]
  write.csv(original_difference, file=paste0("original_difference_",tissue,".csv"), row.names=T)
  
  permutated_difference_list <- list()
  for (idx in 2:length(subject_data2)) {
    #cat(idx, " ")
    x <- subject_data2[[idx]]
    permutated_difference <- purrr::map(1:100, function(j) {
        #cat(j, " ")
        difference <- 1:ncol(x) %>%
          purrr::map(function(i) {
            value1 <- subject_data[i, sample(1:ncol(subject_data),length(subject_data2[[1]][, i]), replace=FALSE,), drop=TRUE] %>% unlist() # 随机打乱
            value2 <- subject_data[i, sample(1:ncol(subject_data),length(subject_data2[[1]][, i]), replace=FALSE,), drop=TRUE] %>% unlist() # 随机打乱
            abs(mean(value1)-mean(value2))
          }) %>%
          unlist()
        difference
      }) %>%
      do.call(cbind, .) %>%
      as.data.frame()
    permutated_difference_list[[as.character(idx)]] <- permutated_difference
  }
  save(permutated_difference_list, file=paste0("permutated_difference_",tissue,".rdata"))
  
  fc_p_value_permutation <- 2:length(subject_data2) %>%
    purrr::map(function(idx) {
      permutated_difference <- permutated_difference_list[[as.character(idx)]]
      p_value <- purrr::map(1:nrow(original_difference), function(i) {
        original_value <- original_difference[i, idx - 1]
        permutation_value <- permutated_difference[i,] %>% as.numeric()
        permutation_value <- sample(permutation_value, 10000, replace=TRUE)
        sum(permutation_value > original_value) / 10000
        }) %>% 
        unlist()
      data.frame(p_value=p_value, variable_id=rownames(original_difference))
    })
  names(fc_p_value_permutation) <- names(original_difference)
  save(fc_p_value_permutation, file=paste0("permutation_fc_p_value_",tissue,".rdata"))
  
  marker_each_point_permutation <- lapply(fc_p_value_permutation, function(x) {
      idx1 <- which(p.adjust(x$p_value, method="fdr") < 0.05)
      gene1 <- try(data.frame(x[idx1, ], stringsAsFactors=FALSE), silent=TRUE)
      if (class(gene1) == "try-error") {gene1 <- NULL}
      gene1
    })
  save(marker_each_point_permutation, file=paste0("permutation_marker_each_point_",tissue,".rdata"))
}

#### Evaluation of the age reflected by omics data ####
library(dplyr)
library(plyr)
library(purrr)
library(tidyverse)
library(mixOmics) 

setwd("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/1.Transcriptomics/2.linear_regression/3.age_evaluation")

for (tissue in c("duodenum","jejunum","ileum","cecum","colon")){
  # load data
  load(paste0("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/1.Transcriptomics/2.linear_regression/1.linear_changing/object_cross_section_",tissue,".rdata"))
  evaluation_data <- object_cross_section[["expression"]]
  
  # Spearman correlation between PC1 and age
  pca_result <- prcomp(t(evaluation_data))
  cor_data <- seq_len(nrow(evaluation_data)) %>%
    purrr::map(function(i) {
      value <- as.numeric(evaluation_data[i, , drop=TRUE])
      cor_result <- cor.test(value, pca_result$x[,1], method="spearman")
      data.frame(variable_id=rownames(evaluation_data)[i], cor_p=cor_result$p.value, spearman_cor=cor_result$estimate)
    }) %>%
    dplyr::bind_rows() %>%
    as.data.frame() %>% 
    {rownames(.) <- NULL; .}
  write.csv(cor_data, file=paste0("correlation_between_PC1_and_age_",tissue,".csv"))
  
  # PLS
  temp_data <- evaluation_data %>% t(.) %>% as.data.frame()
  pls_result <- pls(X=temp_data, Y=object_cross_section[["inf"]]$day, ncomp=10) # assuming 10 components, but adjust as needed
  summary(pls_result)
  pls_performance <- perf(pls_result, validation="Mfold", folds=5) 
  save(pls_performance, file=paste0("pls_performance_",tissue,".rdata")) #pls_performance$measures$R2
}

pls_performance$measures$Q2$summary
cumsum(pls_performance$measures$R2$summary$mean)


#### 3.modified_DEswan ####
rm(list=ls())

#library(DEswan)
library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)
setwd("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/1.Transcriptomics/3.modified_DEswan")

# 定义改良DE-SWAN函数 (离散节点替换滑窗；Wilcoxon检验替换线性回归；随机置换检验)
DE_node_analysis <- function(data, nodes, permutate=FALSE) {
  
  if(permutate) {data$age <- sample(data$age)} # 随机打乱年龄标签
  
  node_results <- list()
  features <- colnames(data)[-1]
  
  for(i in 1:(length(nodes)-1)) {
    low_age <- nodes[i]
    high_age <- nodes[i+1]
    group_low <- data[data$age == low_age, -1]
    group_high <- data[data$age == high_age, -1]
    
    # 对每个 feature 进行 Wilcoxon 检验
    res <- sapply(features, function(f) {
      val_low <- group_low[[f]]
      val_high <- group_high[[f]]
      wt <- wilcox.test(val_high, val_low, exact=FALSE)
      
      #FDR校正
      p_value_adj <- p.adjust(wt$p.value, "fdr")
      
      # 计算中位数差异作为效应量 (Effect Size)
      med_diff <- median(val_high) - median(val_low)
      return(c(Median_Diff=med_diff, p_value=wt$p.value, p_value_adjust=p_value_adj))
    })
    node_results[[paste0(low_age, "_to_", high_age)]] <- t(res)
  }
  return(node_results)
}

for (tissue in c("duodenum","jejunum","ileum","cecum","colon")){
  # load data
  load(paste0("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/1.Transcriptomics/2.linear_regression/1.linear_changing/object_cross_section_",tissue,".rdata"))
  expression_data <- object_cross_section[["expression"]]
  sample_info <- object_cross_section[["inf"]]
  
  # 改良DE-SWAN
  age_nodes <- c(0, 60, 90, 180, 240)
  temp_data <- data.frame(age=sample_info$day, t(expression_data)) %>% {colnames(.) <- gsub("\\.", "|", colnames(.)); .}
  actual_results <- DE_node_analysis(temp_data, age_nodes, permutate=FALSE)
  save(actual_results, file=paste0("actual_results_",tissue,".rdata"))
  
  # 执行置换检验
  n_permutations <- 100
  permutation_results_list <- lapply(1:n_permutations, function(i) {
    perm_res <- DE_node_analysis(temp_data, age_nodes, permutate=TRUE)
    return(perm_res)
  })
  save(permutation_results_list, file=paste0("permutation_results_list_",tissue,".rdata"))
}

validation <- NULL
n_permutations <- 100
p_value_adjust_threshold <- 0.0001
for (tissue in c("duodenum","jejunum","ileum","cecum","colon")){
  load(paste0("actual_results_",tissue,".rdata"))
  load(paste0("permutation_results_list_",tissue,".rdata"))
  
  actual_results_count <- sapply(actual_results, function(x) sum(x[, "p_value_adjust"] < p_value_adjust_threshold))
  perm_results_count <- sapply(permutation_results_list, function(perm_res) {
    sapply(perm_res, function(x) sum(x[, "p_value_adjust"] < p_value_adjust_threshold))
  })
  perm_mean <- rowMeans(perm_results_count)
  perm_sd <- apply(perm_results_count, 1, sd)
  empirical_num <- sapply(1:length(actual_results_count), function(i) {
    num_greater_equal <- sum(perm_results_count[i, ] >= actual_results_count[i]) # 逻辑：统计随机结果中，有多少次产生的显著特征数量 >= 实际观测到的数量
    return(num_greater_equal)
  })

  validation_table <- data.frame(
    Omics="1.Transcriptomics",
    Tissue=tissue,
    Interval=names(actual_results_count),
    Actual_Count=actual_results_count,
    Greater_Count=empirical_num,
    Perm_Mean=perm_mean,
    Perm_Sd=perm_sd,
    Enrich_Fold=actual_results_count / (perm_mean + 0.1), # 加0.1防止除以0
    Empirical_P=empirical_num / n_permutations
  )
  validation <- rbind(validation,validation_table)
}
write.csv(validation,paste0("validation_padj<",p_value_adjust_threshold,".csv"))







