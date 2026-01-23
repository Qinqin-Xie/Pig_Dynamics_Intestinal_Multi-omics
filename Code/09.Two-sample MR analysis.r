setwd("/disk213/xieqq/JINHUA138/Genome_analysis.5.MR_new")

####exposure####
library(data.table)
library(R.utils)
library(dplyr)
library(tidyr)

#eqtl.database <- fread(file="/disk191_2/yupf/GUT_RNA/04_SNPCALLING/eQTL/qtl/result/smallint/smallint.cis_qtl_pairs.significant.txt.gz",header=T)
#ee <- eqtl.database[which(eqtl.database$phenotype_id=="ENSSSCG00000007507" & eqtl.database$pval_nominal<1e-5),]

#eqtl.database <- fread(file="./PigGTEx_v0.eQTL_full_summary/Small_intestine.cis_qtl_pairs.1_18.txt",header=T)
#ee <- eqtl.database %>% filter(grepl("ENSSSCG00000007507", phenotype_id), pval_nominal<1e-5)

eqtl.database <- fread(file="./PigGTEx_v0.eeQTL_full_summary/Small_intestine.cis_qtl_pairs.1_18.txt",header=T)
ee <- eqtl.database %>% filter(grepl("ENSSSCG00000007507", phenotype_id), pval_nominal<1e-5)
#pval_adj_BH = p.adjust(eqtl.database$pval_beta, method = "BH")
ee$phenotype_id <- sapply(strsplit(ee$phenotype_id, ":"), function(x) x[1])
ee$effect_allele <- sapply(strsplit(ee$variant_id, "_"), function(x) x[3])
ee$non_effect_allele <- sapply(strsplit(ee$variant_id, "_"), function(x) x[4])
ee <- ee %>% select("phenotype_id","variant_id","effect_allele","non_effect_allele","tss_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
write.csv(ee,"exposure.csv",row.names=F)



####GWAS####
exp <- read.csv("exposure.csv",check.names=F,header=T)
exp$panel_variant_id <- paste0(sapply(strsplit(unique(exp$variant_id), "_"), function(x) x[1]),"_",sapply(strsplit(unique(exp$variant_id), "_"), function(x) x[2]))
SNP <- paste0("chr",sapply(strsplit(unique(exp$variant_id), "_"), function(x) x[1]),":", sapply(strsplit(unique(exp$variant_id), "_"), function(x) x[2]))
colname_set <- c("chromosome","rs","position","n_miss","effect_allele","non_effect_allele","frequency","effect_size",
                 "standard_error","logl_H1","l_remle","pvalue","panel_variant_id","variant_id")

file_list <- list.files(path="/disk194/gujm/P2P/gemma", pattern="^Meta_", full.names=TRUE)
for (file in file_list) {
  data <- read.table(file, check.names=F, header=T)
  data <- subset(data,rs %in% SNP)
  data$panel_variant_id <- paste0(data$chr,"_",sapply(strsplit(unique(data$rs), ":"), function(x) x[2]))
  if (nrow(data)>0){
    data <- left_join(data,exp[,c("variant_id","panel_variant_id")],by="panel_variant_id")
    colnames(data) <- colname_set
    data <- data %>% select("variant_id","panel_variant_id","chromosome","position","effect_allele","non_effect_allele",
                            "frequency","pvalue","effect_size","standard_error")
    write.table(data, paste0("./output/",gsub("/disk194/gujm/P2P/gemma/","",file)), quote=F, row.names=F)
  }
}


####MR####

library(tidyverse)
library(TwoSampleMR)
library(grid)
library(forestploter)
library(ggplot2)

exp <- read.csv("exposure.csv",check.names=F,header=T)
exp$id <- "PCK1 expression"
exp_dat <- format_data(exp, type="exposure", snp_col="variant_id", beta_col="slope", se_col="slope_se", 
                       effect_allele_col="effect_allele", other_allele_col="non_effect_allele",
                       eaf_col="af", pval_col="pval_nominal", samplesize_col="ma_samples",id_col="id")

file_list <- list.files(path="./output", pattern=".txt", full.names=TRUE)
#for (file in file_list) {
#  data <- read.table(file, check.names=F, header=F)
#  colnames(data) <- c("variant_id","panel_variant_id","chromosome","position","effect_allele","non_effect_allele",
#                      "frequency","pvalue","zscore","effect_size","standard_error","sample_size")
#  write.table(data,file,quote=F,row.names=F)
#}

## step 1
file_list <- list.files(path="./output", pattern="assoc.txt", full.names=TRUE)
out <- NULL
for (file in file_list) {
  trait <- gsub("./output/Meta_|.assoc.txt", "", file)
  data <- read.table(file, check.names=F, header=T)
  data$id <- trait
  t2d_out <- format_data(data, type="outcome", snp_col="variant_id", beta_col="effect_size", se_col="standard_error", 
                         effect_allele_col="effect_allele", other_allele_col="non_effect_allele",
                         eaf_col="frequency", pval_col="pvalue", samplesize_col="sample_size", id_col="id")
  mydata <- harmonise_data(exposure_dat=exp_dat, outcome_dat=t2d_out, action=2)
  heterogeneity <- mr_heterogeneity(mydata,method_list=c("mr_two_sample_ml","mr_ivw","mr_ivw_radial","mr_ivw_radial"))
  data.test <- data.frame(trait=trait,method=heterogeneity$method,heterogeneity_pval=heterogeneity$Q_pval)
  out <- rbind(out,data.test)
}
write.csv(out,"MR.method.heterogeneity.test.csv")

## step 2
file_list <- list.files(path="./output", pattern=".txt", full.names=TRUE)
out <- NULL
for (file in file_list) {
  trait <- gsub("./output/Meta_|.assoc.txt", "", file)
  #mr_method_list: https://blog.csdn.net/qq994327432/article/details/140248005
  method <- c("mr_egger_regression_bootstrap","mr_penalised_weighted_median","mr_weighted_mode_nome","mr_simple_mode_nome","mr_sign","mr_uwr",
              "mr_two_sample_ml",                                   
              "mr_simple_median","mr_weighted_median",             
              "mr_egger_regression",                               
              "mr_simple_mode",                                    
              "mr_weighted_mode",                                   
              "mr_ivw","mr_ivw_radial","mr_ivw_mre","mr_ivw_fe")    
  data <- read.table(file, check.names=F, header=T)
  data$id <- trait
  t2d_out <- format_data(data, type="outcome", snp_col="variant_id", beta_col="effect_size", se_col="standard_error", 
                         effect_allele_col="effect_allele", other_allele_col="non_effect_allele",
                         eaf_col="frequency", pval_col="pvalue", samplesize_col="sample_size", id_col="id")
  mydata <- harmonise_data(exposure_dat=exp_dat, outcome_dat=t2d_out, action=2)
  res <- mr(mydata, method_list=method)
  pleiotropy <- mr_pleiotropy_test(mydata)
  MR <- data.frame(trait=trait,method=res$method,mr_pval=res$pval,pleiotropy_pval=pleiotropy$pval)
  out <- rbind(out,MR)
}
out1 <- out[which(out$mr_pval<0.05&out$pleiotropy_pval>0.05),]
write.csv(out1,"MR.method.csv")

## step 3
file_list <- list.files(path="./output", pattern=".txt", full.names=TRUE)
out <- NULL
for (file in file_list) {
  trait <- gsub("./output/Meta_|.assoc.txt", "", file)
  method <- c("mr_ivw", "mr_ivw_mre", "mr_ivw_fe", "mr_ivw_radial",
              "mr_egger_regression", "mr_weighted_median", "mr_simple_median",
              "mr_two_sample_ml",
              "mr_simple_mode","mr_weighted_mode",
              "mr_sign","mr_penalised_weighted_median")
  data <- read.table(file, check.names=F, header=T)
  data$id <- trait
  t2d_out <- format_data(data, type="outcome", snp_col="variant_id", beta_col="effect_size", se_col="standard_error", 
                         effect_allele_col="effect_allele", other_allele_col="non_effect_allele",
                         eaf_col="frequency", pval_col="pvalue", samplesize_col="sample_size", id_col="id")
  mydata <- harmonise_data(exposure_dat=exp_dat, outcome_dat=t2d_out, action=2)
  res <- mr(mydata, method_list=method)
  result_single <- mr_singlesnp(mydata, all_method=method)
  p2 <- mr_forest_plot(result_single)
  ggsave(filename=paste0("./figure/mr_forest_plot_",trait,".pdf"),plot=p2[[paste0("PCK1 expression.",trait)]],width=8,height=8)
  generate_odds_ratios <- generate_odds_ratios(res)
  out <- rbind(out,generate_odds_ratios)
}
out1 <- out[,c(2,5,9,7,8,10,11)]
write.csv(out1,"MR.generate_odds_ratios.csv",row.names=F)

## step 4
data <- read.csv("MR.generate_odds_ratios.plot.csv")
colnames(data)[3] <- "p_value"
colnames(data)[4] <- "forest_we_need_more_spaces_for_plot_plus_plus_plus_plus_plus"
data[is.na(data)] <- ""
data$b <- as.numeric(data$b)
data$lo_ci <- as.numeric(data$lo_ci)
data$up_ci <- as.numeric(data$up_ci)

max(data$up_ci,na.rm=T)
min(data$lo_ci,na.rm=T)

tm <- forest_theme(base_size=10,
                   # Confidence interval point
                   ci_pch=15, ci_col="#762a83", ci_fill=NULL, ci_alpha=1, ci_lty=1, ci_lwd=1.5, ci_Theight=0.3,
                   # Reference line
                   refline_gp=gpar(lwd=1, lty="dashed", col="grey20"),
                   # Summary
                   summary_fill="yellow", summary_col="#4575b4",
                   # Footnote font
                   footnote_gp=gpar(cex=0.6, fontface="italic", col="red"))

pf <- forest(data[,c(1,3,4)], est=data$b, lower=data$lo_ci, upper=data$up_ci, sizes=0.5, ci_column=3, 
             xlim=c(-0.3,0.2), ticks_at=c(-0.3,-0.2,-0.1,0,0.1,0.2), xlab="Estimates (95% CI)",
             ref_line=0, theme=tm)

row=c(1,7,13,19,25,31,37,42)

pf <- edit_plot(pf, row=row+1, col=3, which="ci", gp=gpar(col="#66C2A5"))
pf <- edit_plot(pf, row=row+2, col=3, which="ci", gp=gpar(col="#FC8D62"))
pf <- edit_plot(pf, row=row+3, col=3, which="ci", gp=gpar(col="#8DA0CB"))
pf <- edit_plot(pf, row=row+4, col=3, which="ci", gp=gpar(col="#E78AC3"))
pf <- edit_plot(pf, row=row+5, col=3, which="ci", gp=gpar(col="#A6D854"))


pf <- edit_plot(pf, row=row, gp=gpar(fontface="bold"))

#"#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F"

ggsave(filename="./figure/mr_forest_plot.pdf",plot=pf,width=8,height=12)
