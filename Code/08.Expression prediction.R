
#### data processed ####
setwd("/disk221/xieqq/JINHUA138/Database.4.TWAS_server")
library("data.table")
library("parallel")
library('GHap')

# 将snp文件按染色体拆分
A <- function(i) {
  dir.create("./data_raw", showWarnings=F, recursive=T)
  bcftools <- '/disk191/zzy/software/bcftools-1.12/bcftools'
  system(paste0(bcftools,' view -Oz -r ',i,' /disk212/liucheng/TWASmethod.haplo_8.30/TWAS_application/SHXX_TB/data_raw/GS.DD.all.keep.snpID.vcf.gz -o ./data_raw/chr',i,'.vcf.gz'))
}

mclapply(1:18, A, mc.cores=10) 


# conform & impute
conform_gt <- '/disk191/chenzt/software/conform-gt.24May16.cee.jar'
beagle <- '/disk191/chenzt/software/beagle.21Apr21.304.jar'
bcftools <- '/disk191/zzy/software/bcftools-1.12/bcftools'

con_imp <- function(tissue,i){
  # ===== 1. 确保目录存在 =====
  dir.create(paste0("./data_processed/",tissue,"/SNP"), showWarnings=F, recursive=T)
  # ===== 2. 拼接路径 =====
  ref_vcf <- paste0("/disk191/liucheng/FarmGTEx_DATA/pigGTEx/genotype/SNP/",tissue,"/chr",i,".vcf.gz")
  raw_vcf <- paste0("./data_raw/chr",i,".vcf.gz")
  cof_vcf_1 <- paste0("./data_processed/",tissue,"/SNP/cof_chr",i)
  cof_vcf_2 <- paste0(cof_vcf_1,".vcf.gz")
  imp_vcf_1 <- paste0("./data_processed/",tissue,"/SNP/imp_chr",i)
  imp_vcf_2 <- paste0(imp_vcf_1,".vcf.gz")
  # ===== 3. conform-gt =====
  cmd1 <- paste0('java -jar ',conform_gt,' ref=',ref_vcf,' match=POS chrom=',i,' gt=',raw_vcf,' out=',cof_vcf_1)
  system(cmd1)
  # ===== 4. Beagle =====
  cmd2 <- paste0('java -Xmx30g -Djava.io.tmpdir=./tmp -jar ',beagle,' ref=',ref_vcf,' ne=100 chrom=',i,' gt=',cof_vcf_2,' seed=4137 nthreads=10 out=',imp_vcf_1)
  system(cmd2)
  # ===== 5. bcftools index =====
  cmd3 <- paste0(bcftools," index ",imp_vcf_2)
  system(cmd3)
}

tissue <- c('Adipose') #'Adipose','Colon','Duodenum','Ileum','Jejunum','Large_intestine','Liver','Muscle','Small_intestine'
for (t in tissue) {
  mclapply(1:18, function(i) con_imp(tissue=t, i), mc.cores=5)
}


# 划分单体型SW
plink <- '/disk191/chenzt/software/plink'
plink2 <- '/disk191/chenzt/software/plink2'

SW <- function(sw,tissue,i){
  # ===== 1. 确保目录存在 =====
  SW <- paste0("SW_",sw)
  dir.create(paste0("./data_processed/",tissue,"/",SW), showWarnings=F, recursive=T)
  # ===== 2. 拼接路径 =====
  imp_vcf <- paste0("./data_processed/",tissue,"/SNP/imp_chr",i,".vcf.gz")
  chr_vcf <- paste0("./data_processed/",tissue,"/",SW,"/chr",i)
  chr_vcf_2 <- paste0(chr_vcf,".vcf")
  # ===== 3. plink2 =====
  cmd1 <- paste0(plink2,' --vcf ',imp_vcf,' --export haps --out ',chr_vcf)
  system(cmd1)
  ghap.oxford2phase(input.files=chr_vcf, out.file=chr_vcf, ncores=10)
  ghap.compress(input.file=chr_vcf, out.file=chr_vcf, ncores=10)
  phase <- ghap.loadphase(chr_vcf)
  blocks.mkr <- ghap.blockgen(phase, windowsize=sw, slide=sw, unit="marker")
  ghap.haplotyping(phase, blocks.mkr, outfile=chr_vcf, binary=T, ncores=10)
  haplo <- ghap.loadhaplo(chr_vcf)
  ghap.hap2plink(haplo, outfile=chr_vcf)
  # ===== 4. plink =====
  cmd2 <- paste0(plink," --bfile ",chr_vcf," --export vcf --out ",chr_vcf)
  system(cmd2)
  cmd3 <- paste0("bgzip -f ",chr_vcf_2)
  system(cmd3)
}

tissue <- c('Adipose','Colon','Duodenum','Ileum','Jejunum','Large_intestine','Liver','Muscle','Small_intestine')
for (i in c(1:18)){SW(sw=10,tissue='Small_intestine',i)}

# 合并染色体
cd /disk221/xieqq/JINHUA138/Database.4.TWAS_server
FILES=""
for i in {1..18}; do
    FILES="$FILES ./data_processed/Small_intestine/SW_10/chr${i}.vcf.gz"
done
/disk191/zzy/software/bcftools-1.12/bcftools concat -Oz -o ./data_processed/Small_intestine/SW_10/chr_all.vcf.gz $FILES
/disk191/zzy/software/bcftools-1.12/bcftools index -t ./data_processed/Small_intestine/SW_10/chr_all.vcf.gz


#### 运行TWAS ####
cd /disk213/xieqq/JINHUA138/Database.4.TWAS_server

python3 /disk212/liucheng/software/MetaXcan/software/Predict.py \
--model_db_path /disk213/xieqq/JINHUA138/Database.4.TWAS_server/SW_10_dbs/PigGTEx_Liver_ElasticNet_models.db \
--model_db_snp_key varID \
--vcf_genotypes /disk221/xieqq/JINHUA138/Database.4.TWAS_server/data_processed/Liver/SW_10/chr_all.vcf.gz \
--vcf_mode genotyped \
--prediction_output /disk213/xieqq/JINHUA138/Database.4.TWAS_server/Liver_sw10_predict.txt \
--prediction_summary_output /disk213/xieqq/JINHUA138/Database.4.TWAS_server/Liver_sw10_summary.txt \
--verbosity 9 \
--throw


#### PCK1 ####
library(data.table)
library(tidyr)
library(dplyr)

setwd("/disk213/xieqq/JINHUA138/Database.4.TWAS_server")

genenames <- c("ENSSSCG00000007371","ENSSSCG00000030921","ENSSSCG00000005423")
#PCK1, ENSSSCG00000007507
#HNF4A, ENSSSCG00000007371
#APOA1, ENSSSCG00000030921
#ABCA1, ENSSSCG00000005423

out1 <- NULL
for (i in c("Adipose","Liver","Muscle","Small_intestine")){
  predict <- fread(paste0(i,"_sw5_predict.txt"), data.table=F) %>% select(any_of(c("FID",genenames)))
  if (ncol(predict)>1){
    predict <- predict %>% pivot_longer(cols=-FID, names_to="gene", values_to="value") %>% mutate(Tissue=i)
    out1 <- rbind(out1,predict)
  }
}

out2 <- fread("Small_intestine_sw5_predict.txt") %>% select(any_of(c("FID","ENSSSCG00000007507"))) %>% 
  pivot_longer(cols=-FID, names_to="gene", values_to="value") %>% mutate(Tissue="SI")

A_group <- out2 %>% mutate(Group=case_when(value<=0.90746 ~ 1, value>0.90746&value<=0.9091849 ~ 2, value>0.9091849 ~ 3))

G_group <- out1 %>% left_join(.,A_group[,c("FID","Group")],by="FID") %>% filter(gene=="ENSSSCG00000005423") %>% filter(Tissue=="Liver") %>%
  arrange(Group,value) 

#相关分析
cor.test(G_group$value, G_group$Group)


#### APOA1 + ABCA1 ####
out1 <- NULL
for (i in c("Adipose","Liver","Muscle","Small_intestine")){
  predict <- fread(paste0(i,"_sw2_predict.txt"), data.table=F) %>% select(any_of(c("FID","ENSSSCG00000005423")))
  if (ncol(predict)>1){
    colnames(predict) <- predict %>% pivot_longer(cols=-FID, names_to="gene", values_to="value") %>% mutate(Tissue=i)
    out1 <- rbind(out1,predict)
  }
}

i="Muscle"
out1 <- fread(paste0(i,"_sw2_predict.txt")) %>% select(any_of(c("FID","ENSSSCG00000005423"))) %>% #ABCA1
  pivot_longer(cols=-FID, names_to="gene", values_to="value") %>% mutate(Tissue="SI")

out2 <- fread("Small_intestine_sw2_predict.txt") %>% select(any_of(c("FID","ENSSSCG00000030921"))) %>% #APOA1
  pivot_longer(cols=-FID, names_to="gene", values_to="value") %>% mutate(Tissue="SI")

out <- left_join(out1[,c("FID","value")],out2[,c("FID","value")],by="FID",suffix = c(".ABCA1", ".APOA1"))

#相关分析
cor.test(out$value.ABCA1, out$value.APOA1)
