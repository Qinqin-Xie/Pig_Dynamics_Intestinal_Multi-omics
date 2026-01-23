

#### data preparation ####

#### 处理基因文件 ####
library(dplyr)
library(tidyr)
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)

setwd("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/5.all_omics/3.functional_analysis_1")

load("final_cluster_info.rdata")
load("count_ID.rdata")

out_gene <- final_cluster_info %>% filter(old_omics!="2.Untargeted_metabolomics") %>% 
  mutate(molecule=str_split_fixed(molecule, "\\|", 2)[,1]) %>% 
  mutate(molecule=ifelse(molecule %in% count_ID$ensembl_gene_id, count_ID$external_gene_name[match(molecule, count_ID$ensembl_gene_id)], molecule))

# KEGG 富集结果
for (i in 1:12) {
  cluster_gene <- out_gene %>% filter(!grepl("^ENSSSCG", molecule)) %>% filter(!grepl("\\(", molecule)) %>% filter(cluster==as.character(i))
  kegg_list <- list()
  for (j in unique(cluster_gene$omics)){
    # 1. 提取对应 cluster 的基因
    cluster_gene1 <- cluster_gene %>% filter(omics==j)
    if (nrow(cluster_gene1) == 0) {kegg_list[[j]]=NULL; next}
  
    # 2. SYMBOL → ENTREZID
    trans <- bitr(cluster_gene1$molecule, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    if (is.null(trans) || nrow(trans) == 0) {kegg_list[[j]]=NULL; next}
  
    # 3. KEGG 富集分析
    ekegg <- enrichKEGG(gene=trans$ENTREZID,organism="hsa",pAdjustMethod="fdr",pvalueCutoff=0.05)
    if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {kegg_list[[j]] <- NULL
    } else {kegg_list[[j]] <- as.data.frame(ekegg) %>% {rownames(.) <- NULL; .}}
  }
  save(kegg_list, file=paste0("kegg_list_cluster_",i,".rdata"))
}

# GO 富集结果
for (i in 1:12) {
  cluster_gene <- out_gene %>% filter(!grepl("^ENSSSCG", molecule)) %>% filter(!grepl("\\(", molecule)) %>% filter(cluster==as.character(i))
  GO_list <- list()
  for (j in unique(cluster_gene$omics)){
    # 1. 提取对应 cluster 的基因
    cluster_gene1 <- cluster_gene %>% filter(omics==j)
    if (nrow(cluster_gene1) == 0) {GO_list[[j]]=NULL; next}
    
    # 2. SYMBOL → ENTREZID
    trans <- bitr(cluster_gene1$molecule, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    if (is.null(trans) || nrow(trans) == 0) {GO_list[[j]]=NULL; next}
    
    # 3. GO 富集分析
    ego <- enrichGO(gene=trans$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", 
                    pAdjustMethod="fdr", pvalueCutoff=0.05, qvalueCutoff=0.05)
    
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {GO_list[[j]] <- NULL
    } else {GO_list[[j]] <- as.data.frame(ego) %>% {rownames(.) <- NULL; .}}
  }
  save(GO_list, file=paste0("GO_list_cluster_",i,".rdata"))
}



#### 处理代谢文件 ####
# source /disk213/xieqq/miniconda/bin/activate /disk213/xieqq/miniconda/envs/r_massdb
library(massdatabase)
library(tidymass)
library(purrr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

setwd("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/5.all_omics/3.functional_analysis_1")

load("final_cluster_info.rdata")
data("kegg_hsa_pathway")

metabolite_count <- sapply(kegg_hsa_pathway@compound_list, function(compounds) {nrow(compounds)})
length(unique(sapply(kegg_hsa_pathway@compound_list, function(x) {x$`KEGG.ID`})))
pathway_metabolites <- data.frame(
  pathway_id = kegg_hsa_pathway@pathway_id,
  pathway_name = kegg_hsa_pathway@pathway_name,
  metabolite_count = metabolite_count
)
save(pathway_metabolites, file="pathway_metabolites.rdata")

metabolomics <- read.csv("/disk213/xieqq/JINHUA138/Metabolome_analysis.1.Day_compare/data/metabolomics_inf.csv") %>%
  select(Metabolomics_name,Kegg_ID) %>%
  filter(Kegg_ID != "")

out_meta <- final_cluster_info %>% filter(old_omics=="2.Untargeted_metabolomics") %>% 
  mutate(molecule=ifelse(molecule %in% metabolomics$Metabolomics_name, metabolomics$Kegg_ID[match(molecule, metabolomics$Metabolomics_name)], molecule))

met_list <- list()
for (i in 1:12) {
  # 1. 提取对应 cluster 的代谢物
  cluster_meta <- out_meta %>% filter(!grepl("^metabolomics_", molecule)) %>% filter(cluster==as.character(i))
  if (nrow(cluster_meta) == 0) {met_list[[paste0("cluster_",i)]]=NULL; next}
    
  # 2. 代谢物的KEGG ID
  query_id_kegg <- cluster_meta %>% select(molecule) %>% unlist()
  if (is.null(query_id_kegg) || length(query_id_kegg) == 0) {met_list[[paste0("cluster_",i)]]=NULL; next}
    
  # 3. KEGG 富集分析
  met_enrich <- enrich_kegg(query_id=query_id_kegg, query_type="compound", pathway_database=kegg_hsa_pathway,
                            method="hypergeometric", p_adjust_method="fdr", p_cutoff=0.05)
  
  if (is.null(met_enrich)) {met_list[[paste0("cluster_",i)]] <- NULL
  } else if (nrow(as.data.frame(met_enrich@result)) == 0) {met_list[[paste0("cluster_",i)]] <- NULL
  } else {met_list[[paste0("cluster_",i)]] <- met_enrich@result %>% {rownames(.) <- NULL; .}}
}
save(met_list, file="met_list.rdata")

#load("met_list.rdata") 
met_enrich <- map_dfr(met_list, ~ as.data.frame(.x)[, c("pathway_name","pathway_id","p_value","mapped_id")], .id="GeneSet") %>%
  filter(mapped_id!="") %>% filter(p_value<0.05) %>% 
  distinct()


#### 去冗余--GO：Wang；KEGG：Jaccard ####
library(simplifyEnrichment)
library(igraph)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(ggforce)
jaccard_similarity <- function(a, b) {length(intersect(a, b)) / length(union(a, b))}

setwd("/disk213/xieqq/JINHUA138/Multi-omics_analysis.1.Dynamics/5.all_omics/3.functional_analysis_1")
dir.create("./figures")
dir.create("./results")

i=9

## Wang相似度 + 构建GO模块
load(paste0("GO_list_cluster_",i,".rdata"))
go_description <- map_dfr(GO_list, ~ as.data.frame(.x)[, c("Description","ID")], .id="GeneSet") %>%
  select(Description,ID) %>%
  distinct()
go_all <- map_dfr(GO_list, ~ as.data.frame(.x)[, c("ID","p.adjust","geneID")], .id="GeneSet")
go_terms <- unique(go_all$ID)
go_sim <- simplifyEnrichment::GO_similarity(go_terms, measure="Wang", ont="BP")

go_cutoff <- 0.7
edges_go <- which(go_sim > go_cutoff, arr.ind=TRUE) %>% 
  as.data.frame() %>% 
  filter(row!=col) %>%
  mutate(from=rownames(go_sim)[row], to=colnames(go_sim)[col], weight=go_sim[cbind(row, col)])
g_go <- graph_from_data_frame(edges_go[, c("from", "to", "weight")], directed=FALSE)
go_modules <- cluster_edge_betweenness(g_go)
go_module_df <- tibble(term=names(membership(go_modules)), module=membership(go_modules))
go_represent <- go_all %>%
  inner_join(go_module_df, by=c("ID"="term")) %>%
  group_by(module) %>% 
  slice_min(p.adjust, n=1) %>% 
  ungroup() %>%
  mutate(DB="GO")

## Jaccard相似度 + 构建KEGG模块
load(paste0("kegg_list_cluster_",i,".rdata"))
kegg_description <- map_dfr(kegg_list, ~ as.data.frame(.x)[, c("Description","ID")], .id="GeneSet") %>%
  select(Description,ID) %>%
  distinct()
kegg_all <- map_dfr(kegg_list, ~ as.data.frame(.x)[, c("ID","p.adjust","geneID")], .id="GeneSet")
kegg_members <- kegg_all %>% mutate(members=strsplit(geneID, "/")) %>% select(ID, members)
kegg_terms <- unique(kegg_members$ID)
jaccard_mat <- outer(
  seq_along(kegg_terms),
  seq_along(kegg_terms),
  Vectorize(function(i, j) {jaccard_similarity(kegg_members$members[[i]],kegg_members$members[[j]])})
)
dimnames(jaccard_mat) <- list(kegg_terms, kegg_terms)
  
kegg_cutoff <- 0.5
edges_kegg <- which(jaccard_mat > kegg_cutoff, arr.ind=TRUE) %>%
  as.data.frame() %>%
  filter(row!=col) %>%
  mutate(from=rownames(jaccard_mat)[row], to=colnames(jaccard_mat)[col], weight=jaccard_mat[cbind(row, col)])
g_kegg <- graph_from_data_frame(edges_kegg[, c("from","to","weight")],directed=FALSE)
kegg_modules <- cluster_edge_betweenness(g_kegg)
kegg_module_df <- tibble(term=names(membership(kegg_modules)), module=membership(kegg_modules))
kegg_represent <- kegg_all %>%
  inner_join(kegg_module_df, by=c("ID"="term")) %>%
  group_by(module) %>%
  slice_min(p.adjust, n=1) %>%
  ungroup() %>%
  mutate(DB="KEGG")


#### 跨数据库整合(Jaccard相似度 + 模块) ####
combined_terms <- bind_rows(go_represent, kegg_represent)
term2gene <- combined_terms %>% mutate(members=strsplit(geneID, "/")) %>% select(term=ID, members)
terms <- term2gene$term
cross_jaccard <- outer(
  seq_along(terms),
  seq_along(terms),
  Vectorize(function(i, j) {jaccard_similarity(term2gene$members[[i]],term2gene$members[[j]])})
)
dimnames(cross_jaccard) <- list(terms, terms)

cutoff <- 0.1 #测试
edges <- which(cross_jaccard > cutoff, arr.ind=TRUE) %>%
  as.data.frame() %>%
  filter(row!=col) %>%
  mutate(from=rownames(cross_jaccard)[row], to=colnames(cross_jaccard)[col], weight=cross_jaccard[cbind(row, col)])
graph <- graph_from_data_frame(edges[, c("from","to","weight")],directed=FALSE)
modules <- cluster_edge_betweenness(graph)

entity_annotation <- bind_rows(
  go_represent %>% left_join(go_description,by="ID") %>%
    mutate(entity=paste0("GO_", module)) %>%
    select(entity, ID, Description, p.adjust, DB),
  kegg_represent %>% left_join(kegg_description,by="ID") %>%
    mutate(entity=paste0("KEGG_", module)) %>%
    select(entity, ID, Description, p.adjust, DB)
) %>% 
  mutate(logP=-log10(p.adjust)) 

set.seed(123)
modules <- cluster_edge_betweenness(graph)
layout <- layout_with_fr(graph)

node_df <- tibble(ID=names(membership(modules)), cluster=membership(modules)) %>%
  left_join(entity_annotation, by="ID") %>%
  mutate(x=layout[,1], y=layout[,2])

node_df2 <- node_df %>%
  group_by(type=case_when(grepl("^GO", ID) ~ "GO", grepl("^hsa", ID) ~ "KEGG")) %>%
  slice_max(logP, n=15, with_ties=FALSE) %>%
  ungroup() %>%
  select(ID)

cluster_hulls <- node_df %>% group_by(cluster) %>% filter(n()>=3)

edge_df <- igraph::as_data_frame(graph, what="edges") %>% 
  distinct() %>%
  left_join(node_df %>% select(ID, x, y), by=c("from"="ID")) %>%
  dplyr::rename(x_from=x, y_from=y) %>%
  left_join(node_df %>% select(ID, x, y), by=c("to"="ID")) %>%
  dplyr::rename(x_to=x, y_to=y) %>% 
  select(from,x_from,y_from,to,x_to,y_to,weight)

result <- list()
result[["go_represent"]] <- go_represent
result[["kegg_represent"]] <- kegg_represent
result[["entity_annotation"]] <- entity_annotation
result[["edge_df"]] <- edge_df
save(result,file=paste0("./results/Pathway_enrichment_cluster_",i,".rdata"))

plot <- ggplot() +
  geom_mark_ellipse(data=cluster_hulls, aes(x=x, y=y, group=cluster), fill="grey90", alpha=0.35, expand=unit(3,"mm")) +
  geom_segment(data=edge_df, aes(x=x_from, y=y_from, xend=x_to, yend=y_to, linewidth=weight), color="grey20", alpha=0.7) +
  geom_point(data=node_df, aes(x=x, y=y, size=logP, fill=DB), shape=21, color="black") +
  geom_text(data=node_df, aes(x=x, y=y, label=Description), size=2, vjust=0, check_overlap=TRUE) +
  scale_linewidth_continuous(name="sim", range=c(0.2, 1)) +
  scale_size_continuous(name="-log10(P.adjust)", range=c(1, 5)) +
  scale_fill_manual(values=c(GO="#4C72B0", KEGG="#DD8452")) +
  theme_void() +
  theme(legend.position="right")
ggsave(paste0("./figures/Pathway_enrichment_cluster_",i,".pdf"),plot=plot,width=7,height=6)
