## using BAML rawcounts; 1) normalize to VST 2) calculate zeng's gene signature scores
## ct.mat not included - appended behind "data/oshu/rnaseq_sampleinfo.csv"

library(tidyverse)
library(R.utils)
library(DESeq2)
source("src/function/functions.R")

## applying onto dataset: baml
sampleinfo <- read.csv("data/baml/baml_sampleinfo_rnaseq.csv") 
wave <- sampleinfo[,c("Sample_ID" ,"dbgap_rnaseq_sample", "Disease.status", "cohort", "cohort2")] %>%
  mutate(cohort = ifelse(is.na(cohort), cohort2, cohort))

rownames(wave) <- wave$Sample_ID

## rawcount
baml_rc <- readRDS("data/baml/baml_full_raw_count.rds")
rownames(baml_rc) <- baml_rc$GeneID
gene <- baml_rc[,1:2]
baml_rc <- baml_rc[3:720]

## VST (for WGCNA)
dds <- DESeqDataSetFromMatrix(countData = baml_rc,
                              colData = wave,
                              design = ~ cohort)

## Perform VST
vst_data <- vst(dds)
normalized_counts <- assay(vst_data)
# write.csv(cbind(gene, normalized_counts), file = "data/baml/baml_FULL_vst_matrix.csv", row.names = F)
# saveRDS(normalized_counts, file = "data/baml/baml_FULL_vst_matrix.rds")

####### gene signature scores ########
data <- normalized_counts

genes <- read.csv("data/baml/fullgene_baml_75vs75vs90_hngc.csv")
genes <- genes[genes$gene_id %in% rownames(data),]
rownames(genes) <- seq(1:nrow(genes))

andyzeng.gs <- load_genesets("data/AMLCellType_Genesets_all.gmt", ignore_cols = 2) 
andyzeng_ <- hgnc.index(unique(unlist(andyzeng.gs)), genes) # 509
andyzeng_ <- cbind(andyzeng_,genes[andyzeng_$index,]) %>% na.omit(.) # 491
andyzeng_ <- setNames(andyzeng_$gene_id, andyzeng_$og)

geneset <-rbindlist(lapply(names(andyzeng.gs), function(x){
  gene = andyzeng.gs[[x]]
  gene = na.omit(as.character(andyzeng_[gene]))
  tmp = data.table(gene)[,rank:=.I]
  tmp[,type:=x] 
  tmp[,.(rank, type, gene)]
}))

x <- plage_like.scores(data, geneset, sampleinfo)
ct.mat <- reshape2::acast(Sample_ID~type, value.var="PC1", data=x$summary)

# write.csv(ct.mat, "data/BAML_plage_like.csv") # 37 healthy samples not included
