
####### functions ########

#' Annotate drugs with class and approval metadata
annotate_drug_classes <- function(dataset, drug_detail, match_col) {
  join_and_suffix <- function(drug_col, suffix) {
    left_join(dataset[drug_col], drug_detail, by = setNames(match_col, drug_col)) %>%
      rename_with(~ paste0(., suffix), c("MOA", "Drug_Sub_Class", "Drug_Class", "Approval_status"))
  }
  
  d1 <- join_and_suffix("Drug1", "1")
  d2 <- join_and_suffix("Drug2", "2")
  
  combined <- bind_cols(d1, d2)
  
  drug_class_pairs <- combined %>%
    transmute(
      Drug_Class1, Drug_Class2,
      Label = paste(pmin(Drug_Class1, Drug_Class2), "-", pmax(Drug_Class1, Drug_Class2))
    )
  
  drug_class_counts <- count(drug_class_pairs, Label, name = "Count") %>%
    arrange(desc(Count))
  
  list(combined, drug_class_counts, drug_class_pairs)
}


#' Count unique drugs and their proportion
count_unique_drugs <- function(data) {
  n_combos <- nrow(data)
  
  col_names <- c("Drug", "MOA", "Drug_Sub_Class", "Drug_Class", "Approval_status")
  drug1 <- data[, 1:5] %>% setNames(col_names)
  drug2 <- data[, 6:10] %>% setNames(col_names)
  
  all_drugs <- rbind(drug1, drug2)
  
  counts <- all_drugs %>%
    group_by(across(everything())) %>%
    summarise(
      Count = n(),
      Proportion = (n() / n_combos) * 100,
      .groups = "drop"
    ) %>%
    arrange(desc(Proportion))
  
  message("Unique drugs: ", nrow(counts))
  counts
}

############### gsea ################
run_msigdb <- function(tmp, ensbl.col="gene", 
                       rank.col="avg_log2FC", 
                       selected_set_list,
                       pvalueCutoff = 0.25,
                       title="",
                       basename="plot/rnaseq/gsea/msigdb/",
                       save.op=NULL){
  require(clusterProfiler)
  
  ## ranked gene list
  original_gene_list <- tmp[[rank.col]] 
  names(original_gene_list) <- tmp[[ensbl.col]]
  gene_list<-na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  ## running gsea
  list_gse = list()
  for (i in names(selected_set_list)) {
    selected_set = selected_set_list[[i]]
    gse <- GSEA(gene_list, TERM2GENE = selected_set[,c("gs_name", "ensembl_gene")], seed = TRUE, verbose = F, minGSSize=10,  pvalueCutoff = pvalueCutoff)
    
    ##saving
    if (nrow(gse@result)!=0){ #str_to_title
      gse@result$Description <- str_replace_all(gse@result$Description, c("^[^_]*_"="", "_"=" "))
      
      ## plotting
      if (!is.null(basename)){
        gsea.dotplot(gse@result, n=10, title=title, savefile=paste0(basename, "_",i,".tiff"))
      }
      
      ## saving in csv
      if (!is.null(save.op)){
        x <- as.data.frame(unique(subset(selected_set, select=-c(ensembl_gene, gene_symbol))) )
        gse@result <- left_join(gse@result, x, by= c("ID"="gs_name"))
        write.csv(gse@result, file=paste0(save.op, "_",i,".csv"), row.names = F)
      }
    }
    
    list_gse[[i]] = gse
  }
  
  return(list_gse)
}


####### statistical tests ########

.run.cor <- function(group.m, x.m, min.n=30, method="spearman", p.adj.m="fdr"){
  require(rstatix)
  
  output = bind_rows(lapply(setNames(colnames(group.m), colnames(group.m)), function(i){

    tmp_ <- rownames(group.m)[!is.na(group.m[[i]])] # remove na
    df <- data.frame(Sample_ID=tmp_, group=group.m[tmp_,i]) %>% 
      cbind(.,x.m[tmp_,]) 
    
    df <- df %>%
      pivot_longer(3:ncol(df), names_to = "Combination", values_to = "DSS") %>%
      filter(!is.na(DSS)) %>% # removes invalid variable
      group_by(Combination) %>%
      mutate(noofgroup = n() ) %>% # at least n TRUE per comparison
      filter(noofgroup >= min.n)
    
    if (nrow(df)!=0) { 
      res <- df %>% 
        group_by(Combination, noofgroup) %>% 
        reframe(cor_test_res = list(suppressWarnings(cor.test(DSS, group, method = method)))) %>%
        mutate(cor = map_dbl(cor_test_res, ~ .x$estimate),
               statistic = map_dbl(cor_test_res, ~ .x$statistic),
               p = map_dbl(cor_test_res, ~ .x$p.value) ) %>%
        dplyr::select(-cor_test_res) %>%
        ungroup() %>%
        adjust_pvalue("p", "p.adj", method = p.adj.m) %>%
        add_significance("p.adj") 
      
    } else { # when all groups are not of enough size for comparison
      res = data.frame()
    }
    res
  }), .id = "Group"
  )
  
  return(output)
}

.run.pwc.wilcox <- function(group.m, x.m, min.n=5, p.adj.m="fdr"){
  
  output = bind_rows(lapply(setNames(colnames(group.m), colnames(group.m)), function(i){
    print(i)
    
    tmp_ <- rownames(group.m)[!is.na(group.m[[i]])] # remove na
    df <- data.frame(Sample_ID=tmp_, group=group.m[tmp_,i]) %>% 
      cbind(.,x.m[tmp_,]) 
    
    df <- df %>%
      pivot_longer(3:ncol(df), names_to = "Combination", values_to = "DSS") %>%
      filter(!is.na(DSS)) %>% # removes invalid variable
      mutate(group=factor(group)) %>%
      
      group_by(Combination, group) %>% # remove those combinations w insufficient sample size for groups
      mutate(noofsubgroup = n()) %>%
      filter(noofsubgroup >= min.n) %>%
      
      group_by(Combination) %>%
      mutate(noofgroup = length(unique(group))) %>% # remove those combinations w 1 factor
      filter(noofgroup > 1) 
    
    if (nrow(df)!=0) { 
      res_ <- df %>% group_by(Combination) %>%
        wilcox_test(DSS~group, detailed = T, p.adjust.method = p.adj.m)  %>%
        adjust_pvalue("p", "global.p.adj", method = p.adj.m) %>% #for each group
        add_significance("global.p.adj") %>%
        add_significance(p.col="p", output.col="p.signif") 
      
      res <- df %>% group_by(Combination) %>% # calculates effect size
        wilcox_effsize(DSS~group) %>%
        left_join(res_, ., by=c("Combination", ".y.", "group1", "group2", "n1", "n2")) 
    } else { # when all groups are not of enough size for comparison
      res = data.frame()
    }
    
    res
  }), .id = "Group"
  )
  
  return(output)
}


# # Override rstatix's as_tidy_cor
# as_tidy_cor <- function(x){
#   estimate <- cor <- statistic <- p <- conf.low <- conf.high <- method <- NULL
#   
#   res <- x %>%
#     as_tidy_stat() %>%  # still use rstatix::as_tidy_stat
#     rename(cor = estimate) %>%
#     mutate(cor = signif(cor, 5))  # <- custom precision here
#   
#   if (res$method == "Pearson") {
#     res %>% select(cor, statistic, p, conf.low, conf.high, method)
#   } else {
#     res %>% select(cor, statistic, p, method)
#   }
# }
# 
# # assign them to the factoextra namespace
# environment(as_tidy_cor) <- asNamespace("rstatix")
# assignInNamespace("as_tidy_cor",as_tidy_cor,"rstatix")

r_squared <- function(y_val, y_pred) {
  ss_res <- sum((y_val - y_pred)^2)
  ss_tot <- sum((y_val - mean(y_val))^2)
  1 - (ss_res / ss_tot)
}


####### gene signature ########

## Loading gmt genesets
load_genesets <- function(path, ignore_cols = 1){
  x <- scan(path, what="", sep="\n")
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  for(i in 1:ignore_cols){
    y <- lapply(y, `[`, -1) 
  }
  return(y)
}

#' @param genes vector of genes considered - remove gs with <= 5 genes
load_gs <- function(category = NULL, subcategory = NULL, genes=NULL) {
  require(msigdbr)
  gs <- msigdbr(category = category, subcategory = subcategory) %>% 
    dplyr::select(gs_name, ensembl_gene, gene_symbol, gs_description, gs_exact_source)
  
  if (!is.null(genes)) {
    # more than 5 genes are present in the gene set
    # subset genes in the gene set to what is available
    gs <- gs %>%
      group_by(gs_name) %>%
      mutate(gene.in = ensembl_gene %in% genes,
             filter.out = sum(gene.in)) %>%
      filter(filter.out > 5, gene.in) %>% 
      select(-c(gene.in, filter.out))
  }
  
  return(gs)
}

## finding indexes for gene symbol
# Alias symbol - the input has been added by an HGNC editor as aa alias of the approved symbol.
# Previous symbol - the input was previously an approved gene symbol, but the gene has since been updated with the approved symbol shown.

hgnc.index <- function(df.gene, hgnc.gene = NULL){
  require(data.table)
  require(tidyverse)
  
  if (is.null(hgnc.gene)) {
    hgnc.gene <- read.csv("data/hgnc_complete_set.txt", sep="\t")
  }
  
  df <- data.table(index=character(), og=character())
  
  find_name <- toupper(df.gene)
  temp2 <- gsub("\\|$", "", paste(hgnc.gene$alias_symbol, hgnc.gene$prev_symbol, sep="|"))
  temp2 <- strsplit(toupper(temp2),"\\|")
  temp3 <- toupper(hgnc.gene$symbol)
  
  pb <- txtProgressBar(min = 0, max = length(find_name), style = 3, width = 50, char = "=")
  start_time <- Sys.time()
  for (i in seq_along(find_name)){
    index <- match(find_name[i], temp3) # checks the standardised symbol
    if (is.na(index)) { 
      index <- which(sapply(temp2, function(x) find_name[i] %in% x)) # checks if in alias+prvs symbol
    }
    df <- rbind(df, data.table(index=paste(index, collapse=","), og=find_name[i]))
    setTxtProgressBar(pb, i)
  }
  end_time <- Sys.time()
  close(pb)
  
  time_taken <- end_time - start_time
  print(paste("Time taken: ", time_taken))
  
  # Check for rows with multiple detected outputs and remove duplicates for those that alr hv exact match
  df[, index := sapply(index, function(x) {
    if (grepl(",", x)) {
      indices <- unlist(strsplit(x, ","))
      unique_indices <- unique(indices)
      # Remove duplicates seen in other rows
      unique_indices <- unique_indices[!unique_indices %in% unlist(df$index[!grepl(",", df$index)], ",")]
      paste(unique_indices, collapse = ",")
    } else {
      x
    }
  })]
  
  return(df)
}

## calculating gene sig as described in Bottomly et al. Cancer Cell 2022
## source from https://github.com/biodev/beataml2_manuscript
plage_like.scores <- function(rna_seq, geneset, rna.clin, save=NULL){
  rna.clin <- data.table(rna.clin)
  rna_seq <- t(rna_seq)
  geneset <- geneset[gene %in% colnames(rna_seq)] #& rank %in% 1:30]
  # rna.clin <- split(clin.dat, by="cohort")
  
  #center and scale separately between BM and PB to mitigate distributional differences
  
  #scaling
  bm.scale.exprs <- t(scale(rna_seq[rna.clin[specimenType == "Bone Marrow Aspirate",Sample_ID],]))
  pb.scale.exprs <- t(scale(rna_seq[rna.clin[specimenType != "Bone Marrow Aspirate",Sample_ID],]))
  scale.exprs <- cbind(bm.scale.exprs, pb.scale.exprs[rownames(bm.scale.exprs),])
  
  split.vg <- split(geneset, by="type")
  
  eigct.list <- lapply(split.vg, function(x){
    
    tmp <- prcomp(t(scale.exprs[x$gene,]), center=F, scale=F)
    
    #align with average expression
    av.expr <- colMeans(scale.exprs[x$gene,])
    
    if(sign(cor(tmp$x[,"PC1"], av.expr)) < 0){ # higher expression with higher PC
      
      tmp.dt <- data.table(Sample_ID=rownames(tmp$x), PC1=-tmp$x[,"PC1"], PC2=tmp$x[,"PC2"])
      
    }else{
      tmp.dt <- data.table(Sample_ID=rownames(tmp$x), PC1=tmp$x[,"PC1"], PC2=tmp$x[,"PC2"])
    }
    
    #roughly group
    tmp.dt[,groups:="intermediate"]
    tmp.dt[PC1 < quantile(PC1, .25), groups:="low"]
    tmp.dt[PC1 > quantile(PC1, .75), groups:="high"]
    
    #cors with eigengene
    
    kme.mat <- cor(t(scale.exprs[x$gene,])[tmp.dt$Sample_ID,], tmp.dt$PC1)
    kme.dt <- data.table(gene=rownames(kme.mat), kme=kme.mat[,1])
    kme.dt <- merge(x, kme.dt, by="gene")
    
    #also record eigenvector
    
    kme.dt <- cbind(kme.dt, ev=tmp$rotation[kme.dt$gene,"PC1"])
    
    list(summary=tmp.dt, pcs=tmp, kme=kme.dt)
  })
  
  
  eigct.dt <- cbind(rbindlist(lapply(eigct.list, "[[", "summary"), idcol="type"))
  
  eigct.dt <- merge(eigct.dt, rna.clin[,.(Sample_ID, cohort, specimenType)], by=c("Sample_ID"))
  
  if (!is.null(save)){
    write.csv(reshape2::acast(Sample_ID~type, value.var="PC1", data=eigct.dt), save)
  }
  
  return(list(summary=eigct.dt, results=eigct.list, sexprs=scale.exprs, 
              bm_scale=list(center=attr(bm.scale.exprs,"scaled:center"), scale=attr(bm.scale.exprs,"scaled:scale")), 
              pb_scale=list(center=attr(pb.scale.exprs,"scaled:center"), scale=attr(pb.scale.exprs,"scaled:scale"))))
  
}


