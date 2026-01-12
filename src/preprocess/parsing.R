library(tidyverse)
library(data.table)

#' HEAVY MANUAL EDITS required, not reproducible by code 
#' code not updated to be used for GITHUB
#' 
#' For OSHU clinical summary: parsing immunophenotype, mutations and cytogenetics
#' manual edits are required for immunophenotypes and cytogenetics
############################# mutation ##################################
new.2 <- read.csv("../../baml_combi/clinical_summary2.csv")

# 162 unique genes
unique_gene = new.2$mutationsSummary %>%
  unique(.) %>%
  gsub("[[:space:]]", "", .) %>%
  strsplit(., split="\\|") %>%
  unlist(.) %>%
  gsub("\\(.*\\)", "", .) %>%
  unique(.) %>%
  sort(.)

x <- new.2$mutationsSummary %>%
  gsub("[[:space:]]", "", .) %>%
  strsplit(., split="\\|") 

x <- lapply(x, gsub,pattern="\\(.*\\)", replacement="") 

tp <- bind_rows(lapply(setNames(seq(1,600),seq(1,600)), function(i){
  if(is_empty(x[[i]]) | sum(is.na(x[[i]]))){
    res <- setNames(rep(NA, length(unique_gene)), unique_gene) # not everything has values
  } else {
    res <- sapply(unique_gene, function(y) y %in% x[[i]]) 
    # if u wanna get specific VAF etc - use grep("ASXL1", x[[1]]) - then extract the values from the bracket
  }
  res
}))

tp <- data.frame(tp, row.names = new.2$Sample_ID )

# using the alr annotated columns as a reference to change any positive to T
y = case_when(is.na(new.2$FLT3_ITDCall)|(new.2$FLT3_ITDCall == "")  ~ NA,
              grepl("Positive", new.2$FLT3_ITDCall) ~ T,
              .default = F)
tp$FLT3 = tp$FLT3|tp$`FLT3.ITD`|tp$`FLT3.TKD`|tp$`FLT3.D835`|tp$`FLT3.V491L` | grepl("Positive", new.2$FLT3_ITDCall)
tp$FLT3[is.na(tp$FLT3)] = y[is.na(tp$FLT3)]

tp$`FLT3.ITD` = tp$`FLT3.ITD`| grepl("Positive", new.2$FLT3_ITDCall)
tp$`FLT3.ITD`[is.na(tp$`FLT3.ITD`)] = y[is.na(tp$`FLT3.ITD`)]

y = case_when(is.na(new.2$NPM1Call)|(new.2$NPM1Call == "") ~ NA,
              grepl("Positive", new.2$NPM1Call) ~ T,
              .default = F)
tp$NPM1 = tp$NPM1| grepl("Positive", new.2$NPM1Call)
tp$NPM1[is.na(tp$NPM1)] = y[is.na(tp$NPM1)]

# this column doesnt have "Negative"
y = case_when(is.na(new.2$CEBPA_Biallelic)|(new.2$CEBPA_Biallelic == "")|(new.2$CEBPA_Biallelic == "N/A") ~ F, .default = T) 
tp$CEBPA = tp$CEBPA| y 


# other mutation - 370 x 1264
wes <- read.csv("../../baml_combi/oshu_wes_targeted_mutation.csv") %>%
  mutate(value = 1) %>%
  reshape2::acast(Sample_ID ~ symbol, value.var = "value", fun.aggregate = function(x) as.integer(sum(x) > 0), fill = 0, drop = FALSE)

# overlapping
wes_overlap <- wes[, colnames(wes) %in% colnames(tp)] # 68 overlaps
tp_overlap <- tp[rownames(wes_overlap), colnames(wes_overlap)]
tp_overlap[is.na(tp_overlap)] <- FALSE
tp_overlap <- tp_overlap | wes_overlap

tp[rownames(wes_overlap), colnames(wes_overlap)] <- tp_overlap
tp[rownames(wes_overlap), ][is.na(tp[rownames(wes_overlap), ])] <- FALSE # assuming that all the other mutations are negative

# not overlapping - ignoring
wes_non_overlap <- wes[, !colnames(wes) %in% colnames(tp)]
wes_non_overlap <- wes_non_overlap[, colSums(wes_non_overlap, na.rm = TRUE) >= 5] # "CELSR2" "CPS1" "NFATC4" "PDS5B" "SMC3" "TRIO" "TRPM3" "WRNIP1"

# removing those genes with less than 5 samples
tp <- tp[, colSums(tp, na.rm = TRUE) >= 5] # 162 -> 72
tp <- rownames_to_column(tp, var = "Sample_ID")

write.csv(tp, "data/oshu/mutation.csv", row.names = F)

############################# immunophenotype ###############################
##  requires heavy manual curation
data.og <- read.csv("../../baml_combi/clinical_summary2.csv") %>%
  select(Sample_ID, surfaceAntigensImmunohistochemicalStains) %>%
  dplyr::rename(immunotype = surfaceAntigensImmunohistochemicalStains) %>%
  mutate(immunotype=ifelse(immunotype %in% c("", "not done", "Not Run", "negative (not done)", NA), NA, immunotype))

data <- data.og %>%
  filter(!is.na(immunotype))


columns = c("Sample_ID", "manual", "CD1a", "CD2","CD3","CD4","CD5","CD7","CD9",
            "CD10","CD11b","CD11c","CD13","CD14","CD15","CD16","CD19",
            "CD20","CD22","CD23","CD25",
            "CD30","CD33","CD34","CD36","CD38",
            "CD41","CD42b","CD45",
            "CD56","CD57","CD58",
            "CD61","CD64","CD71","CD79a","CD117","CD123",
            "Glyco A","Lysozyme","Myeloperoxidase","HLA-DR", "TdT","MPO","PAX5")
df = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(df) = columns


result <- rbindlist(fill=T,lapply(setNames(1:length(data$immunotype), data$Sample_ID), function(x){
  op = c(data$Sample_ID[x])
  
  tmp <- toupper(data$immunotype[x])
  tmp <- ifelse(is.na(tmp), tmp, gsub("HLA( |)DR", "HLA-DR", tmp))
  op = c(op, ifelse(grepl("negative", tmp, ignore.case = T), "TRUE", "FALSE") ) # whether manual curation is required
  tmp <- strsplit(tmp, split = c("\\, |and | \\| |:"))[[1]] #split by comma, 'and', |, colon
  
  search <- sapply(columns[3:length(columns)], function(antigen) {
    res <- paste(tmp[grep(paste0(antigen,"\\b"), tmp, ignore.case = T)], collapse = "")
    # print(res)
    if(is_empty(res) | res==""){
      res = "-"
    } else if (op[2]) { # when manual curation is required
      res <- "Check"
    } else if (sum(grepl("\\(.*\\)", res))>0) { # whether it contains a bracket
      res <- gsub(".*\\(([^)]+)\\).*", "\\1", res)
    } else {
      res <- gsub("[[:space:]]", "", gsub(antigen, "", res))
      res <- ifelse(res=="", "+", res)
    }
  } 
  ) 
  
  op = c(op, search)
  
  # print(x)
  tmp <- data.frame(t(op))
  colnames(tmp) <- columns
  tmp
}))

result <- data.frame(result)
# result[grepl("partial",result, ignore.case = T)] = "partial"

data <- left_join(data.og, result)

write.csv(data, "data/oshu/immuno2.csv", row.names = F)
##  requires heavy manual curation


############################# cytogenetics ###############################
#' functions for the common AML fusions
#' @param x dataframe with a column named 'karyotype'
cyto_search <- function(x){
  require(tidyverse)
  
  x$karyotype <- tolower(x$karyotype) %>%
    gsub("[[:space:]]", "", .)
  
  x$`PML-RARA` = grepl("t\\(15;17\\)\\(q2[24];q21.*\\)", x$karyotype)
  
  x$`RARA_re` = (grepl("t\\(.*;17\\)\\(.*;q21.*\\)", x$karyotype) | grepl("t\\(17;.*\\)\\(q21.*;.*\\)", x$karyotype) | grepl("[invdel]\\(17\\)\\(q21.*\\)", x$karyotype) ) & (!x$`PML-RARA`)
  
  
  x$`RUNX1-RUNX1T1` = grepl("t\\(.*8;21\\)\\(.*q2[12].*;q22.*\\)", x$karyotype)
  
  x$`CBFB-MYH11` = grepl("inv\\(16\\)\\(p13.*q22.*\\)", x$karyotype) | grepl("t\\(16;16\\)\\(p13.*;q22.*\\)", x$karyotype)
  
  # x$`MLLT3-KMT2A` = grepl("t\\(.*9;11\\)\\(.*p2.*;q23\\)", x$karyotype)
  
  x$`KMT2A_re` = (grepl("t\\(.*;11\\)\\(.*;q23.*\\)", x$karyotype)| grepl("t\\(11;.*\\)\\(q23.*;.*\\)", x$karyotype)) # & (!x$`MLLT3-KMT2A`)
  
  
  x$`DEK-NUP214` = grepl("t\\(6;9\\)\\(p2[123].*;q34\\)", x$karyotype)
  
  # x$`GATA2-MECOM` = grepl("inv\\(3\\)\\(q21.*q26.*\\)", x$karyotype) | grepl("t\\(3;3\\)\\(q21.*;q26.*\\)", x$karyotype) 
  
  x$`MECOM_re` = grepl("inv\\(3\\)\\(q21.*q26.*\\)", x$karyotype) | (grepl("t\\(.*;3\\)\\(.*;q26.*\\)", x$karyotype) | grepl("t\\(3;.*\\)\\(q26.*;.*\\)", x$karyotype)) # & (!x$`GATA2-MECOM`)
  
  x$`NUP98_re` = grepl("t\\(.*;11\\)\\(.*;p15.*\\)", x$karyotype) | grepl("t\\(11;.*\\)\\(p15.*;.*\\)", x$karyotype)
  
  x$`BCR-ABL1` = grepl("t\\(9;22\\)\\(q34.*;q11.*\\)", x$karyotype) 
  
  x$`Other_re` = grepl("t\\(1;3\\)\\(p36.*;q21.*\\)", x$karyotype) | grepl("t\\(1;22\\)\\(p13.*;q13.*\\)", x$karyotype) | grepl("t\\(3;5\\)\\(q25.*;q35.*\\)", x$karyotype) | grepl("t\\(8;16\\)\\(p11.*;p13.*\\)", x$karyotype) | grepl("t\\(10;11\\)\\(p13;(q14|q21)\\)", x$karyotype)| grepl("t\\(7;12\\)\\(p36.*;q13.*\\)", x$karyotype) | grepl("t\\(16;21\\)\\((p11|q24).*;q22.*\\)", x$karyotype) | grepl("inv\\(16\\)\\(p13.*;q24.*\\)", x$karyotype)
  
  x$`Other_abn` = grepl("((,|^)-5|(,|^)-7|(,|^)\\+8|(,|^)-17)", x$karyotype) | grepl("[adelt]\\((|.*;)5(|;.*)\\)\\((|.*;)q(.*|.*;.*)\\)", x$karyotype) | grepl("del\\(7\\)\\(q", x$karyotype) | grepl("[adelt]\\((|.*;)12(|;.*)\\)\\((|.*;)p(.*|.*;.*)\\)", x$karyotype) | grepl("[adel]\\(17\\)\\(p", x$karyotype) | grepl("i\\(17\\)\\(q", x$karyotype) | grepl("del\\(20\\)\\(q", x$karyotype)
  
  return(x)
}

sampleinfo <- read.csv("../../baml_combi/clinical_summary2.csv")
x = sampleinfo["karyotype"]


checking = sampleinfo[c("Sample_ID","karyotype","otherCytogenetics")] %>%
  mutate_all(.funs=~trimws(tolower(.x), which = "both")) %>%
  mutate(karyo_mask = case_when(karyotype %in% c("not available", "not done", "not performed", "", "n/a", NA) ~ "na",
                                grepl("^46,x[xy]\\[20\\](?:$|\\s)", karyotype) ~ "negative", #need atleast 20
                                grepl("^normal(?:$|\\s)", karyotype) ~ "negative",
                                .default = "present"),
         other_mask = case_when(otherCytogenetics %in% c("","not available", "n/a", NA) ~ "na",
                                grepl("^normal(?:$|\\s)", otherCytogenetics) ~ "negative",
                                grepl("negative", otherCytogenetics) ~ "negative",
                                .default = "present")
  ) %>%
  mutate(na_mask = ifelse(karyo_mask=="na" & other_mask=="na", T, F),
         normal_mask = case_when(karyo_mask=="negative" & other_mask=="negative" ~ T,
                                 karyo_mask=="na" & other_mask=="negative" ~ T,
                                 karyo_mask=="negative" & other_mask=="na" ~ T,
                                 .default = F )
  )

# write.csv(checking, "data/oshu/oshu_cyto_temp.csv") ## requires manual curation

##  requires heavy manual curation
mask <- read.csv("data/oshu/oshu_cyto_temp.csv")[["na_mask"]] # combines outcome of karyotype and FISH
mask_normal <- read.csv("data/oshu/oshu_cyto_temp.csv")[["not_detect_mask"]] # combines outcome of karyotype and FISH
othercyto <- read.csv("data/oshu/oshu_cyto_temp.csv")[["other_results"]] # outcomes of FISH
othercyto[othercyto=="NotDetected"] = NA 

x <- cyto_search(x)
x[mask,2:12] = rep(NA, 11)
# write.csv(x, "data/oshu/oshu_karyo_temp.csv") ## requires manual curation

### summarising both outcome ###
fusion = c("BCR-ABL1","CBFB-MYH11","DEK-NUP214","KMT2A_re","MECOM_re", "NUP98_re","PML-RARA","RARA_re","RUNX1-RUNX1T1")
other = c("Other_re","Other_abn(,|$)","Other_abn2(,|$)")
other2 = c("Other_re","Other_abn","Other_abn2")

x <- read.csv("data/oshu/oshu_karyo_temp.csv", check.names = F) 
x$consensusAMLFusion = sampleinfo$consensusAMLFusions
x$otherCytogenetics <- sampleinfo$otherCytogenetics
x = x[c(1:2,18,3:17)]

y <- read.csv("data/oshu/oshu_cyto_temp.csv")[,c(1:3,6:8)]
y$other_results[is.na(y$other_results)] = ""

x$AMLFusion=""
x$Other=""

for (i in 1:600) {
  tp <- str_detect(y$other_results[i], fusion)
  tp2 <- str_detect(x$consensusAMLFusion[i], fusion)
  x[i,fusion] <- x[i,fusion]|tp|tp2
  x$AMLFusion[i] = ifelse(sum(x[i,fusion]) > 0, paste0(fusion[t(x[i,fusion])], collapse = ","), NA)
  
  tp <- str_detect(y$other_results[i], other)
  x[i,other2] <- x[i,other2]|tp
  x$Other[i] = ifelse(sum(x[i,other2]) > 0, paste0(other2[t(x[i,other2])], collapse = ","), NA)
}

x <- mutate(x, outcome = case_when(!is.na(AMLFusion) ~ AMLFusion,
                                   consensusAMLFusion != "" ~ consensusAMLFusion,
                                   !is.na(Other) ~ "Other",
                                   not_detect_mask == TRUE ~ "NotDetected",
                                   na_mask == TRUE ~ NA,
                                   .default = consensusAMLFusion))

## output file to save
write.csv(x,"data/oshu/oshu_cyto2.csv", row.names = F) ## requires manual edits 
