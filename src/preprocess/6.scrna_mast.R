library(tidyverse)
library(Seurat)

source("src/colors.R")
set.seed(10)

all_ <- readRDS(file="output/all_samples/all_.rds")

all_[["RNA"]] <- JoinLayers(all_[["RNA"]])

all_$CellType <- factor(all_$CellType, levels=hierarchy_order[hierarchy_order %in% all_$CellType])
all_$celltype <- case_when(all_$CellType == 'HSC-like' ~ "HSC",
                           all_$CellType == 'Prog-like' ~ "Prog", 
                           all_$CellType == 'GMP-like' ~ "GMP", 
                           all_$CellType == 'ProMono-like' ~ "ProMono", 
                           all_$CellType == 'Mono-like' ~ "Mono", 
                           all_$CellType == 'cDC-like' ~ "cDC",
                           .default = all_$CellType)
all_$celltype <- factor(all_$celltype, levels=hierarchy_order[hierarchy_order %in% all_$celltype])

all_$cell_status <- paste0(ifelse(grepl("AML",all_$orig.ident) ,"AML","Healthy"), "_", all_$PredictionRefined)
all_$cell_status <- ifelse(all_$cell_status == "AML_normal", "AML_nonmalignant", all_$cell_status)

# to remove NA-783 & unclear -1133 (all belongs to BM5-34p38n) -> 36,494 cells
all_ <- subset(all_, cells = colnames(all_)[!all_$PredictionRefined %in% c("", "unclear")])

########### lineage specific: DEG between primitive and mature in each cell status ###############
selected <- subset(all_, subset = celltype %in% c("Prog", "ProMono", "HSC", "cDC", "Mono")) 

selected$group <- ifelse(selected$celltype %in% c("ProMono","Mono", "cDC"), "mature", "primitive")

#######
selected_ <- subset(selected, subset = PredictionRefined == "malignant") 

df <- FindMarkers(selected_, group.by = "group", ident.1 = "primitive", ident.2 = "mature", test.use = "MAST")

write.csv(df, "output/all_samples/pseudo2/findmarkers_malignant_primitive_mature.csv")

#######
selected_ <- subset(selected, subset = cell_status == "AML_nonmalignant") 

df <- FindMarkers(selected_, group.by = "group", ident.1 = "primitive", ident.2 = "mature", test.use = "MAST")

write.csv(df, "output/all_samples/pseudo2/findmarkers_nonmalignant_primitive_mature.csv")

#######
selected_ <- subset(selected, subset = cell_status == "Healthy_normal") 

df <- FindMarkers(selected_, group.by = "group", ident.1 = "primitive", ident.2 = "mature", test.use = "MAST")

write.csv(df, "output/all_samples/pseudo2/findmarkers_normal_primitive_mature.csv")

#######
df <- FindMarkers(selected, group.by = "group", ident.1 = "primitive", ident.2 = "mature", test.use = "MAST")

write.csv(df, "output/all_samples/pseudo2/findmarkers_primitive_mature.csv")


############# lineage specific: DEG between malignant, nonmalignant, normal in mature myeloid compartment #############
selected <- subset(all_, subset = celltype %in% c("GMP", "Prog", "ProMono", "HSC", "cDC", "Mono")) #, "pDC"
df <- FindMarkers(selected, group.by = "PredictionRefined", ident.1 = "malignant", ident.2 = "normal", test.use = "MAST")
write.csv(df, "output/all_samples/pseudo2/findmarkers.csv")

#######

selected <- subset(all_, subset = celltype %in% c("ProMono","Mono", "cDC"))
df <- FindMarkers(selected, group.by = "PredictionRefined", ident.1 = "malignant", ident.2 = "normal", test.use = "MAST")
write.csv(df, "output/all_samples/pseudo2/findmarkers_mono_promono_cDC.csv")

#######

selected <- subset(all_, subset = celltype %in% c("ProMono","Mono", "cDC"))

df <- FindMarkers(selected, group.by = "cell_status", ident.1 = "AML_nonmalignant", ident.2 = "Healthy_normal", test.use = "MAST")
write.csv(df, "output/all_samples/pseudo2/findmarkers_mono_promono_cDC_nonmalignant_normal.csv")

#######

selected <- subset(all_, subset = celltype %in% c("ProMono","Mono", "cDC"))

df <- FindMarkers(selected, group.by = "cell_status", ident.1 = "AML_malignant", ident.2 = "AML_nonmalignant", test.use = "MAST")
write.csv(df, "output/all_samples/pseudo2/findmarkers_mono_promono_cDC_malignant_nonmalignant.csv")

#######

selected <- subset(all_, subset = celltype %in% c("ProMono","Mono", "cDC"))

df <- FindMarkers(selected, group.by = "cell_status", ident.1 = "AML_malignant", ident.2 = "Healthy_normal", test.use = "MAST")
write.csv(df, "output/all_samples/pseudo2/findmarkers_mono_promono_cDC_malignant_normal.csv")

