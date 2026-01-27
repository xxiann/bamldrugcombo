## initialising renv
renv::init(bioconductor = "3.19")

## installing packages
renv::install("rstatix@0.7.2", "R.utils@2.12.3", "reshape2@1.4.4", "coin@1.4-3", "caret@6.0-94", "pROC@1.18.5") # installs new versions Rcpp@1.0.13

renv::install("tidyverse@2.0.0", "ggplot2@3.5.1", "cowplot@1.1.3", "patchwork@1.2.0", "svglite@2.1.3") # installs new versions; if break install "data.table@1.15.4", "ggplot2@3.5.1", "RColorBrewer@1.1-3"

renv::install("bioc::ComplexHeatmap@2.20.0", "dendsort@0.3.4", "magick") # installs new versions if break install "circlize@0.4.16", "matrixStats@1.5.0"

renv::install("VennDiagram@1.7.3", "gtsummary@2.0.1", "fastDummies@1.7.4")

renv::install("viridis@0.6.5", "pals@1.9",
              "survminer@0.4.9", "survival@3.6-4",
              "ggrepel@0.9.5", "ggpubr@0.6.0", 
              "ggExtra@0.10.1", "ggbeeswarm@0.7.2", "ggfortify@0.4.17")

renv::install("WGCNA@1.72-5", "CEMiTool@1.28.0", "bioc::ROTS@1.32.0", "glmnet@4.1-8",
              "limma@3.60.6", "bioc::DESeq2@1.44.0", "bioc::edgeR@4.2.1", 
              "bioc::msigdbr@7.5.1", "bioc::enrichR@3.2", 
              "bioc::clusterProfiler@4.12.6", "enrichplot@1.24.2")

renv::install("extrafont", "Cairo", "showtext", "svglite")
extrafont::font_import() # only once

## for DSS
# packages <- c("reshape","scales", "drc", "caTools", "MESS", "svMisc", "egg", "pheatmap", "bioc::sva", "bioc::pcaMethods"
# )
# packages <- c("drc", "e1071", "furrr", "MESS", "mice", "minpack.lm", "svMisc")
# renv::install(packages)


## older package version
renv::install("package/Matrix_1.6-1.1.tar.gz", rebuild = T)
renv::install("Rcpp@1.0.13") 

renv::install("package/seurat-4.3.0.tar.gz", rebuild = T)
renv::install("package/seurat-object-4.1.3.tar.gz", rebuild = T)
renv::install("package/irlba_2.3.5.1.tar.gz", rebuild = T)
renv::install("package/harmony_1.0.1.tar.gz", rebuild = T)
renv::install("uwot@0.2.2") 

## for scrna
renv::install("slingshot@2.12.0", "tradeSeq@1.18.0", "huayc09/SeuratExtend", "DelayedMatrixStats@1.26.0")

## updating
# renv::update()

## saving
renv::snapshot()

## restoring
# renv::activate()
# renv::restore()