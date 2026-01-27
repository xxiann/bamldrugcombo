# BeatAML Drug Combo Manuscript

Citation 

# Setup

R version: 4.4.1, Bioconductor 3.19, renv 1.0.11

To reproduce fully, refer to `src/init.R` for specific versions of packages used

```r
## restoring renv
renv::activate()
renv::restore()
```

Note: Seurat versions (4.3.0 vs 5.1.0) used to process single cell is different from the loaded renv and not included.

# Data

Abbreviations for datasets involved

- Beat AML drug combo - bamlcombo; bamlcombi; oshu
- Beat AML 2.0 ([Tyner et al 2018](https://doi.org/10.1038/s41586-018-0623-z), [Bottomly et al 2022](https://doi.org/10.1016/j.ccell.2022.07.002)) - baml
- Beat AML Ven Combinations ([Eide et al 2023](https://doi.org/10.1158/2643-3230.bcd-23-0014)) - bamlvencombo
- Helsinki ([Malani et al 2022](https://doi.org/10.1158/2159-8290.cd-21-0410)) - fpmtb

## Source data of bamlcombo

| Table | Files | Remarks |
| --- | --- | --- |
| Sample mapping | `data/oshu/oshu_sample_mapping.csv` | Sample identifier mapping across all datasets involved in this study, including Beat AML and Beat AML-Ven Combo. Information pertaining to the data sources and source-specific identifiers were included. |
| Clinical summary | `data/oshu/clinical_summary.csv` |  |
| WES targeted sequencing | `data/oshu/mutation.csv` | Retrieved from [BeatAML2 Github page](https://biodev.github.io/BeatAML2/) <br/> Parsed - refer to `src/preprocess/parsing.R` |
| Surface antigen | `data/oshu/clinical_immuno.csv` | Parsed and manually curated - refer to `src/preprocess/parsing.R` |
| Consensus AML fusion | `data/oshu/clinical_cyto.csv` | Parsed and manually curated - refer to `src/preprocess/parsing.R` |
| RNA sequencing | (raw count) `data/oshu/rnaseq_rawcount.csv` <br/> (VST) `data/oshu/rnaseq_vst.rds` <br/> (sampleinfo) `data/oshu/rnaseq_sampleinfo.csv` | Retrieved from [BeatAML2 Github page](https://biodev.github.io/BeatAML2/) and [Vizome](http://www.vizome.org/additional_figures_BeatAML.html) <br/> Refer to `src/preprocess/4.rnaseq_preprocess.Rmd`  for VST |
| Monotherapy details | `data/oshu/mono_detail.csv` | Internal reference for drug names |
| Combination details | `data/oshu/combi_detail.csv` | Internal reference for drug names |
| RawDrugResponse | (raw viability)`data/oshu/raw/viability_oshu.csv` <br/> (DSS calculated)`data/oshu/raw/DSS_oshu_raw.csv` | Refer to `src/preprocess/running_dss.R` for DSS calculation from viability |
| DSS monotherapy | `data/oshu/DSS_mono.csv` |  |
| DSS combination | `data/oshu/DSS_combi.csv` |  |
| Global proteomics | `data/oshu/global_proteomics.csv` <br/> (sampleinfo) `data/oshu/global_proteomics_metadata.csv` | Retrieved from [Synapse repository](http://synapse.org/ptrc) from Pacific Northwest National Laboratory |
| DSS of healthy samples | (raw viability) `data/oshu/raw/viability_oshu_healthy.csv` <br/> (DSS calculated) `data/oshu/raw/DSS_oshu_healthy_raw.csv` |  |
| Sample information + DSS (wide format) | `data/oshu/sampleinfo_dss.csv` |  |

## External data

| Data | Files | Remarks |
| --- | --- | --- |
| 246-gene inflammation gene signature | `data/inflammatory_sig_246.csv` <br/> (rds with ENSG-ID and name) `data/inflammatory_sig_genes.rds` | From Lasry et al 2023 supplemetary table 6 |
| scRNAseq data | (raw) `data/scrna/merged_seurat.rds` <br/> (clustered) `data/scrna/all_.rds` <br/> `data/scrna/all_v3.rds` <br/> (metadata) `data/scrna/vangalen_sampleinfo.csv` | From van Galen et al 2019 and [GSE116256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116256) <br/> Refer to `src/preprocess/6.scrna_vg_preprocess.Rmd` for preprocessing and clustering of dataset |
| AML cell hierarchy gene signatures | `data/AMLCellType_Genesets_all.gmt` | From Zeng et al 2022 and corresponding [github](https://github.com/andygxzeng/AMLHierarchies) |
| BeatAML RNAseq data | (raw count) `data/baml/baml_full_raw_count.rds` <br/> (VST) `data/baml/baml_FULL_vst_matrix.rds` <br/> (ensembl gene list - GRCh 37.75) `data/baml/filtered_gene_baml_rnaseq.csv` <br/> (sampleinfo) `data/baml/baml_sampleinfo_rnaseq.csv` | From Tyner et al 2018 and Bottomly et al 2022 <br/> Retrieved from [BeatAML2 Github page](https://biodev.github.io/BeatAML2/) and [Vizome](http://www.vizome.org/additional_figures_BeatAML.html) |
| BeatAML drug response | (drug details) `data/baml/mono_detail.csv` <br/> (viability + DSS - full dataset) `data/baml/DSS_baml_indiv_raw.csv` <br/> (DSS - full dataset, post filtering) `data/baml/DSS_BAML_postfilt.csv` <br/> (DSS - excluding overlapping in bamlcombo) `data/baml/DSS_BAML_woOSHU.csv` <br/> (sampleinfo + DSS) `data/baml/sampleinfo_dss.csv`  | From Tyner et al 2018 and Bottomly et al 2022 <br/> Retrieved from [BeatAML2 Github page](https://biodev.github.io/BeatAML2/)  <br/> Refer to `src/preprocess/running_dss.R` for DSS calculation from viability |
| BeatAML-VenCombo  | (combination details) `data/ven_combi_detail.csv` <br/> (Combination response) `data/ven_combi_response.csv` | From Eide et al 2023 supplementary table |
| FPMTB drug response | (drug details) `data/fpmtb/mono_detail.csv` <br/> (DSS - full MCM dataset) `data/fpmtb/DSS_FPMTB_MCM_conc_full.csv` <br/> (overlapping drugs with BAML) `data/fpmtb/DSS_FPMTB_MCM_conc.csv` <br/> (DSS - common drugs with BAML) `data/fpmtb/conc_common_fp_baml.csv` | From Malani et al 2022 <br/> Retrieved from [zenodo public repository](https://zenodo.org/records/7370747) |
| Clinical response to VenAza (from NUS and oshu cohort) | `data/all_clinical_response.csv` | From Sahib et al 2024 and clinical summary; including samples of patients collected at diagnosis stage and treated with VenAza, with their individual DSS of ven and aza, combination DSS of VenAza, and clinical outcomes  |
| HGNC complete set gene names | `data/hgnc_complete_set.txt` | Obtained from [https://www.genenames.org/download/](https://www.genenames.org/download/); used for gene symbol standardizing |

# Workflow

## Functions

| Functions | Remarks |
| --- | --- |
| `src/function/colour.R` <br/> `src/function/scrna_colour.R` | Color codes used for plots |
| `src/function/DSS_HelperFunctions.R` <br/> `src/function/DSS.R` | Codes taken from [yingjchen/DSS-v2.0](https://github.com/yingjchen/DSS-v2.0) <br/> Modified version that return QC |
| `src/function/ida_predict.R` | Modified from [Alexander-Ling/IDACombo](https://github.com/Alexander-Ling/IDACombo) for use with DSS |
| `src/function/ida_valid.R` | runs ida prediction, but also calculates the actual DSS, HR, IDAComboscore |
| `src/function/loading_data.R` | Standard parameters used for plots |
| `src/function/plot.R` | Plot functions |
| `src/function/functions.R` | Self-defined functions used in analysis |
| `src/function/supp/calculate_sensitivity_score.R` <br/> `src/function/supp/calculate_synergy_score.R` <br/> `src/function/supp/curvefit.R` <br/> `src/function/supp/fit_dose_response_model.R` <br/> `src/function/supp/loewe_fn.R` <br/> `src/function/supp/Reshape_data.R` | Codes taken from [shuyuzheng/synergyfinder](https://github.com/shuyuzheng/synergyfinder) <br/> Used for modelling loewe synergy |

## Analysis

R objects, python objects

| Analysis/Code | Remarks | Files | Figures |
| --- | --- | --- | --- |
| `src/analysis/1.overview_combi.Rmd`  | sample demographics and cohort sample overlap are calculated not exported | Supp Tables <br/> `output/oshu_drug_mapping.csv` <br/> `output/oshu_drug_sens_correlation.csv` <br/> `output/oshu_combinationindex.csv` <br/> <br/> Other tables <br/> `output/supp/oshu_HSA_eval.csv` | Fig 1 <br/> SuppFig 1 <br/> Fig 2 <br/> SuppFig 2 |
| `src/analysis/2.predict_IDA.Rmd`  | Refer to `src/analysis/2a.synergy_calculation.Rmd` for calculation of Bliss and Loewe | Supp Tables <br/> `output/baml_IDA_prediction.csv` <br/> `output/combine_baml_fp_IDA_prediction.csv` <br/> Other tables <br/> `output/supp/baml_woOSHU_IDA_prediction.csv` <br/> `output/supp/fp_full_IDA_prediction.csv` <br/> `output/supp/fp_commondrug_IDA_prediction.csv` <br/> Other rds <br/> (folder) `output/supp/synergy/` | Figure 3 <br/> SuppFig 3 <br/> SuppFig 4 |
| `src/analysis/3.catboost.Rmd`  | Refer to `src/analysis/catboost_py.ipynb` for codes modeling | Supp Tables <br/> `output/all_feature_summary_df.csv` <br/> `output/model_perf_summary_perfold.csv` <br/> `output/shap_summary_df.csv` <br/> `output/top_feature_summary_df.csv` <br/> `output/wilcox_pwc_combi.csv` <br/> `output/spear_numeric_combi.csv` <br/> Other tables <br/> (folder) `output/supp/ML/` <br/> <br/> Other rds <br/> `output/rds/cv_results_all170.pkl` | Figure 4 <br/> SuppFig 5 |
| `src/analysis/4.RNA_clustering.Rmd`  | Refer to `src/preprocess/4.rnaseq_preprocess.Rmd` for VST normalization and cellular hierarchy gene signature calculation using PLAGE-like method | Supp Tables <br/> `output/oshu_GSEA_hall_cluster.csv` <br/> `output/oshu_GSEA_hall_gobp_cluster.csv` <br/> `output/oshu_GSEA_gobp_cluster.csv` <br/> Other rds <br/> `output/rds/oshu_cluster_info.rds` <br/> `output/rds/oshu_seurat_obj.rds` | Figure 5a-d <br/> SuppFig 6a-d |
| `src/analysis/5.WGCNA.Rmd`  |  | Supp Tables  <br/> `output/oshu_wgcna_module_genes.csv` <br/> `output/oshu_wgcna_ora_hall_q25_go_q05.csv` <br/> `output/oshu_corr_ME_trait_results.csv` <br/> Other rds <br/> `output/rds/oshu_wgcna.rds` | Figure 5e <br/> SuppFig e-f |
| `src/analysis/6.scrna.Rmd`  | Refer to `src/preprocess/6.scrna_preprocessing.Rmd` for preprocessing <br/> Refer to `src/preprocess/6.scrna_mast.R` for DEG | Supp Tables <br/> `output/scrna_deg_all.csv` <br/> Other Tables <br/> (folder) `output/supp/scrna/` <br/> Other rds <br/> `output/rds/scrna_selected_obj.rds` <br/> `output/rds/scrna_selected_pseudotime_fitGAM.rds` <br/> `output/rds/scrna_selected_pseudotime_slingshot.rds` | Figure 6 <br/> SuppFig 7 |
| `src/analysis/7.clusterassoc.Rmd`  |  | Supp Tables <br/> `output/oshu_feature_cluster_fisher.csv` <br/> `output/oshu_feature_cluster_pwc.csv` <br/> `output/oshu_drug_cluster_kw.csv` <br/> `output/oshu_drug_cluster_wilcox_inflam.csv` <br/> `output/oshu_drug_cluster_posthoc_pwc_ref.csv` <br/> `output/oshu_genesig_drugresp.csv` | Figure 7a-c |
| `src/analysis/8.integrating_rna_prot.Rmd` |  | Supp Tables <br/> `output/oshu_corr_global_prot_rna_per_gene.csv` <br/> `output/oshu_rna_prot_overlap.csv` <br/> Other Tables <br/> `output/supp/corr_prot_combi_all.csv` <br/> `output/supp/corr_prot_rna_combi.csv` <br/> `output/supp/corr_rna_combi_all.csv` <br/> `output/supp/inflam_DEG_protein.csv` <br/> `output/supp/inflam_DEG_transcript.csv` | Figure 7d |
| `src/analysis/9.biomarker_elastic.Rmd` |  | Supp Tables <br/> `output/oshu_biomarkers.csv` <br/> `output/biomarkers_elasticnet_summary.csv` <br/> Other Tables <br/> `output/supp/biomarkers_elasticnet_long.csv` <br/> `output/supp/elasticnet_perf_summary.csv` <br/> Other rds <br/> `output/rds/biomarkers_elasticnet.rds` | Figure 7e-g <br/> SuppFig 8 |