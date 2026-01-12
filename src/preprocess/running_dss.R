## Calculating DSS values from viability (OSHU)
## date: "18 March"

source("src/function/DSS.R")
source("src/function/HelperFunctions.R")
lapply(c("matrixStats","reshape","reshape2", "scales", "drc", "caTools", "data.table", "MESS", "tidyverse"), library, character.only = T)

input.dir <- "data/oshu/raw"
output.dir <- "data/oshu/raw"

df_dose.responses <- read.csv(paste0(input.dir, "viability_oshu.csv")) %>%
  filter(drug != "Control")

sampleinfo1 <- unique(df_dose.responses[,c(1:2,6,11)]) 
sampleinfo1$sampleid <- seq(1:nrow(sampleinfo1))
sampleinfo1 <- sampleinfo1 %>%
  mutate(sampleid = paste0("P", sampleid))

df_dose.responses <- left_join(df_dose.responses, sampleinfo1) %>%
  transmute(DRUG = drug,
            CONCENTRATION_nM = 1000*concentration,
            SCREEN_NAME = sampleid,
            CELL_VIABILITY = normalized_probit_result/100)

df_dose.responses.list <- DOSE_RESPONSE_PROCESS(df_dose.responses, viability = T)

df.metrics <- CALC_METRICS.2(dose_responses = df_dose.responses.list[[1]], dose_responses_grouped = df_dose.responses.list[[2]], graph = F, failed = F)

df.metrics <- left_join(df.metrics, sampleinfo1, by=c("Patient.num" = "sampleid")) 

write.csv(df.metrics, file = paste0(output.dir,"DSS_oshu_raw.csv"), row.names = F)

# -------------------------------------------------------------------------------

# check if there is replicates
df_dose.responses <- read.csv(paste0(input.dir, "viability_oshu_healthy.csv")) %>%
  filter(drug != "Control")

sampleinfo1 <- unique(df_dose.responses[,c(1:2,6,11)]) 
sampleinfo1$sampleid <- seq(1:nrow(sampleinfo1))
sampleinfo1 <- sampleinfo1 %>%
  mutate(sampleid = paste0("P", sampleid))

df_dose.responses <- left_join(df_dose.responses, sampleinfo1) %>%
  transmute(DRUG = drug,
            CONCENTRATION_nM = 1000*concentration,
            SCREEN_NAME = sampleid,
            CELL_VIABILITY = normalized_probit_result/100)

df_dose.responses.list <- DOSE_RESPONSE_PROCESS(df_dose.responses, viability = T)

df.metrics <- CALC_METRICS.2(dose_responses = df_dose.responses.list[[1]], dose_responses_grouped = df_dose.responses.list[[2]], graph = F, failed = F)

df.metrics <- left_join(df.metrics, sampleinfo1, by=c("Patient.num" = "sampleid")) 

write.csv(df.metrics, file = paste0(output.dir,"DSS_oshu_healthy_raw.csv"), row.names = F)

#------------------------------------------------------------------------------------------
## Processing the DSS values (OSHU)
library(tidyverse)

combidet <- read.csv("data/oshu/combi_detail.csv")
rownames(combidet) <- combidet$Combination
drugdet <- read.csv("data/oshu/mono_detail.csv")
rownames(drugdet) <- drugdet$Chemical_compound
sampleinfo <- read.csv("data/oshu/clinical_summary.csv")[,c("patientId","labId","sampleid","Sample_ID","dbgap_subject_id","dbgap_rnaseq_sample","dbgap_dnaseq_sample","Status")] %>%
  dplyr::rename(lab_id = labId, patient=patientId, Disease.status=Status)

# changing the labels for the samples
df <- read.csv(paste0(output.dir,"DSS_oshu_raw.csv")) %>% # 35456
  mutate(Concentration = paste0(as.character(Min.Conc.tested), "_",sprintf("%.f", Max.Conc.tested)),
         s_dss = DSS/50) %>%
  dplyr::select(3,5,20,24:28,30:31) %>%
  inner_join(., sampleinfo)

## EDA
table(df$believe_DSS)
table(df$cond[df$believe_DSS==FALSE])
# FALSE  TRUE
# 587 34869
# 2
# 587

## duplicated
# check <- df_mono[df_mono$believe_DSS,] %>% mutate(label=paste(Sample_ID, Chemical_compound))
# check <- check[check$label %in% check$label[duplicated(check$label)],]

## for single drugs
df_mono <- df %>%
  filter(drug %in% drugdet$Chemical_compound) %>%
  filter(drug != "Quizartinib - NVP-TAE684")
df_mono$Chemical_compound <- drugdet[df_mono$drug, "Drug"]
df_mono <- df_mono %>% #12705
  filter(believe_DSS) %>%
  group_by(Sample_ID, Disease.status, Chemical_compound, Concentration,
           dbgap_subject_id, dbgap_dnaseq_sample, dbgap_rnaseq_sample) %>% #12582
  reframe(DSS=mean(DSS), s_dss=mean(s_dss)) %>% #10065
  dplyr::select(1:4, 8, 9, 5:7)

write.csv(df_mono, "data/oshu/DSS_mono.csv", row.names = F)

## for combinations
combi <- df %>%
  filter(drug %in% combidet$Combination) %>%
  filter(drug != "Quizartinib - NVP-TAE684 - Arsenic Trioxide")
combi <- cbind(combi, combidet[combi$drug, c(1,2,4)])
combi <- combi %>% #22717
  filter(believe_DSS) %>%
  group_by(Sample_ID, Disease.status, Drug1, Drug2, combo, Concentration,
           dbgap_subject_id, dbgap_dnaseq_sample, dbgap_rnaseq_sample) %>% #22253
  reframe(DSS=mean(DSS), s_dss=mean(s_dss)) %>% #22150
  dplyr::select(1:6, 10, 11, 7:9) %>%
  dplyr::rename(Combination = combo)

write.csv(combi, "data/oshu/DSS_combi.csv", row.names = F)

#-------------------------------------------------
drugs <- data.frame(OGname = combidet$Combination, Drug = combidet$combo) %>%
  rbind(data.frame(OGname = drugdet$Chemical_compound, Drug = drugdet$Drug))
drugs <- drugs[-227,]
rownames(drugs) <- drugs$OGname

## for healthy
df <- read.csv(paste0(output.dir,"DSS_oshu_healthy.csv")) %>% # 272
  mutate(Concentration = paste0(as.character(Min.Conc.tested), "_",sprintf("%.f", Max.Conc.tested)),
         s_dss = DSS/50) %>%
  dplyr::select(3,20,24:28,30:31) %>%
  inner_join(., sampleinfo)

df$Chemical_compound <- drugs[df$drug, "Drug"]
df <- df %>%
  filter(believe_DSS) %>%
  group_by(Sample_ID, Disease.status, Chemical_compound, Concentration,
           dbgap_subject_id, dbgap_dnaseq_sample, dbgap_rnaseq_sample) %>% #264
  reframe(DSS=mean(DSS), s_dss=mean(s_dss)) %>% #264
  dplyr::select(1:4, 8, 9, 5:7)

write.csv(df, "data/oshu/DSS_healthy.csv", row.names = F)