#' runs ida prediction, but also calculates the actual DSS, HR, IDAComboscore
#' note: HARD CODED

ida_valid <- function(mono, combi, listofcombi){
  IDA_main <- data.frame()
  monodata <- data.frame()
  combidata <- data.frame()
  for (i in 1:nrow(listofcombi)){
    combi_ = combi[combi$Chemical_compound %in% listofcombi[i, 3],] 
    rownames(combi_) <- combi_$Sample_ID
    
    drug1_ = mono[mono$Chemical_compound %in% listofcombi[i, 1],] 
    rownames(drug1_) <- drug1_$Sample_ID
    
    drug2_ = mono[mono$Chemical_compound %in% listofcombi[i, 2],] 
    rownames(drug2_) <- drug2_$Sample_ID
    
    common <- Reduce(intersect, list(rownames(combi_), rownames(drug1_), rownames(drug2_)))
    
    drug1_ = drug1_[common,]; drug2_ = drug2_[common,]; combi_ = combi_[common,]
    res <- data.frame(Drug1 = drug1_$Chemical_compound[1],
                      Drug2 = drug2_$Chemical_compound[1],
                      Drug1Dose = drug1_$Concentration[1],
                      Drug2Dose = drug2_$Concentration[1],
                      Mean_Drug1_s_dss = mean(drug1_$s_dss),
                      Mean_Drug2_s_dss = mean(drug2_$s_dss),
                      Mean_Combo_s_dss = mean(pmax(drug1_$s_dss, drug2_$s_dss))) %>%
      mutate(HR_vs_Drug1 = (1-Mean_Combo_s_dss)/(1-Mean_Drug1_s_dss),
             HR_vs_Drug2 = (1-Mean_Combo_s_dss)/(1-Mean_Drug2_s_dss),
             delta_hazard = pmin((1-Mean_Drug1_s_dss),(1-Mean_Drug2_s_dss)) - (1-Mean_Combo_s_dss),
             HR_C_over_Mbest = pmax(HR_vs_Drug1, HR_vs_Drug2),
             IDA_Comboscore = delta_hazard - delta_hazard * HR_C_over_Mbest,
             Actual_s_dss = mean(combi_$s_dss),
             Actual_CombiDose = combi_$Concentration[1],
             Actual_HR_vs_Drug1 = (1-Actual_s_dss)/(1-Mean_Drug1_s_dss),
             Actual_HR_vs_Drug2 = (1-Actual_s_dss)/(1-Mean_Drug2_s_dss),
             Actual_delta_hazard = pmin((1-Mean_Drug1_s_dss),(1-Mean_Drug2_s_dss)) - (1-Actual_s_dss),
             Actual_HR_C_over_Mbest = pmax(Actual_HR_vs_Drug1, Actual_HR_vs_Drug2),
             Actual_IDA_Comboscore = ifelse(Actual_HR_C_over_Mbest > 1, NA, Actual_delta_hazard - Actual_delta_hazard * Actual_HR_C_over_Mbest),
             n = length(common)) %>%
      dplyr::select(-HR_C_over_Mbest, -Actual_HR_C_over_Mbest)
    
    IDA_main <- rbind(IDA_main, res)
    monodata <- rbind(monodata, rbind(drug1_, drug2_))
    combidata <- rbind(combidata, combi_)
  }
  
  # collects the values that were involved in the prediction
  monodata <- monodata[!duplicated(monodata[,1:3]),]
  combidata <- combidata[!duplicated(combidata[,1:3]),]
  
  ## returning outputs
  return(list(IDA_main=IDA_main, monodata=monodata, combidata=combidata))
}