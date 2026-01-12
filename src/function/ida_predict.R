#' faster prediction for all possible pairwise combinations
#' edited from Alexander-Ling/IDACombo
#' date: 5 April

#' predict all possible pairwise combinations from the data
ida_predict <-  function(Monotherapy_Data, 
                         output="IDA_prediction.csv", 
                         Efficacy_Column = "s_dss", 
                         Cell_Line_Name_Column = "Sample_ID",
                         Drug_Name_Column = "Chemical_compound",
                         Drug_Concentration_Column = "Concentration", 
                         min.sample = 3, 
                         retun.rank = TRUE,
                         rm.idacomboscore = TRUE) {
  require(tidyverse)
  
  if (typeof(Monotherapy_Data)=="character") {
    Monotherapy_Data <- read.csv(Monotherapy_Data)
  }
  
  Data <- Monotherapy_Data[,c(Cell_Line_Name_Column, Drug_Name_Column, Drug_Concentration_Column, Efficacy_Column)]
  
  # rm(Monotherapy_Data)
  
  Data <- Data[complete.cases(Data),] # ensure no NA
  colnames(Data) <- c("CellLine", "Drug", "Conc", "Efficacy")
  Data$Label <- paste0(Data$Drug, "_", Data$Conc)
  Data <- arrange(Data, Label)
  
  min.sample = max(min.sample,3)
  
  list.drug.conc <- unique(Data[,c(2,3,5)])
  list.drug <- unique(Data[,2])
  
  # averages if duplicates exists
  # dataframe for all unique Drug_Conc column
  Data <-reshape2::dcast(Data, CellLine ~ Label, value.var = "Efficacy", fun.aggregate = function(x){if (length(x)==0) {NaN} else {mean(x, na.rm = T)} })
  
  Data <- Data[,c("CellLine", list.drug.conc$Label)] # sorting the column in order
  
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0, max = length(list.drug), style = 3, width = 50, char = "=")
  
  res=c()
  i=1 # tracks which unique drug
  count=1 # tracks the column (Drug_Conc) that has been used 
  
  # loops through each drug - sorted in order
  while (i < length(list.drug)) {
    drug <- list.drug.conc[list.drug.conc$Drug == list.drug[i],] # extracts all Drug_Conc of the drug
    count = count + nrow(drug) 
    b = as.matrix(Data[,(count+1):ncol(Data)]) # matrix of drug after the selected drug column
    
    drugs = list.drug.conc[count:nrow(list.drug.conc), ] # details of b
    
    # loops through each Drug_Conc of the drug
    for (j in 1:nrow(drug)) {
      a = drug[j,"Label"]; a_drug = drug[j,"Drug"]; a_conc = drug[j,"Conc"]
      c = c()
      a_col <- Data[,a]
      n = ncol(b)
      
      # create matrix same ncol as b
      if(n==1){ # base data
        drug1_m <- as.matrix(a_col)
      } else {
        drug1_m <- matrix(replicate(n, a_col), nrow = length(a_col))
      }
      
      b_combi_efficacy <- pmax(drug1_m, b)
      mask <- is.nan(b_combi_efficacy) # to remove uncommon sample
      drug1_m[mask] = NA
      b[mask] = NA
      b_combi_efficacy[mask] = NA
      
      c$Drug1 <- rep(a_drug, n)
      c$Drug2 <- drugs$Drug
      c$Drug1Dose <- rep(a_conc, n)
      c$Drug2Dose <- drugs$Conc
      c$Mean_Drug1_Efficacy <- as.numeric(colMeans(drug1_m, na.rm = T))
      c$Mean_Drug2_Efficacy <- as.numeric(colMeans(b, na.rm = T))
      c$Mean_Combo_Efficacy <- as.numeric(colMeans(b_combi_efficacy, na.rm = T))
      c$HR_vs_Drug1 <- (1-c$Mean_Combo_Efficacy)/(1-c$Mean_Drug1_Efficacy)
      c$HR_vs_Drug2 <- (1-c$Mean_Combo_Efficacy)/(1-c$Mean_Drug2_Efficacy)
      c$delta_hazard <- pmin((1-c$Mean_Drug1_Efficacy),(1-c$Mean_Drug2_Efficacy)) - (1-c$Mean_Combo_Efficacy)
      c$HR_C_over_Mbest <- pmax(c$HR_vs_Drug1, c$HR_vs_Drug2)
      c$IDA_Comboscore <- c$delta_hazard - c$delta_hazard * c$HR_C_over_Mbest
      c$n <- as.numeric(colSums(!is.na(b_combi_efficacy)))
      
      c <- as.data.frame(c) %>%
        filter(n>=min.sample) # at least samples
      
      res <- rbind(res,c)
    }
    i=i+1
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  res <- res %>%
    rowwise() %>%
    mutate(Drug1_temp=pmin(Drug1, Drug2), Drug2_temp=pmax(Drug1, Drug2),
           Combination = paste(Drug1_temp,"-",Drug2_temp)) %>%
    dplyr::select(-Drug1_temp, -Drug2_temp) %>%
    arrange(desc(IDA_Comboscore)) %>%
    ungroup()
    
  
  if (retun.rank) {
    q75 = quantile.n(Monotherapy_Data[[Efficacy_Column]], n=.75)
    res <- res %>%
      mutate(rank = ifelse(Mean_Combo_Efficacy > q75 & IDA_Comboscore > 0, -1, NA)) %>% 
      group_by(rank) %>%
      mutate(rank = as.integer(ifelse(rank == -1, rank(-IDA_Comboscore), NA))) %>%
      arrange(rank, desc(IDA_Comboscore)) 
  }
  
  if (rm.idacomboscore) {
    res <- res[res$IDA_Comboscore > 0,]
  }
  
  colnames(res) <- gsub("Efficacy", Efficacy_Column, colnames(res))
  print(paste("Saving...", output))
  print(paste0("Using unique drug: ", ncol(Data)-1, "; sample: ", nrow(Data)))
  
  write.csv(res, output, row.names = F)
  return(res)
}

#' outlier removal + quantile
quantile.n <-  function(data, n=0.75){
  q = quantile(data, probs = c(.05, .95), na.rm = T)
  new = data[which(data > q[[1]] & data < q[[2]])]
  q = quantile(new, probs = n, na.rm = T)
  print(paste("Threshold at",n, ":", q))
  return(q)
}