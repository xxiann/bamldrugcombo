packages.required <- packages.required <- c("matrixStats","reshape","reshape2","drc","data.table", "tidyverse")

lapply(packages.required, library, character.only = T)

source("src/function/DSS.R")
source("src/function/HelperFunctions.R")

#' return DSS 
#' undone
#' 
CALC_DSS_single <- compiler::cmpfun(function(inhibition, dose, drug_name, DSS_typ=2, return_DSS_only = FALSE, path = "", graph=FALSE, failed=FALSE){
  
  tryCatch({
    
    inhibition2 = inhibition
    viability2 = 100 - inhibition2; # with 2 used for ploting of real values.
    # if there are identical values in inhibition, add a bit noise
    if(all(inhibition <= 0)) inhibition <- rep(0, length(inhibition))
    if(any(duplicated(inhibition))) inhibition <- seq(from = 0, length.out = length(inhibition), by = 0.01) + inhibition;
    
    viability = 100-inhibition; believe_ = T;
    
    #combine the data and sort by dose.
    mat_tbl <- data.frame(inhibition,dose,logconc = log10(dose),viability, inhibition2, viability2)
    mat_tbl <- mat_tbl[order(mat_tbl[,"dose"]),]

    
    if(nrow(mat_tbl) <= 3 || (length(unique(mat_tbl$dose)) <= 3)){
      
      print("Less than 3 rows... skipping...")
      NULL
    } else {
      
      #############################
      #############    IC50
      
      estimate_param <- tryCatch({drm(inhibition ~ logconc, data = mat_tbl, fct = LL.4(fixed = c(NA, NA, NA,NA),names = c("SLOPE","MIN","MAX","IC50")),logDose=10,control = drmc(errorm = F))},
                                 warning=function(w){drm(inhibition ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)},
                                 error=function(e){drm(inhibition ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)})
      # (extract and name coefficients)
      coef_estim <- coef(estimate_param); names(coef_estim) <- c("SLOPE","MIN","MAX","IC50")
      # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696819/
      coef_estim["SLOPE"] <- coef_estim["SLOPE"]*-1
      
      # if curve decreases or IC50 is higher than max (i.e. IC50 is "outlier"), set IC50 to max conc.
      coef_estim["IC50"] <- ifelse(coef_estim["MAX"]<=coef_estim["MIN"] | coef_estim["IC50"]>max(mat_tbl$dose,na.rm=T), max(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      # if IC50 is less than 0 set it to min. conc. and if even min. conc. < 0, then set IC50 to mean of all conc.
      coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,min(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,mean(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      # similar to previous step but now compare log10(IC50) with log(min. conc.).
      coef_estim["IC50"] <- log10(coef_estim["IC50"])
      coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<min(mat_tbl$logconc),max(mat_tbl$logconc),coef_estim["IC50"])
      # if all inhib. < 0 set IC50 to max. log. conc !!!!! not obvious why!
      coef_estim["IC50"] <- ifelse(all(mat_tbl$inhibition<0),max(mat_tbl$logconc,na.rm=T),coef_estim["IC50"])
      #(Trying to fix curves that need outlier kickout)
      coef_estim["MIN"] <- 0; coef_estim["MAX"] <- max(mat_tbl$inhibition,na.rm=T)
      #(Fix off minimums) Find lowest inhibition value. If it is not in (0:100), fix it whether to 0 or 99.
      min_lower <- ifelse(min(mat_tbl$inhibition,na.rm=T) > 0,min(mat_tbl$inhibition,na.rm=T),0)
      min_lower <- ifelse(min_lower >= 100,99,min_lower)
      #similar to previous step but for MAX
      coef_estim["MAX"] <- ifelse(coef_estim["MAX"]>100,100,coef_estim["MAX"])
      coef_estim["MAX"] <- ifelse(coef_estim["MAX"]<0,100,coef_estim["MAX"])
      #max_lower and max_upper - lower and upper bounds for 'nl2sol' algorithm in nonlinear least-squares
      max_lower <- ifelse(max(mat_tbl$inhibition,na.rm=T)>100,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T))
      max_lower <- ifelse(max_lower < 0,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T))
      max_lower <- ifelse(max_lower < 0,0,max_lower)
      max_lower <- ifelse(max_lower > 100,100,max_lower)
      #(Fix upper maximum for negative slopes)
      run_avg <- caTools::runmean(mat_tbl$inhibition, 10)
      max_upper <- ifelse(any(run_avg[-nrow(mat_tbl)]>run_avg[nrow(mat_tbl)]),max(mat_tbl$inhibition[run_avg>run_avg[nrow(mat_tbl)]]),coef_estim["MAX"])
      max_upper <- ifelse(any(mat_tbl$inhibition > max_upper),mean(mat_tbl$inhibition[mat_tbl$inhibition > max_upper])+5,max_upper)
      max_upper <- ifelse(max_upper < 0,coef_estim["MAX"],max_upper)
      max_upper <- ifelse(max_upper > 100,100,max_upper) #coef_estim["MAX"]
      max_upper <- ifelse(max_lower > max_upper,coef_estim["MAX"],max_upper)
      # left it as it was, just rewritten a bit (ALEKS). not clear how values 25, 60 and 5 are chosen.
      mean_inh_last = mean(tail(mat_tbl$inhibition,2),na.rm=T)
      if(mean_inh_last < 60) {
        if(mean_inh_last > 25) coef_estim["IC50"] <- mean(mat_tbl$logconc,na.rm=T)
        else if(mean_inh_last < 25) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)}
      if(mean(mat_tbl$inhibition[1:3],na.rm=T)<5) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)
      #add a bit of positive noise to MAX if it is the same as MIN.
      if(unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) coef_estim["MAX"] <- coef_estim["MAX"] + 0.001
      
      #adaptive nonlinear Least-Squares algorithm NL2SOL to estimate parameters.
      nls_result_ic50_old <- function(){
        tryCatch({
          nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]), lower=list(SLOPE=0.1,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)), upper=list(SLOPE=2.5,MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
        }, error = function(e) {
          
          # allows higher residual sum-of-squares
          minpack.lm::nlsLM(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl,
                            start=list(SLOPE=1, MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]),
                            lower=c(SLOPE=0.1, MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),
                            upper=c(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)))
        })
      }
      
      # IC50 first
      nls_result_ic50 <- nls_result_ic50_old();
      
      # IC50 second
      nls_result_ic50_2 <- tryCatch({
        # allows higher residual sum-of-squares
        nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",  start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]], IC50=median(mat_tbl$logconc)),lower=list(SLOPE=0.1,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),upper=list(SLOPE=2.5,MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
      },warning = function(w) {
        nls_result_ic50_old()
      },error = function(e) {
        nls_result_ic50_old()
      })
      
      #element (4, 4) is zero, so the inverse cannot be computed
      nls_result_ic50 = tryCatch({summary(nls_result_ic50); nls_result_ic50},error=function(e){nls_result_ic50_2})
      
      #Calculate the standard error scores
      sumIC50 = list(summary(nls_result_ic50), summary(nls_result_ic50_2))
      
      ic50std_resid <- round(sqrt(sum((sumIC50[[1]]$residuals)^2)/(length(sumIC50[[1]]$residuals)-1)),1);
      ic50std_resid2 <- round(sqrt(sum((sumIC50[[2]]$residuals)^2)/(length(sumIC50[[2]]$residuals)-1)),1);
      
      # continue with the best
      switch_ = which.min(c(ic50std_resid, ic50std_resid2))
      nls_result_ic50 = list(nls_result_ic50, nls_result_ic50_2)[[switch_]]
      
      
      #if SLOPE <= 0.2, decrease IC50, change lower bound for SLOPE to 0.1 and repeat.
      if(coef(nls_result_ic50)["SLOPE"] <= 0.2)
      {
        if(mean_inh_last > 60)
          coef_estim["IC50"] <- min(mat_tbl$logconc,na.rm=T)
        nls_result_ic50 <- nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),lower=list(SLOPE=0.1,MIN=min_lower,MAX=max_lower,IC50=min(mat_tbl$logconc)),upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
      }
      
      #Calculate the standard error scores
      sumIC50 = summary(nls_result_ic50);
      ic50std_Error <- sumIC50$coefficients["IC50","Std. Error"]; #tec50std_Error <- sumTEC50$coefficients["TEC50","Std. Error"]
      ic50std_resid <- round(sqrt(sum((sumIC50$residuals)^2)/(length(sumIC50$residuals)-1)),1);
      max_signal <- max(mat_tbl$dose,na.rm=T); min_signal <- min(mat_tbl$dose,na.rm=T)
      
      #############################
      #############   Final modification & STD error
      
      #prepare final data and convert IC50 back from log scale (inverse)
      coef_ic50 <- coef(nls_result_ic50)[c("IC50", "SLOPE","MAX","MIN")]; coef_ic50["IC50"] <- 10^coef_ic50["IC50"]
      #(Fix ic50 for curves in wrong direction)
      coef_ic50["IC50"] <- ifelse(coef_ic50["SLOPE"]<0,max_signal,coef_ic50["IC50"])
      #(Fix based on MAX)
      coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<0,max_signal,coef_ic50["IC50"])
      coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<10,max_signal,coef_ic50["IC50"])
      coef_ic50["MAX"] <- ifelse(coef_ic50["MAX"]<0,0,coef_ic50["MAX"])
      #(Fix over sensitive drugs)
      coef_ic50["IC50"] <- ifelse(all(c(max(mat_tbl$inhibition,na.rm=T),min(mat_tbl$inhibition,na.rm=T))>50),min_signal,coef_ic50["IC50"])
      
      # for ploting
      x <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc), length=100)
      yic <- predict(nls_result_ic50, data.frame(logconc=x))
      auc <- MESS::auc(x,yic)
      
      ##average replicates
      mat_tblCp <- mat_tbl[, c("inhibition", "dose")]
      cols_ <- colnames(mat_tblCp)[!grepl("inhibition", colnames(mat_tblCp))] # columns which should be equal to average PI
      X <- as.data.table(mat_tblCp)
      mat_tblCp <- as.data.frame(X[,list(inhibition = mean(inhibition)),cols_], stringAsFactors = !1)
      
      
      perInh <- t(matrix(mat_tblCp[,"inhibition"],dimnames=
                           list(paste0(rep("D", length(mat_tblCp[,"inhibition"])), 1:length(mat_tblCp[,"inhibition"])))))
      
      
      #############################
      #############    DSS
      
      dss_cond <- dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=as.integer(DSS_typ));
      
      dss_score <- round(as.numeric(dss_cond[[1]]),1);
      
      coef_ic50 <- c(coef_ic50,Min.Conc.tested=min_signal,Max.Conc.tested=max_signal,IC50_std_error=ic50std_Error)
      
      
      #dataframe for IC50
      IC50_dataframe <- data.frame(ID=drug_name,DRUG_NAME=drug_name,ANALYSIS_NAME="IC50", t(as.matrix(coef_ic50)), perInh,
                                   GRAPH=NA, DSS = as.numeric(dss_score), sDSS = paste(round(sumIC50$residuals,2), collapse = ",", sep = ","), SE_of_estimate = as.numeric(ic50std_resid),AUC=auc)
      
      #round by 2 dex. all the numeric colums
      numeric_cols <- sapply(IC50_dataframe, is.numeric); IC50_dataframe[,numeric_cols] <- round(IC50_dataframe[,numeric_cols],1)
      

      ####################################
      ## check believe (QC for dose response curve fit)
      resid_ = as.numeric(sumIC50$residuals); resp_ = as.numeric(perInh); cond = 0; outlier = floor(length(resid_ )/2) # min tolerable number
      
      # if(ic50std_resid>70) {believe_ = !1; cond = 1} # when SE  is very large
      # if(sum(abs(resid_)>30)> outlier) {believe_ = !1; cond = 2};  # when ard half the points are outlying
      # if( sum((abs(resid_)>30) == (resp_ < 2)) == length(resid_ ) ) {believe_ = !0; cond = 2.1};  # if the outliers are below 0
      # if(sum(resp_ < 2) > outlier) {believe_ = !0; cond = 3} # when majority of the points are below 0
      # if(dss_score == 0) {believe_ = !0; cond = 3.1} # when DSS score is zero 
      
      if(ic50std_resid>50) {believe_ = !1; cond = 1} # set FALSE when SE  is very large
      if(sum(abs(resid_)>25)> outlier) {believe_ = !1; cond = 2} # set FALSE when ard half the points are outlying
      if( sum((abs(resid_)>25) == (resp_ < 2)) == length(resid_ ) ) {believe_ = !0; cond = 2.1}  # set TRUE if the outliers are below 0
      if(sum(abs(resid_)>15)<= outlier) {believe_ = !0; cond = 2.2} # even when resid is high, allowing less than half outliers
      
      if(sum(resp_ < 2) > outlier) {believe_ = !0; cond = 3} # when majority of the points are below 0
      if(dss_score == 0) {
        believe_ = !0
        cond = ifelse(dss_cond[[2]]!=-1, dss_cond[[2]], 3.1)
      } # when DSS score is zero 
      if(coef_ic50["IC50"]==min_signal){cond = 4} # when IC50 == Min.conc.tested
      
      if (graph | (failed & !believe_)){
        #plot IC50
        mat_tbl$inhibition = xpr_tbl$inhibition_percent; # if we have all values < 0, they will be replaced
        mat_tbl$viability = 100 - mat_tbl$inhibition;  # we are replacing them back here.
        icpl <- ggplot2::ggplot(mat_tbl, aes(logconc, inhibition2)) + 
          scale_x_continuous(breaks=mat_tbl$logconc,labels=round(mat_tbl$dose,2)) +
          geom_point(color = "blue", size = 2.8) + 
          geom_line(data = data.frame(x = x, y = yic), aes(x, yic), color="blue", size = 0.8) +
          geom_vline(xintercept = log10(coef_ic50["IC50"]), colour="grey", size = 0.8) + 
          ggtitle(paste0(strtrim(drug_name, 30)," (dss:",dss_score,")\n"))+
          theme_bw() + labs(y = "% inhibition", x = "conc(nM)")  +  ylim(-25, 125) +
          geom_text(mapping=aes(x2,y2,label = text2), data=data.frame(x2=log10(coef_ic50["IC50"])*0.95, y2=115, text2="IC50"), color="grey", parse=T) +
          theme(plot.background = element_rect(fill = ifelse(failed, "#ffecf2", "transparent"),colour = NA),
                panel.background =element_rect(fill = "transparent",colour = NA), 
                plot.title = element_text(hjust = 0.5, size = 12.5),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
          )
        graphics.off()
        filename_ = file.path(path, paste0(drug_name,"_IC50_curve_drug.png"))
        ggsave(filename_, icpl, height = 3, width = 3.7, bg= 'white', dpi = 300, create.dir = TRUE)

      } 
      
      if(return_DSS_only){
        return(dss_score)
      } else {return(data.frame(DSS=dss_score/50, believe=believe_, condition=cond))}
      
    }
  }, error = function(e) {
    print(e);
  })
})



#' calculates for dose response function
#' modified from CALC_IC50_EC50_DSS
#' can return either the nls function or the parameter
CALC_DR <- compiler::cmpfun(function(xpr_tbl, return_param=T, graph=F){
  
  tryCatch({
    
    inhibition = inhibition2 <- xpr_tbl$inhibition_percent; viability2 = 100 - inhibition2; # with 2 used for ploting of real values.
    # if there are identical values in inhibition, add a bit noise
    if(all(inhibition <= 0)) inhibition <- rep(0, length(inhibition))
    if(any(duplicated(inhibition))) inhibition <- seq(from = 0, length.out = length(inhibition), by = 0.01) + inhibition;
    
    viability = 100-inhibition; believe_ = T;
    
    # extract concentrations, unique drug names and product ids for wells with drugs in current plate
    dose <- as.numeric(xpr_tbl$Concentration)
    drug_name <- unique(as.character(xpr_tbl$ID))
    product_id <- unique(as.character(xpr_tbl$ID))
    
    #combine the data and sort by dose.
    mat_tbl <- data.frame(inhibition,dose,logconc = log10(dose),viability, inhibition2, viability2)
    mat_tbl <- mat_tbl[order(mat_tbl[,"dose"]),]

    
    if(nrow(mat_tbl) <= 3 || (length(unique(mat_tbl$dose)) <= 3)){
      print("Less than 3 rows... skipping...")
      NULL
    } else {
      
      #############################
      #############    IC50
      
      estimate_param <- tryCatch({drm(inhibition ~ logconc, data = mat_tbl, fct = LL.4(fixed = c(NA, NA, NA,NA),names = c("SLOPE","MIN","MAX","IC50")),logDose=10,control = drmc(errorm = F))},
                                 warning=function(w){drm(inhibition ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)},
                                 error=function(e){drm(inhibition ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)})
      # (extract and name coefficients)
      coef_estim <- coef(estimate_param); names(coef_estim) <- c("SLOPE","MIN","MAX","IC50")
      # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696819/
      coef_estim["SLOPE"] <- coef_estim["SLOPE"]*-1
      
      # if curve decreases or IC50 is higher than max (i.e. IC50 is "outlier"), set IC50 to max conc.
      coef_estim["IC50"] <- ifelse(coef_estim["MAX"]<=coef_estim["MIN"] | coef_estim["IC50"]>max(mat_tbl$dose,na.rm=T), max(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      # if IC50 is less than 0 set it to min. conc. and if even min. conc. < 0, then set IC50 to mean of all conc.
      coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,min(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,mean(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      
      
      # similar to previous step but now compare log10(IC50) with log(min. conc.).
      coef_estim["IC50"] <- log10(coef_estim["IC50"])
      coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<min(mat_tbl$logconc),max(mat_tbl$logconc),coef_estim["IC50"])
      # if all inhib. < 0 set IC50 to max. log. conc !!!!! not obvious why!
      coef_estim["IC50"] <- ifelse(all(mat_tbl$inhibition<0),max(mat_tbl$logconc,na.rm=T),coef_estim["IC50"])
      
      
      #(Trying to fix curves that need outlier kickout)
      coef_estim["MIN"] <- 0; coef_estim["MAX"] <- max(mat_tbl$inhibition,na.rm=T)
      #(Fix off minimums) Find lowest inhibition value. If it is not in (0:100), fix it whether to 0 or 99.
      min_lower <- ifelse(min(mat_tbl$inhibition,na.rm=T) > 0,min(mat_tbl$inhibition,na.rm=T),0)
      min_lower <- ifelse(min_lower >= 100,99,min_lower)
      #similar to previous step but for MAX
      coef_estim["MAX"] <- ifelse(coef_estim["MAX"]>100,100,coef_estim["MAX"])
      coef_estim["MAX"] <- ifelse(coef_estim["MAX"]<0,100,coef_estim["MAX"])
      #max_lower and max_upper - lower and upper bounds for 'nl2sol' algorithm in nonlinear least-squares
      max_lower <- ifelse(max(mat_tbl$inhibition,na.rm=T)>100,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T))
      max_lower <- ifelse(max_lower < 0,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T))
      max_lower <- ifelse(max_lower < 0,0,max_lower)
      max_lower <- ifelse(max_lower > 100,100,max_lower)
      #(Fix upper maximum for negative slopes)
      run_avg <- caTools::runmean(mat_tbl$inhibition, 10)
      max_upper <- ifelse(any(run_avg[-nrow(mat_tbl)]>run_avg[nrow(mat_tbl)]),max(mat_tbl$inhibition[run_avg>run_avg[nrow(mat_tbl)]]),coef_estim["MAX"])
      max_upper <- ifelse(any(mat_tbl$inhibition > max_upper),mean(mat_tbl$inhibition[mat_tbl$inhibition > max_upper])+5,max_upper)
      max_upper <- ifelse(max_upper < 0,coef_estim["MAX"],max_upper)
      max_upper <- ifelse(max_upper > 100,100,max_upper) #coef_estim["MAX"]
      max_upper <- ifelse(max_lower > max_upper,coef_estim["MAX"],max_upper)
      # left it as it was, just rewritten a bit (ALEKS). not clear how values 25, 60 and 5 are chosen.
      mean_inh_last = mean(tail(mat_tbl$inhibition,2),na.rm=T)
      if(mean_inh_last < 60) {
        if(mean_inh_last > 25) coef_estim["IC50"] <- mean(mat_tbl$logconc,na.rm=T)
        else if(mean_inh_last < 25) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)}
      if(mean(mat_tbl$inhibition[1:3],na.rm=T)<5) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)
      #add a bit of positive noise to MAX if it is the same as MIN.
      if(unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) coef_estim["MAX"] <- coef_estim["MAX"] + 0.001
      
      #adaptive nonlinear Least-Squares algorithm NL2SOL to estimate parameters.
      nls_result_ic50_old <- function(){
        tryCatch({
          nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]), lower=list(SLOPE=0.1,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)), upper=list(SLOPE=2.5,MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
        }, error = function(e) {
          
          # allows higher residual sum-of-squares
          minpack.lm::nlsLM(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl,
                            start=list(SLOPE=1, MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]),
                            lower=c(SLOPE=0.1, MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),
                            upper=c(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)))
        })
      }
      
      # IC50 first
      nls_result_ic50 <- nls_result_ic50_old();
      
      # IC50 second
      nls_result_ic50_2 <- tryCatch({
        # allows higher residual sum-of-squares
        nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",  start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]], IC50=median(mat_tbl$logconc)),lower=list(SLOPE=0.1,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),upper=list(SLOPE=2.5,MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
      },warning = function(w) {
        nls_result_ic50_old()
      },error = function(e) {
        nls_result_ic50_old()
      })
      
      #element (4, 4) is zero, so the inverse cannot be computed
      nls_result_ic50 = tryCatch({summary(nls_result_ic50); nls_result_ic50},error=function(e){nls_result_ic50_2})
      
      #Calculate the standard error scores
      sumIC50 = list(summary(nls_result_ic50), summary(nls_result_ic50_2))
      
      ic50std_resid <- round(sqrt(sum((sumIC50[[1]]$residuals)^2)/(length(sumIC50[[1]]$residuals)-1)),1);
      ic50std_resid2 <- round(sqrt(sum((sumIC50[[2]]$residuals)^2)/(length(sumIC50[[2]]$residuals)-1)),1);
      
      # continue with the best
      switch_ = which.min(c(ic50std_resid, ic50std_resid2))
      nls_result_ic50 = list(nls_result_ic50, nls_result_ic50_2)[[switch_]]
      
      
      #if SLOPE <= 0.2, decrease IC50, change lower bound for SLOPE to 0.1 and repeat.
      if(coef(nls_result_ic50)["SLOPE"] <= 0.2)
      {
        if(mean_inh_last > 60)
          coef_estim["IC50"] <- min(mat_tbl$logconc,na.rm=T)
        nls_result_ic50 <- nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),lower=list(SLOPE=0.1,MIN=min_lower,MAX=max_lower,IC50=min(mat_tbl$logconc)),upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
      }
      
      #Calculate the standard error scores
      sumIC50 = summary(nls_result_ic50);
      ic50std_Error <- sumIC50$coefficients["IC50","Std. Error"]; 
      ic50std_resid <- round(sqrt(sum((sumIC50$residuals)^2)/(length(sumIC50$residuals)-1)),1);
      max_signal <- max(mat_tbl$dose,na.rm=T); min_signal <- min(mat_tbl$dose,na.rm=T)
      
      #############################
      #############   Final modification & STD error
      
      #prepare final data and convert IC50 back from log scale (inverse)
      coef_ic50 <- coef(nls_result_ic50)[c("IC50", "SLOPE","MAX","MIN")]; coef_ic50["IC50"] <- 10^coef_ic50["IC50"]
      #(Fix ic50 for curves in wrong direction)
      coef_ic50["IC50"] <- ifelse(coef_ic50["SLOPE"]<0,max_signal,coef_ic50["IC50"])
      #(Fix based on MAX)
      coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<0,max_signal,coef_ic50["IC50"])
      coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<10,max_signal,coef_ic50["IC50"])
      coef_ic50["MAX"] <- ifelse(coef_ic50["MAX"]<0,0,coef_ic50["MAX"])
      #(Fix over sensitive drugs)
      coef_ic50["IC50"] <- ifelse(all(c(max(mat_tbl$inhibition,na.rm=T),min(mat_tbl$inhibition,na.rm=T))>50),min_signal,coef_ic50["IC50"])
      
      # for ploting
      x <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc), length=100)
      yic <- predict(nls_result_ic50, data.frame(logconc=x))
      auc <- MESS::auc(x,yic)
      
      ##average replicates
      mat_tblCp <- mat_tbl[, c("inhibition", "dose")]
      cols_ <- colnames(mat_tblCp)[!grepl("inhibition", colnames(mat_tblCp))] # columns which should be equal to average PI
      X <- as.data.table(mat_tblCp)
      mat_tblCp <- as.data.frame(X[,list(inhibition = mean(inhibition)),cols_], stringAsFactors = !1)
      
      
      perInh <- t(matrix(mat_tblCp[,"inhibition"],dimnames=
                           list(paste0(rep("D", length(mat_tblCp[,"inhibition"])), 1:length(mat_tblCp[,"inhibition"])))))
      
      coef_ic50 <- c(coef_ic50,Min.Conc.tested=min_signal,Max.Conc.tested=max_signal,IC50_std_error=ic50std_Error)
      
      if (graph){
        #plot IC50
        mat_tbl$inhibition = xpr_tbl$inhibition_percent; # if we have all values < 0, they will be replaced
        mat_tbl$viability = 100 - mat_tbl$inhibition;  # we are replacing them back here.
        icpl <- ggplot2::ggplot(mat_tbl, aes(logconc, inhibition2)) + scale_x_continuous(breaks=mat_tbl$logconc,labels=round(mat_tbl$dose,2)) +
          geom_point(color = "blue", size = 2.8) + 
          geom_line(data = data.frame(x = x, y = yic), aes(x, yic), color="blue", size = 0.8) +
          geom_vline(xintercept = log10(coef_ic50["IC50"]), colour="grey", size = 0.8) + 
          theme_bw() + labs(y = "% inhibition", x = "conc(nM)")  +  ylim(-25, 125) +
          geom_text(mapping=aes(x2,y2,label = text2), data=data.frame(x2=log10(coef_ic50["IC50"])*0.95, y2=115, text2="IC50"), color="grey", parse=T) +
          theme(plot.background = element_rect(fill = "transparent",colour = NA),
                panel.background =element_rect(fill = "transparent",colour = NA), 
                plot.title = element_text(hjust = 0.5, size = 12.5),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
          )
        
        graphics.off()
        filename_ = paste0(product_id,"_IC50_curve_drug.png")
        ggsave(filename_, icpl, height = 3, width = 3.7, bg= 'white', dpi = 300, create.dir = TRUE)

      }
      
      #dataframe for IC50
      IC50_dataframe <- data.frame(ID=product_id,DRUG_NAME=drug_name,t(as.matrix(coef_ic50)), perInh, SE_of_estimate = as.numeric(ic50std_resid),AUC=auc)

      # check believe
      resid_ = as.numeric(sumIC50$residuals); resp_ = as.numeric(perInh); cond = 0; outlier = floor(length(resid_ )/2) # min tolerable number
      
      if(ic50std_resid>70) {believe_ = !1; cond = 1} # when SE  is very large
      if(sum(abs(resid_)>30)> outlier) {believe_ = !1; cond = 2};  # when ard half the points are outlying
      if( sum((abs(resid_)>30) == (resp_ < 2)) == length(resid_ ) ) {believe_ = !0; cond = 2.1};  # if the outliers are below 0
      
      if(sum(resp_ < 2) > outlier) {believe_ = !0; cond = 3} # when majority of the points are below 0

      list(IC50_dataframe,
           believe_,
           coef_ic50["IC50"],
           cond
      )
      
      if (return_param){
        return (as.data.frame(t(coef_ic50)))
      } else {
        return (nls_result_ic50) 
      }
      
      
    }
  }, error = function(e) {
    print(e);
  })
})
