# referenced from https://github.com/shuyuzheng/synergyfinder

packages.required <- c("nleqslv","reshape2","drc", "magrittr", "data.table", "tidyverse")

lapply(packages.required, library, character.only = T)

source("src/function/Reshape_data.R")
source("src/function/calculate_sensitivity_score.R")
source("src/function/calculate_synergy_score.R")
source("src/function/fit_dose_response_model.R")

# modified functions

#' @param response should contain 3 rows: conc1 conc2 response (\%inhibition)
#' @param single_drug_data list of the 2 single drugs used; contains dataframe with 2 columns - conc & response
#' @description
#' doesnt include single drug as part of calculations
#' 
Loewe_2 <- function(single_drug_data, response,
                    Emin = NA,
                    Emax = NA,
                    quiet = TRUE) {
  if (quiet) {
    options(warn = -1)
  }
  
  single_drug_model <- lapply(
    single_drug_data,
    function(x)  suppressWarnings(FitDoseResponse(x, Emin = Emin, Emax = Emax))
  )
  
  single_par <- lapply(single_drug_model, FindModelPar)
  single_type <- lapply(single_drug_model, FindModelType)
  
  y_loewe <- c()
  dist_loewe <- c()
  ci_loewe <- c()
  for (i in 1:nrow(response)) {
    x <- response[i, 1:2] # concentrations of drugs
    y <- response$response[i] # the observed combination response
    
    # find the dose of single drugs that achieve the observed combination response
    x_cap <- mapply(function(par, type) .SolveExpDose(y, par, type),
                    par = single_par, type = single_type
    )
    
    if (all(!is.finite(x_cap))) { # if none of drugs achieve the combination response
      # max of the single drug response
      y_loewe[i] <- max(mapply(function(model) {
        suppressWarnings(.PredictResponseFromModel(model, sum(x)))
      },
      model = single_drug_model
      ))
      dist_loewe[i] <- NA
    } else {
      # determine the minimal distance
      tmp <- .SolveLoewe(x, single_par, single_type, nsteps = 100)
      y_loewe[i] <- tmp$y_loewe
      dist_loewe[i] <- tmp$distance
    }
    ci_loewe[i] <- sum(x/x_cap)
    
  }
  
  output <- list(
    Loewe_ref = y_loewe,
    Excess_loewe = response$response - y_loewe,
    CI = ci_loewe
  )
  
  # Output results as a list
  return(output)
  options(warn = 0)
}