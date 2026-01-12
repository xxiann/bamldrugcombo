require(tidyverse)

## named vector
init_drugnames <- function(){
  ## changing column names
  combi <- read.csv("data/oshu/combi_detail.csv") %>% select(4)  %>%
    distinct() %>% mutate(newname = make.names(combo)) %>% dplyr::rename(Drug=combo)
  mono <- read.csv("data/oshu/mono_detail.csv") %>% select(2) %>% distinct() %>% mutate(newname = make.names(Drug))
  listofnames <- rbind(mono, combi)
  listofnames <- setNames(listofnames$Drug, listofnames$newname)
  return(listofnames)
}