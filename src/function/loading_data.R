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

listofnames <- init_drugnames()


## saving pdf
extrafont::loadfonts(device = "pdf", quiet = T)

library(showtext)
font_add("Arial", regular = "C:/Windows/Fonts/arial.ttf")
# showtext_auto()
# showtext_opts(dpi = 600)

theme_pdf <- theme(
  panel.background = element_rect(fill = "transparent", colour = NA), 
  plot.background = element_rect(fill = "transparent", colour = NA),
  legend.background = element_rect(fill = "transparent", colour = NA),
  text = element_text(family = "Arial")
)
