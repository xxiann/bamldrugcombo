require(circlize)
colour = c("aquamarine", "blue", "blue4", "blueviolet", "brown", "brown1", 
           "chartreuse", "chartreuse4", "chocolate", "coral", 
           "cornflowerblue", "cyan4", "darkgoldenrod", "darkgreen", 
           "darkmagenta", "darkolivegreen1", "darkorchid1", "darksalmon",
           "darkseagreen1", "darkslategray3", "darkturquoise", 
           "deeppink", "deepskyblue", "dodgerblue1", "indianred1", 
           "khaki", "lavender", "lavenderblush", "lightblue", 
           "lightsteelblue4", "maroon2", "mediumorchid1", 
           "mediumpurple", "mediumseagreen", "plum", "royalblue", 
           "tomato", "turquoise4", "yellowgreen", "gold")

##### colour ramp #####

col_fun = colorRamp2(c(0, 0.108, 1), c("blue", "white", "red"))
# col_fun2 = colorRamp2(c(0, 60, 100), c("blue", "white", "red"))
color.1 <- colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("darkblue","blue","white","pink","red"))

col_fun3 = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
##### factor #####

disease.col3 = c("Diagnosis"="mediumpurple", "Relapse"="firebrick", "Residual"="darkorange", "Remission"="gold","PostChemo"="gold", "R/R"="pink", "Unknown"="grey65", "Healthy"="seagreen" )

disease.order3 = c("Diagnosis", "Relapse", "Residual", "Remission", "Unknown")

specimen.col = c("Bone Marrow Aspirate"="beige", "Leukapheresis"="yellow", "Peripheral Blood"="darkred")

subtype.col1 = c("M1"="darkolivegreen1",
                 # "M1; M0" = "deepskyblue",
                 "M0; M1" = "darkslategray3",
                 "M2" = "cyan4",
                 "M4" = "deepskyblue",
                 "M4 eos" = "royalblue",
                 "M4; M5" = "plum",
                 "M5" = "darkmagenta")

subtype.col2 = c("M0"="brown1",
                 "M1"="darkolivegreen1",
                 "M0/M1" = "deepskyblue",
                 "M2" = "cyan4",
                 "M3" = "khaki",
                 "M4" = "royalblue",
                 "M4eo" = "darkslategray3",
                 "M5" = "darkmagenta",
                 "M6" = "maroon2",
                 "M6a" = "firebrick",
                 "M7" = "darksalmon",
                 "NOS" = "grey65")

subtype.col3 = c("M0"="brown1",
                 "M1"="darkolivegreen3",
                 "M2" = "cyan4",
                 "M3" = "khaki3",
                 "M4" = "royalblue",
                 "M5" = "darkmagenta",
                 "M6" = "maroon2",
                 "M6a" = "firebrick",
                 "M7" = "darksalmon",
                 "NOS" = "grey65")


cohort.col = c("Waves1+2" = "orange",
               "Waves3+4" = "royalblue",
               "Both" = "seagreen",
               "NotAvail" = "grey65"
)

eln.col = c("Favorable" = "seagreen",
            "Intermediate" = "orange",
            "Adverse" = "#E41A1C",
            "NonInitial" = "lavender",
            "Unknown" = "grey65"
)

eln.col.order = c("Favorable", "Intermediate", "Adverse", "NonInitial", "Unknown")

fusion <- c("BCR-ABL1" = "#F8766D",
            "CBFB-MYH11" = "#CD9600",
            "DEK-NUP214" = "#7CAE00",
            "KMT2A_re" = "#00BA38", #"RPN1-EVI1" = "#00BA38",
            "MECOM_re" = "#00C08B" ,
            "NUP98_re" = "#00BFC4", #"MLLT3-MLL" = "#00BFC4",
            "PML-RARA" = "#00A9FF",
            "RARA_re" = "#C77CFF",
            "RUNX1-RUNX1T1" = "#FF61CC",
            "Other" = "#FFD92F",
            "NotDetected" = "grey90",
            "Unknown" = "grey65"
)

fusion.order = c("BCR-ABL1","CBFB-MYH11","DEK-NUP214","KMT2A_re","MECOM_re", "NUP98_re","PML-RARA","RARA_re","RUNX1-RUNX1T1","Other", "NotDetected", NA)

##### drugs #####
drugclass = c("Apoptotic modulator" = "#E41A1C", 
              "Chemotherapeutics" ="#FFD92F",
              "Differentiating/epigenetic modifier" = "#4DAF4A",
              "Hormone" = "#984EA3",
              "HSP inhibitor" = "#8DA0CB",
              "Immunomodulator" ="#FF7F00",
              "Kinase inhibitor" = "royalblue",
              "Kinesin inhibitor" = "darksalmon",
              "Metabolic modifier" = "darkolivegreen2",
              "mTOR inhibitor" = "turquoise4",
              "Other" = "plum"
)

combi.partner = c("Apoptotic modulator - Kinase inhibitor",
                  "Apoptotic modulator - Metabolic modifier",
                  "Apoptotic modulator - Chemotherapeutics",
                  "Apoptotic modulator - HSP inhibitor",
                  "Apoptotic modulator - mTOR inhibitor",
                  "Immunomodulator - Kinase inhibitor",
                  "Kinase inhibitor - Kinase inhibitor",
                  "Chemotherapeutics - Kinase inhibitor",
                  "Kinase inhibitor - Metabolic modifier",
                  "Apoptotic modulator - Immunomodulator",
                  "Chemotherapeutics - Metabolic modifier",
                  "Differentiating/epigenetic modifier - Kinase inhibitor",
                  "Apoptotic modulator - Differentiating/epigenetic modifier",
                  "Kinase inhibitor - Other",
                  "Chemotherapeutics - Immunomodulator",
                  "Chemotherapeutics - Chemotherapeutics",
                  "Chemotherapeutics - Differentiating/epigenetic modifier" ,
                  "HSP inhibitor - Kinase inhibitor",
                  "Chemotherapeutics - mTOR inhibitor",
                  "Apoptotic modulator - Apoptotic modulator",
                  "Differentiating/epigenetic modifier - Immunomodulator",
                  "Differentiating/epigenetic modifier - Other",
                  "Differentiating/epigenetic modifier - Differentiating/epigenetic modifier"
)

##### other ######
cell.type.name = c("LSPC-Quiescent", "LSPC-Primed", "LSPC-Cycle", "GMP-like", "ProMono-like", "Mono-like", "cDC-like")

cell.type.order = c("LSPC.Quiescent", "LSPC.Primed.Top100", "LSPC.Cycle.Top100", "GMP.like.Top100", "ProMono.like.Top100", "Mono.like.Top100", "cDC.like.Top100")

hierarchy_colors = c(LSPC.Quiescent = '#A7303099', 
                     LSPC.Primed.Top100 = '#CD534C99', 
                     LSPC.Cycle.Top100 = '#8F770099', 
                     GMP.like.Top100 = '#3B3B3B99', 
                     ProMono.like.Top100 = '#0073C299', 
                     Mono.like.Top100 = '#003C6799', 
                     cDC.like.Top100 = '#EFC00099',
                     `LSPC-Quiescent` = '#A7303099', 
                     `LSPC-Primed` = '#CD534C99', 
                     `LSPC-Cycle` = '#8F770099', 
                     `GMP-like` = '#3B3B3B99', 
                     `ProMono-like` = '#0073C299', 
                     `Mono-like` = '#003C6799', 
                     `cDC-like` = '#EFC00099')
##### genes #####
goi = c("TP53", "DNMT3A", "ASXL1", "TET2", "SF3B1", "SRSF2", 
        "NPM1", "CEBPA", "RUNX1", "FLT3", "FLT3-ITD", "KIT", 
        "WT1", "PTPN11", "CEBPA", "U2AF1", "ZRSR2", "EZH2", 
        "BCOR", "STAG2", "ETV6", "GATA2", "KRAS", "NRAS", 
        "NF1", "RAD21", "SMC1A", "SMC3", "IDH1", "IDH2", 
        "KMT2A", "JAK2")

##### brewer pal ######
# brewer.pal.info
# display.brewer.all()
dark2 <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E",
           "#E6AB02","#A6761D","#666666")

Set3 <- c(                    
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
  "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F" )

Set2 <- c(
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"
)

Set1 <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
  "#A65628", "#F781BF", "#999999"
)

Custom <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
  "#A65628", "darksalmon", "#999999","#7570B3",
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#66C2A5", "#80B1D3", "#FDB462", "#B3DE69",
  "#FCCDE5", "#D9D9D9", "#BC80BD", "#FB8072", "#F781BF", "#FFED6F",
  "darkolivegreen1")

combi.partner.col = setNames(Custom[1:23], combi.partner)



## seurat RNA seq colours
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=6)

