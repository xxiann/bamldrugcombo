## scRNAseq colours
cell.type.order = c("LSPC-Quiescent", "LSPC-Primed-Top100", "LSPC-Cycle-Top100", "GMP-like-Top100", "ProMono-like-Top100", "Mono-like-Top100", "cDC-like-Top100")

cell.type.name = c("LSPC-Quiescent", "LSPC-Primed", "LSPC-Cycle", "GMP-like", "ProMono-like", "Mono-like", "cDC-like")

old_hierarchy_order = c('HSC-like', 'Prog-like', 'GMP-like', 'ProMono-like', 
                        'Mono-like', 'cDC-like',
                        'HSC', 'Prog', 'GMP', 'ProMono', 'Monocyte',
                        'cDC', 'pDC', 'earlyEry', 'lateEry', 'ProB', 'B', 'Plasma', 'T', 'CTL', 'NK')

new_hierarchy_order = c(cell.type.name, 'Monocyte', 'cDC', 'B', 'Plasma', 'T', 'CTL', 'NK')

hierarchy_colors = c(`HSC-like` = "#BD3030", `Prog-like` = "#ce9866", 
                     `LSPC-Quiescent` = '#A7303099', 
                     `LSPC-Primed` = '#CD534C99', `LSPC-Cycle` = '#8F770099', 
                     `GMP-like` = '#3B3B3B', `ProMono-like` = '#0073C2', 
                     `Mono-like` = '#003C67', `cDC-like` = '#EFC000',
                     
                     HSC = "#A65628", Prog = "#FF7F00", GMP = "#FDB462", 
                     ProMono = "#66A61E", Mono ="#B3DE69", Monocyte ="#B3DE69",
                     cDC = "#66C2A5", pDC = "#8DD3C7", 
                     earlyEry = "#E7298A", lateEry = "#F781BF", 
                     ProB = "#984EA3", B = "#BEBADA", Plasma = "#FCCDE5", 
                     `T` = "royalblue", CTL = "dodgerblue1", NK = "darkturquoise")

hierarchy_order = c('HSC-like', 'Prog-like', 'LSPC-Quiescent', 'LSPC-Primed', 'LSPC-Cycle', 
                    'GMP-like', 'ProMono-like', 'Mono-like', 'cDC-like',
                    'HSC', 'Prog', 'GMP', 'ProMono', 'Mono', 'Monocyte',
                    'cDC', 'pDC', 'earlyEry', 'lateEry', 'ProB', 'B', 'Plasma', 'T', 'CTL', 'NK')


color_celltype <- c("AML_nonmalignant" = "#8DA0CB", "AML_malignant" = "#FC8D62", "Healthy_normal" = "#66C2A5")