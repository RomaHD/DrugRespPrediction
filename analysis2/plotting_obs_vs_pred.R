# 11.10.2017 plotting obs vs. pred. with lattice, edited on 28.05.2018

install.packages("lattice")
library(lattice)

path <- "/data/kurilov/genestack/phd/work_2018/DrugRespPrediction/"
setwd(file.path(path, "analysis2"))
source("functions.R")

# reading data
dt_gcsi <- read.table("perform_table_dt_gcsi.txt")
dt_nibr <- read.table("perform_table_dt_nibr.txt")
erl_gcsi <- read.table("perform_table_erl_gcsi.txt")
gem_gcsi <- read.table("perform_table_gem_gcsi.txt")
pac_breast_gcsi <- read.table("perform_table_pac_breast_gcsi.txt")
pac_lung_gcsi <- read.table("perform_table_pac_lung_gcsi.txt")
erl_nibr <- read.table("perform_table_erl_nibr.txt")
gem_nibr <- read.table("perform_table_gem_nibr.txt")
pac_brca_nibr <- read.table("perform_table_pac_brca_nibr.txt")
pac_nsclc_nibr <- read.table("perform_table_pac_nsclc_nibr.txt")

# preprocessing

d1 <-  dt_gcsi
d2 <- dt_nibr
d3 <- erl_gcsi
d4 <- gem_gcsi
d5 <- pac_breast_gcsi
d6 <- pac_lung_gcsi
d7 <- erl_nibr[,1:2]
d8 <- erl_nibr[,3:4]
d9 <- gem_nibr[,1:2]
d10 <- gem_nibr[,3:4]
d11 <- pac_brca_nibr[,1:2]
d12 <- pac_brca_nibr[,3:4]
d13 <- pac_nsclc_nibr[,1:2]
d14 <- pac_nsclc_nibr[,3:4]
nn <- c("observed","predicted")
colnames(d1) <- nn
colnames(d2) <- nn
colnames(d3) <- nn
colnames(d4) <- nn
colnames(d5) <- nn
colnames(d6) <- nn
colnames(d7) <- nn
colnames(d8) <- nn
colnames(d9) <- nn
colnames(d10) <- nn
colnames(d11) <- nn
colnames(d12) <- nn
colnames(d13) <- nn
colnames(d14) <- nn

titles <- c("Doubling time (cell lines)",
            "Slope (xenografts)",
            "Erlotinib (cell lines)",
            "Gemcitabine (cell lines)",
            "Paclitaxel, bresat (cell lines)",
            "Paclitaxel, lung (cell lines)",
            "Erlotinib, volume (xenografts)",
            "Erlotinib, slope (xenografts)",
            "Gemcitabine, volume (xenografts)",
            "Gemcitabine, slope (xenografts)",
            "Paclitaxel, BRCA, volume (xenografts)",
            "Paclitaxel, BRCA, slope (xenografts)",
            "Paclitaxel, NSCLC, volume (xenografts)",
            "Paclitaxel, NSCLC, slope (xenografts)")

table_lat <- matrix(NA, ncol=3,nrow=0)

for (i in 1:14)
{
  ds <- eval(parse(text = paste0("d", i)))
  r2 <- signif((cor(ds[,1],ds[,2]))^2, digits = 2)
  ci_val <- signif(ci(ds[,1],ds[,2]), digits = 2)
  title <- paste0(titles[i], "\nR2=",r2, "   conc. index=",ci_val)
  ds <- cbind(ds, title)
  
  table_lat <- rbind(table_lat, ds)
}
colnames(table_lat) <- c("observed","predicted","label")

# plotting

pdf(file="fig_5a.pdf", width=11.7, height=8.3)
xyplot(predicted ~ observed | label, table_lat,
       grid = TRUE,
       scales=list(relation="free"),
       layout=c(5,3),
       par.strip.text = list(cex=0.7, lines=2),
       skip=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
       index.cond=list(c(8,10,12,14,2,7,9,11,13,1,3,4,5,6)),
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)  ## First call the default panel function for 'xyplot'
         panel.abline(a=c(0,1), lty=2)  ## Add a horizontal line at the median
       })
dev.off()
