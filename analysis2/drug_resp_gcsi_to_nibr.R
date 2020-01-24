# 22.11.2017, modified on 28.05.2018 
# applying models trained on cell lines to xenograft data

source("http://bioconductor.org/biocLite.R")
biocLite("survcomp")
biocLite("PharmacoGx")
install.packages("caret")
install.packages("gam")
install.packages("lattice")
library(survcomp)
library(PharmacoGx)
library(caret)
library(gam)
library(lattice)

path <- "/data/kurilov/genestack/phd/work_2018/DrugRespPrediction/"
setwd(file.path(path, "analysis2"))
source("functions.R")

###  0:loading gCSI and NIBR PDXE data  ###

## gCSI
load("data/gcsi.genomics.rda")
load("data/gcsi.line.info.rda")
load("data/gcsi.genomics.feature.info.rda")
#let's remove all cl. features
f_cl <- grep("cln.", colnames(gcsi.genomics), v=F)
gcsi_table <- gcsi.genomics[,-f_cl]

gCSI <- downloadPSet("gCSI")
gcsi_auc <- summarizeSensitivityProfiles(gCSI,drugs=drugNames(gCSI) ,sensitivity.measure = "auc_recomputed")
gcsi_auc <- t(gcsi_auc)
gcsi_auc <- 1-gcsi_auc
rownames(gcsi_auc) <- rownames(gcsi.line.info)

## NIBR PDXE
load("data/tissue_nibr.RData")
rownames(tissue_nibr) <- tissue_nibr$sample
load("data/mol_nibr.RData")
# removing tissue variables from mol_nibr
mol_nibr <- mol_nibr[,1:72468]
mol_nibr <- mol_nibr[-152,] # sample X-4215 with missing data in copy numbers
load("data/slope_table_nibr.RData")
load("data/resp_nibr.RData")

int <- intersect(rownames(gcsi.line.info), rownames(gcsi.genomics))

int_erl <- intersect(int, rownames(gcsi.line.info)[which(gcsi.line.info$TissueMetaclass=="Lung")])
int_gem <- intersect(int, rownames(gcsi.line.info)[which(gcsi.line.info$TissueMetaclass=="Pancreas")])
int_pac_breast <- intersect(int, rownames(gcsi.line.info)[which(gcsi.line.info$TissueMetaclass=="Breast")])
int_pac_lung <- intersect(int, rownames(gcsi.line.info)[which(gcsi.line.info$TissueMetaclass=="Lung")])

# NIBR PDXE response data
resp_erlotinib <- resp_nibr$erlotinib
resp_gemcitabine <- resp_nibr$gemcitabine
resp_paclitaxel <- resp_nibr$paclitaxel

# converting gCSI column names for molecular features

columns <- unlist(strsplit(colnames(gcsi_table), split="\\."))
col_table <- matrix(columns, ncol=2, byrow=T)
ord <- match(col_table[,2], gcsi.genomics.feature.info$GeneID)
col_table <- cbind(col_table, gcsi.genomics.feature.info$Symbol[ord])
empty_symbol <- which(col_table[,3]=="")
col_table[empty_symbol,3] <- 1:length(empty_symbol)

intersect(col_table[which(col_table[,1]=="mutp"),3], col_table[which(col_table[,1]=="hot"),3])

columns2 <- apply(col_table,1, function(x) {
  paste0(x[1], "_", x[3])
})
colnames(gcsi_table) <- columns2
columns3 <- gsub("vsd", "exp", columns2, fixed = T)
colnames(gcsi_table) <- columns3

f_int <- intersect(colnames(gcsi_table), colnames(mol_nibr))

# model training

f1_erl <- names(sort(gam_fs(x=gcsi_table[int,f_int], y=gcsi_auc[int,"Erlotinib"])))[1:400]
f1_gem <- names(sort(gam_fs(x=gcsi_table[int,f_int], y=gcsi_auc[int,"Gemcitabine"])))[1:400]
f1_pac_breast <- names(sort(gam_fs(x=gcsi_table[int,f_int], y=gcsi_auc[int,"paclitaxel"])))[1:400]
f1_pac_lung <- names(sort(gam_fs(x=gcsi_table[int,f_int], y=gcsi_auc[int,"paclitaxel"])))[1:400]

f2_erl <- names(sort(gam_fs(x=gcsi_table[int_erl,f_int], y=gcsi_auc[int_erl,"Erlotinib"])))[1:100]
f2_gem <- names(sort(gam_fs(x=gcsi_table[int_gem,f_int], y=gcsi_auc[int_gem,"Gemcitabine"])))[1:100]
f2_pac_breast <- names(sort(gam_fs(x=gcsi_table[int_pac_breast,f_int], y=gcsi_auc[int_pac_breast,"paclitaxel"])))[1:100]
f2_pac_lung <- names(sort(gam_fs(x=gcsi_table[int_pac_lung,f_int], y=gcsi_auc[int_pac_lung,"paclitaxel"])))[1:100]

m1_erl <- simple_regression(y=gcsi_auc[int,"Erlotinib"], x=gcsi_table[int,f1_erl], method="rf")
m1_gem <- simple_regression(y=gcsi_auc[int,"Gemcitabine"], x=gcsi_table[int,f1_gem], method="rf")
m1_pac_breast <- simple_regression(y=gcsi_auc[int,"paclitaxel"], x=gcsi_table[int,f1_pac_breast], method="rf")
m1_pac_lung <- simple_regression(y=gcsi_auc[int,"paclitaxel"], x=gcsi_table[int,f1_pac_lung], method="rf")

m2_erl <- simple_regression(y=gcsi_auc[int_erl,"Erlotinib"], x=gcsi_table[int_erl,f2_erl], method="rf")
m2_gem <- simple_regression(y=gcsi_auc[int_gem,"Gemcitabine"], x=gcsi_table[int_gem,f2_gem], method="rf")
m2_pac_breast <- simple_regression(y=gcsi_auc[int_pac_breast,"paclitaxel"], x=gcsi_table[int_pac_breast,f2_pac_breast], method="rf")
m2_pac_lung <- simple_regression(y=gcsi_auc[int_pac_lung,"paclitaxel"], x=gcsi_table[int_pac_lung,f2_pac_lung], method="rf")


# predictions for NIBR samples

int_erl <- intersect(rownames(resp_erlotinib), rownames(mol_nibr))
int_gem <- intersect(rownames(resp_gemcitabine), rownames(mol_nibr))
int_pac <- intersect(rownames(resp_paclitaxel), rownames(mol_nibr))
int_pac_brca <- intersect(int_pac, tissue_nibr$sample[which(tissue_nibr$tissue=="BRCA")])
int_pac_nsclc <- intersect(int_pac, tissue_nibr$sample[which(tissue_nibr$tissue=="NSCLC")])

p1_erl <- predict.train(m1_erl, newdata=mol_nibr[int_erl,f1_erl])
p1_gem <- predict.train(m1_gem, newdata=mol_nibr[int_gem,f1_gem])
p1_pac_brca <- predict.train(m1_pac_breast, newdata=mol_nibr[int_pac_brca,f1_pac_breast])
p1_pac_nsclc <- predict.train(m1_pac_lung, newdata=mol_nibr[int_pac_nsclc,f1_pac_lung])

p2_erl <- predict.train(m2_erl, newdata=mol_nibr[int_erl,f2_erl])
p2_gem <- predict.train(m2_gem, newdata=mol_nibr[int_gem,f2_gem])
p2_pac_brca <- predict.train(m2_pac_breast, newdata=mol_nibr[int_pac_brca,f2_pac_breast])
p2_pac_nsclc <- predict.train(m2_pac_lung, newdata=mol_nibr[int_pac_nsclc,f2_pac_lung])

cor1 <- c(cor(p1_erl, resp_erlotinib[int_erl,1]), cor(p1_erl, resp_erlotinib[int_erl,2]),
          cor(p1_erl, resp_erlotinib[int_erl,3]), cor(p1_erl, resp_erlotinib[int_erl,4]),
          cor(p1_gem, resp_gemcitabine[int_gem,1]), cor(p1_gem, resp_gemcitabine[int_gem,2]),
          cor(p1_gem, resp_gemcitabine[int_gem,3]), cor(p1_gem, resp_gemcitabine[int_gem,4]),
          cor(p1_pac_brca, resp_paclitaxel[int_pac_brca,1]), cor(p1_pac_brca, resp_paclitaxel[int_pac_brca,2]),
          cor(p1_pac_brca, resp_paclitaxel[int_pac_brca,3]), cor(p1_pac_brca, resp_paclitaxel[int_pac_brca,4]),
          cor(p1_pac_nsclc, resp_paclitaxel[int_pac_nsclc,1]), cor(p1_pac_nsclc, resp_paclitaxel[int_pac_nsclc,2]),
          cor(p1_pac_nsclc, resp_paclitaxel[int_pac_nsclc,3]), cor(p1_pac_nsclc, resp_paclitaxel[int_pac_nsclc,4]))
cor1 <- matrix(cor1, ncol=4,nrow=4, byrow=T)
colnames(cor1) <- colnames(resp_erlotinib)
rownames(cor1) <- c("erlotinib", "gemcitabine", "pacl_brca", "pacl_nsclc")

cor2 <- c(cor(p2_erl, resp_erlotinib[int_erl,1]), cor(p2_erl, resp_erlotinib[int_erl,2]),
          cor(p2_erl, resp_erlotinib[int_erl,3]), cor(p2_erl, resp_erlotinib[int_erl,4]),
          cor(p2_gem, resp_gemcitabine[int_gem,1]), cor(p2_gem, resp_gemcitabine[int_gem,2]),
          cor(p2_gem, resp_gemcitabine[int_gem,3]), cor(p2_gem, resp_gemcitabine[int_gem,4]),
          cor(p2_pac_brca, resp_paclitaxel[int_pac_brca,1]), cor(p2_pac_brca, resp_paclitaxel[int_pac_brca,2]),
          cor(p2_pac_brca, resp_paclitaxel[int_pac_brca,3]), cor(p2_pac_brca, resp_paclitaxel[int_pac_brca,4]),
          cor(p2_pac_nsclc, resp_paclitaxel[int_pac_nsclc,1]), cor(p2_pac_nsclc, resp_paclitaxel[int_pac_nsclc,2]),
          cor(p2_pac_nsclc, resp_paclitaxel[int_pac_nsclc,3]), cor(p2_pac_nsclc, resp_paclitaxel[int_pac_nsclc,4]))
cor2 <- matrix(cor2, ncol=4,nrow=4, byrow=T)
colnames(cor2) <- colnames(resp_erlotinib)
rownames(cor2) <- c("erlotinib", "gemcitabine", "pacl_brca", "pacl_nsclc")

erl_table <- cbind(p1_erl, p2_erl, resp_erlotinib[int_erl,])
gem_table <- cbind(p1_gem, p2_gem, resp_gemcitabine[int_gem,])
pac_brca_table <- cbind(p1_pac_brca, p2_pac_brca, resp_paclitaxel[int_pac_brca,])
pac_nsclc_table <- cbind(p1_pac_nsclc, p2_pac_nsclc, resp_paclitaxel[int_pac_nsclc,])

# saving results
write.table(cor1, file="cor_gcsi_to_nibr_pred_all_lines.txt")
write.table(cor2, file="cor_gcsi_to_nibr_pred_tissue_spec.txt")
write.table(erl_table, file="erl_table_gcsi_to_nibr_pred.txt")
write.table(gem_table, file="gem_table_gcsi_to_nibr_pred.txt")
write.table(pac_brca_table, file="pac_brca_table_gcsi_to_nibr_pred.txt")
write.table(pac_nsclc_table, file="pac_nsclc_table_gcsi_to_nibr_pred.txt")


# plotting with lattice

erl_table <- read.table("erl_table_gcsi_to_nibr_pred.txt")
gem_table <- read.table("gem_table_gcsi_to_nibr_pred.txt")
pac_brca_table <- read.table("pac_brca_table_gcsi_to_nibr_pred.txt")
pac_nsclc_table <- read.table("pac_nsclc_table_gcsi_to_nibr_pred.txt")

d1 <-  erl_table[,c(1,3)]
d2 <- erl_table[,c(1,4)]
d3 <- erl_table[,c(1,5)]
d4 <- erl_table[,c(1,6)]
d5 <-  gem_table[,c(1,3)]
d6 <- gem_table[,c(1,4)]
d7 <- gem_table[,c(1,5)]
d8 <- gem_table[,c(1,6)]
d9 <-  pac_brca_table[,c(1,3)]
d10 <- pac_brca_table[,c(1,4)]
d11 <- pac_brca_table[,c(1,5)]
d12 <- pac_brca_table[,c(1,6)]
d13 <-  pac_nsclc_table[,c(1,3)]
d14 <- pac_nsclc_table[,c(1,4)]
d15 <- pac_nsclc_table[,c(1,5)]
d16 <- pac_nsclc_table[,c(1,6)]
nn <- c("predicted", "observed")
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
colnames(d15) <- nn
colnames(d16) <- nn

titles <- c("Erlotinib, volume", "Erlotinib, AUC", "Erlotinib, slope", "Erlotinib, diff.slope",
            "Gemcitabine, volume", "Gemcitabine, AUC", "Gemcitabine, slope", "Gemcitabine, diff.slope",
            "Paclitaxel (BRCA), volume", "Paclitaxel (BRCA), AUC", "Paclitaxel (BRCA), slope", "Paclitaxel (BRCA), diff.slope",
            "Paclitaxel (NSCLC), volume", "Paclitaxel (NSCLC), AUC", "Paclitaxel (NSCLC), slope", "Paclitaxel (NSCLC), diff.slope")

table_lat <- matrix(NA, ncol=3,nrow=0)

for (i in 1:16)
{
  ds <- eval(parse(text = paste0("d", i)))
  cor <- signif(cor(ds[,1],ds[,2]), digits = 2)
  title <- paste0(titles[i], "\ncor=",cor)
  ds <- cbind(ds, title)
  
  table_lat <- rbind(table_lat, ds)
}
colnames(table_lat) <- c("predicted","observed","label")
table_lat <- as.data.frame(table_lat)
table_lat$observed <- as.numeric(as.character(table_lat$observed))
table_lat$predicted <- as.numeric(as.character(table_lat$predicted))

pdf(file="fig_s10.pdf", width=11.7, height=8.3)
xyplot(predicted ~ observed | label, table_lat,
       grid = TRUE,
       layout=c(4,4),
       scales=list(relation="free"),
       par.strip.text = list(cex=0.7, lines=2),
       index.cond=list(c(16,13,15,14,12,9,11,10,8,5,7,6,4,1,3,2)))
dev.off()

# plotting only volume and slope:
num <- c(grep("AUC",table_lat[,3]), grep("diff.", table_lat[,3]))
table_lat2 <- table_lat[-num,]

pdf(file="fig_5b.pdf", width=11.7, height=7)
xyplot(predicted ~ observed | label, table_lat2,
       grid = TRUE,
       as.table=T,
       layout=c(4,2),
       scales=list(relation="free"),
       par.strip.text = list(cex=1.2, lines=2),
       xlab="observed drug response values (volume or slope)",
       ylab="predicted drug response values (AUC)",
       index.cond=list(c(1,3,5,7,2,4,6,8)))
dev.off()

