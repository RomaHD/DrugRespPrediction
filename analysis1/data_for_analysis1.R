# 23.01.18 getting data from pharmacoGx for rerunning Analysis I (CCLE, CTRP, GDSC comparison)
# edited on 19.04.18
source("http://bioconductor.org/biocLite.R")
biocLite("PharmacoGx")
biocLite("org.Hs.eg.db")
biocLite("AnnotationDbi")
library(PharmacoGx)
library(org.Hs.eg.db) 
library(AnnotationDbi)
library(Biobase)

path <- "/data/kurilov/genestack/phd/work_2018/DrugRespPrediction/"
setwd(file.path(path, "analysis1"))

CCLE <- downloadPSet("CCLE")
CTRP <- downloadPSet("CTRPv2")
GDSC <- downloadPSet("GDSC")

#CCLE/CTRP
ccle_exp <- summarizeMolecularProfiles(pSet=CCLE, mDataType="rna", cell.lines=cellNames(CCLE), summary.stat = 'mean', fill.missing = TRUE, verbose=TRUE)
ccle_cnv <- summarizeMolecularProfiles(pSet=CCLE, mDataType="cnv", cell.lines=cellNames(CCLE), summary.stat = 'mean', fill.missing = TRUE, verbose=TRUE)
ccle_mut <- summarizeMolecularProfiles(pSet=CCLE, mDataType="mutation", cell.lines=cellNames(CCLE), summary.stat = 'or', fill.missing = TRUE, verbose=TRUE)

ccle_exp <- t(exprs(ccle_exp))
ccle_cnv <- t(exprs(ccle_cnv))
ccle_mut <- t(exprs(ccle_mut))

#GDSC
gdsc_exp <- summarizeMolecularProfiles(pSet=GDSC, mDataType="rna", cell.lines=cellNames(GDSC), summary.stat = 'mean', fill.missing = TRUE, verbose=TRUE)
gdsc_cnv <- summarizeMolecularProfiles(pSet=GDSC, mDataType="cnv", cell.lines=cellNames(GDSC), summary.stat = 'mean', fill.missing = TRUE, verbose=TRUE)
gdsc_mut <- summarizeMolecularProfiles(pSet=GDSC, mDataType="mutation", cell.lines=cellNames(GDSC), summary.stat = 'or', fill.missing = TRUE, verbose=TRUE)

gdsc_exp <- t(exprs(gdsc_exp))
gdsc_cnv <- t(exprs(gdsc_cnv))
gdsc_mut <- t(exprs(gdsc_mut))

# ensemble -> gene symbol conversion for ccle
gene_order <- match(colnames(ccle_exp), featureInfo(CCLE, "rna")[,"Probe"])
colnames(ccle_exp) <- featureInfo(CCLE, "rna")[gene_order,"Symbol"]
na_names <- which(duplicated(colnames(ccle_exp)))
ccle_exp <- ccle_exp[,-na_names]

# ensemble -> gene symbol conversion for gdsc
gene_order <- match(colnames(gdsc_exp), featureInfo(GDSC, "rna")[,"Probe"])
colnames(gdsc_exp) <- featureInfo(GDSC, "rna")[gene_order,"Symbol"]
na_names <- which(duplicated(colnames(gdsc_exp)))
gdsc_exp <- gdsc_exp[,-na_names]

#prefixes
colnames(ccle_exp) <- paste0("exp_", colnames(ccle_exp))
colnames(gdsc_exp) <- paste0("exp_", colnames(gdsc_exp))
colnames(ccle_cnv) <- paste0("cnv_", colnames(ccle_cnv))
colnames(gdsc_cnv) <- paste0("cnv_", colnames(gdsc_cnv))
colnames(ccle_mut) <- paste0("mut_", colnames(ccle_mut))
colnames(gdsc_mut) <- paste0("mut_", colnames(gdsc_mut))

# drug_response
drugNames(CTRP) <- CTRP@drug$cpd_name

drugs_ccle <- c("Erlotinib", "paclitaxel", "lapatinib", "Nilotinib", "Nutlin-3", "PLX4720", "Sorafenib")
drugs_ctrp <- c("erlotinib", "paclitaxel", "lapatinib", "nilotinib", "nutlin-3", "PLX-4720", "sorafenib")
drugs_gdsc <- c("Erlotinib", "paclitaxel", "lapatinib", "Nilotinib", "Nutlin-3", "PLX4720", "Sorafenib")

ccle_ic50 <- t(summarizeSensitivityProfiles(CCLE, sensitivity.measure="ic50_recomputed", drugs=drugs_ccle))
ccle_auc <- t(summarizeSensitivityProfiles(CCLE, sensitivity.measure="auc_recomputed", drugs=drugs_ccle))
ctrp_ic50 <- t(summarizeSensitivityProfiles(CTRP, sensitivity.measure="ic50_recomputed", drugs=drugs_ctrp))
ctrp_auc <- t(summarizeSensitivityProfiles(CTRP, sensitivity.measure="auc_recomputed", drugs=drugs_ctrp))
gdsc_ic50 <- t(summarizeSensitivityProfiles(GDSC, sensitivity.measure="ic50_recomputed", drugs=drugs_gdsc))
gdsc_auc <- t(summarizeSensitivityProfiles(GDSC, sensitivity.measure="auc_recomputed", drugs=drugs_gdsc))

#load viability_1uM data

load("viability_data/viab_ccle.RData")
load("viability_data/viab_ctrp.RData")
load("viability_data/viab_gdsc.RData")

ccle_viab <- viab_ccle[,c("Erlotinib", "Paclitaxel", "Lapatinib", "Nilotinib", "Nutlin-3", "PLX4720", "Sorafenib")]
ctrp_viab <- viab_ctrp[,c("erlotinib", "paclitaxel", "lapatinib", "nilotinib", "nutlin-3", "PLX-4720", "sorafenib")]
gdsc_viab <- viab_gdsc[,c("Erlotinib", "Paclitaxel", "Lapatinib", "Nilotinib", "Nutlin-3a", "PLX4720", "Sorafenib")]

# renaming rownames for viability data
temp <- read.csv("viability_data/viab_names_diff.txt")
temp <- temp[,2:3]
rownames(ccle_viab) <- cellInfo(CCLE)$cellid[match(rownames(ccle_viab), cellInfo(CCLE)$CCLE.cellid)]
ord_ctrp <- match(temp[,1], rownames(ctrp_viab))
rownames(ctrp_viab)[ord_ctrp] <- as.character(temp[,2])
rownames(gdsc_viab) <- cellInfo(GDSC)$cellid[match(rownames(gdsc_viab), cellInfo(GDSC)$Sample.name)]

#drug names and prefixes
drug_list <- c("Erlotinib",  "Paclitaxel", "Lapatinib",  "Nilotinib", "Nutlin-3","PLX4720", "Sorafenib")
colnames(ccle_ic50) <- paste0("ic50_", drug_list)
colnames(ccle_auc) <- paste0("auc_", drug_list)
colnames(ccle_viab) <- paste0("viab_", drug_list)
colnames(ctrp_ic50) <- paste0("ic50_", drug_list)
colnames(ctrp_auc) <- paste0("auc_", drug_list)
colnames(ctrp_viab) <- paste0("viab_", drug_list)
colnames(gdsc_ic50) <- paste0("ic50_", drug_list)
colnames(gdsc_auc) <- paste0("auc_", drug_list)
colnames(gdsc_viab) <- paste0("viab_", drug_list)

# 1-AUC correction for AUC data
ccle_auc <- 1-ccle_auc
ctrp_auc <- 1-ctrp_auc
gdsc_auc <- 1-gdsc_auc

#combining data
ccle_table <- cbind(ccle_exp, ccle_cnv, ccle_mut, ccle_ic50, ccle_auc, ccle_viab[rownames(ccle_exp),])
lines_ctrp <- intersect(cellInfo(CCLE)$cellid, cellInfo(CTRP)$cellid)
ctrp_table <- cbind(ccle_exp[lines_ctrp,], ccle_cnv[lines_ctrp,], ccle_mut[lines_ctrp,], ctrp_ic50[lines_ctrp,], ctrp_auc[lines_ctrp,], ctrp_viab[lines_ctrp,])
gdsc_table <- cbind(gdsc_exp, gdsc_cnv, gdsc_mut, gdsc_ic50, gdsc_auc, gdsc_viab[rownames(gdsc_exp),])

#removing redundant cols 
x <- apply(ccle_table,2,function(x){sd(x, na.rm=T)})
sdz <- which(x==0)
ccle_table <- ccle_table[,-sdz]

x <- apply(ctrp_table,2,function(x){sd(x, na.rm=T)})
sdz <- which(x==0)
ctrp_table <- ctrp_table[,-sdz]

x <- apply(gdsc_table,2,function(x){sd(x, na.rm=T)})
sdz <- which(x==0)
sdn <- which(is.na(x))
gdsc_table <- gdsc_table[,-c(sdz,sdn[1:16])]

##
y1 <- apply(ccle_table[,1:45989],1,function(x){length(which(is.na(x)))==length(x)})
y2 <- apply(ccle_table[,45990:46010],1,function(x){length(which(is.na(x)))==length(x)})
sdn <- which(y1)
sdn2 <- which(y2)
ccle_table <- ccle_table[-c(sdn,sdn2),]

y1 <- apply(ctrp_table[,1:45987],1,function(x){length(which(is.na(x)))==length(x)})
y2 <- apply(ctrp_table[,45988:46008],1,function(x){length(which(is.na(x)))==length(x)})
sdn <- which(y1)
sdn2 <- which(y2)
ctrp_table <- ctrp_table[-c(sdn, sdn2),]

y1 <- apply(gdsc_table[,1:36721],1,function(x){length(which(is.na(x)))==length(x)})
y2 <- apply(gdsc_table[,36722:36742],1,function(x){length(which(is.na(x)))==length(x)})
sdn <- which(y1)
sdn2 <- which(y2)
gdsc_table <- gdsc_table[-c(sdn, sdn2),]

ccle_table <- as.matrix(ccle_table)
ctrp_table <- as.matrix(ctrp_table)
gdsc_table <- as.matrix(gdsc_table)
mode(ccle_table) <- "numeric"
mode(ctrp_table) <- "numeric"
mode(gdsc_table) <- "numeric"

# substituting outlier ic50 values (those which are higher than 85 percintile) with 85percintile
outlier_capping <- function(table){
  for (i in grep("ic50", colnames(table))){
    fin <- which(is.finite(table[,i]))
    qv <- quantile(table[fin,i], probs=0.85, na.rm=T)
    out <- which(table[,i]>qv)
    table[out,i] <- qv
  }
  return(table)
}

ccle_table <- outlier_capping(ccle_table)
ctrp_table <- outlier_capping(ctrp_table)
gdsc_table <- outlier_capping(gdsc_table)


# saving tables
save(ccle_table, file="ccle_table.RData")
save(ctrp_table, file="ctrp_table.RData")
save(gdsc_table, file="gdsc_table.RData")
