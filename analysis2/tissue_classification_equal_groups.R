# 28.05.2018, tissue classification with equal groups

#source("http://bioconductor.org/biocLite.R")

install.packages("caret")
install.packages("gam")
install.packages("dplyr")
install.packages("gplots")

library(caret)
library(gam)
library(dplyr)
library(gplots)

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

## NIBR PDXE
load("data/tissue_nibr.RData")
rownames(tissue_nibr) <- tissue_nibr$sample
load("data/mol_nibr.RData")
# removing tissue variables from mol_nibr
mol_nibr <- mol_nibr[,1:72468]
mol_nibr <- mol_nibr[-152,] # sample X-4215 with missing data in copy numbers

###  tissue type predictions with number of samples per tissue being equal  ###

## gCSI
int <- intersect(rownames(gcsi.line.info), rownames(gcsi.genomics))

# using only 6 biggest tissue types: Colon, Skin,  Pancreas, Breast, Blood, Lung
major <- rownames(gcsi.line.info)[which(gcsi.line.info$TissueMetaclass %in% c("Colon", "Skin",  "Pancreas", "Breast", "Blood", "Lung"))]
int3 <- intersect(int, major)

tissues <- c("Lung", "Blood", "Breast", "Pancreas", "Colon", "Skin")

t_table <- cbind(rownames(gcsi.line.info[int3,]),gcsi.line.info[int3,2])
res <- matrix(NA, nrow=10, ncol=length(tissues)+1)
colnames(res) <- c("accuracy", paste0("acc_",tissues))

pred_obs <- matrix(NA, nrow=0, ncol=2)
colnames(pred_obs) <- c("pred", "obs")

for (i in 1:10){
  print(i)
  lung_sampl <- sample(t_table[which(t_table[,2]=="Lung"),1], 23)
  blood_sampl <- sample(t_table[which(t_table[,2]=="Blood"),1], 23)
  breast_sampl <- sample(t_table[which(t_table[,2]=="Breast"),1], 23)
  pancreas_sampl <- sample(t_table[which(t_table[,2]=="Pancreas"),1], 23)
  colon_sampl <- sample(t_table[which(t_table[,2]=="Colon"),1], 23)
  skin_sampl <- sample(t_table[which(t_table[,2]=="Skin"),1], 23)
  
  train_sampl <- c(sample(lung_sampl,16), sample(blood_sampl,16), sample(breast_sampl,16),
                   sample(pancreas_sampl,16), sample(colon_sampl,16), sample(skin_sampl,16))
  test_sampl <- setdiff(c(lung_sampl, blood_sampl, breast_sampl, pancreas_sampl, colon_sampl, skin_sampl), train_sampl)
  
  f_anova <- names(sort(anova_fs(x=gcsi_table[train_sampl,], y=gcsi.line.info[train_sampl,2])))[1:400]
  
  model_xgbt <- simple_classification(y=gcsi.line.info[train_sampl,2], x=gcsi_table[train_sampl,f_anova], method="xgbTree")
  temp_res <- accuracy_check_class2(model_xgbt, test_sampl, f_anova, gcsi_table, gcsi.line.info, 2) 
  res[i,] <- temp_res[[1]]
  pred_obs <- rbind(pred_obs, temp_res[[2]])
} 

write.table(res, file="res_tissue_gcsi_equal_groups.txt")
write.table(pred_obs, file="pred_obs_tissue_gcsi_equal_groups.txt")

# plotting

conf_table <- matrix(NA, nrow=length(tissues), ncol=length(tissues))
rownames(conf_table) <- tissues
colnames(conf_table) <- tissues
pred_obs <- as.data.frame(pred_obs)

for (i in 1:length(tissues)){
  for (j in 1:length(tissues)){
    conf_table[i,j] <- pred_obs %>% filter(obs==tissues[i] & pred==tissues[j]) %>% nrow()
  }
}
pdf(file="fig_4.pdf")
heatmap.2(conf_table, Rowv=F, Colv = F, density.info = "none", trace="none", col=bluered, margins = c(7,7) )
dev.off()

## NIBR PDXE

int <- intersect(tissue_nibr$sample, rownames(mol_nibr))
tissues <- c("BRCA", "CM", "CRC", "NSCLC", "PDAC")
t_table <- tissue_nibr[int,]
t_table <- as.matrix(t_table)
res <- matrix(NA, nrow=10, ncol=length(tissues)+1)
colnames(res) <- c("accuracy", paste0("acc_",tissues))

pred_obs <- matrix(NA, nrow=0, ncol=2)
colnames(pred_obs) <- c("pred", "obs")

for (i in 1:10){
  print(i)
  brca_sampl <- sample(t_table[which(t_table[,2]=="BRCA"),1], 23)
  cm_sampl <- sample(t_table[which(t_table[,2]=="CM"),1], 23)
  crc_sampl <- sample(t_table[which(t_table[,2]=="CRC"),1], 23)
  nsclc_sampl <- sample(t_table[which(t_table[,2]=="NSCLC"),1], 23)
  pdac_sampl <- sample(t_table[which(t_table[,2]=="PDAC"),1], 23)
  
  
  train_sampl <- c(sample(brca_sampl,16), sample(cm_sampl,16), sample(crc_sampl,16),
                   sample(nsclc_sampl,16), sample(pdac_sampl,16))
  test_sampl <- setdiff(c(brca_sampl, cm_sampl, crc_sampl, nsclc_sampl, pdac_sampl), train_sampl)
  
  #feature selection
  f_anova <- names(sort(anova_fs(x=mol_nibr[train_sampl,], y=tissue_nibr[train_sampl,"tissue"])))[1:400]
  
  model_xgbt <- simple_classification(y=tissue_nibr[train_sampl,2], x=mol_nibr[train_sampl,f_anova], method="xgbTree")
  temp_res <- accuracy_check_class2(model_xgbt, test_sampl, f_anova, mol_nibr, tissue_nibr, 2) 
  res[i,] <- temp_res[[1]]
  pred_obs <- rbind(pred_obs, temp_res[[2]])
} 

write.table(res, file="res_tissue_nibr_equal_groups.txt")
write.table(pred_obs, file="pred_obs_tissue_nibr_equal_groups.txt")

# plotting
conf_table <- matrix(NA, nrow=length(tissues), ncol=length(tissues))
rownames(conf_table) <- tissues
colnames(conf_table) <- tissues
pred_obs <- as.data.frame(pred_obs)

for (i in 1:length(tissues)){
  for (j in 1:length(tissues)){
    conf_table[i,j] <- pred_obs %>% filter(obs==tissues[i] & pred==tissues[j]) %>% nrow()
  }
}

heatmap.2(conf_table, Rowv=F, Colv = F, density.info = "none", trace="none", col=bluered, margins = c(7,7) )

