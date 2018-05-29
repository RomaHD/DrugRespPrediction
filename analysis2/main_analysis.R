# 27.05.2018 script for main analysis
# I: tissue type prediction, II: doubling time/slope prediction, III: drug response prediction

source("http://bioconductor.org/biocLite.R")
biocLite("survcomp")
biocLite("PharmacoGx")
install.packages("caret")
install.packages("gam")
install.packages("ggfortify")
library(survcomp)
library(PharmacoGx)
library(caret)
library(gam)
library(ggfortify)


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
save(gcsi_auc, file="data/gcsi_auc.RData")

## NIBR PDXE
load("data/tissue_nibr.RData")
rownames(tissue_nibr) <- tissue_nibr$sample
load("data/mol_nibr.RData")
# removing tissue variables from mol_nibr
mol_nibr <- mol_nibr[,1:72468]
mol_nibr <- mol_nibr[-152,] # sample X-4215 with missing data in copy numbers
load("data/slope_table_nibr.RData")
load("data/resp_nibr.RData")

###  I: tissue type prediction  ###

## I: gCSI ##
int <- intersect(rownames(gcsi.line.info), rownames(gcsi.genomics))

# let's exclude tissues with only 5 or less samples: Bone, Bone Marrow, Cervix, Head and Neck, 
# Muscle, Prostate, Urinary, Uterus
excl_t <- c("Bone", "Bone Marrow", "Cervix", "Head and Neck", "Muscle", "Prostate", "Urinary", "Uterus", "Unclassified")
excl <- rownames(gcsi.line.info)[which(gcsi.line.info$TissueMetaclass %in% excl_t)]
int2 <- setdiff(int, excl)

# test pca plot
library(ggfortify)
prcomponents <- prcomp(gcsi.genomics[int2,])
p <- autoplot(prcomponents, data = gcsi.line.info[int2,] ,colour =  "TissueMetaclass")
p +geom_text(aes(label=gcsi.line.info[int2,2], size=8, colour =  gcsi.line.info[int2,2]))
ggsave("fig_s11_1.pdf")  

# for variation with only 6 biggest tissue types: Colon, Skin,  Pancreas, Breast, Blood, Lung
major <- rownames(gcsi.line.info)[which(gcsi.line.info$TissueMetaclass %in% c("Colon", "Skin",  "Pancreas", "Breast", "Blood", "Lung"))]
int3 <- intersect(int2, major)

res <- rep(NA, 10)
for (i in 1:10){
  print(i)
  train_sampl <- sample(int2, round(0.7*length(int2)))
  test_sampl <- setdiff(int2, train_sampl)
  
  f_anova <- names(sort(anova_fs(x=gcsi_table[train_sampl,], y=gcsi.line.info[train_sampl,2])))[1:400]
  model_xgbt <- simple_classification(y=gcsi.line.info[train_sampl,2], x=gcsi_table[train_sampl,f_anova], method="xgbTree")
  res[i] <- accuracy_check_class(model_xgbt, test_sampl, f_anova, gcsi_table, gcsi.line.info, 2)
} 

#mean accuracy
mean(res)

write.table(res, file="res_tissue_gcsi.txt")


## I: NIBR PDXE ##

int <- intersect(tissue_nibr$sample, rownames(mol_nibr))


#removing outlier sample X-4145
int <- int[-150]

# test pca plot
prcomponents <- prcomp(mol_nibr[int,])
p <- autoplot(prcomponents, data = tissue_nibr[int,] ,colour =  "tissue")
p +geom_text(aes(label=tissue_nibr[int,"tissue"], size=8, colour =  tissue_nibr[int,"tissue"]))
ggsave("fig_s11_2.pdf") 

res <- rep(NA, 10)
for (i in 1:10){
  print(i)
  train_sampl <- sample(int, round(0.7*191))
  test_sampl <- setdiff(int, train_sampl)

  f_anova <- names(sort(anova_fs(x=mol_nibr[train_sampl,], y=tissue_nibr[train_sampl,"tissue"])))[1:400]
  model_xgbt <- simple_classification(y=tissue_nibr[train_sampl,2], x=mol_nibr[train_sampl,f_anova], method="xgbTree")
  res[i] <- accuracy_check_class(model_xgbt, test_sampl, f_anova, mol_nibr, tissue_nibr,2)
} 

#mean accuracy
mean(res)

write.table(res, file="res_tissue_nibr.txt")

###  II: doubling time/slope prediction  ###

## II: gCSI ##
int <- intersect(rownames(gcsi.line.info), rownames(gcsi.genomics))
int <- int[which(!(is.na(gcsi.line.info[int,3])))]

res <- matrix(NA, 10,3)
colnames(res) <- c("RMSE", "R2", "concordance_ind")
for (i in 1:10){
  print(i)
  train_sampl <- sample(int, round(0.7*length(int)))
  test_sampl <- setdiff(int, train_sampl)
  
  f_gam <- names(sort(gam_fs(x=gcsi_table[train_sampl,], y=gcsi.line.info[train_sampl,3])))[1:400]
  model_rf <- simple_regression(y=gcsi.line.info[train_sampl,3], x=gcsi_table[train_sampl,f_gam], method="rf")
  res[i,] <- accuracy_check(model_rf, test_sampl, f_gam, gcsi_table, gcsi.line.info,3)
} 

#mean RMSE, R2, conc. index
apply(res,2, mean)

write.table(res, file="res_doubl_time_gcsi.txt")

# saving data for "predicted vs. observed" plots
pred <- predict.train(model_rf, newdata=gcsi_table[test_sampl,f_gam])
obs <- gcsi.line.info[test_sampl,3]
perform_table <- cbind(obs, pred)
write.table(perform_table, "perform_table_dt_gcsi.txt")

## II: NIBR PDXE ##

int <- intersect(slope_table[,1], rownames(mol_nibr))

res <- matrix(NA, 10,3)
colnames(res) <- c("RMSE", "R2", "concordance_ind")
for (i in 1:10){
  print(i)
  train_sampl <- sample(int, round(0.7*length(int)))
  test_sampl <- setdiff(int, train_sampl)
  
  f_gam <- names(sort(gam_fs(x=mol_nibr[train_sampl,], y=slope_table[train_sampl,2])))[1:400]
  model_rf <- simple_regression(y=slope_table[train_sampl,2], x=mol_nibr[train_sampl,f_gam], method="rf")
  res[i,] <- accuracy_check(model_rf, test_sampl, f_gam, mol_nibr, slope_table,2)
} 

#mean RMSE, R2, conc. index
apply(res,2, mean)

write.table(res, file="res_slope_nibr.txt")

# saving data for "predicted vs. observed" plots
pred <- predict.train(model_rf, newdata=mol_nibr[test_sampl,f_gam])
obs <- slope_table[test_sampl,2]
perform_table <- cbind(obs, pred)
write.table(perform_table, "perform_table_dt_nibr.txt")

###  III: drug response prediction  ###

## III: gCSI ##
int <- intersect(rownames(gcsi.line.info), rownames(gcsi.genomics))

int_erl <- intersect(int, rownames(gcsi.line.info)[which(gcsi.line.info$TissueMetaclass=="Lung")])
int_gem <- intersect(int, rownames(gcsi.line.info)[which(gcsi.line.info$TissueMetaclass=="Pancreas")])
int_pac_breast <- intersect(int, rownames(gcsi.line.info)[which(gcsi.line.info$TissueMetaclass=="Breast")])
int_pac_lung <- intersect(int, rownames(gcsi.line.info)[which(gcsi.line.info$TissueMetaclass=="Lung")])

res_erl <- matrix(NA, 10,3)
res_gem <- matrix(NA, 10,3)
res_pac_breast <- matrix(NA, 10,3)
res_pac_lung <- matrix(NA, 10,3)
colnames(res_erl) <- c("RMSE", "R2", "concordance_ind")
colnames(res_gem) <- c("RMSE", "R2", "concordance_ind")
colnames(res_pac_breast) <- c("RMSE", "R2", "concordance_ind")
colnames(res_pac_lung) <- c("RMSE", "R2", "concordance_ind")

for (i in 1:10){
  print(i)

  train_erl <- sample(int_erl, round(0.7*length(int_erl)))
  test_erl <- setdiff(int_erl, train_erl)
  train_gem <- sample(int_gem, round(0.7*length(int_gem)))
  test_gem <- setdiff(int_gem, train_gem)
  train_pac_breast <- sample(int_pac_breast, round(0.7*length(int_pac_breast)))
  test_pac_breast <- setdiff(int_pac_breast, train_pac_breast)
  train_pac_lung <- sample(int_pac_lung, round(0.7*length(int_pac_lung)))
  test_pac_lung <- setdiff(int_pac_lung, train_pac_lung)
  
  #feature selection
  f_erl <- names(sort(gam_fs(x=gcsi_table[train_erl,], y=gcsi_auc[train_erl,"Erlotinib"])))[1:100]
  f_gem <- names(sort(gam_fs(x=gcsi_table[train_gem,], y=gcsi_auc[train_gem,"Gemcitabine"])))[1:100]
  f_pac_breast <- names(sort(gam_fs(x=gcsi_table[train_pac_breast,], y=gcsi_auc[train_pac_breast,"paclitaxel"])))[1:100]
  f_pac_lung <- names(sort(gam_fs(x=gcsi_table[train_pac_lung,], y=gcsi_auc[train_pac_lung,"paclitaxel"])))[1:100]
  
  model_erl <- simple_regression(y=gcsi_auc[train_erl,"Erlotinib"], x=gcsi_table[train_erl,f_erl], method="rf")
  model_gem <- simple_regression(y=gcsi_auc[train_gem,"Gemcitabine"], x=gcsi_table[train_gem,f_gem], method="rf")
  model_pac_breast <- simple_regression(y=gcsi_auc[train_pac_breast,"paclitaxel"], x=gcsi_table[train_pac_breast,f_pac_breast], method="rf")
  model_pac_lung <- simple_regression(y=gcsi_auc[train_pac_lung,"paclitaxel"], x=gcsi_table[train_pac_lung,f_pac_lung], method="rf")
  
  res_erl[i,] <- accuracy_check(model_erl, test_erl, f_erl, gcsi_table, gcsi_auc, "Erlotinib")
  res_gem[i,] <- accuracy_check(model_gem, test_gem, f_gem, gcsi_table, gcsi_auc, "Gemcitabine")
  res_pac_breast[i,] <- accuracy_check(model_pac_breast, test_pac_breast, f_pac_breast, gcsi_table, gcsi_auc, "paclitaxel")
  res_pac_lung[i,] <- accuracy_check(model_pac_lung, test_pac_lung, f_pac_lung, gcsi_table, gcsi_auc, "paclitaxel")
} 

#mean RMSE, R2, conc. index
apply(res_erl,2, mean)
apply(res_gem,2, mean)
apply(res_pac_breast,2, mean)
apply(res_pac_lung,2, mean)

#saving the results
write.table(res_erl, file="res_erlotinib_gcsi.txt")
write.table(res_gem, file="res_gemcitabine_gcsi.txt")
write.table(res_pac_breast, file="res_paclitaxel_breast_gcsi.txt")
write.table(res_pac_lung, file="res_paclitaxel_lung_gcsi.txt")

# saving data for "predicted vs. observed" plots
pred_erl <- predict.train(model_erl, newdata=gcsi_table[test_erl,f_erl])
obs_erl <- gcsi_auc[test_erl,"Erlotinib"]
pred_gem <- predict.train(model_gem, newdata=gcsi_table[test_gem,f_gem])
obs_gem <- gcsi_auc[test_gem,"Gemcitabine"]
pred_pac_breast <- predict.train(model_pac_breast, newdata=gcsi_table[test_pac_breast,f_pac_breast])
obs_pac_breast <- gcsi_auc[test_pac_breast,"paclitaxel"]
pred_pac_lung <- predict.train(model_pac_lung, newdata=gcsi_table[test_pac_lung,f_pac_lung])
obs_pac_lung <- gcsi_auc[test_pac_lung,"paclitaxel"]

perform_table_erl <- cbind(obs_erl, pred_erl)
write.table(perform_table_erl, "perform_table_erl_gcsi.txt")
perform_table_gem <- cbind(obs_gem, pred_gem)
write.table(perform_table_gem, "perform_table_gem_gcsi.txt")
perform_table_pac_breast <- cbind(obs_pac_breast, pred_pac_breast)
write.table(perform_table_pac_breast, "perform_table_pac_breast_gcsi.txt")
perform_table_pac_lung <- cbind(obs_pac_lung, pred_pac_lung)
write.table(perform_table_pac_lung, "perform_table_pac_lung_gcsi.txt")

## III: NIBR PDXE ##
resp_erlotinib <- resp_nibr$erlotinib
resp_gemcitabine <- resp_nibr$gemcitabine
resp_paclitaxel <- resp_nibr$paclitaxel

int_erl <- intersect(rownames(resp_erlotinib), rownames(mol_nibr))
int_gem <- intersect(rownames(resp_gemcitabine), rownames(mol_nibr))
int_pac <- intersect(rownames(resp_paclitaxel), rownames(mol_nibr))
int_pac_brca <- intersect(int_pac, tissue_nibr$sample[which(tissue_nibr$tissue=="BRCA")])
int_pac_nsclc <- intersect(int_pac, tissue_nibr$sample[which(tissue_nibr$tissue=="NSCLC")])

res_erl <- matrix(NA, 10,12)
res_gem <- matrix(NA, 10,12)
res_pac_brca <- matrix(NA, 10,12)
res_pac_nsclc <- matrix(NA, 10,12)
columns <- c("vol_RMSE", "vol_R2", "vol_concordance_ind",
             "auc_RMSE", "auc_R2", "auc_concordance_ind",
             "slope_RMSE", "slope_R2", "slope_concordance_ind",
             "ds_RMSE", "ds_R2", "ds_concordance_ind")
colnames(res_erl) <- columns
colnames(res_gem) <- columns
colnames(res_pac_brca) <- columns
colnames(res_pac_nsclc) <- columns

fn <- 100
method <- "rf"
for (i in 1:10){
  print(i)
  train_erl <- sample(int_erl, round(0.7*length(int_erl)))
  test_erl <- setdiff(int_erl, train_erl)
  train_gem <- sample(int_gem, round(0.7*length(int_gem)))
  test_gem <- setdiff(int_gem, train_gem)
  train_pac_brca <- sample(int_pac_brca, round(0.7*length(int_pac_brca)))
  test_pac_brca <- setdiff(int_pac_brca, train_pac_brca)
  train_pac_nsclc <- sample(int_pac_nsclc, round(0.7*length(int_pac_nsclc)))
  test_pac_nsclc <- setdiff(int_pac_nsclc, train_pac_nsclc)
  
  #feature selection
  f_erl_vol <- names(sort(gam_fs(x=mol_nibr[train_erl,], y=resp_erlotinib[train_erl,1])))[1:400]
  f_erl_auc <- names(sort(gam_fs(x=mol_nibr[train_erl,], y=resp_erlotinib[train_erl,2])))[1:400]
  f_erl_slope <- names(sort(gam_fs(x=mol_nibr[train_erl,], y=resp_erlotinib[train_erl,3])))[1:400]
  f_erl_ds <- names(sort(gam_fs(x=mol_nibr[train_erl,], y=resp_erlotinib[train_erl,4])))[1:400]
  
  f_gem_vol <- names(sort(gam_fs(x=mol_nibr[train_gem,], y=resp_gemcitabine[train_gem,1])))[1:400]
  f_gem_auc <- names(sort(gam_fs(x=mol_nibr[train_gem,], y=resp_gemcitabine[train_gem,2])))[1:400]
  f_gem_slope <- names(sort(gam_fs(x=mol_nibr[train_gem,], y=resp_gemcitabine[train_gem,3])))[1:400]
  f_gem_ds <- names(sort(gam_fs(x=mol_nibr[train_gem,], y=resp_gemcitabine[train_gem,4])))[1:400]
  
  f_pac_brca_vol <- names(sort(gam_fs(x=mol_nibr[train_pac_brca,], y=resp_paclitaxel[train_pac_brca,1])))[1:400]
  f_pac_brca_auc <- names(sort(gam_fs(x=mol_nibr[train_pac_brca,], y=resp_paclitaxel[train_pac_brca,2])))[1:400]
  f_pac_brca_slope <- names(sort(gam_fs(x=mol_nibr[train_pac_brca,], y=resp_paclitaxel[train_pac_brca,3])))[1:400]
  f_pac_brca_ds <- names(sort(gam_fs(x=mol_nibr[train_pac_brca,], y=resp_paclitaxel[train_pac_brca,4])))[1:400]
  
  f_pac_nsclc_vol <- names(sort(gam_fs(x=mol_nibr[train_pac_nsclc,], y=resp_paclitaxel[train_pac_nsclc,1])))[1:400]
  f_pac_nsclc_auc <- names(sort(gam_fs(x=mol_nibr[train_pac_nsclc,], y=resp_paclitaxel[train_pac_nsclc,2])))[1:400]
  f_pac_nsclc_slope <- names(sort(gam_fs(x=mol_nibr[train_pac_nsclc,], y=resp_paclitaxel[train_pac_nsclc,3])))[1:400]
  f_pac_nsclc_ds <- names(sort(gam_fs(x=mol_nibr[train_pac_nsclc,], y=resp_paclitaxel[train_pac_nsclc,4])))[1:400]
  
  #model fitting
  m_erl_vol <- simple_regression(y=resp_erlotinib[train_erl,1], x=mol_nibr[train_erl,f_erl_vol[1:fn]] , method=method)
  m_erl_auc <- simple_regression(y=resp_erlotinib[train_erl,2], x=mol_nibr[train_erl,f_erl_auc[1:fn]] , method=method)
  m_erl_slope <- simple_regression(y=resp_erlotinib[train_erl,3], x=mol_nibr[train_erl,f_erl_slope[1:fn]] , method=method)
  m_erl_ds <- simple_regression(y=resp_erlotinib[train_erl,4], x=mol_nibr[train_erl,f_erl_ds[1:fn]] , method="rf")
  
  m_gem_vol <- simple_regression(y=resp_gemcitabine[train_gem,1], x=mol_nibr[train_gem,f_gem_vol[1:fn]] , method=method)
  m_gem_auc <- simple_regression(y=resp_gemcitabine[train_gem,2], x=mol_nibr[train_gem,f_gem_auc[1:fn]] , method=method)
  m_gem_slope <- simple_regression(y=resp_gemcitabine[train_gem,3], x=mol_nibr[train_gem,f_gem_slope[1:fn]] , method=method)
  m_gem_ds <- simple_regression(y=resp_gemcitabine[train_gem,4], x=mol_nibr[train_gem,f_gem_ds[1:fn]] , method=method)
  
  m_pac_brca_vol <- simple_regression(y=resp_paclitaxel[train_pac_brca,1], x=mol_nibr[train_pac_brca,f_pac_brca_vol[1:fn]] , method=method)
  m_pac_brca_auc <- simple_regression(y=resp_paclitaxel[train_pac_brca,2], x=mol_nibr[train_pac_brca,f_pac_brca_auc[1:fn]] , method=method)
  m_pac_brca_slope <- simple_regression(y=resp_paclitaxel[train_pac_brca,3], x=mol_nibr[train_pac_brca,f_pac_brca_slope[1:fn]] , method=method)
  m_pac_brca_ds <- simple_regression(y=resp_paclitaxel[train_pac_brca,4], x=mol_nibr[train_pac_brca,f_pac_brca_ds[1:fn]] , method=method)
  
  m_pac_nsclc_vol <- simple_regression(y=resp_paclitaxel[train_pac_nsclc,1], x=mol_nibr[train_pac_nsclc,f_pac_nsclc_vol[1:fn]] , method=method)
  m_pac_nsclc_auc <- simple_regression(y=resp_paclitaxel[train_pac_nsclc,2], x=mol_nibr[train_pac_nsclc,f_pac_nsclc_auc[1:fn]] , method=method)
  m_pac_nsclc_slope <- simple_regression(y=resp_paclitaxel[train_pac_nsclc,3], x=mol_nibr[train_pac_nsclc,f_pac_nsclc_slope[1:fn]] , method=method)
  m_pac_nsclc_ds <- simple_regression(y=resp_paclitaxel[train_pac_nsclc,4], x=mol_nibr[train_pac_nsclc,f_pac_nsclc_ds[1:fn]] , method=method)
  
  #results
  res_erl[i,] <- c(accuracy_check(m_erl_vol, test_erl, f_erl_vol, mol_nibr, resp_erlotinib,1),
                   accuracy_check(m_erl_auc, test_erl, f_erl_auc, mol_nibr, resp_erlotinib,2),
                   accuracy_check(m_erl_slope, test_erl, f_erl_slope, mol_nibr, resp_erlotinib,3),
                   accuracy_check(m_erl_ds, test_erl, f_erl_ds, mol_nibr, resp_erlotinib,4))
  # remove problematic X-4215 sample
  s4215 <- which(test_gem=="X-4215")
  if (length(s4215)>0){test_gem <- test_gem[-s4215]} 
  
  res_gem[i,] <- c(accuracy_check(m_gem_vol, test_gem, f_gem_vol, mol_nibr, resp_gemcitabine,1),
                   accuracy_check(m_gem_auc, test_gem, f_gem_auc, mol_nibr, resp_gemcitabine,2),
                   accuracy_check(m_gem_slope, test_gem, f_gem_slope, mol_nibr, resp_gemcitabine,3),
                   accuracy_check(m_gem_ds, test_gem, f_gem_ds, mol_nibr, resp_gemcitabine,4))
  
  res_pac_brca[i,] <- c(accuracy_check(m_pac_brca_vol, test_pac_brca, f_pac_brca_vol, mol_nibr, resp_paclitaxel,1),
                        accuracy_check(m_pac_brca_auc, test_pac_brca, f_pac_brca_auc, mol_nibr, resp_paclitaxel,2),
                        accuracy_check(m_pac_brca_slope, test_pac_brca, f_pac_brca_slope, mol_nibr, resp_paclitaxel,3),
                        accuracy_check(m_pac_brca_ds, test_pac_brca, f_pac_brca_ds, mol_nibr, resp_paclitaxel,4))
  
  res_pac_nsclc[i,] <- c(accuracy_check(m_pac_nsclc_vol, test_pac_nsclc, f_pac_nsclc_vol, mol_nibr, resp_paclitaxel,1),
                         accuracy_check(m_pac_nsclc_auc, test_pac_nsclc, f_pac_nsclc_auc, mol_nibr, resp_paclitaxel,2),
                         accuracy_check(m_pac_nsclc_slope, test_pac_nsclc, f_pac_nsclc_slope, mol_nibr, resp_paclitaxel,3),
                         accuracy_check(m_pac_nsclc_ds, test_pac_nsclc, f_pac_nsclc_ds, mol_nibr, resp_paclitaxel,4))
} 

#mean RMSE, R2, conc. index
apply(res_erl,2, mean)
apply(res_gem,2, mean)
apply(res_pac_brca,2, mean)
apply(res_pac_nsclc,2, mean)

#saving the results
write.table(res_erl, file="res_erlotinib_nibr.txt")
write.table(res_gem, file="res_gemcitabine_nibr.txt")
write.table(res_pac_brca, file="res_paclitaxel_brca_nibr.txt")
write.table(res_pac_nsclc, file="res_paclitaxel_nsclc_nibr.txt")

# saving data for "predicted vs. observed" plots

pred_erl_vol <- predict.train(m_erl_vol, newdata=mol_nibr[test_erl,f_erl_vol])
obs_erl_vol <- resp_erlotinib[test_erl,1]
pred_erl_slope <- predict.train(m_erl_slope, newdata=mol_nibr[test_erl,f_erl_slope])
obs_erl_slope <- resp_erlotinib[test_erl,3]

pred_gem_vol <- predict.train(m_gem_vol, newdata=mol_nibr[test_gem,f_gem_vol])
obs_gem_vol <- resp_gemcitabine[test_gem,1]
pred_gem_slope <- predict.train(m_gem_slope, newdata=mol_nibr[test_gem,f_gem_slope])
obs_gem_slope <- resp_gemcitabine[test_gem,3]

pred_pac_brca_vol <- predict.train(m_pac_brca_vol, newdata=mol_nibr[test_pac_brca,f_pac_brca_vol])
obs_pac_brca_vol <- resp_paclitaxel[test_pac_brca,1]
pred_pac_brca_slope <- predict.train(m_pac_brca_slope, newdata=mol_nibr[test_pac_brca,f_pac_brca_slope])
obs_pac_brca_slope <- resp_paclitaxel[test_pac_brca,3]

pred_pac_nsclc_vol <- predict.train(m_pac_nsclc_vol, newdata=mol_nibr[test_pac_nsclc,f_pac_nsclc_vol])
obs_pac_nsclc_vol <- resp_paclitaxel[test_pac_nsclc,1]
pred_pac_nsclc_slope <- predict.train(m_pac_nsclc_slope, newdata=mol_nibr[test_pac_nsclc,f_pac_nsclc_slope])
obs_pac_nsclc_slope <- resp_paclitaxel[test_pac_nsclc,3]

perform_table_erl <- cbind(obs_erl_vol, pred_erl_vol, obs_erl_slope, pred_erl_slope)
write.table(perform_table_erl, "perform_table_erl_nibr.txt")

perform_table_gem <- cbind(obs_gem_vol, pred_gem_vol, obs_gem_slope, pred_gem_slope)
write.table(perform_table_gem, "perform_table_gem_nibr.txt")

perform_table_pac_brca <- cbind(obs_pac_brca_vol, pred_pac_brca_vol, obs_pac_brca_slope, pred_pac_brca_slope)
write.table(perform_table_pac_brca, "perform_table_pac_brca_nibr.txt")

perform_table_pac_nsclc <- cbind(obs_pac_nsclc_vol, pred_pac_nsclc_vol, obs_pac_nsclc_slope, pred_pac_nsclc_slope)
write.table(perform_table_pac_nsclc, "perform_table_pac_nsclc_nibr.txt")
