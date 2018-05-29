# 26.01.18 analysis script
# edited on 19.04.18

source("http://bioconductor.org/biocLite.R")
biocLite("survcomp")
install.packages("caret")
install.packages("gam")
library(survcomp)
library(caret)
library(gam)


path <- "/data/kurilov/genestack/phd/work_2018/DrugRespPrediction/"
setwd(file.path(path, "analysis1"))

# get the processed data for modelling
load("ccle_table.RData")
load("ctrp_table.RData")
load("gdsc_table.RData")

## funtions for model fitting and calculation the accuracy metrics
# model fitting via caret
simple_regression <- function(x, y, method) {
  
  ctrl <- trainControl(method = "cv", number = 10, search = "random")
  na <- which(is.na(y) | is.infinite(y))
  if(length(na)>0) {
    x <- as.matrix(x[-na,])
    y <- y[-na]
  } else {  x <- as.matrix(x)}
  
  print(dim(x))
  Tune <- train(x, y,
                method = method,
                tuneLength = 30,
                trControl = ctrl,
                preProc = c("medianImpute","center","scale"))
  
  return(Tune)
}

# feature selection with gamScores
gam_fs <- function(x, y) {
  inf_y <- which(is.infinite(y))
  if (length(inf_y)>0){
    y <- y[-inf_y]
    x <- x[-inf_y,]
  }
  p_vect <- apply(x,2, function(z) {gamScores(z, y)})
  return(p_vect)
}

# concordance index
ci <- function(v1, v2) {
  value <- concordance.index(x=v1, surv.time=v2, surv.event=rep(1,length(v1)))
  return(1-value$c.index)
}

# function for saving predicted values and accuracy calculation (R2, conc. index) overall and for each tissue separately
accuracy_check <- function(model, test, f_set) {
  raw_pred_index <<- raw_pred_index +1
  pred <- predict.train(model, newdata=table[test,f_set])
  names(pred) <- test
  raw_predictions[test,raw_pred_index] <<- pred
  
  obs <- y[test]
  inf_obs <- which(is.infinite(obs) | is.na(obs))
  if (length(inf_obs)>0){
    obs <- obs[-inf_obs]
    pred <- pred[-inf_obs]
  }
  ci <- ci(pred,obs)
  overall_acc <- c(postResample(pred=pred, obs=obs)[c("RMSE","Rsquared")], ci)
  r2_all_tissues <- rep(NA,18)
  rmse_all_tissues <- rep(NA,18)
  ci_all_tissues <- rep(NA,18)
  for (i in 1:18)
  {
    tissue <- tissue_list[i]
    lines_t <- intersect(names(pred), tissue_table$cellid[which(tissue_table$tissueid==tissue)])
    pred_t <- pred[lines_t]
    obs_t <- obs[lines_t]
    if (length(pred_t>1)){
      res_t <- postResample(pred=pred_t, obs=obs_t)[c("RMSE","Rsquared")]
      r2_t <- res_t[2]
      rmse_t <- res_t[1]
      ci_t <- ci(pred_t,obs_t)
    } else {
      r2_t <- NA
      rmse_t <- NA
      ci_t <- NA}
    r2_all_tissues[i] <- r2_t
    rmse_all_tissues[i] <- rmse_t
    ci_all_tissues[i] <- ci_t }
    
  
  return(c(raw_pred_index, overall_acc,rmse_all_tissues,r2_all_tissues, ci_all_tissues))
}

# overview of the analysis' loops:
# (3)dataset -> (7)drug -> (3)metric -> FS -> (2)var_type -> (4)var_number 

#creating a list with major tissues (by number of cell lines)
tissue_list <- c("lung", "haematopoietic_and_lymphoid_tissue", "skin", "central_nervous_system",
                 "breast", "large_intestine", "pancreas", "ovary", "stomach", "oesophagus", 
                 "upper_aerodigestive_tract", "bone", "urinary_tract", "endometrium", 
                 "autonomic_ganglia", "liver", "soft_tissue", "kidney")

# table for results and table for storing predicted values
results_table <- matrix(NA, ncol=73, nrow=504)
colnames(results_table) <- c("s_size","dataset", "drug", "metric", "var_type", "var_number", "rmse", "r2","conc_index" ,
                             paste0("rmse_", tissue_list),paste0("r2_",tissue_list), paste0("ci_",tissue_list),
                             paste0("raw",1:10))
raw_predictions <- matrix(NA, ncol=5040, nrow=1119)
rownames(raw_predictions) <- unique(c(rownames(ccle_table), rownames(ctrp_table), rownames(gdsc_table)))
raw_pred_index <- 0

# lists for loops
datasets <- c("ccle", "ctrp", "gdsc")
drugs <- c("Erlotinib",  "Paclitaxel", "Lapatinib",  "Nilotinib", "Nutlin-3","PLX4720", "Sorafenib")
metrics <- c("ic50", "auc", "viab")

# for parallel calculation
#library(doMC)
#registerDoMC(cores = 4)

# main analysis loop
for (i in 1:3){
  dataset <- datasets[i]
  for (j in 1:7){
    drug <- drugs[j]
    
    for (k in 1:3){
      metric <- metrics[k]
      
      print(paste0(dataset," ",drug," ",metric))
      
      table <- eval(parse(text = paste0(dataset, "_table")))
      table <- table[,c(1:(min(grep("ic50", colnames(table)))-1),grep(paste0(metric,"_",drug), colnames(table)))]
      fin <- which(is.finite(table[,dim(table)[2]]))
      table <- table[fin,]
      x <- apply(table,2,function(x){sd(x, na.rm=T)})
      sdz <- which(x==0)
      if (length(sdz)>0){table <- table[,-sdz]}
      
      y <- table[,dim(table)[2]]
      names(y) <- rownames(table)
      table <- table[,-dim(table)[2]]
      
      inter_results <- matrix(NA, nrow=10, ncol=456)
      raw_data_inf <- matrix(NA, nrow=10, ncol=8)
      for (l in 1:10){
        print(l)
        train <- sample(rownames(table), round(0.7*length(rownames(table))))
        test <- setdiff(rownames(table), train)
       
        #gam
        f_set <- gam_fs(x=table[train,], y=y[train])
        f_set_exp <- f_set[1:max(grep("exp_",colnames(table)))]
        f_set_exp <- sort(f_set_exp)
        fs1 <- names(f_set_exp)[1:10]
        fs2 <- names(f_set_exp)[1:50]
        fs3 <- names(f_set_exp)[1:200]
        fs4 <- names(f_set_exp)[1:500]
        f_set <- sort(f_set)
        fs5 <- names(f_set)[1:10]
        fs6 <- names(f_set)[1:50]
        fs7 <- names(f_set)[1:200]
        fs8 <- names(f_set)[1:500]
        
        m1 <- simple_regression(x=table[train,fs1], y=y[train], method="svmRadial")
        m2 <- simple_regression(x=table[train,fs2], y=y[train], method="svmRadial")
        m3 <- simple_regression(x=table[train,fs3], y=y[train], method="svmRadial")
        m4 <- simple_regression(x=table[train,fs4], y=y[train], method="svmRadial")
        m5 <- simple_regression(x=table[train,fs5], y=y[train], method="svmRadial")
        m6 <- simple_regression(x=table[train,fs6], y=y[train], method="svmRadial")
        m7 <- simple_regression(x=table[train,fs7], y=y[train], method="svmRadial")
        m8 <- simple_regression(x=table[train,fs8], y=y[train], method="svmRadial")
        
        print(paste0("raw_pred_index = ",raw_pred_index))
        res1 <- accuracy_check(m1, test, fs1)
        res2 <- accuracy_check(m2, test, fs2)
        res3 <- accuracy_check(m3, test, fs3)
        res4 <- accuracy_check(m4, test, fs4)
        res5 <- accuracy_check(m5, test, fs5)
        res6 <- accuracy_check(m6, test, fs6)
        res7 <- accuracy_check(m7, test, fs7)
        res8 <- accuracy_check(m8, test, fs8)
        print(paste0("raw_pred_index = ",raw_pred_index))
        
        inter_results[l,] <- c(res1[-1],res2[-1],res3[-1],res4[-1],res5[-1],res6[-1],res7[-1],res8[-1])
        raw_data_inf[l,] <- c(res1[1],res2[1],res3[1],res4[1],res5[1],res6[1],res7[1],res8[1])
      }
    inter_results2 <- apply(inter_results,2,function(x){mean(x, na.rm=T)})
    inter_results2 <- matrix(inter_results2, nrow=8, byrow=T)
    
    num <- (i-1)*168 + (j-1)*24 + (k-1)*8 +1
    results_table[num:(num+7),] <- cbind(rep(length(train),8),rep(dataset,8), rep(drug,8),rep(metric,8),c(rep("exp",4),rep("all",4)),rep(c(10,50,200,500),2),
                                                             inter_results2,t(raw_data_inf))
    print(results_table[num:(num+7),])
    }
  }
}

# saving results and predicted values
save(results_table, file="results_table.RData")
save(raw_predictions, file="raw_predictions.RData")
