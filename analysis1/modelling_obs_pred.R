# rerunning some modelling in order to produce observed/predicted plots
# edited on 20.04.18

path <- "/data/kurilov/genestack/phd/work_2018/DrugRespPrediction/"
setwd(file.path(path, "analysis1"))

library(survcomp)
library(caret)

load("ccle_table.RData")
load("ctrp_table.RData")
load("gdsc_table.RData")

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

gam_fs <- function(x, y) {
  inf_y <- which(is.infinite(y))
  if (length(inf_y)>0){
    y <- y[-inf_y]
    x <- x[-inf_y,]
  }
  p_vect <- apply(x,2, function(z) {gamScores(z, y)})
  return(p_vect)
}

drugs <- c("Erlotinib",  "Paclitaxel", "Lapatinib",  "Nilotinib", "Nutlin-3","PLX4720", "Sorafenib")
datasets <- c("ccle", "ctrp", "gdsc")

raw_predictions <- matrix(NA, ncol=63, nrow=1119)
rownames(raw_predictions) <- unique(c(rownames(ccle_table), rownames(ctrp_table), rownames(gdsc_table)))
tasks <- cbind(rep(c("ic50","auc", "viab"),times=1,each=21), rep(drugs,9), rep(datasets, times=3, each=7))
colnames(tasks) <- c("metric","drug","dataset")
tasks <- as.data.frame(tasks)

for (i in 1:63){
print(i)
  dataset <- tasks$dataset[i]
  metric <- tasks$metric[i]
  drug <- tasks$drug[i]
  
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
  
    
    train <- sample(rownames(table), round(0.7*length(rownames(table))))
    test <- setdiff(rownames(table), train)
    
    #gam
    f_set <- gam_fs(x=table[train,], y=y[train])
    fs8 <- names(sort(f_set))[1:500]
    m8 <- simple_regression(x=table[train,fs8], y=y[train], method="svmRadial")
    
    pred <- predict.train(m8, newdata=table[test,fs8])
    names(pred) <- test
    raw_predictions[test,i] <- pred

}

save(raw_predictions, file="raw_predictions2.RData")

