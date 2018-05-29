# 27.05.2018 script with functions

# for feature selection
gam_fs <- function(x, y) {
  p_vect <- apply(x,2, function(z) {gamScores(z, y)})
  return(p_vect)
}

anova_fs <- function(x, y) {
  p_vect <- apply(x,2, function(z) {anovaScores(z, y)})
  return(p_vect)
}

# for model fitting
simple_classification <- function(y, x, method) {
  
  ctrl <- trainControl(method = "cv", 
                       number = 10, 
                       search = "random", 
                       classProbs = T,
                       summaryFunction = multiClassSummary)
  
  na <- which(is.na(y) | y=="NA")
  if(length(na)>0) {
    x <- as.matrix(x[-na,])
    y <- y[-na]
  } else {  x <- as.matrix(x)}
  
  print(dim(x))
  mode(x) <- "numeric"
  
  print(summary(as.factor(y)))
  y <- factor(y, levels=unique(y))
  
  Tune <- train(x, y,
                method = method,
                tuneLength = 30,
                trControl = ctrl,
                metric="Accuracy",
                maximize=T,
                preProc = c("center","scale","medianImpute"))
  
  return(Tune)
}

simple_regression <- function(y, x, method) {
  ctrl <- trainControl(method = "cv", number = 10, search = "random")
  na <- which(is.na(y))
  if(length(na)>0) {
    x <- as.matrix(x[-na,])
    y <- y[-na]
  } else {  x <- as.matrix(x)}
  
  print(dim(x))
  Tune <- train(x, y,
                method = method,
                tuneLength = 30,
                trControl = ctrl,
                preProc = c("center","scale","medianImpute"))
  
  return(Tune)
}

# concordance.index 
ci <- function(v1, v2) {
  value <- concordance.index(x=v1, surv.time=v2, surv.event=rep(1,length(v1)))
  return(1-value$c.index)
}

# for accuracy calculations

#calssification
accuracy_check_class <- function(model, test_sampl, f_set, mol_data, class_data, class_column){
  pred <- predict.train(object=model, newdata=mol_data[test_sampl,f_set], type="raw")
  obs <- class_data[test_sampl,class_column]
  df <- cbind(as.character(pred),as.character(obs))
  compare <- apply(df,1, function(x){x[1]==x[2]})
  return(length(which(compare))/length(compare))
}

# classification for the script "tissue_classification_equal_greoups.R"
accuracy_check_class2 <- function(model, test_sampl, f_set, mol_data, class_data, class_column){
 
  pred <- predict.train(object=model, newdata=mol_data[test_sampl,f_set], type="raw")
  obs <- class_data[test_sampl,class_column]
  df <- cbind(as.character(pred),as.character(obs))
  compare <- apply(df,1, function(x){x[1]==x[2]})
  acc <- length(which(compare))/length(compare)
  acc_t <- vector()
  for (i in 1:length(tissues)){
    compare_i <- compare[which(obs==tissues[i])]
    acc_t[i] <- length(which(compare_i))/length(compare_i)
  }
  return(list(c(acc, acc_t), df))
}

#regression
accuracy_check <- function(model, test_sampl, f_gam, mol_data, obs_data, obs_column) {
  pred <- predict.train(model, newdata=mol_data[test_sampl,f_gam])
  obs <- obs_data[test_sampl,obs_column]
  ci <- ci(pred,obs)
  return(c(postResample(pred=pred, obs=obs)[c("RMSE","Rsquared")], ci))
}