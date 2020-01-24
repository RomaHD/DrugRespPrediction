# 13.02.2018 collecting top 500 features + p-values
# edited on 20.04.18

path <- "/abi/data/kurilov/work_2018/DrugRespPrediction/"
setwd(file.path(path, "analysis1"))

library(survcomp)
library(caret)
library(ggplot2)

load("ccle_table.RData")
load("ctrp_table.RData")
load("gdsc_table.RData")

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

drugs <- c("Erlotinib",  "Paclitaxel", "Lapatinib",  "Nilotinib", "Nutlin-3","PLX4720", "Sorafenib")
datasets <- c("ccle", "ctrp", "gdsc")
metrics <- c("ic50","auc", "viab")

features <- matrix(NA, ncol=63, nrow=500)
pv <- matrix(NA, ncol=63, nrow=500)

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
  
  #gam
  f_set <- gam_fs(x=table, y=y)
  fs <- sort(f_set)[1:500]

  features[,i] <- names(fs)
  pv[,i] <- fs
  
}

colnames(features) <- apply(tasks,1, function(x){paste(x, collapse="_")})
features_inf <- list("tasks"=tasks, "features"=features, "p.values"=pv)
save(features_inf, file="features_inf.RData")

# number of top-500 features in overlop between CCLE, CTRP and GDSC for each metric-drug combination
f_number <- vector()
for (i in c(1:7, 22:28, 43:49)){
f_number_i <- length(intersect(intersect(features[,i], features[,i+7]),features[,i+14]))
f_number <- c(f_number, f_number_i)
}

table <- cbind(f_number, rep(drugs,3), rep(metrics, each=7))
table <- as.data.frame(table)
table$f_number <- as.numeric(as.character(table$f_number))
colnames(table) <- c("number", "drug", "metric")
table$drug <- factor(table$drug, levels=drugs)
table$metric <- factor(table$metric, levels=metrics)

theme_set(theme_bw(base_size = 18))
g1 <- ggplot(data=table, aes(x=drug, y=number)) + geom_bar(aes(fill=metric), stat="identity", position="dodge") + labs(y="# of common features",x=NULL) + 
      geom_hline(yintercept=42, colour="red", linetype="dotted")+ theme(legend.text=element_text(size=12), legend.position = "bottom")
ggsave("fig_2c.pdf", plot=g1)
