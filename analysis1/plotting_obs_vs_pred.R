#4.02.2018, observed vs. predicted values: analysis I
# edited on 20.04.18

#path <- "/data/kurilov/genestack/phd/work_2018/DrugRespPrediction/"
path <- "/abi/data/kurilov/work_2018/DrugRespPrediction/"
setwd(file.path(path, "analysis1"))

library(lattice)

load("ccle_table.RData")
load("ctrp_table.RData")
load("gdsc_table.RData")

load("raw_predictions2.RData")

drugs <- c("Erlotinib",  "Paclitaxel", "Lapatinib",  "Nilotinib", "Nutlin-3","PLX4720", "Sorafenib")
datasets <- c("ccle", "ctrp", "gdsc")
tasks0 <- cbind(rep(c("ic50","auc", "viab"),times=1,each=21), rep(drugs,9), rep(datasets, times=3, each=7))
colnames(tasks0) <- c("metric","drug","dataset")
tasks0 <- as.data.frame(tasks0)

colnames(raw_predictions) <- apply(tasks0, 1, function(x){paste(x, collapse="_")})

# 1. Plotting Figure 3. with lattice
tasks <- cbind(rep(c("ic50","auc", "viab"),times=1,each=3),
               rep(c("Paclitaxel","Nutlin-3","Sorafenib"), 3),
               rep("ctrp",9))
colnames(tasks) <- c("metric","drug","dataset")
tasks <- as.data.frame(tasks)

num <- match(apply(tasks, 1, function(x){paste(x, collapse="_")}), apply(tasks0, 1, function(x){paste(x, collapse="_")})) 
raw_predictions_sel <- raw_predictions[,num]

#for concordance index calculation
library(survcomp)
ci <- function(v1, v2) {
  value <- concordance.index(x=v1, surv.time=v2, surv.event=rep(1,length(v1)))
  return(1-value$c.index)
}




table_lat <- matrix(NA, ncol=3,nrow=0) 
for (i in 1:9)
{
  n <- which(!(is.na(raw_predictions_sel[,i])))
  table <- eval(parse(text = paste0(tasks$dataset[i], "_table")))
  
  ds <- cbind(table[rownames(raw_predictions_sel)[n],paste(tasks$metric[i],tasks$drug[i], sep="_")],raw_predictions_sel[n,i])
  r2 <- signif((cor(ds[,1],ds[,2]))^2, digits=2)
  ci_val <- signif(ci(ds[,1],ds[,2]), digits = 2)
  title <- paste(unlist(tasks[i,1:2]), collapse=", ")
  title <- paste0(title, "\nR2=",r2, "   conc. index=",ci_val)
  ds <- cbind(ds, title)
  
  table_lat <- rbind(table_lat, ds)
}

colnames(table_lat) <- c("observed","predicted","label")
table_lat <- as.data.frame(table_lat)
table_lat[,1] <- as.numeric(as.character(table_lat[,1]))
table_lat[,2] <- as.numeric(as.character(table_lat[,2]))

#changing titles
table_lat$label <- gsub("ic50", "IC50", table_lat$label)
table_lat$label <- gsub("auc", "AUC", table_lat$label)
table_lat$label <- gsub("viab", "Viability_1uM", table_lat$label)
#table_lat$label <- gsub("ccle", "CCLE", table_lat$label)
#table_lat$label <- gsub("ctrp", "CTRP", table_lat$label)
#table_lat$label <- gsub("gdsc", "GDSC", table_lat$label)
levels(table_lat$label) <- unique(table_lat$label)

pdf(file="fig_3.pdf", height=8.3, width=11.7)
xyplot(predicted ~ observed | label, table_lat,
       grid = TRUE,
       as.table=T,
       scales=list(relation="free"),
       layout=c(3,3),
       par.strip.text = list(cex=1.3, lines=2),
       #skip=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
       xlab="observed drug response values",
       ylab="predicted drug response values",
       index.cond=list(c(4,5,6,1,2,3,7,8,9)),
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)  ## First call the default panel function for 'xyplot'
         panel.abline(a=c(0,1), lty=2)  ## Add a line
       })
dev.off()

# # 1. Plotting Figure 3. with lattice
# tasks <- cbind(rep(c("ic50","auc", "viab"),times=1,each=3),
#                c("Paclitaxel", "Paclitaxel","Nilotinib","Nutlin-3","Nutlin-3","Sorafenib","Paclitaxel","Paclitaxel","Sorafenib"),
#                c("ctrp","gdsc","ctrp","ccle","ctrp","gdsc","gdsc","ccle","gdsc"))
# colnames(tasks) <- c("metric","drug","dataset")
# tasks <- as.data.frame(tasks)
# 
# raw_predictions_sel <- raw_predictions[,c(9,16,11,26,33,42,58,44,63)]
# 
# #for concordance index calculation
# library(survcomp)
# ci <- function(v1, v2) {
#   value <- concordance.index(x=v1, surv.time=v2, surv.event=rep(1,length(v1)))
#   return(1-value$c.index)
# }
# 
# 
# 
# 
# table_lat <- matrix(NA, ncol=3,nrow=0) 
# for (i in 1:9)
# {
#   n <- which(!(is.na(raw_predictions_sel[,i])))
#   table <- eval(parse(text = paste0(tasks$dataset[i], "_table")))
#   
#   ds <- cbind(table[rownames(raw_predictions_sel)[n],paste(tasks$metric[i],tasks$drug[i], sep="_")],raw_predictions_sel[n,i])
#   r2 <- signif((cor(ds[,1],ds[,2]))^2, digits=2)
#   ci_val <- signif(ci(ds[,1],ds[,2]), digits = 2)
#   title <- paste(unlist(tasks[i,]), collapse=", ")
#   title <- paste0(title, "\nR2=",r2, "   conc. index=",ci_val)
#   ds <- cbind(ds, title)
#   
#   table_lat <- rbind(table_lat, ds)
# }
# 
# colnames(table_lat) <- c("observed","predicted","label")
# table_lat <- as.data.frame(table_lat)
# table_lat[,1] <- as.numeric(as.character(table_lat[,1]))
# table_lat[,2] <- as.numeric(as.character(table_lat[,2]))
# 
# #changing titles
# table_lat$label <- gsub("ic50", "IC50", table_lat$label)
# table_lat$label <- gsub("auc", "AUC", table_lat$label)
# table_lat$label <- gsub("viab", "Viability_1uM", table_lat$label)
# table_lat$label <- gsub("ccle", "CCLE", table_lat$label)
# table_lat$label <- gsub("ctrp", "CTRP", table_lat$label)
# table_lat$label <- gsub("gdsc", "GDSC", table_lat$label)
# 
# pdf(file="fig_3.pdf", height=8.3, width=11.7)
# xyplot(predicted ~ observed | label, table_lat,
#        grid = TRUE,
#        as.table=T,
#        scales=list(relation="free"),
#        layout=c(3,3),
#        par.strip.text = list(cex=0.9, lines=2),
#        #skip=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
#        index.cond=list(c(5,6,4,1,2,3,7,8,9)),
#        panel = function(x, y, ...) {
#          panel.xyplot(x, y, ...)  ## First call the default panel function for 'xyplot'
#          panel.abline(a=c(0,1), lty=2)  ## Add a line
#        })
# dev.off()


# 2. Plotting Figures s5. with lattice
drugs <- c("Erlotinib",  "Paclitaxel", "Lapatinib",  "Nilotinib", "Nutlin-3","PLX4720", "Sorafenib")
datasets <- c("ccle", "ctrp", "gdsc")

tasks <- cbind(rep(c("ic50","auc", "viab"),times=1,each=21), rep(drugs,9), rep(datasets, times=3, each=7))
colnames(tasks) <- c("metric","drug","dataset")
tasks <- as.data.frame(tasks)

table_lat <- matrix(NA, ncol=3,nrow=0) 
for (i in 1:63)
{
  n <- which(!(is.na(raw_predictions[,i])))
  table <- eval(parse(text = paste0(tasks$dataset[i], "_table")))
  
  ds <- cbind(table[rownames(raw_predictions)[n],paste(tasks$metric[i],tasks$drug[i], sep="_")],raw_predictions[n,i])
  r2 <- signif((cor(ds[,1],ds[,2]))^2, digits=2)
  title <- paste(unlist(tasks[i,]), collapse=", ")
  title <- paste0(title, "\nR2=",r2)
  ds <- cbind(ds, title)
  
  table_lat <- rbind(table_lat, ds)
}

colnames(table_lat) <- c("observed","predicted","label")
table_lat <- as.data.frame(table_lat)
table_lat[,1] <- as.numeric(as.character(table_lat[,1]))
table_lat[,2] <- as.numeric(as.character(table_lat[,2]))

max_ic50 <- max(grep("ic50", table_lat$label))
max_auc <- max(grep("auc", table_lat$label))
max_viab <- max(grep("viab", table_lat$label))

table_lat_ic50 <- table_lat[1:max_ic50,]
table_lat_auc <- table_lat[(1+max_ic50):max_auc,]
table_lat_viab <- table_lat[(1+max_auc):max_viab,]

table_lat_ic50$label <- substr(table_lat_ic50$label, 7,30)
table_lat_auc$label <- substr(table_lat_auc$label, 6,30)
table_lat_viab$label <- substr(table_lat_viab$label, 7,30)

pdf(file="fig_s5_1.pdf", height=8.3, width=11.7)
xyplot(predicted ~ observed | label, table_lat_ic50,
       main="IC50",
       grid = TRUE,
       as.table=T,
       scales=list(relation="free"),
       layout=c(7,3),
       par.strip.text = list(cex=0.9, lines=2),
       #skip=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
       index.cond=list(c(1,4,7,10,13,16,19,2,5,8,11,14,17,20,3,6,9,12,15,18,21)),
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)  ## First call the default panel function for 'xyplot'
         panel.abline(a=c(0,1), lty=2)  ## Add a horizontal line at the median
       })
dev.off()

pdf(file="fig_s5_2.pdf", height=8.3, width=11.7)
xyplot(predicted ~ observed | label, table_lat_auc,
       main="AUC",
       grid = TRUE,
       as.table=T,
       #scales=list(relation="free"),
       layout=c(7,3),
       par.strip.text = list(cex=0.9, lines=2),
       #skip=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
       index.cond=list(c(1,4,7,10,13,16,19,2,5,8,11,14,17,20,3,6,9,12,15,18,21)),
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)  ## First call the default panel function for 'xyplot'
         panel.abline(a=c(0,1), lty=2)  ## Add a horizontal line at the median
       })
dev.off()

pdf(file="fig_s5_3.pdf", height=8.3, width=11.7)
xyplot(predicted ~ observed | label, table_lat_viab,
       main="Viabilitu_1_uM",
       grid = TRUE,
       as.table=T,
       #scales=list(relation="free"),
       layout=c(7,3),
       par.strip.text = list(cex=0.9, lines=2),
       #skip=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
       index.cond=list(c(1,4,7,10,13,16,19,2,5,8,11,14,17,20,3,6,9,12,15,18,21)),
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)  ## First call the default panel function for 'xyplot'
         panel.abline(a=c(0,1), lty=2)  ## Add a horizontal line at the median
       })
dev.off()

