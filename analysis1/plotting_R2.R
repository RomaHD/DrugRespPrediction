# 2.02.18 plotting analysis I results
# edited on 19.04.18

library(ggplot2)
library(plyr)
library(gridExtra)

path <- "/data/kurilov/genestack/phd/work_2018/DrugRespPrediction/"
setwd(file.path(path, "analysis1"))
load("results_table.RData")

results_table <- as.data.frame(results_table)
results_table$s_size <- as.numeric(as.character(results_table$s_size))
results_table$var_number <- as.numeric(as.character(results_table$var_number))
results_table$rmse <- as.numeric(as.character(results_table$rmse))
results_table$r2 <- as.numeric(as.character(results_table$r2))
results_table$s_size <- as.numeric(as.character(results_table$s_size))
results_table$conc_index <- as.numeric(as.character(results_table$conc_index))
results_table$dataset <- revalue(results_table$dataset, c(ccle="CCLE", ctrp="CTRP", gdsc="GDSC"))
results_table$metric <- revalue(results_table$metric, c(ic50="IC50", auc="AUC", viab="viability_1uM"))

drugs <- c("Erlotinib",  "Paclitaxel", "Lapatinib",  "Nilotinib", "Nutlin-3","PLX4720", "Sorafenib")
metrics <- c("IC50", "AUC", "viability_1uM")
datasets <- c("CCLE", "CTRP", "GDSC")
tissue_list <- c("lung", "haematopoietic_and_lymphoid_tissue", "skin", "central_nervous_system",
                 "breast", "large_intestine", "pancreas", "ovary", "stomach", "oesophagus", 
                 "upper_aerodigestive_tract", "bone", "urinary_tract", "endometrium", 
                 "autonomic_ganglia", "liver", "soft_tissue", "kidney")

results_table$drug <- factor(results_table$drug, levels=drugs)
results_table$metric <- factor(results_table$metric, levels=metrics)
results_table$dataset <- factor(results_table$dataset, levels=datasets)

# function for extracting a legend from a ggplot (for grid.arrange plotting)
g_legend <- function(a.gplot){
  
  a.gplot <- a.gplot + theme(legend.direction = "horizontal") + guides(fill=guide_legend(nrow=1))
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

theme_set(theme_bw(base_size = 18))
# 1. Plot with R2 differences (all features - exp features)
diff_all_exp <- rep(NA,7)
for (i in 1:7){
  diff_all_exp[i] <- mean(results_table[which(results_table$var_type=="all" & results_table$drug==drugs[i]),"r2"])-mean(results_table[which(results_table$var_type=="exp" & results_table$drug==drugs[i]),"r2"])
}
table1 <- cbind(drugs, diff_all_exp)
table1 <- as.data.frame(table1)
colnames(table1)[2] <- "difference in R2"
table1$drugs <- factor(table1$drugs, levels=drugs)
table1$`difference in R2` <- as.numeric(as.character(table1$`difference in R2`))

ggplot(data=table1, aes(y=`difference in R2`,x=drugs)) +geom_bar(stat="identity")
ggsave("fig_s6.pdf")  

# 2. Plot with all results for "all features" models

table2 <- results_table[which(results_table$var_type=="all"),]

for (j in 1:7) {
  drug <- drugs[j]
  
  for (i in 1:3)
  {
    n <- which(table2$metric==metrics[i] & table2$drug==drug)
    table2_2 <- table2[n,c("dataset", "var_number", "r2")]
    table2_2$var_number <- factor(table2_2$var_number, levels=c(10,50,200,500))
    theme_set(theme_grey(base_size = 10)) 
    
    g <- ggplot(table2_2, aes(dataset, r2))
    g2 <- g + geom_bar(aes(fill=var_number),position = "dodge", stat="identity") + labs(title = paste0(metrics[i],", ",drug))  + ylim(0,0.5) 
    if (i==1 & j==1) {legend <- g_legend(g2)}
    g2 <- g2 + theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),legend.position="none")
    name <-paste0("p",j,i)
    assign(name, g2)
  }
}

fig2 <- grid.arrange(arrangeGrob(p11,p12,p13,p21,p22,p23,p31,p32,p33,p41,p42,p43,p51,p52,p53,p61,p62,p63,p71,p72,p73, ncol=3, nrow =7), legend, nrow=2, heights=c(10, 1), left="R2")
ggsave("fig_s4.pdf",plot=fig2, width = 20, height = 30, units = "cm")

# 3. plot with all results where var_number=500
#R2
table5 <- table2[which(table2$var_number==500),]
for (i in 1:3)
{
  n <- which(table5$metric==metrics[i])
  table6 <- table5[n,c("dataset","metric","drug","r2","conc_index")]
  
  g <- ggplot(table6, aes(dataset, r2))
  g2 <- g + geom_bar(aes(fill = drug), position = "dodge", stat="identity") + labs(title = metrics[i], y="R2") + ylim(0,0.5) +theme_bw(base_size = 20) 
  if (i==1) {legend <- g_legend(g2)}
  g2 <- g2 + theme(axis.title.x=element_blank(),axis.title.y =element_blank(),  legend.position="none")  
  
  g3 <- ggplot(table6, aes(dataset, conc_index))
  g4 <- g3 + geom_bar(aes(fill = drug), position = "dodge", stat="identity") + labs(title = metrics[i], y="concordance index") + ylim(0,0.8) +theme_bw(base_size = 20) 
  #if (i==1) {legend <- g_legend(g2)}
  g4 <- g4 + theme(axis.title.x=element_blank(),axis.title.y =element_blank(),  legend.position="none") 
  
  name1 <-paste0("p",i)
  name2 <-paste0("c",i)
  assign(name1, g2)
  assign(name2, g4)
}

fig4a <- grid.arrange(arrangeGrob(p1,p2,p3, ncol=3, nrow =1),legend,ncol=1,nrow=2,heights=c(10,1))
figs3 <- grid.arrange(arrangeGrob(c1,c2,c3, ncol=3, nrow =1),legend,ncol=1,nrow=2,heights=c(10,1))
ggsave("fig_2a.pdf",plot=fig4a)
ggsave("fig_s3.pdf",plot=figs3)


# 4. results for different tissues

table_t2 <- cbind(rep(NA,70), rep(drugs,1,each=10), rep(tissue_list[1:10],7))
for (i in 1:7){
  x <-apply(table2[which(table2$metric=="AUC" & table2$drug==drugs[i]),28:37],2,function(x){mean(as.numeric(as.character(x)), na.rm=T)})
  #print(x)
  table_t2[((i-1)*10+1):((i-1)*10+10),1] <- x
}

table_t2 <- as.data.frame(table_t2)
colnames(table_t2) <- c("R2", "drug", "tissue")
table_t2$tissue <- revalue(table_t2$tissue, replace=c("haematopoietic_and_lymphoid_tissue"="blood"))
table_t2$R2 <- as.numeric(as.character(table_t2$R2))
table_t2$drug <- factor(table_t2$drug, levels = drugs)
table_t2$tissue <- factor(table_t2$tissue, levels=c("lung", "blood", "skin", "central_nervous_system",
                                                    "breast", "large_intestine", "pancreas", "ovary", "stomach", "oesophagus"))

g <- ggplot(table_t2, aes(drug, R2))
g2 <- g + geom_bar(aes(fill = tissue), position = "dodge", stat="identity") +labs(x=NULL)+ theme(legend.text=element_text(size=10), legend.position = "bottom")
ggsave("fig_s8.pdf",plot=g2)

# 5. results (R2, RMSE, concordance index) for different drugs
rmse_r2_ci_drugs_auc <- matrix(NA,ncol=3,nrow=7)
rownames(rmse_r2_ci_drugs_auc) <- drugs
colnames(rmse_r2_ci_drugs_auc) <- c("rmse","r2","conc_index")

for (i in 1:7){
  rmse_r2_ci_drugs_auc[i,] <- apply(table2[which(table2$metric=="AUC" & table2$drug==drugs[i]),c("rmse","r2","conc_index")],2,mean)
}

sd_r2_drugs <- rep(NA, 7)
for (i in 1:7){
  sd_r2_drugs[i] <- sd(table2[which(table2$metric=="AUC" & table2$drug==drugs[i]),"r2"])
}

sd_rmse_drugs <- rep(NA, 7)
for (i in 1:7){
  sd_rmse_drugs[i] <- sd(table2[which(table2$metric=="AUC" & table2$drug==drugs[i]),"rmse"])
}

sd_ci_drugs <- rep(NA, 7)
for (i in 1:7){
  sd_ci_drugs[i] <- sd(table2[which(table2$metric=="AUC" & table2$drug==drugs[i]),"conc_index"])
}

r2_rmse_drugs_auc <- rmse_r2_ci_drugs_auc
r2_rmse_drugs_auc <- cbind(rownames(r2_rmse_drugs_auc),r2_rmse_drugs_auc, sd_r2_drugs, sd_rmse_drugs, sd_ci_drugs)
colnames(r2_rmse_drugs_auc) <- c("drug", "RMSE", "R2","conc_index", "sd", "sd_rmse","sd_ci")
table7 <- as.data.frame(r2_rmse_drugs_auc)
table7$R2 <- as.numeric(as.character(table7$R2))
table7$RMSE <- as.numeric(as.character(table7$RMSE))
table7$conc_index<- as.numeric(as.character(table7$conc_index))
table7$sd <- as.numeric(as.character(table7$sd))
table7$sd_rmse <- as.numeric(as.character(table7$sd_rmse))
table7$sd_ci <- as.numeric(as.character(table7$sd_ci))
table7$drug <- factor(table7$drug, levels=drugs)

theme_set(theme_bw(base_size = 18))
p1 <- ggplot(data=table7, aes(y=R2,x=drugs)) +geom_bar(stat="identity") + geom_errorbar(aes(x=drug,ymin=R2-sd, ymax=R2+sd)) +labs(x=NULL)
p2 <- ggplot(data=table7, aes(y=RMSE,x=drugs)) +geom_bar(stat="identity") + geom_errorbar(aes(x=drug,ymin=RMSE-sd_rmse, ymax=RMSE+sd_rmse)) +labs(x=NULL)
p3 <- ggplot(data=table7, aes(y=conc_index,x=drugs)) +geom_bar(stat="identity") + geom_errorbar(aes(x=drug,ymin=conc_index-sd_ci, ymax=conc_index+sd_ci)) +labs(x=NULL)
ggsave("fig_2b.pdf",plot=p1)
ggsave("fig_s7.pdf",plot=p3)

# 6. results for comparison with DREAM challenge's 2nd method (used in Table 2)
table8 <- table2[which(table2$dataset=="CCLE" & table2$metric=="AUC" & table2$var_number==500),c("drug","r2")]


# 7. R2 (sample size)
table9 <- table2[which(table2$var_number==500),]

for (i in 1:3){
   mt <- metrics[i]
   table9_2 <- table9[which(table9$metric==mt),]
   g <- ggplot(data=table9_2, aes(x=s_size,y=r2, group=drug, colour=drug)) +geom_point() +geom_line() + labs(title = mt) +xlab("sample size") 
   if (i==1) {leg <- g_legend(g)}
   g2 <- g + theme(legend.position="none") 
   name <- paste0("p",i)
   assign(name, g2)
  }

fig_s3 <- grid.arrange(arrangeGrob(p1,p2,p3, nrow=1), leg, nrow=2, heights=c(10, 1))
ggsave("fig_s9.pdf",plot=fig_s3)

#.8 /2019/ results with GR_AOC metric and results with higher number of varibles(500-5000)  
# 8.1 correlation results
cor_res <- cbind(drugs, R2=c(0.046, 0.101, 0.099, -0.0167, -0.0338, -0.0843, 0.0248))
cor_res <- data.frame(cor_res)
cor_res$R2 <- as.numeric(as.character(cor_res$R2))
ggplot(cor_res, aes(x=drugs, y=R2)) + geom_bar(stat="identity") + ylim(-1,1) 

# 8.2 results for different metrics including GR_AOC metric
metrics_res <- cbind(metrics=c("AUC", "GR_AOC", "IC50", "viability_1uM"), R2=c(0.176, 0.143, 0.091, 0.185))
metrics_res <- data.frame(metrics_res)
metrics_res$R2 <- as.numeric(as.character(metrics_res$R2))
metrics_res$metrics <- factor(metrics_res$metrics, levels=c("IC50","AUC", "viability_1uM","GR_AOC"))
ggplot(metrics_res, aes(x=metrics, y=R2)) + geom_bar(stat="identity")
