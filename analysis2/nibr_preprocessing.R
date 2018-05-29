# 27.05.2018, NIBR PDXE data preprocessing

install.packages("flux")
library(flux)

path <- "/data/kurilov/genestack/phd/work_2018/DrugRespPrediction/"
setwd(file.path(path, "analysis2"))

# reading csv files extracted from NIBR PDXE paper supplementary data 
# https://media.nature.com/original/nature-assets/nm/journal/v21/n11/extref/nm.3954-S2.xlsx
exp_data <- read.csv("data/exp.csv", check.names=FALSE, header=TRUE )
cn_data <- read.csv("data/cn.csv", check.names=FALSE)
mut_ind_data <- read.csv("data/mut_ind.csv", check.names=FALSE)
drug_resp_raw <- read.csv("data/drug_resp_raw.csv", check.names=FALSE)

#tissue type
tissue_nibr <- unique(drug_resp_raw[,1:2])
colnames(tissue_nibr) <- c("sample","tissue")
# delete two duplicated sample records (with no tissue data)
tissue_nibr <- tissue_nibr[-c(94,219),]
save(tissue_nibr, file="data/tissue_nibr.RData")

#expression
rownames(exp_data) <- exp_data[,1]

exp_data <- exp_data[,-1]
exp_nibr <- exp_data
exp_nibr <- t(exp_nibr)
colnames(exp_nibr) <- paste0("exp_", colnames(exp_nibr))

#copy number
cn_data <- cn_data[-c(1,2),]
# there are multiple entries for some genes in cn data, let's calculate average for these cases (MIR1244-1  MIR1244-2  MIR1244-3  LOC100133331 MIR548AR etc.)
a <- which(cn_data[,1]=="MIR1244-3")
cn_data[a[1],-1]<- (cn_data[a[1],-1]+cn_data[a[2],-1]+cn_data[a[3],-1]+cn_data[a[4],-1])/4
cn_data <- cn_data[-a[-1],]

b <- which(cn_data[,1]=="MIR548AR")
cn_data[b[1],-1]<- (cn_data[b[1],-1]+cn_data[b[2],-1]+cn_data[b[3],-1])/3
cn_data <- cn_data[-b[-1],]

test <- which(duplicated(cn_data[,1]))
for (i in 1:length(test))
{
  cn_data[(test[i]-1),-1]<- (cn_data[(test[i]-1),-1]+cn_data[test[i],-1])/2
}
cn_data <- cn_data[-test,]

rownames(cn_data) <- cn_data[,1]
cn_nibr <- cn_data[,-1]
cn_nibr <- round(cn_nibr)
cn_nibr <- t(cn_nibr)
colnames(cn_nibr) <- paste0("cn_", colnames(cn_nibr))

#mutation, indels
sample <- unique(mut_ind_data$Sample)
gene_name <- unique(mut_ind_data$Gene)

mut <- matrix(0, length(sample), length(gene_name))
rownames(mut) <- sample
colnames(mut) <- gene_name

ind <- mut

for (i in 1:(dim(mut_ind_data)[1])) {
  row <- which(rownames(mut)==mut_ind_data$Sample[i])
  col <- which(colnames(mut)==mut_ind_data$Gene[i])
  
  if (mut_ind_data$Category[i] %in% c("MutKnownFunctional", "MutLikelyFunctional", "MutNovel"))
  {
    mut[row,col]<- mut[row,col]+1
  }
  
  if (mut_ind_data$Category[i] %in% c("Amp5", "Amp8", "Del0.8"))
  {
    ind[row,col]<- ind[row,col]+1
  }
}


colnames(mut) <- paste("mut", colnames(mut), sep="_")
colnames(ind) <- paste("ind", colnames(ind), sep="_")
mut_ind_nibr <- cbind(mut, ind)
sd_vect <- apply(mut_ind_nibr, 2, sd)
mut_ind_nibr <- mut_ind_nibr[,-(which(sd_vect==0))]

# combining all molecular data

lines <- intersect(tissue_nibr$sample, rownames(exp_nibr))
dummies_tissue = model.matrix(~tissue_nibr$tissue)
rownames(dummies_tissue) <- tissue_nibr$sample

order_exp <- match(lines, rownames(exp_nibr))
order_cn <- match(lines, rownames(cn_nibr))
order_mut <- match(lines, rownames(mut_ind_nibr))
order_tissue <- match(lines, rownames(dummies_tissue))

mol_nibr <- cbind(exp_nibr[order_exp,],cn_nibr[order_cn,],mut_ind_nibr[order_mut,],dummies_tissue[order_tissue,])
null_var <- which(apply(mol_nibr,2,sd)==0)
mol_nibr <- mol_nibr[,-null_var]
colnames(mol_nibr)[72469:72473] <- c("tissue_BRCA", "tissue_CM", "tissue_CRC", "tissue_NSCLC", "tissue_PDAC")
save(mol_nibr, file="data/mol_nibr.RData")

## slope of untreated tumor growth curve ##

drug_resp_raw0 <- drug_resp_raw[which(drug_resp_raw$Treatment=="untreated"),]
sampl <- unique(drug_resp_raw0$Model)
slope <- rep(NA, length(sampl))
for (i in 1:length(sampl)){
  s <- sampl[i]
  df <- drug_resp_raw0[which(drug_resp_raw0$Model==s),c(4,6)]
  reg <- lm(df[,1]~df[,2])
  slope[i] <- coef(reg)[2]
}
slope_table <- cbind(as.character(sampl), slope)
rownames(slope_table) <- slope_table[,1]
slope_table <- as.data.frame(slope_table)
slope_table[,2] <- as.numeric(as.character(slope_table[,2]))
save(slope_table, file="data/slope_table_nibr.RData")

## drug response

drug_resp_raw <- as.data.frame(drug_resp_raw)

# drugs of interest
drug_int <- c("untreated", "erlotinib", "gemcitabine-50mpk","paclitaxel")

# function for interpolating drug response data (volume) for the days on which there were no measurement.
fill_gap <- function(vect, position)
{
  left=position-1
  right=position+min(which(!(is.na(vect[(position+1):length(vect)]))))
  if (is.na(right)) {
    right <- left
    left <- left-1}
  
  y <- vect[c(left, right)]
  x <- c(left-1, right-1)
  new <- data.frame(x=position-1)
  gap <- predict(lm(y ~ x), newdata=new)
  return(gap)
}

# volume (volume at each day, normalised by the vol. at 0 day (day0 volume=1))
for (k in 1:length(drug_int)) {
  drug_resp_raw1 <- drug_resp_raw[which(drug_resp_raw$Treatment ==drug_int[k]),]
  
  
  model <- unique(drug_resp_raw1$Model)
  days <- sort(unique(drug_resp_raw1$`Days Post T0`))
  resp_table <- matrix(NA, nrow=length(model), ncol=length(days))
  colnames(resp_table) <- days
  rownames(resp_table) <- model
  for (i in 1:nrow(drug_resp_raw1)) {
    col <- as.character(drug_resp_raw1$`Days Post T0`[i])
    row <- as.character(drug_resp_raw1$Model[i])
    resp_table[row,col] <- drug_resp_raw1$`Volume (mm3)`[i]
  }
  
  colnames(resp_table) <- paste0("day", days)
  
  if(drug_int[k]=="paclitaxel") {
    num <- which(rownames(resp_table)=="X-2042")
    resp_table <- resp_table[-num,] }
  
  resp_table2 <- resp_table
  for (i in 1:nrow(resp_table2)) {
    for (j in 2:22){
      if (is.na(resp_table2[i,j])) {
        resp_table2[i,j] <- fill_gap(resp_table2[i,],j)
      }
    }}
  
  resp_table3 <- t(apply(resp_table2,1, function(z){
    z=z/z[1]
    z
  }))
  colnames(resp_table3) <- colnames(resp_table2)
  
  assign(paste0("vol_",drug_int[k]), resp_table3)
}

vol_nibr <- list("erlotinib"=vol_erlotinib,"gemcitabine-50mpk"=`vol_gemcitabine-50mpk`, "paclitaxel"=vol_paclitaxel)

# integral response (control AUC - drug treatment AUC)

# function for subtracting a drug treatment AUC from the control AUC
control_adj <- function(auc_table){
  int <- intersect(rownames(auc_table), rownames(auc_untreated))
  auc_table_new <- auc_untreated[int,]-auc_table[int,]
  return(auc_table_new)
}

for (k in 1:length(drug_int)) {
  print(k)
  drug_resp_raw1 <- drug_resp_raw[which(drug_resp_raw$Treatment ==drug_int[k]),]
  
  
  model <- unique(drug_resp_raw1$Model)
  days <- sort(unique(drug_resp_raw1$`Days Post T0`))
  resp_table <- matrix(NA, nrow=length(model), ncol=length(days))
  colnames(resp_table) <- days
  rownames(resp_table) <- model
  for (i in 1:nrow(drug_resp_raw1)) {
    col <- as.character(drug_resp_raw1$`Days Post T0`[i])
    row <- as.character(drug_resp_raw1$Model[i])
    resp_table[row,col] <- drug_resp_raw1$`Volume (mm3)`[i]
  }
  
  colnames(resp_table) <- paste0("day", days)
  
  if(drug_int[k]=="paclitaxel") {
    num <- which(rownames(resp_table)=="X-2042")
    resp_table <- resp_table[-num,] }
  
  resp_table2 <- resp_table
  for (i in 1:nrow(resp_table2)) {
    for (j in 2:22){
      if (is.na(resp_table2[i,j])) {
        resp_table2[i,j] <- fill_gap(resp_table2[i,],j)
      }
    }}
  
  auc_table <- t(apply(resp_table2,1, function(z){
    z=z/z[1]
    for (i in 1:21) {
      auc_n <- flux:::auc(x=0:i, y=z[1:(i+1)])
      assign(paste0("auc",i), auc_n)
    }
    c(auc1,auc2,auc3,auc4,auc5,auc6,auc7, auc8, auc9,auc10,
      auc11,auc12,auc13,auc14,auc15,auc16,auc17,auc18,auc19,auc20,auc21)
  }))
  colnames(auc_table) <- colnames(resp_table2)[2:22]
  
  if (drug_int[k]!="untreated") {auc_table <- control_adj(auc_table)}
  
  assign(paste0("auc_",drug_int[k]), auc_table)
}

auc_nibr <- list("erlotinib"=auc_erlotinib,"gemcitabine-50mpk"=`auc_gemcitabine-50mpk`, "paclitaxel"=auc_paclitaxel)

# slope of tumor growth curve

for (k in 1:length(drug_int)) {
  print(k)
  drug_resp_raw1 <- drug_resp_raw[which(drug_resp_raw$Treatment ==drug_int[k]),]
  
  sampl <- unique(drug_resp_raw1$Model)
  slope <- rep(NA, length(sampl))
  for (i in 1:length(sampl)){
    s <- sampl[i]
    df <- drug_resp_raw1[which(drug_resp_raw1$Model==s),c(4,6)]
    reg <- lm(df[,1]~df[,2])
    slope[i] <- coef(reg)[2]
  }
  slope_table <- cbind(as.character(sampl), slope)
  rownames(slope_table) <- slope_table[,1]
  slope_table <- as.data.frame(slope_table)
  slope_table[,2] <- as.numeric(as.character(slope_table[,2]))
  
  assign(paste0("slope_",drug_int[k]), slope_table)
}

#renaming
slope_gemcitabine <- `slope_gemcitabine-50mpk`

# differential slope (control slope - drug treatment slope)
diff_slope <- function(slope_treated){
  int <- intersect(rownames(slope_untreated), rownames(slope_treated))
  new_table <- cbind(int, slope_untreated[int,2]-slope_treated[int,2])
  
  rownames(new_table) <- new_table[,1]
  new_table <- as.data.frame(new_table)
  new_table[,2] <- as.numeric(as.character(new_table[,2]))
  return(new_table)
}
diff_slope_erlotinib <- diff_slope(slope_erlotinib)
diff_slope_gemcitabine <- diff_slope(slope_gemcitabine)
diff_slope_paclitaxel <- diff_slope(slope_paclitaxel)

## saving all drug response data
samp_erl <- intersect(intersect(rownames(vol_nibr$erlotinib), rownames(auc_nibr$erlotinib)), rownames(slope_erlotinib))
samp_gem <- intersect(intersect(rownames(vol_nibr$`gemcitabine-50mpk`), rownames(auc_nibr$`gemcitabine-50mpk`)), rownames(slope_gemcitabine))
samp_pac <- intersect(intersect(rownames(vol_nibr$paclitaxel), rownames(auc_nibr$paclitaxel)), rownames(slope_paclitaxel))
resp_erlotinib <- cbind(vol_nibr$erlotinib[samp_erl,"day21"], auc_nibr$erlotinib[samp_erl,"day21"], slope_erlotinib[samp_erl,2], diff_slope_erlotinib[samp_erl,2])
resp_gemcitabine <- cbind(vol_nibr$`gemcitabine-50mpk`[samp_gem,"day21"], auc_nibr$`gemcitabine-50mpk`[samp_gem,"day21"], slope_gemcitabine[samp_gem,2], diff_slope_gemcitabine[samp_gem,2])
resp_paclitaxel <- cbind(vol_nibr$paclitaxel[samp_pac,"day21"], auc_nibr$paclitaxel[samp_pac,"day21"], slope_paclitaxel[samp_pac,2], diff_slope_paclitaxel[samp_pac,2])
resp_cols <- c("vol", "auc", "slope", "diff_slope")
colnames(resp_erlotinib) <- resp_cols
colnames(resp_gemcitabine) <- resp_cols
colnames(resp_paclitaxel) <- resp_cols
resp_nibr <- list(erlotinib=resp_erlotinib, gemcitabine=resp_gemcitabine, paclitaxel=resp_paclitaxel)
save(resp_nibr, file="data/resp_nibr.RData")
