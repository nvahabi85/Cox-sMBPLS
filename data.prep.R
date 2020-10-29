### Data preparation for the main supervised Cox-sMBPLS function ###

### Dir & Lib ----
library(ggplot2)
library(rlang)
library(rstatix)
library(survminer)
library(ggpubr)
library(biomaRt) # QTL annotation
library(glmnet) # Cindex()
library(ranger)
library(randomForestSRC)
library(mlr)
library(vctrs)
library(tuneRanger)
library(data.table)
library(parallelMap)
library(survival) # survConcordance()
library(survAUC)
library(coxed)
library(stats) # lm()
library(Hmisc)
library(plsgenomics)
library(pec) # predictSurvProb()
library(survAUC) # AUC.cd()
library(magrittr) # %>%
#library(mefa4) # %notin%()
library(dplyr) # %>%
library(purrr) # map()
library(ade4) # mbpls()
library(MASS) # ginv()
library(base) # svd()
library(survcomp) # concordance.index()
library(rlist) #list.cbind
library(c060) # predictProb()
library(peperr) # predictProb()
library(irace)
library(earth)
library(coxed) #Y simulation
library(risksetROC) #AUC
library(MLmetrics) #R2score()


# DATA: Omics
GENE = fread("gene_sample.csv")
GENE = data.frame(GENE, row.names=GENE$ID)[, -1]
SNP = fread("snp_sample.csv")
SNP = data.frame(SNP, row.names=SNP$ID)[, -1]
CPG = fread("cpg_sample.csv")
CPG = data.frame(CPG, row.names=CPG$ID)[, -1]
time_status = fread("time_status.csv")
time_status = data.frame(time_status, row.names=time_status$ID)[, -1]
time = y = time_status[, 1, drop=F]
status = time_status[, 2, drop=F]
Y = cbind(y, status)


## DATA -QT-
eQTL_sample = fread("Whole_Blood.allpairs.maf0.05.p0.1.rs#.gene-name.csv")
eQTM_sample = fread("2015_09_02_cis_eQTMsFDR0.05-CpGLevel.csv")
meQTL_sampel = fread("sic-meQTL.csv")
meQTL_sampel$aaa = paste("X", meQTL_sampel$SNPChr, sep="") #Make X1.1234567 to be able to match with my own data which do not have the rs#!
meQTL_sampel$SNP_id_pos = paste(meQTL_sampel$aaa, meQTL_sampel$SNPChrPos, sep=".")

## QT Genes/SNPs/CpGs
QT_genes = c(eQTL_sample$gene_symbol, eQTM_sample$HGNCName) # List of all "Genes" which are in (eQTL+eQTM)
QT_genes = QT_genes[!duplicated(QT_genes)] # Remove duplicated genes

QTL_snps = c(eQTL_sample$rs_id, meQTL_sampel$SNP_id_pos) # List of all "SNPs" which are in (eQTL+meQTL)
QTL_snps = QTL_snps[!duplicated(QTL_snps)]

QT_cpgs = c(eQTM_sample$SNPName, meQTL_sampel$ProbeName) # List of all "CpGs" which are in eQT (eQTM+meQTL)
QT_cpgs = QT_cpgs[!duplicated(QT_cpgs)]

## Gene-SNP
z1 = eQTL_sample$gene_symbol
z2 = eQTL_sample$Xchr.pos #X1.1234567
z2_2 = eQTL_sample$rs_id_dbSNP147_GRCh37p13 #RS
eQTL_sample_z = data.frame(z1, z2, z2_2, stringsAsFactors = FALSE)
eQTL_sample_z = eQTL_sample_z[order(z1), ]

# Gene-CpG
zz1 = eQTM_sample$HGNCName
zz2 = eQTM_sample$SNPName
eQTM_sample_zz = data.frame(zz1, zz2, stringsAsFactors = FALSE)
eQTM_sample_zz = eQTM_sample_zz[order(zz1), ]

# SNP-CpG
zzz1 = meQTL_sampel$SNP_id_pos
zzz2 = meQTL_sampel$ProbeName
meQTL_sample_zzz = data.frame(zzz1, zzz2, stringsAsFactors = FALSE)
meQTL_sample_zzz = meQTL_sample_zzz[order(meQTL_sample_zzz$zzz1), ]

eQTL_sample_z = data.frame(z1, z2, stringsAsFactors = FALSE) #For the rest of the code (to be matched with my SNPs' names)

##################################################################################################################################################

X1 = GENE[ , order(colnames(GENE))]
X1 = as.matrix(X1)
X2 = SNP[ , order(colnames(SNP))]
X2 = as.matrix(X2)
X3 = CPG[ , order(colnames(CPG))]
X3 = as.matrix(X3)

X.list = list(X1, X2, X3) #combine
names(X.list) = c("GENE", "SNP", "CPG")
X = cbind.data.frame(X.list) #this is just to get the gene-name and snp-name like: GENE.XXX & SNP.XXX to be able to match the beta.T (resulted from MBsPLS) with it
X.col.names = colnames(X)

# 1. Each GENE with its ALL associated SNPs (to split the GENE data)
eQTL_sample_z2 = eQTL_sample_z[which( eQTL_sample_z$z1  %in%  colnames(X1) ), ]
eQTL_sample_z22 = eQTL_sample_z2[which( eQTL_sample_z2$z2  %in%  colnames(X2) ), ]

if (dim(eQTL_sample_z22)[1] != 0) { # In case there is no eQTL-gene in the sampled X1!
  print("Yes1_eQTL>0")
  eQTL_sample_splited_by_gene = eQTL_sample_z22 %>% split(.$z1) %>% map(pull, z2)
  table(duplicated(names(eQTL_sample_splited_by_gene))) #0
  
  # Split X1
  qtl_names1 = names(eQTL_sample_splited_by_gene)
  X1_qtl_names = as.matrix(X1[, which(colnames(X1) %in% qtl_names1), drop=FALSE])
  X1_Non_qtl_names = as.matrix(X1[ , -which(colnames(X1) %in% colnames(X1_qtl_names)), drop=FALSE])
  ncol(X1_qtl_names)
  
  # lm1
  depVarList1 = colnames(X1_qtl_names)
  my_lms1 <- lapply(depVarList1, function(x) 
    lm( X1_qtl_names[, 1] ~ X2[ , which(colnames(X2) %in% eQTL_sample_splited_by_gene[[x]])] ) )
  names(my_lms1) = depVarList1
  
  # u1 (u1 | Gene_non.eqtl)
  X1_qtl_names_u = matrix(0, nrow = nrow(X1_qtl_names), ncol = ncol(X1_qtl_names))
  X1_qtl_names_u = sapply(my_lms1, residuals)
  
  # Update data X1 (GENE) with u
  X1_new = cbind(X1_qtl_names_u, X1_Non_qtl_names) #residuals (u) instead of the gene expression for the eQTL part.
  X1_new = X1_new[order(rownames(X1_new)), order(colnames(X1_new))] #order
  
} else {
  
  print("No1_eQTL=0")
  X1_new = X1
}
X1_new = X1_new[order(rownames(X1_new)), order(colnames(X1_new))] #order


# 2. Each SNP with its ALL associated CPGs (to split the SNP data)
eQTL_sample_zzz2 = meQTL_sample_zzz[which( meQTL_sample_zzz$zzz1 %in% colnames(X2) ), ]
eQTL_sample_zzz22 = eQTL_sample_zzz2[which( eQTL_sample_zzz2$zzz2  %in%  colnames(X3) ), ]

if (dim(eQTL_sample_zzz22)[1] != 0) { # In case there is no eQTL-gene in the sampled X1!
  print("Yes2_eQTL>0")
  eQTL_sample_splited_by_snp = eQTL_sample_zzz22 %>% split(.$zzz1) %>% map(pull, zzz2)
  table(duplicated(names(eQTL_sample_splited_by_snp))) #738 (# of the unique SNPs in meQTL)
  
  # Split X2
  qtl_names2 = names(eQTL_sample_splited_by_snp)
  X2_qtl_names = as.matrix(X2[, which(colnames(X2) %in% qtl_names2), drop=FALSE])
  X2_Non_qtl_names = as.matrix(X2[ , -which(colnames(X2) %in% colnames(X2_qtl_names)), drop=FALSE])
  ncol(X2_qtl_names)
  
  # lm2
  depVarList2 = colnames(X2_qtl_names)
  my_lms2 <- lapply(depVarList2, function(x) 
    lm( X2_qtl_names[, x] ~ X3[ , which(colnames(X3) %in% eQTL_sample_splited_by_snp[[x]])] ) )
  names(my_lms2) = depVarList2
  
  # u2 (u2 | SNP_non.meqtl)
  X2_qtl_names_u = matrix(0, nrow = nrow(X2_qtl_names), ncol = ncol(X2_qtl_names))
  X2_qtl_names_u = sapply(my_lms2, residuals)
  X2_qtl_names_u[1:5,]
  table(duplicated(colnames(X2_qtl_names_u))) #0
  
  # Update X2
  X2_new = cbind(X2_qtl_names_u, X2_Non_qtl_names) #residuals (u) instead of the gene expression for the eQTL part.
  X2_new = X2_new[order(rownames(X2_new)), ] #order
  table(duplicated(colnames(X2_new))) #0
  
} else {
  print("No2_eQTL=0!")
  X2_new = X2
  X2_new = X2_new[order(rownames(X2_new)), ] #order
}


# 3. Each CPG with its ALL associated GENEs (to split the CPG data)
eQTL_sample_zz2 = eQTM_sample_zz[which( eQTM_sample_zz$zz2  %in%  colnames(X3) ), ]
eQTL_sample_zz22 = eQTL_sample_zz2[which( eQTL_sample_zz2$zz1 %in% colnames(X1) ), ]

if (dim(eQTL_sample_zz22)[1] != 0) { # In case there is no eQTL-gene in the sampled X1!
  print("Yes3_eQTL>0")
  eQTL_sample_splited_by_cpg = eQTL_sample_zz22 %>% split(.$zz2) %>% map(pull, zz1)
  table(duplicated(names(eQTL_sample_splited_by_cpg))) #11,797 (# of the unique CpGs in eQTM)
  
  # Split X3
  qtl_names3 = names(eQTL_sample_splited_by_cpg)
  X3_qtl_names = as.matrix(X3[, which(colnames(X3) %in% qtl_names3), drop=FALSE])
  X3_Non_qtl_names = as.matrix(X3[ , -which(colnames(X3) %in% colnames(X3_qtl_names)), drop=FALSE])
  ncol(X3_qtl_names) #11,265***
  
  # lm3
  depVarList3 = colnames(X3_qtl_names)
  my_lms3 <- lapply(depVarList3, function(x) 
    lm( X3_qtl_names[, x] ~ X1[ , which(colnames(X1) %in% eQTL_sample_splited_by_cpg[[x]])] ) )
  names(my_lms3) = depVarList3
  
  # u3 (u3 | CpG_non.eqtm)
  X3_qtl_names_u = matrix(0, nrow = nrow(X3_qtl_names), ncol = ncol(X3_qtl_names))
  X3_qtl_names_u = sapply(my_lms3, residuals)
  X3_qtl_names_u[1:5, ]
  table(duplicated(colnames(X3_qtl_names_u))) #0
  
  # Update X3
  X3_new = cbind(X3_qtl_names_u, X3_Non_qtl_names) #residuals (u) instead of the gene expression for the eQTL part.
  X3_new = X3_new[order(rownames(X3_new)), order(colnames(X3_new))] #order
  table(duplicated(colnames(X3_new))) #0
  
} else {
  print("No3_eQTL=0!")
  X3_new = X3
}
X3_new = X3_new[order(rownames(X3_new)), ] #order


## [Scale]
# X1
mean_X1 = apply(X1_new, 2, mean)
sigma_X1 = apply(X1_new, 2, sd)

for (i in 1:ncol(X1_new)) {
  if(any( sigma_X1[i] < .Machine$double.eps )){
    sigma_X1[i] = 0.0001}
}
sX1_new <- scale( X1_new, center = mean_X1, scale = sigma_X1)
sX1_new = sX1_new[order(rownames(sX1_new)), order(colnames(sX1_new))]

# X2
mean_X2 = apply(X2_new, 2, mean)
sigma_X2 = apply(X2_new, 2, sd)

for (i in 1:ncol(X2_new)) {
  if(any( sigma_X2[i] < .Machine$double.eps )){
    sigma_X2[i] = 0.0001}
}
sX2_new <- scale( X2_new, center = mean_X2, scale = sigma_X2 )
sX2_new = sX2_new[ , order(colnames(sX2_new))]

# X3
mean_X3 = apply(X3_new, 2, mean)
sigma_X3 = apply(X3_new, 2, sd)

for (i in 1:ncol(X3_new)) {
  if(any( sigma_X3[i] < .Machine$double.eps )){
    sigma_X3[i] = 0.0001}
}
sX3_new = scale( X3_new, center = mean_X3, scale = sigma_X3 )
sX3_new = sX3_new[ , order(colnames(sX3_new))]

# Final scaled data 
sX = cbind(sX1_new, sX2_new, sX3_new)


## [Split] 30/70
# X  
Ntrain = floor(0.70 * nrow(X1_new))
train.index = sample(nrow(X1), Ntrain)

X1.train = X1[train.index, ]
X1.train = X1.train[order(rownames(X1.train)), order(colnames(X1.train))]

X1.test = X1[-train.index, ]
X1.test = X1.test[order(rownames(X1.test)), order(colnames(X1.test))]

X2.train = X2[rownames(X1.train), ]
X2.train = X2.train[order(rownames(X2.train)), order(colnames(X2.train))]

X2.test = X2[rownames(X1.test), ]
X2.test = X2.test[order(rownames(X2.test)), order(colnames(X2.test))]

X3.train = X3[train.index, ]
X3.train = X3.train[order(rownames(X3.train)), order(colnames(X3.train))]

X3.test = X3[-train.index, ]
X3.test = X3.test[order(rownames(X3.test)), order(colnames(X3.test))]

X.train.dataframe = cbind(X1.train, X2.train, X3.train) # Combined (data.frame) 
X.test.dataframe = cbind(X1.test, X2.test, X3.test)
X.col.names = colnames(X.train.dataframe)

X.train.list = list(X1.train, X2.train, X1.train) # combined (list)
names(X.train.list) = c("GENE", "SNP", "CPG")

# Y
Y.train = time_status[rownames(X1.train), ]
Y.test = time_status[rownames(X1.test), ]


# Final splited data::**MAIN DATA FOR ALL OTHER MODELS (except MBsPLS)***
dat.train = cbind(Y.train, X1.train, X2.train, X3.train) # data.frame
dat.test = cbind(Y.test, X1.test, X2.test, X3.test)

surv.train = with(dat.train, Surv(y,status))
surv.test = with(dat.test, Surv(y,status))

###################################################################################################
#### Supervised Cox-sMBPLS mODEL ####
# Fit the supervised Cox-sMBPLS model
# For the proposed model (supervised Cox-sMBPLS),we run a Train/Test/Valid experiments, to also validate the model.
# therefore, we use 85% of the train-set as a valid-set to validate the model. 

## Data Preparation & Tuning
# Train-set
X1_new_train = X1_new[which(rownames(X1_new) %in% rownames(X1.train)), ]
X2_new_train = X2_new[which(rownames(X2_new) %in% rownames(X2.train)), ]
X3_new_train = X3_new[which(rownames(X3_new) %in% rownames(X3.train)), ]

Y.rew.train = Y.train[order(Y.train$y), ] 
sc = survfit(Surv(y, 1-status) ~ 1, data = Y.rew.train)$surv 
sc = c(1, sc[-nrow(Y.rew.train)])
sc[sc==0] = 0.0001
y = Y.rew.train$y * Y.rew.train$status/sc 
Y.rew.train = cbind(y, Y.rew.train[, -1, drop=FALSE])
Y.rew.train = Y.rew.train[order(rownames(Y.rew.train)), ] 

# Test-set
X1_new_test = X1_new[which(rownames(X1_new) %in% rownames(X1.test)), ] 
X2_new_test = X2_new[which(rownames(X2_new) %in% rownames(X2.test)), ]
X3_new_test = X3_new[which(rownames(X3_new) %in% rownames(X3.test)), ]

Y.rew.test = Y.test[order(Y.test$y), ] 
sc = survfit(Surv(y, 1-status) ~ 1, data = Y.rew.test)$surv 
sc = c(1, sc[-nrow(Y.rew.test)]) 
sc[sc==0] = 0.0001
y = Y.rew.test$y * Y.rew.test$status/sc 
Y.rew.test = cbind(y, Y.rew.test[, -1, drop=FALSE])
Y.rew.test = Y.rew.test[order(rownames(Y.rew.test)), ] 

# Split to train/validate sets (85/15)
num_train = floor(0.85 * nrow(X1_new_train))
index_train = sample(rownames(X1_new_train), num_train)

XX1.train = X1_new_train[index_train, ]
XX1.train = XX1.train[order(rownames(XX1.train)), ]
XX1.valid = X1_new_train[ !( rownames(X1_new_train) %in% rownames(XX1.train) ), , drop=F]
XX1.valid = XX1.valid[order(rownames(XX1.valid)), , drop=F]
XX1.test = X1_new_test[which(rownames(X1_new_test) %in% rownames(X1.test)), ]
XX1.test = XX1.test[order(rownames(XX1.test)), ]


XX2.train = X2_new_train[index_train, ]
XX2.train = XX2.train[order(rownames(XX2.train)), ]
XX2.valid = X2_new_train[ !(rownames(X2_new_train) %in% rownames(XX2.train)), , drop=F]
XX2.valid = XX2.valid[order(rownames(XX2.valid)), , drop=F]
XX2.test = X2_new_test[which(rownames(X2_new_test) %in% rownames(X2.test)), ]
XX2.test = XX2.test[order(rownames(XX2.test)), ]

XX3.train = X3_new_train[index_train, ]
XX3.train = XX3.train[order(rownames(XX3.train)), ]
XX3.valid = X3_new_train[ !(rownames(X3_new_train) %in% rownames(XX3.train)), , drop=F]
XX3.valid = XX3.valid[order(rownames(XX3.valid)), , drop=F]
XX3.test = X3_new_test[which(rownames(X3_new_test) %in% rownames(X3.test)), ]
XX3.test = XX3.test[order(rownames(XX3.test)), ]

YY.train = y.train.finalCox[index_train, ]
YY.train = YY.train[order(rownames(YY.train)), ]
YY.valid = y.train.finalCox[ !(rownames(y.train.finalCox) %in% rownames(YY.train)), ]
YY.valid = YY.valid[order(rownames(YY.valid)), ]
YY.test = y.test.finalCox[which(rownames(y.test.finalCox) %in% rownames(Y.test)), ]
YY.test = YY.test[order(rownames(YY.test)), ]

y.train.finalCox = cbind.data.frame(Y.rew.train, Y.train[ , 1, drop=FALSE])
y.test.finalCox = cbind.data.frame(Y.rew.test, Y.test[ , 1, drop=FALSE])

# Final Data for the supervised Cox-sMBPLS model 
dat.train.list = list(rbind(YY.train, YY.valid), rbind(XX1.train, XX1.valid), rbind(XX2.train, XX2.valid), rbind(XX3.train, XX3.valid)) #list
names(dat.train.list) = c("Y", "GENE", "SNP", "CPG")

dat.test.list = list(YY.test, XX1.test, XX2.test, XX3.test) #list
names(dat.test.list) = c("Y.test", "GENE.test", "SNP.test", "CPG.test")
surv.test = with(dat.test.list[[1]], Surv(dat.test.list[[1]][,3], dat.test.list[[1]][,2]))

## Parameter tuning
ptm.mbspls = proc.time()
source("fun.cox.smbpls.R")
ncomp_tune_set = c(2,3,4)
lambda_tune_set = c(0.0001, 0.1, 0.5, 0.95)
tuning_result = expand.grid("ncomp_i" = ncomp_tune_set, # Components no
                            "lambda_i" = lambda_tune_set, # Sparsity
                            "cindex.train_i" = NA, 
                            "cindex.test_i" = NA,
                            stringsAsFactors = FALSE)
for (lambda_i in lambda_tune_set) {
  for (ncomp_i in ncomp_tune_set) {
    MBsPLS.result.tune = MB_spls3_Cox(Xtrain = dat.train.list[2:4], 
                                      Ytrain = dat.train.list[[1]], 
                                      lambda.l1 = lambda_i, 
                                      ncomp = ncomp_i,
                                      weight.mat=NULL, 
                                      Xtest = dat.test.list[2:4], # Turn "Xtest=TRUE" in MBsPLS function
                                      Ytest = dat.test.list[[1]], #Turn "Ytest=TRUE" in MBsPLS function
                                      center.X=TRUE, 
                                      center.Y=TRUE,
                                      scale.X=TRUE, 
                                      scale.Y=TRUE,
                                      weighted.center=FALSE)
    # write the params used
    tune_lambda_result_lambda_ncomp_i = MBsPLS.result.tune$lambda.l1
    tune_ncomp_result_lambda_ncomp_i = MBsPLS.result.tune$ncomp
    train_cindex_result_lambda_ncomp_i = MBsPLS.result.tune$Cindex.train1
    test_cindex_result_lambda_ncomp_i = MBsPLS.result.tune$Cindex.test2
    # Save the results
    tuning_result[tuning_result$ncomp==ncomp_i & tuning_result$lambda_i==lambda_i, 3] = MBsPLS.result.tune$Cindex.train1
    tuning_result[tuning_result$ncomp==ncomp_i & tuning_result$lambda_i==lambda_i, 4] = MBsPLS.result.tune$Cindex.test2
    
  } 
}
# Extract the tuned parameters
tuning_result = tuning_result[order(-tuning_result$cindex.test_i, -tuning_result$cindex.train_i, tuning_result$ncomp_i), ] # Sort
lambda.tuned = tuning_result$lambda_i[1]
ncomp.tuned = tuning_result$ncomp_i[1]
cindex.MBsPLS.traine.tuneFunc = tuning_result$cindex.train_i[1]
cindex.MBsPLS.valid.tuneFunc = tuning_result$cindex.test_i[1]

## Cox-sMBPLS model fit
source("fun.cox.smbpls.R")
dat.train.valid.list = list(rbind(YY.train, YY.valid), rbind(XX1.train, XX1.valid), rbind(XX2.train, XX2.valid), rbind(XX3.train, XX3.valid))
names(dat.train.valid.list) = c("Y.trainTotal", "GENE.TrainTotal", "SNP.TrainTotal", "CPG.TrainTotal")
MBsPLS.result = MB_spls3_Cox(Xtrain = dat.train.valid.list[2:4], 
                             Ytrain = dat.train.valid.list[[1]], 
                             lambda.l1 = lambda.tuned, 
                             ncomp = ncomp.tuned,
                             weight.mat=NULL, 
                             Xtest = dat.test.list[2:4], # Turn "Xtest=TRUE" in MBsPLS function
                             Ytest = dat.test.list[[1]], #Turn "Ytest=TRUE" in MBsPLS function
                             center.X=TRUE, 
                             center.Y=TRUE,
                             scale.X=TRUE, 
                             scale.Y=TRUE,
                             weighted.center=FALSE)
cindex.MBsPLS.train.tuned = MBsPLS.result$Cindex.train1
cindex.MBsPLS.test.tuned = MBsPLS.result$Cindex.test2_2


## AUC
fit.mbspls.cox = MBsPLS.result$Cox.mbspls
data.test.final.cox = MBsPLS.result$data.test.final.cox
surv.train.valid = with(dat.train.valid.list[[1]], Surv(dat.train.valid.list[[1]][,3], dat.train.valid.list[[1]][,2]))

# Chambless & Diao AUC
max_time = max(Y$y)
times = seq(10, max_time, 20) 
auc.mbspls = AUC.cd(surv.train.valid, 
                    surv.test, 
                    predict(fit.mbspls.cox), 
                    predict(fit.mbspls.cox, data.test.final.cox), 
                    times)
MBsPLS.iauc1 = auc.mbspls$iauc 

# Heagerty D/I AUC
w.ROC = risksetROC(Stime = Y.test$y,  
                   status = Y.test$status, 
                   marker = predict(fit.mbspls.cox, data.test.final.cox), 
                   predict.time = median(Y.test$y), 
                   method = "Cox", 
                   main = paste("OOB Survival ROC Curve at median(t)"), 
                   lwd = 3, 
                   col = "red" )
MBsPLS.auc2 = w.ROC$AUC

# Uno's AUC
AUC_Uno.mbspla = AUC.uno(surv.train.valid, 
                         surv.test,
                         predict(fit.mbspls.cox, data.test.final.cox),
                         times)
MBsPLS.uno = AUC_Uno.mbspla$auc[1]
MBsPLS.iuno = AUC_Uno.mbspla$iauc

# Run-time
run.time.mbspls = proc.time() - ptm.mbspls
run.time.mbspls = run.time.mbspls[1]
# End
