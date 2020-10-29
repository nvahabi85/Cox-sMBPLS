### Data preparation for the main supervised Cox-sMBPLS function ###

#### Dir & Lib ####
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

# DATA: QT
eQTL_sample = fread("Whole_Blood.allpairs.maf0.05.p0.1.rs#.gene-name.csv")
eQTM_sample = fread("2015_09_02_cis_eQTMsFDR0.05-CpGLevel.csv")
meQTL_sampel = fread("sic-meQTL.csv")
meQTL_sampel$chr.pos = paste("X", meQTL_sampel$SNPChr, sep="") # Making a new column: "XChr#.position"
meQTL_sampel$SNP_id_pos = paste(meQTL_sampel$chr.pos, meQTL_sampel$SNPChrPos, sep=".")

## QT Genes/SNPs/CpGs
QT_genes = c(eQTL_sample$gene_symbol, eQTM_sample$HGNCName) # List of all "Genes" which are in (eQTL+eQTM)
QT_genes = QT_genes[!duplicated(QT_genes)] # Remove duplicated genes

QTL_snps = c(eQTL_sample$rs_id, meQTL_sampel$SNP_id_pos) # List of all "SNPs" which are in (eQTL+meQTL)
QTL_snps = QTL_snps[!duplicated(QTL_snps)]

QT_cpgs = c(eQTM_sample$SNPName, meQTL_sampel$ProbeName) # List of all "CpGs" which are in eQT (eQTM+meQTL)
QT_cpgs = QT_cpgs[!duplicated(QT_cpgs)]

## Gene-SNP
z1 = eQTL_sample$gene_symbol
z2 = eQTL_sample$Xchr.pos #chr.pos
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
eQTL_sample_z = data.frame(z1, z2, stringsAsFactors = FALSE)


#### Updating the Omics data using the QT-residuals ####
# Define:
X1 = as.matrix(GENE)
X2 = as.matrix(SNP)
X3 = as.matrix(CPG)

# 1. Each GENE with its ALL associated SNPs (to split the GENE data)
eQTL_sample_z2 = eQTL_sample_z[which( eQTL_sample_z$z1  %in%  colnames(X1) ), ]
eQTL_sample_z22 = eQTL_sample_z2[which( eQTL_sample_z2$z2  %in%  colnames(X2) ), ]
eQTL_sample_splited_by_gene = eQTL_sample_z22 %>% split(.$z1) %>% map(pull, z2)

# 2. Each SNP with its ALL associated CPGs (to split the SNP data)
eQTL_sample_zzz2 = meQTL_sample_zzz[which( meQTL_sample_zzz$zzz1 %in% colnames(X2) ), ]
eQTL_sample_zzz22 = eQTL_sample_zzz2[which( eQTL_sample_zzz2$zzz2  %in%  colnames(X3) ), ]
eQTL_sample_splited_by_snp = eQTL_sample_zzz22 %>% split(.$zzz1) %>% map(pull, zzz2)

# 3. Each CPG with its ALL associated GENEs (to split the CPG data)
eQTL_sample_zz2 = eQTM_sample_zz[which( eQTM_sample_zz$zz2  %in%  colnames(X3) ), ]
eQTL_sample_zz22 = eQTL_sample_zz2[which( eQTL_sample_zz2$zz1 %in% colnames(X1) ), ]
eQTL_sample_splited_by_cpg = eQTL_sample_zz22 %>% split(.$zz2) %>% map(pull, zz1)

# Split X1
qtl_names1 = names(eQTL_sample_splited_by_gene)
X1_qtl_names = as.matrix(X1[, which(colnames(X1) %in% qtl_names1), drop=FALSE])
X1_Non_qtl_names = as.matrix(X1[ , -which(colnames(X1) %in% colnames(X1_qtl_names)), drop=FALSE])

# Split X2
qtl_names2 = names(eQTL_sample_splited_by_snp)
X2_qtl_names = as.matrix(X2[, which(colnames(X2) %in% qtl_names2), drop=FALSE])
X2_Non_qtl_names = as.matrix(X2[ , -which(colnames(X2) %in% colnames(X2_qtl_names)), drop=FALSE])

# Split X3
qtl_names3 = names(eQTL_sample_splited_by_cpg)
X3_qtl_names = as.matrix(X3[, which(colnames(X3) %in% qtl_names3), drop=FALSE])
X3_Non_qtl_names = as.matrix(X3[ , -which(colnames(X3) %in% colnames(X3_qtl_names)), drop=FALSE])

# lm1 (Extract residuals u1)
depVarList1 = colnames(X1_qtl_names)
my_lms1 = lapply(depVarList1, function(x) 
  lm( X1_qtl_names[, 1] ~ X2[ , which(colnames(X2) %in% eQTL_sample_splited_by_gene[[x]])] ) )

# lm2 (Extract residuals u2)
depVarList2 = colnames(X2_qtl_names)
my_lms2 <- lapply(depVarList2, function(x) 
  lm( X2_qtl_names[, x] ~ X3[ , which(colnames(X3) %in% eQTL_sample_splited_by_snp[[x]])] ) )
names(my_lms2) = depVarList2

# lm3 (Extract residuals u3)
depVarList3 = colnames(X3_qtl_names)
my_lms3 <- lapply(depVarList3, function(x) 
  lm( X3_qtl_names[, x] ~ X1[ , which(colnames(X1) %in% eQTL_sample_splited_by_cpg[[x]])] ) )
names(my_lms3) = depVarList3


# Update QT-X1 with u1
X1_qtl_names_u = matrix(0, nrow = nrow(X1_qtl_names), ncol = ncol(X1_qtl_names))
X1_qtl_names_u = sapply(my_lms1, residuals)

# Update QT-X2 with u2
X2_qtl_names_u = matrix(0, nrow = nrow(X2_qtl_names), ncol = ncol(X2_qtl_names))
X2_qtl_names_u = sapply(my_lms2, residuals)

# Update QT-X3 with u3
X3_qtl_names_u = matrix(0, nrow = nrow(X3_qtl_names), ncol = ncol(X3_qtl_names))
X3_qtl_names_u = sapply(my_lms3, residuals)

# Split X1 [u1 | Gene_non.eqtl]
X1_new = cbind(X1_qtl_names_u, X1_Non_qtl_names) #residuals (u) instead of the gene expression for the eQTL part.
X1_new = X1_new[order(rownames(X1_new)), order(colnames(X1_new))] #order

# Split X2: [u2 | SNP_non.meqtl]
X2_new = cbind(X2_qtl_names_u, X2_Non_qtl_names) #residuals (u) instead of the gene expression for the eQTL part.
X2_new = X2_new[order(rownames(X2_new)), order(colnames(X2_new))] #order

# Split X3: [u3 | CpG_non.eqtm]
X3_new = cbind(X3_qtl_names_u, X3_Non_qtl_names) #residuals (u) instead of the gene expression for the eQTL part.
X3_new = X3_new[order(rownames(X3_new)), order(colnames(X3_new))] #order

### Y Re-weighting ###
# Y.train
Y = as.data.frame(Y)
Y.rew = Y[order(Y$y), ] #sort by y=time
sc = survfit(Surv(y, 1-status) ~ 1, data = Y.rew)$surv #survival func of censoring 
sc = c(1, sc[-nrow(Y.rew)]) #shift by 1
y = (Y.rew$y * Y.rew$status)/(sc+0.0000000001) #New-y
Y.rew = cbind(y, Y.rew[, -1, drop=FALSE])
Y.rew = Y.rew[order(rownames(Y.rew)), ] #order-back by ID
is.nan.data.frame = function(x) do.call(cbind, lapply(x, is.nan))
Y.rew[is.nan(Y.rew)] = 0


#### Final Data (Input to Cox-sMBPLS function) ####
# Train/Test (70/30)
y_finalCox = cbind.data.frame(Y.rew, Y)
num_train = floor(0.70 * nrow(X1_new))
index_train = sample(nrow(X1_new), num_train)

XX1.train = X1_new[index_train, ]
XX1.train = XX1.train[order(rownames(XX1.train)), order(colnames(XX1.train))]
XX1.test = X1_new[-index_train, ]
XX1.test = XX1.test[order(rownames(XX1.test)), order(colnames(XX1.test))]

XX2.train = X2_new[rownames(XX1.train), ]
XX2.train = XX2.train[order(rownames(XX2.train)), order(colnames(XX2.train))]
XX2.test = X2_new[rownames(XX1.test), ]
XX2.test = XX2.test[order(rownames(XX2.test)), order(colnames(XX2.test))]

XX3.train = X3_new[index_train, ]
XX3.train = XX3.train[order(rownames(XX3.train)), order(colnames(XX3.train))]
XX3.test = X3_new[-index_train, ]
XX3.test = XX3.test[order(rownames(XX3.test)), order(colnames(XX3.test))]

YY.train = y_finalCox[rownames(XX1.train), ]
YY.test = y_finalCox[rownames(XX1.test), ]


dat.train.list = list(YY.train, XX1.train, XX2.train, XX3.train) #list
names(dat.train.list) = c("Y", "GENE", "SNP", "CPG")
dat.test.list = list(YY.test, XX1.test, XX2.test, XX3.test) #list
names(dat.test.list) = c("Y.test", "GENE.test", "SNP.test", "CPG.test")


#### Parameter Tuning ####
# Tuning (1) sparsity parameter, (2) number of the latent components
source("fun.cox.smbpls.R")
ncomp_tune_set = c(2, 5, 10, 15, 25, 56, 100)
lambda_tune_set = c(0.00001, 0.0001, 0.01, 0.25, 0.63, 0.95)
tuning_result = expand.grid("ncomp_i" = ncomp_tune_set, # Components no
                            "lambda_i" = lambda_tune_set, # Sparsity
                            "cindex.train_i" = NA, 
                            "cindex.test_i" = NA,
                            stringsAsFactors = FALSE)
for (lambda_i in lambda_tune_set) {
  for (ncomp_i in ncomp_tune_set) {
    MBsPLS.result.tune = cox.smbpls(Xtrain = dat.train.list[2:4],
                                    Ytrain = dat.train.list[[1]],
                                    lambda.l1 = lambda_i,
                                    ncomp = ncomp_i,
                                    weight.mat=NULL,
                                    Xtest = dat.test.list[2:4],
                                    Ytest = dat.test.list[[1]], 
                                    center.X=TRUE, 
                                    center.Y=TRUE,
                                    scale.X=TRUE,
                                    scale.Y=TRUE,
                                    weighted.center=FALSE)
    # write the params used
    tune_lambda_result_lambda_ncomp_i = MBsPLS.result.tune$lambda.l1
    tune_ncomp_result_lambda_ncomp_i = MBsPLS.result.tune$ncomp
    train_cindex_result_lambda_ncomp_i = MBsPLS.result.tune$Cindex.train1
    test_cindex_result_lambda_ncomp_i = MBsPLS.result.tune$Cindex.test1
    # Save the results
    tuning_result[tuning_result$ncomp==ncomp_i & tuning_result$lambda_i==lambda_i, 3] = MBsPLS.result.tune$Cindex.train1
    tuning_result[tuning_result$ncomp==ncomp_i & tuning_result$lambda_i==lambda_i, 4] = MBsPLS.result.tune$Cindex.test1
  } 
}
lambda.tuned = tuning_result$lambda_i[1]
ncomp.tuned = tuning_result$ncomp_i[1]

# Final Cox-sMBPLS model using the tuned parameters
source("fun.cox.smbpls.R")
MBsPLS.result = cox.smbpls(Xtrain = dat.train.list[2:4],
                           Ytrain = dat.train.list[[1]], 
                           lambda.l1 = lambda.tuned, 
                           ncomp = ncomp.tuned,
                           weight.mat=NULL, 
                           Xtest = dat.test.list[2:4], 
                           Ytest = dat.test.list[[1]], 
                           center.X=TRUE,
                           center.Y=TRUE,
                           scale.X=TRUE, 
                           scale.Y=TRUE,
                           weighted.center=FALSE)
