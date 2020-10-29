### supervised Cox-sMBPLS main function (Input: processed -Omics matrices as a list)
# Xtrain: List of Omics-matrices. Here it is a list with three elements (1) GENE: gene-expression matrix. (2) SNP: genotypes-matrix, and (3) CPG: DNA methylation-matrix.
# Ytrain: Re-weighted survival time (Section 2.1).
# lambda.l1: Sparsity parameter. This parameter vary between 0 and 1 and is tuned using CV.
# ncomp: Number of latent components.
# Xtest: TRUE if tuning the function using a train-set and testing the performance using the test-set.
# Ytest: TRUE if tuning the function using a train-set and testing the performance using the test-set.
# center.X: TRUE if want to center the Xtrain using the mean value.
# center.Y: TRUE if want to center the Ytrain using the mean value.
# scale.X: TRUE if want to scale the Xtrain using the standard deviation value.
# scale.Y: TRUE if want to scale the Ytrain using the standard deviation value.
# nblo: Number of the blocks (=number of the Xtrain’s elements, here is 3).
# ntrain: Number of the subjects.

### Supervised Cox-sMBPLS Function ###
MB_spls3_Cox = function(Xtrain,
                        Ytrain,
                        lambda.l1,
                        ncomp,
                        weight.mat=NULL, 
                        Xtest=TRUE,
                        Ytest=TRUE,
                        center.X=TRUE,
                        center.Y=TRUE,
                        scale.X=TRUE, 
                        scale.Y=TRUE,
                        weighted.center=FALSE) {
  #####################################################################
  #### Initialization
  #####################################################################
  list_name = c("GENE", "SNP", "CPG")
  status = Ytrain[, 2, drop=F]
  y.orig = Ytrain[, 3, drop=F] #original (raw) response variable, before re-weighting
  Ytrain = Ytrain[, 1, drop=F] #re-weighted response variable
  
  if (!is.null(Xtest)) {
    status.test = Ytest[, 2, drop=F]
    y.orig.test = Ytest[, 3, drop=F] 
    Ytest = Ytest[, 1, drop=F]} 
  
  nblo = length(Xtrain)
  ntrain = nrow(Ytrain) #number of observations (subject)
  p = lapply(1:nblo, function(x) ncol(Xtrain[[x]])) #number of features per block 
  index.p = lapply(1:nblo, function(x) c(1:ncol(Xtrain[[x]])))
  names(index.p) = list_name 
  q = 1 #number of response variable/s (here for survival analysis, it is q=1)
  blo = c(p[[1]], p[[2]], p[[3]]) # needed for VIMP calculation
  ncomp.max = ncomp #number of the latent components
  
  if (!is.null(Xtest)) {
    ntest = nrow(Ytest) #number of the observations in the test test
  }
  
  cnames = list()
  cnames = lapply(1:nblo, function(x) colnames(Xtrain[[x]]))
  names(cnames) = list_name
  
  #####################################################################
  #### Tests on Input Data
  #####################################################################
  # Xtrain 
  if (!inherits(Xtrain, "list"))
    stop("object 'list' expected")
  
  # nblo (number of the blocks)
  if (nblo <2)
    stop("Minimum two blocks are required")
  if (ncomp < 0)
    ncomp = 2

  # ntrain (number of the subjects)
  if(any(row.names(Xtrain[[1]]) != row.names(Ytrain))) 
    stop("Gene-block and Ytrain must have the same rows")
  if(any(row.names(Xtrain[[2]]) != row.names(Ytrain))) 
    stop("SNP-block and Ytrain must have the same rows")
  if(any(row.names(Xtrain[[3]]) != row.names(Ytrain))) 
    stop("CpG-block and Ytrain must have the same rows")
  
  # number of the features in each block
  if(any(ncol(Xtrain[[1]]) < 2)) 
    stop("Minimum two variables for Gene-block are required")
  if(any(ncol(Xtrain[[2]]) < 2)) 
    stop("Minimum two variables for SNP-block  are required")
  if(any(ncol(Xtrain[[3]]) < 2)) 
    stop("Minimum two variables for CpG-block  are required")
  
  # Weighting matrix V
  if(!is.null(weight.mat)) { # weighting in scalar product (in observation space of dimension n)
    V = as.matrix(weight.mat) 
    
    if ((!is.matrix(V)) || (!is.numeric(V))) {
      stop("V is not of valid type")}
    
    if ((ntrain != ncol(V)) || (ntrain != nrow(V))) {
      stop("Wrong dimension for V, must be a square matrix of size the number of observations in Xtrain")
    }
  } else { # no weighting in scalar product
    V = diag(rep(1, ntrain), nrow=ntrain, ncol=ntrain)
  }
  row.names(V) = row.names(Ytrain)
  
  # lambda.l1
  if ((!is.numeric(lambda.l1)) || (lambda.l1<0) || (lambda.l1>1)) {
    stop("lambda is not of valid type")
  }
  
  # ncomp
  if ((!is.numeric(ncomp)) || (round(ncomp)-ncomp!=0) || (ncomp<1) || (ncomp>p)) {
    stop("ncomp is not of valid type")
  }
  
  # weighted.center
  if ( (weighted.center) && (is.null(weight.mat))) {
    stop("If the centering is weighted, the weighting matrix V should be provided")
  }
  
  #####################################################################
  #### Centering and Scaling
  #####################################################################
  if (!weighted.center) {
    
    # Xtrain: mean & sd
    meanXtrain = lapply(1:nblo, function (x) 
      meanXtrain = apply(Xtrain[[x]], 2, mean))
    
    sigmaXtrain = lapply(1:nblo, function (x) 
      sigmaXtrain = apply(Xtrain[[x]], 2, sd))
    
    # Xtrain: centering & scaling
    sXtrain = lapply(1:nblo, function(x)
      if(center.X && scale.X) {
        sXtrain = scale( Xtrain[[x]], center=meanXtrain[[x]], scale=sigmaXtrain[[x]]+0.000000001)
      } else if(center.X && !scale.X) {
        sXtrain = scale( Xtrain[[x]], center=meanXtrain[[x]], scale=FALSE)
      } else {
        sXtrain = Xtrain[[x]]
      })
    names(sXtrain) = list_name
    
    # Ytrain: mean & sd
    meanYtrain = apply(Ytrain, 2, mean)
    sigmaYtrain = apply(Ytrain, 2, sd)
    
    # Y-train: centering & scaling
    sYtrain = if(center.Y && scale.Y) {
      sYtrain = scale( Ytrain, center=meanYtrain, scale=sigmaYtrain+0.0000000001)
    } else if(center.Y && !scale.Y) {
      sYtrain = scale( Ytrain, center=meanYtrain, scale=FALSE )
    } else {
      sYtrain = Ytrain
    }

    # Xtest
    if (!is.null(Xtest)) {
      ## X-test: centering and scaling
      sXtest = lapply(1:nblo, function(x)
        if(center.X && scale.X) {
          sXtest = scale( Xtest[[x]], center=meanXtrain[[x]], scale=sigmaXtrain[[x]]+0.000000001)
        } else if(center.X && !scale.X) {
          sXtest = scale( Xtest[[x]], center=meanXtrain[[x]], scale=FALSE)
        } else {
          sXtest = Xtest[[x]]
        })
      names(sXtest) = list_name
    }

    # Ytest: centering & scaling
    if (!is.null(Xtest)) {
      sYtest = if(center.Y && scale.Y) {
        sYtest = scale( Ytest, center=meanYtrain, scale=sigmaYtrain+0.0000000001)
      } else if(center.Y && !scale.Y) {
        sYtest = scale( Ytest, center=meanYtrain, scale=FALSE )
      } else {
        sYtest = Ytest
      }}
    
  } else { # Weighted scaling
    sumV = sum(diag(V))
    
    # Xtrain: mean & sd
    meanXtrain = lapply(1:nblo, function (x) 
      meanXtrain = matrix(diag(V), nrow=1) %*% Xtrain[[x]] / sumV)
    
    sigmaXtrain = lapply(1:nblo, function (x) 
      sigmaXtrain = apply(Xtrain[[x]], 2, sd))
    
    # Xtrain: centering & scaling 
    sXtrain = lapply(1:nblo, function(x)
      sXtrain = scale(Xtrain[[x]], center=meanXtrain[[x]], scale=FALSE))
    
    # Ytrain: mean & sd
    meanYtrain = matrix(diag(V), nrow=1) %*% as.matrix(Ytrain) / sumV
    sigmaYtrain = apply(Ytrain, 2, sd)
    
    # Ytrain
    sYtrain = scale(Ytrain, center=meanYtrain, scale=FALSE)
    
    # Xtest:
    if (!is.null(Xtest)) {
      lapply(1:nblo, function(x)
        sXtest = scale(Xtest[[x]], center=meanXtrain[[x]], scale=FALSE))
    }
    if (!is.null(Xtest)) {
      # Ytest: 
      sYtest = scale(Ytest, center=meanYtrain, scale=FALSE)
    }
  } 
  
  #####################################################################
  #### Result objects
  #####################################################################
  X1 = sXtrain
  Y1 = sYtrain
  X = cbind.data.frame(sXtrain) 
  ncolX = ncol(X)
  if (!is.null(Xtest)) {
    X.test = cbind.data.frame(sXtest) # scaled & Residual-updated Xtest
    ncolX.test = ncol(X.test)
    }
  
  dimlab = paste("Ax", 1:ncomp.max, sep = "")
  betahat = lapply(1:nblo, function(b) matrix(0, nrow = ncol(Xtrain[[b]]), ncol = 1, dimnames = list(colnames(Xtrain[[b]]), paste("betahat")))); names(betahat) = list_name
  betamat = lapply(1:nblo, function(b) list());  names(betamat) = list_name

  W = matrix(data=NA, nrow=ncolX, ncol=ncomp.max, dimnames = list(colnames(X), dimlab)) # spls weights; p is the total number of features 
  T = matrix(data=NA, nrow=ntrain, ncol=ncomp.max, dimnames = list(row.names(Ytrain), dimlab)) # spls latent components components
  P = matrix(data=NA, nrow=ncomp.max, ncol=ncolX, dimnames = list(dimlab, colnames(X))) # regression of X over latent components
  Q = matrix(data=NA, nrow=ncomp.max, ncol=q, dimnames = list(dimlab, paste("Loading.y"))) # regression of Y over latent components
  
  lX  = matrix(0, nrow=ntrain, ncol=ncomp.max, dimnames=list(row.names(Ytrain), dimlab))
  lX1 = matrix(0, nrow=ntrain, ncol=ncomp.max, dimnames=list(row.names(Ytrain), dimlab)) 

  Ak = matrix(0, nrow = nblo, ncol = ncomp.max, dimnames=list(list_name, dimlab))
  cov2 = matrix(0, nrow = nblo, ncol = ncomp.max, dimnames=list(list_name, dimlab)) 
  
  what = lapply(1:nblo, function(b)  matrix(0, nrow = ncol(Xtrain[[b]]), ncol = ncomp.max, dimnames = list(colnames(Xtrain[[b]]), dimlab)))
  T_p = rep(list(matrix(0, nrow = ntrain, ncol = ncomp.max, dimnames = list(row.names(Ytrain), dimlab))), nblo)
  
  ust = function(w, lambda.l1) {
    w.th = matrix(0, nrow=length(w), ncol=1)
    val.w = abs(w) - lambda.l1 * max( abs(w) )
    w.th[ val.w >=0 ] = val.w[ val.w>=0 ] * (sign(w))[ val.w>=0 ]
    return(w.th)
  } # soft thresholding function
  
  #####################################################################
  #### Main iteration
  #####################################################################
  
  lapply(1:nblo, function(x)
    if (is.null(colnames(Xtrain[[x]]))) {
      Xnames = index.p[[x]]
    } else { 
      Xnames = colnames(Xtrain[[x]])})

    
  new2As = lapply(1:nblo, function(b) list()); names(new2As) = list_name
  for (k in 1:ncomp.max) {
    M = list()
    Mnorm1 = list()
    # Soft thresholding
    M[[1]] = t(X1[[1]]) %*% (V %*% Y1)
    M[[2]] = t(X1[[2]]) %*% (V %*% Y1)
    M[[3]] = t(X1[[3]]) %*% (V %*% Y1) 
    
    library(base)
    Mnorm1[[1]] = norm(M[[1]], "2")
    Mnorm1[[2]] = norm(M[[2]], "2")
    Mnorm1[[3]] = norm(M[[3]], "2")
    
    M[[1]] = M[[1]] / Mnorm1[[1]]
    M[[2]] = M[[2]] / Mnorm1[[2]]
    M[[3]] = M[[3]] / Mnorm1[[3]]
    
    covutk = rep(0, nblo)
    plsfit = list()
    A = list()
    new2A = list()
    X.A = list()
    what[[1]] = ust(M[[1]], lambda.l1)  # w-partial (TCL)
    what[[2]] = ust(M[[2]], lambda.l1)
    what[[3]] = ust(M[[3]], lambda.l1)
    
    T_p[[1]] = X1[[1]] %*% what[[1]] #t-partial (TLX)
    T_p[[2]] = X1[[2]] %*% what[[2]]
    T_p[[3]] = X1[[3]] %*% what[[3]]
    
    covutk[1] = crossprod(as.matrix(Y1), T_p[[1]])
    covutk[2] = crossprod(as.matrix(Y1), T_p[[2]])
    covutk[3] = crossprod(as.matrix(Y1), T_p[[3]])
    
    cov2[1, k] = covutk[1]^2 # cov2 is "b=Tu" in MBPLS paper # We need cov2 only to calculate the the final T (== lx), so we don't need to have it outside the b-for-loop.
    cov2[2, k] = covutk[2]^2 
    cov2[3, k] = covutk[3]^2 
    
    # Active set A
    A[[1]] = unique( index.p[[1]][what[[1]] !=0 | betahat[[1]][, 1] !=0] )
    A[[2]] = unique( index.p[[2]][what[[2]] !=0 | betahat[[2]][, 1] !=0] )
    A[[3]] = unique( index.p[[3]][what[[3]] !=0 | betahat[[3]][, 1] !=0] )
    
    new2A[[1]] = index.p[[1]][ what[[1]] !=0 & betahat[[1]][, 1] == 0 ] # Can I change the B-hat cut-off!!! (Like: <<0.001) :)
    new2A[[2]] = index.p[[2]][ what[[2]] !=0 & betahat[[2]][, 1] == 0 ] 
    new2A[[3]] = index.p[[3]][ what[[3]] !=0 & betahat[[3]][, 1] == 0 ] 
    
    
    # PLS with selected features (included in A) for each block
    X.A[[1]] = sXtrain[[1]][ , A[[1]], drop=FALSE]
    X.A[[2]] = sXtrain[[2]][ , A[[2]], drop=FALSE]
    X.A[[3]] = sXtrain[[3]][ , A[[3]], drop=FALSE]

    plsfit[[1]] = wpls(Xtrain=X.A[[1]], 
                       Ytrain=sYtrain, 
                       weight.mat=V, 
                       ncomp=min(k,length(A[[1]])), 
                       type="pls1", 
                       center.X=FALSE, 
                       scale.X=FALSE, 
                       center.Y=FALSE, 
                       scale.Y=FALSE, 
                       weighted.center=FALSE )
    plsfit[[2]] = wpls(Xtrain=X.A[[2]], 
                       Ytrain=sYtrain, 
                       weight.mat=V, 
                       ncomp=min(k,length(A[[2]])), 
                       type="pls1", 
                       center.X=FALSE, 
                       scale.X=FALSE, 
                       center.Y=FALSE, 
                       scale.Y=FALSE, 
                       weighted.center=FALSE )
    plsfit[[3]] = wpls(Xtrain=X.A[[3]], 
                       Ytrain=sYtrain, 
                       weight.mat=V, 
                       ncomp=min(k,length(A[[3]])), 
                       type="pls1", 
                       center.X=FALSE, 
                       scale.X=FALSE, 
                       center.Y=FALSE, 
                       scale.Y=FALSE, 
                       weighted.center=FALSE )

    Ak[1, k] = covutk[1] / sqrt(sum(cov2))
    Ak[2, k] = covutk[2] / sqrt(sum(cov2))
    Ak[3, k] = covutk[3] / sqrt(sum(cov2))

    lX[, k] = lX[, k] + Ak[1, k] * T_p[[1]] 
    lX[, k] = lX[, k] + Ak[2, k] * T_p[[2]]
    lX[, k] = lX[, k] + Ak[3, k] * T_p[[3]]

    
    ## Output Storage
    # Weights (w)
    w.global.k = do.call("rbind", what)
    w.global.k = matrix(w.global.k, ncol=1)
    w.global.k = w.global.k / sqrt(as.numeric(t(w.global.k) %*% w.global.k))
    W[,k] = w.global.k
    
    # Latent components (t)
    t.global.k = lX[, k]
    T[, k] = t.global.k
    
    # P
    p.k = (t(X) %*% (V %*% t.global.k)) / as.numeric(t(t.global.k) %*% (V %*% t.global.k)) 
    P[k, ] = t(p.k)
    
    # Q
    q.k = (t(Y1) %*% (V %*% t.global.k)) / as.numeric(t(t.global.k) %*% (V %*% t.global.k))
    Q[k,] = t(q.k)

    
    ## Update
    X1 = sXtrain
    X1[[1]][, A[[1]]] = (sXtrain[[1]][, A[[1]]]) - (plsfit[[1]]$T %*% plsfit[[1]]$P)
    X1[[2]][, A[[2]]] = (sXtrain[[2]][, A[[2]]]) - (plsfit[[2]]$T %*% plsfit[[2]]$P)
    X1[[3]][, A[[3]]] = (sXtrain[[3]][, A[[3]]]) - (plsfit[[3]]$T %*% plsfit[[3]]$P)

    betahat[[1]] = matrix(0, ncol(Xtrain[[1]]), q, dimnames = list(colnames(Xtrain[[1]]), paste("beta_hat(PLS.Coeff)")));  names(betahat) = list_name
    betahat[[2]] = matrix(0, ncol(Xtrain[[2]]), q, dimnames = list(colnames(Xtrain[[2]]), paste("beta_hat(PLS.Coeff)")));  names(betahat) = list_name
    betahat[[3]] = matrix(0, ncol(Xtrain[[3]]), q, dimnames = list(colnames(Xtrain[[3]]), paste("beta_hat(PLS.Coeff)")));  names(betahat) = list_name
    
    betahat[[1]][A[[1]], ] = matrix(plsfit[[1]]$coeff, length(A[[1]]), q, dimnames = list(colnames(X.A[[1]]), paste("beta_hat(PLS.Coeff)")))
    betahat[[2]][A[[2]], ] = matrix(plsfit[[2]]$coeff, length(A[[2]]), q, dimnames = list(colnames(X.A[[2]]), paste("beta_hat(PLS.Coeff)")))
    betahat[[3]][A[[3]], ] = matrix(plsfit[[3]]$coeff, length(A[[3]]), q, dimnames = list(colnames(X.A[[3]]), paste("beta_hat(PLS.Coeff)")))
    betahat.bind = do.call(rbind, betahat)
    
    betamat[[1]][[k]] = betahat[[1]] # for cv.spls
    betamat[[2]][[k]] = betahat[[2]]
    betamat[[3]][[k]] = betahat[[3]]

    new2As[[1]][[k]] = new2A[[1]]
    new2As[[2]][[k]] = new2A[[2]]
    new2As[[3]][[k]] = new2A[[3]]
    
    #Y1 <- sYtrain - as.matrix(X) %*% betahat.bind #Alter!
    plsfit.T.bind = as.matrix(cbind.data.frame(plsfit[[1]]$T, plsfit[[2]]$T, plsfit[[3]]$T))
    plsfit.Q.bind = as.matrix(rbind.data.frame(plsfit[[1]]$Q, plsfit[[2]]$Q, plsfit[[3]]$Q))
    Y1 = sYtrain - plsfit.T.bind %*% plsfit.Q.bind
    
  } # End of the main loop
  
  
  #####################################################################
  #### Return Final Result
  #####################################################################
  hatY = numeric(ntrain)
  bip = Ak^2
  Anames = list()
  for (b in 1:nblo) {
    Anames[[b]] = cnames[[b]][A[[b]]]
  }

  # Final data (train) to fit the Cox-PH model with latent components as the covariates
  data.final.cox = cbind(y.orig, status, T)
  data.final.cox = as.data.frame(data.final.cox)
  Cox.mbspls = coxph(Surv(y, status) ~ ., data = data.final.cox, x=TRUE, y=TRUE)
  # C-index train
  cindex.train1 = summary(Cox.mbspls)$concordance[1]
  
  
  
  # Final data (test) preparation
  ww = W
  bb2 = cov2
  
  # Weights (w)
  w_gene = ww[rownames(ww) %like% "GENE.", ]
  w_gene = as.matrix(w_gene)
  w_snp = ww[rownames(ww) %like% "SNP.X", ]
  w_snp = as.matrix(w_snp)
  w_cpg = ww[rownames(ww) %like% "CPG.", ]
  w_cpg = as.matrix(w_cpg)
  
  # block weights (b)
  b_gene = bb2[1, ] #GENE
  b_snp = bb2[2, ] #SNP
  b_cpg = bb2[3, ] #SNP

  # bw
  dimlab = paste("Ax", 1:ncomp, sep = "")
  bw_gene = matrix(0, nrow=nrow(w_gene), ncol = ncomp, dimnames = list(row.names(w_gene), dimlab))
  bw_snp = matrix(0, nrow=nrow(w_snp), ncol = ncomp, dimnames = list(row.names(w_snp), dimlab))
  bw_cpg = matrix(0, nrow=nrow(w_cpg), ncol = ncomp, dimnames = list(row.names(w_cpg), dimlab))

  for(k in 1:ncomp) {
    bw_gene[, k] = b_gene[k] * w_gene[, k] 
    bw_snp[, k] = b_snp[k] * w_snp [, k]
    bw_cpg[, k]  = b_cpg[k] * w_cpg[, k]
  }
  bw = rbind(bw_gene, bw_snp, bw_cpg)

  # Latent components (T = original.X %*% bw)
  orig.x = X.test
  Ts = as.matrix(orig.x) %*% as.matrix(bw)
  
  # Final test data
  data.final.cox.test = cbind.data.frame(y.orig.test, status.test, as.data.frame(Ts)) #Orig.Y
  # C-index test
  pred.test = predict(Cox.mbspls, data.final.cox.test)
  cindex.test2 = concordance.index(x=pred.test, surv.time=data.final.cox.test$y, surv.event=data.final.cox.test$status, method="noether")
  cindex.test2 = cindex.test2$c.index

  # RESULTS # 
  result = list(Xtrain = Xtrain,
                Ytrain = Ytrain,
                sXtrain = sXtrain,
                sYtrain = sYtrain,
                
                meanXtrain = meanXtrain,
                meanYtrain = meanYtrain, 
                sigmaXtrain = sigmaXtrain,
                sigmaYtrain = sigmaYtrain,
                
                betahat = betahat.bind,
                X.score = T[, 1:ncomp, drop=F],
                X.loading = t(t(P)[, 1:ncomp, drop=F]),
                Y.loading = t(t(Q)[, 1:ncomp, drop=F]), 
                X.weight = W[, 1:ncomp, drop=F],
                
                A = A, 
                Anames = Anames,
                betamat = betamat,
                new2As = new2As,
                lambda.l1 = lambda.l1,
                ncomp = ncomp,
                V = V,
                
                block.weights = cov2[, 1:ncomp, drop=F],
                block.weights2 = Ak[, 1:ncomp, drop=F],
                block.importance = bip[, 1:ncomp, drop=F],
                
                data.train.final.cox = data.final.cox,
                Cox.mbspls = Cox.mbspls,
                Cindex.train1 = cindex.train1,
                data.test.final.cox = data.final.cox.test,
                pred.test = pred.test,
                Cindex.test2 = cindex.test2)
  
  class(result) = c("Multiblock", "MB_spls_Cox", "cox-PH")
  return(result)
  
} # End of the Cox-sMBPLS function
