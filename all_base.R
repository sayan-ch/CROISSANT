library(parallel)
library(dplyr)
library(RSpectra)
library(cluster)
library(Rfast)
library(data.table)
library(readr)
library(irlba)
library(softImpute)
library(Matrix)
library(pROC)
library(IMIFA)
library(rlist)
library(latentnet)

## CROISSANT

##SBM or DCBM generator

blockmodel.gen.fast <- function(n, K, B, psi = rep(1,n), 
                                PI = rep(1/K, K), ncore = 1)
{
  on.exit(gc())
  
  g <- sample.int(K, n, T, PI)
  
  #max.psi = 1 for each community
  psi.scale <- psi
  for(kk in 1:K)
    psi.scale[g==kk] <- psi[g==kk]/max(psi[g==kk])
  
  stor <- do.call('rbind',
                  mclapply(1:(n-1), function(i) {
                    tmp <- which(rbinom(n-i, 1, 
                                        B[g[i],g[(i+1):n]]*psi.scale[i]*psi.scale[(i+1):n]) == 1)
                    
                    if(length(tmp) == 0)
                      return(NULL)
                    else
                      return(cbind(rep(i, length(tmp)), i + tmp))
                  }, mc.cores = ncore))
  
  A <- sparseMatrix(stor[,1], stor[,2], dims = c(n,n), symmetric = T)
  
  return(list(A = A, member = g, psi = psi.scale))
}

## Required functions
sumFast <- function(X){
  if(is.vector(X))  return(sum(X))
  return(sum(rowSums(X)))
}

modal <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

best.perm.label.match <- function(lab, fixed, 
                                  n = length(lab), K = max(lab, fixed)){
  
  if(identical(lab, fixed))
    return(diag(1, K))
  
  if(K == 1)
    return(matrix(1,1,1))
  
  if(K == 2){
    if(sum(lab!=fixed) <= n/2)
      return(diag(1,2))
    else
      return(matrix(c(0,1,1,0),2,2,T))
  }
  
  E <- matrix(0, K, K)
  
  C.lab <- as(sparseMatrix(i = 1:n, j = lab, dims = c(n, K)), 'dMatrix')
  C.fixed <- as(sparseMatrix(i = 1:n, j = fixed, dims = c(n, K)), 'dMatrix')
  M <- crossprod(C.lab, C.fixed)
  
  # M <- Rfast::Table(lab, fixed)
  while(max(M) != -1)
  {
    ind <- which(M == max(M), T)[1,]
    E[ind[2],ind[1]] <- 1
    M[ind[1],] <- rep(-1,K)
    M[,ind[2]]  <- rep(-1,K)
  }
  return(E)
}




matched.lab <- function(lab, fixed, 
                        n = length(lab), K = max(lab, fixed)){
  
  E <- best.perm.label.match(lab, fixed, n = n, K = K)
  
  lmat <- sparseMatrix(i = 1:n, j = lab, dims = c(n,K))
  
  as.vector(tcrossprod(tcrossprod(lmat, E), rbind(1:K)))
}


fast.SBM.est <- function(A, g, n = nrow(A), K = max(g)){
  
  B <- matrix(0, K, K)
  if(K == 1){
    B[K,K] <- sumFast(A)/(n^2-n)
    return(B)
  }
  
  G <- lapply(1:K, function(k) which(g == k))
  nk <- sapply(G, 'length')
  
  for(k in 1:K){
    for(l in k:K){
      B[k,l] <- B[l,k] <- sumFast(A[G[[k]], G[[l]]])/(nk[k]*nk[l])
    }
  }
  
  diag(B) <- diag(B)*nk/(nk-1)
  B[!is.finite(B)] <- 1e-6
  
  return(B)        
}

fast.DCBM.est <- function(A, g, n = nrow(A), K = max(g),
                          psi.omit = 0, p.sample = 1){
  
  B.sum <- matrix(0, K, K)
  if(K == 1){
    B.sum[K,K] <- sumFast(A) + 0.001
    
    if(psi.omit > 0){
      psi <- as.numeric(rowSums(A[-(1:psi.omit),])/B.sum[K,K])
      return(list(Bsum = B.sum/p.sample, psi = psi))
    }
    
    psi <- as.numeric(rowSums(A)/B.sum[K,K]) + 0.0001
    
    return(list(Bsum = B.sum/p.sample, psi = psi))
  }
  
  G <- lapply(1:K, function(k) which(g == k))
  
  for(k in 1:K){
    for(l in k:K){
      B.sum[k,l] <- B.sum[l,k] <- sumFast(A[G[[k]], G[[l]]]) + 0.0001
    }
  }
  
  if(psi.omit > 0){
    psi <- as.numeric(rowSums(A[-(1:psi.omit),])/
                        rowSums(B.sum)[g[-(1:psi.omit)]])
    return(list(Bsum = B.sum/p.sample, psi = psi))
  }
  
  psi <- as.numeric(rowSums(A)/rowSums(B.sum)[g])
  
  return(list(Bsum = B.sum/p.sample, psi = psi))
}

eigen.DCBM.est <- function(A, g, n = nrow(A), K = max(g),
                           psi.omit = 0, p.sample = 1){
  
  U.hat <- irlba::irlba(A, nu = K, nv = K)$v
  
  psi.hat <- rowSums(U.hat^2)^0.5
  
  psi.outer <- tcrossprod(psi.hat)
  
  B.sum <- matrix(0, K, K)
  if(K == 1){
    B.sum[K,K] <- sumFast(A)/sumFast(psi.outer)
    
    if(psi.omit > 0){
      #psi <- as.numeric(rowSums(A[-(1:psi.omit),])/B.sum[K,K])
      return(list(Bsum = B.sum/p.sample, psi = psi.hat[-(1:psi.omit)]))
    }
    
    #psi <- as.numeric(rowSums(A)/B.sum[K,K]) + 0.001
    
    return(list(Bsum = B.sum/p.sample, psi = psi.hat))
  }
  
  G <- lapply(1:K, function(k) which(g == k))
  
  for(k in 1:K){
    for(l in k:K){
      B.sum[k,l] <- B.sum[l,k] <- 
        sumFast(A[G[[k]], G[[l]]])/sumFast(psi.outer[G[[k]], G[[l]]])
    }
  }
  
  if(psi.omit > 0){
    #psi <- as.numeric(rowSums(A[-(1:psi.omit),])/
    #                    rowSums(B.sum)[g[-(1:psi.omit)]])
    return(list(Bsum = B.sum/p.sample, psi = psi.hat[-(1:psi.omit)]))
  }
  
  #psi <- as.numeric(rowSums(A)/rowSums(B.sum)[g])
  
  return(list(Bsum = B.sum/p.sample, psi = psi.hat))
}

l2 <- function(x,y){
  sqrt(sum(rowSums((x-y)^2)))
}

bin.dev <- function(x,y){
  tmp <- -x*log(y) - (1-x)*log(1-y)
  
  tmp[!is.finite(tmp)] <- 0
  
  return(sum(rowSums(tmp)))
}

auroc <- function(score, bool) {
  n1 <- sum(!bool)
  #n2 <- sum(bool)
  n2 <- length(score) - n1
  U  <- sum(rank(score)[!bool]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}

AUC <- function(A, P){
  -Rfast::auc(as(A, "vector"), as(P, "vector"))
  #-auroc(as(P, 'vector'), as(A, 'vector'))
}



################################################################################

croissant.blockmodel <- function(A, K.CAND,
                           s, o, R, 
                           tau = 0,
                           laplace = F,
                           dc.est = 2,
                           loss = c("l2", "bin.dev", "AUC"),
                           ncore = 1){
  
  if(length(K.CAND) == 1) K.CAND <- 1:K.CAND
  
  K.max <- max(K.CAND)
  
  n <- nrow(A)
  m <- (n-o)/s
  
  L <- list()
  
  mod <- c("SBM", "DCBM")
  #mod <- mod.cand
  
  over <- lapply(1:R, function(ii) sample.int(n, o, F))
  non.over <- lapply(1:R, function(ii) sample((1:n)[-over[[ii]]], n-o, replace = F))
  
  raw.ind <- cbind(rep(1:R, each = s), rep(1:s, R))
  
  raw.out <- mclapply(1:nrow(raw.ind), function(ii){
    
    q <- raw.ind[ii, 2]
    r <- raw.ind[ii, 1]
    
    sonn <- c(over[[r]], non.over[[r]][((q-1)*m+1):(q*m)])
    A.sonn <- A[sonn, sonn]
    
    deg <- rowSums(A.sonn)
    avg.deg <- mean(deg)
    
    A.sonn.tau <- A.sonn + tau*avg.deg/n
    d.sonn.tau <- sparseMatrix( i = 1:(o+m), j = 1:(o+m),
                                x = 1/sqrt(deg + tau*avg.deg))
    
    
    L.sonn <- A.sonn.tau
    if(laplace){
      L.sonn <- tcrossprod(crossprod(d.sonn.tau, A.sonn.tau), d.sonn.tau)
      L.sonn[is.na(L.sonn)] <- 0
      }
    
    
    #eig.max <- eigs_sym(L.sonn, K.max, "LM")$vectors
    eig.max <- irlba::partial_eigen(x = L.sonn, n = K.max, 
                                    symmetric = T)$vectors
    
    out.SBM <- list()
    out.DCBM <- list()
    #psi.hat <- list()
    
    for(k.cand in seq_along(K.CAND)){
      if(K.CAND[k.cand] == 1){
        out.SBM[[k.cand]] <- out.DCBM[[k.cand]] <- rep(1,o+m)
        next
      }
      
      work.K <- K.CAND[[k.cand]]
      # out.SBM[[k.cand]] <- as.integer(kmeans(eig.max[,1:work.K],
      #                                        work.K,
      #                                        nstart = 100,
      #                                        iter.max = 10000)$cluster)
      out.SBM[[k.cand]] <- as.integer(pam(eig.max[,1:work.K], work.K,
                                          metric = "euclidean",
                                          do.swap = F,
                                          cluster.only = T,
                                          pamonce = 6))
      
      #psi.hat[[k.cand]] <- sqrt(rowSums(eig.max[,1:work.K]^2))
      rownorm <- sqrt(rowSums(eig.max[, 1:work.K]^2))
      rownorm[rownorm == 0] <- 1
      
      rn.eig <- eig.max[,1:work.K]/rownorm
      
      out.DCBM[[k.cand]] <- as.integer(pam(rn.eig, work.K,
                                          metric = "euclidean",
                                          do.swap = F,
                                          cluster.only = T,
                                          pamonce = 6))
      
      # out.DCBM[[k.cand]] <- as.integer(kmeans(rn.eig, work.K,
      #                                         nstart = 100,
      #                                         iter.max = 10^7)$cluster)
    }
    # if(dc.est == 2){
    #   return(list('SBM' = out.SBM, 'DCBM' = out.DCBM))
    # }
    
    return(list('SBM' = out.SBM, 'DCBM' = out.DCBM))
    # 'psi' = psi.hat))
  },
  mc.cores = ncore)
  
  K.size <- length(K.CAND)
  
  est.out <- mclapply(1:(K.size*nrow(raw.ind)), function(ii){
    
    k.cand <- ii %% K.size
    k.cand <- ifelse(k.cand == 0, K.size, k.cand)
    
    rot <- ceiling(ii / K.size)
    
    q <- raw.ind[rot, 2]
    r <- raw.ind[rot, 1]
    
    #message(paste0("Est. at s=",q, " started"))
    sonn <- c(over[[r]], non.over[[r]][((q-1)*m+1):(q*m)])
    A.sonn <- A[sonn, sonn]
    
    out.SBM.std <- raw.out[raw.ind[,1] == r][[1]]$SBM[[k.cand]]
    out.DCBM.std <- raw.out[raw.ind[,1] == r][[1]]$DCBM[[k.cand]]
    
    out.SBM <- raw.out[raw.ind[,1] == r][[q]]$SBM[[k.cand]]
    out.DCBM <- raw.out[raw.ind[,1] == r][[q]]$DCBM[[k.cand]]
    
    
    work.K <- K.CAND[[k.cand]]
    
    if(work.K == 1){
      mat.SBM <- mat.DCBM <- rep(1, m)
      
      B.SBM <- fast.SBM.est(A.sonn, rep(1,o+m), o+m, 1)
      mat.SBM <- rep(1, m)
      
      if(dc.est == 2){
        tmp <- fast.DCBM.est(A.sonn, rep(1,o+m), o+m, 1, o,
                             p.sample = 1)
      }else{
        #psi.hat <- raw.out[raw.ind[,1] == r][[q]]$psi[[k.cand]]
        tmp <- eigen.DCBM.est(A.sonn, rep(1,o+m), o+m, 1, o,
                              p.sample = 1)
      }
      
      B.DCBM <- tmp$Bsum
      psi.DCBM <- tmp$psi
      mat.DCBM <- rep(1, m)
      
      return(list('gSBM' = mat.SBM,
                  'BSBM' = B.SBM,
                  'gDCBM' = mat.DCBM,
                  'BDCBM' = B.DCBM,
                  'psiDCBM' = psi.DCBM))
    }
    
    E.SBM.kc <- best.perm.label.match(out.SBM[1:o], 
                                      out.SBM.std[1:o],
                                      o, work.K)
    
    E.DCBM.kc <- best.perm.label.match(out.DCBM[1:o], 
                                       out.DCBM.std[1:o],
                                       o, work.K)
    
    tmp.SBM <- sparseMatrix(i = 1:(o+m), j = out.SBM, 
                            dims = c((o+m),work.K))
    
    tmp.DCBM <- sparseMatrix(i = 1:(o+m), j = out.DCBM, 
                             dims = c((o+m),work.K))
    
    mat.SBM <- as.vector(tcrossprod(tcrossprod(tmp.SBM, E.SBM.kc),
                                    rbind(1:work.K)))
    
    mat.DCBM <- as.vector(tcrossprod(tcrossprod(tmp.DCBM, E.DCBM.kc),
                                     rbind(1:work.K)))
    
    
    B.SBM <- fast.SBM.est(A.sonn, mat.SBM, o+m, work.K)
    mat.SBM <- mat.SBM[-(1:o)]
    
    if(dc.est == 2){
      tmp <- fast.DCBM.est(A.sonn, mat.DCBM, o+m, work.K, o, 
                           p.sample = 1)
    }else{
      #psi.hat <- raw.out[raw.ind[,1] == r][[q]]$psi[[k.cand]]
      tmp <- eigen.DCBM.est(A.sonn, mat.DCBM, o+m, work.K, o, 
                            p.sample = 1)
    }
    
    B.DCBM <- tmp$Bsum
    psi.DCBM <- tmp$psi
    mat.DCBM <- mat.DCBM[-(1:o)]
    
    message(paste0("Est. at s=",q, " finished"))
    
    return(list('gSBM' = mat.SBM,
                'BSBM' = B.SBM,
                'gDCBM' = mat.DCBM,
                'BDCBM' = B.DCBM,
                'psiDCBM' = psi.DCBM))
  },
  mc.cores = ncore)
  
  g.SBM <- list()
  B.SBM <- list()
  g.DCBM <- list()
  B.DCBM <- list()
  psi.DCBM <- list()
  
  raw.mat <- cbind(raw.ind[rep(1:nrow(raw.ind), each = K.size), ], 
                   rep(1:K.size, nrow(raw.ind)))
  
  for(r in 1:R){
    g.SBM[[r]] <- list()
    B.SBM[[r]] <- list()
    g.DCBM[[r]] <- list()
    B.DCBM[[r]] <- list()
    psi.DCBM[[r]] <- list()
    for(k.cand in seq_along(K.CAND)){
      tmp.est <- est.out[which(raw.mat[,3] == k.cand & raw.mat[,1] == r)]
      B.SBM[[r]][[k.cand]] <- 0
      B.DCBM[[r]][[k.cand]] <- 0
      
      g.SBM[[r]][[k.cand]] <- list()
      g.DCBM[[r]][[k.cand]] <- list()
      psi.DCBM[[r]][[k.cand]] <- list()
      
      for(q in 1:s){
        B.SBM[[r]][[k.cand]] <- B.SBM[[r]][[k.cand]] + 
          tmp.est[[q]]$BSBM/s
        B.DCBM[[r]][[k.cand]] <- B.DCBM[[r]][[k.cand]] + 
          tmp.est[[q]]$BDCBM/s
        
        g.SBM[[r]][[k.cand]][[q]] <- tmp.est[[q]]$gSBM
        g.DCBM[[r]][[k.cand]][[q]] <- tmp.est[[q]]$gDCBM
        psi.DCBM[[r]][[k.cand]][[q]] <- tmp.est[[q]]$psiDCBM
      }
    }
  }
  
  non.size <- s*(s-1)/2
  non.mat <- matrix(nrow = R*non.size*K.size, ncol = 4)
  cc <- 1
  for(r in 1:R)
    for(k.cand in seq_along(K.CAND))
      for(p in 1:(s-1))
        for(q in (p+1):s){
          non.mat[cc, ] <- c(r, k.cand, p, q)
          cc <- cc + 1
        }
  
  L.all <- mclapply(1:nrow(non.mat), function(ii){
    r <- non.mat[ii, 1]
    k.cand <- non.mat[ii, 2]
    p <- non.mat[ii, 3]
    q <- non.mat[ii, 4]
    
    p.non <- non.over[[r]][((p-1)*m+1):(p*m)]
    q.non <- non.over[[r]][((q-1)*m+1):(q*m)]
    
    A.non <- A[p.non, q.non]
    
    L.temp <- matrix(0, nrow = 2*length(loss), ncol = 1)
    row.names(L.temp) <- paste(rep(mod, each = length(loss)), rep(loss, 2), 
                               sep = "_")
    colnames(L.temp) <- as.character(K.CAND[k.cand])
    
    
    P.SBM <- B.SBM[[r]][[k.cand]][g.SBM[[r]][[k.cand]][[p]],
                                  g.SBM[[r]][[k.cand]][[q]] ]
    
    
    P.DCBM <- B.DCBM[[r]][[k.cand]][g.DCBM[[r]][[k.cand]][[p]],
                                    g.DCBM[[r]][[k.cand]][[q]] ] *
      tcrossprod(psi.DCBM[[r]][[k.cand]][[p]],
                 psi.DCBM[[r]][[k.cand]][[q]])
    
    P.DCBM[P.DCBM < 1e-6] <- 1e-6
    P.DCBM[P.DCBM > 1- 1e-6] <- 1 - 1e-6
    
    for(mq in seq_along(mod)){
      if(mod[mq] == "SBM"){
        for(lq in seq_along(loss)){
          tmp.nm <- paste(mod[mq], loss[lq], sep = "_")
          L.temp[tmp.nm, 1] <-  L.temp[tmp.nm, 1] +
            do.call(loss[lq], list(as.numeric(A.non), P.SBM))/(s*(s-1)*0.5) 
        }
        next}
      else{
        for(lq in seq_along(loss)){
          tmp.nm <- paste(mod[mq], loss[lq], sep = "_")
          L.temp[tmp.nm, 1] <-  L.temp[tmp.nm, 1] +
            do.call(loss[lq], list(as.numeric(A.non), P.DCBM))/(s*(s-1)*0.5)
        }
      }
    }
    #}
    return(L.temp)},
    mc.cores = ncore)
  
  for(r in 1:R){
    L[[r]] <- do.call('cbind', lapply(1:K.size, function(kk){
      Reduce('+', L.all[which(non.mat[,2] == kk & non.mat[,1] == r)])
    }))
    
    row.names(L[[r]]) <- paste(rep(mod, each = length(loss)), rep(loss, 2), 
                               sep = "_")
    colnames(L[[r]]) <- as.character(K.CAND)
  }
  
  obj <- data.table(`Candidate Model` = rep(mod, length(K.CAND)),
                    `Candidate Value` = rep(K.CAND, each = 2)
  )
  
  
  for(lq in seq_along(loss))
    for(r in 1:R){
      obj[[paste0(loss[lq], "-Rep=", r)]] <- 
        as(rbind(L[[r]][paste0("SBM_",loss[lq]), ],
                 L[[r]][paste0("DCBM_",loss[lq]), ]), "vector")
    }
  
  obj2 <- list()
  
  obj2[["Candidate Models"]] <- "SBM and DCBM"
  
  obj2[["Candidate Values"]] <- K.CAND
  
  for(lq in seq_along(loss)){
    
    obj2[[paste0("Mod.K.hat.each.rep (", loss[lq], ")")]] <- sapply(1:R, function(r){
      l.sbm <- min(L[[r]][paste0("SBM_", loss[lq]),])
      l.dcbm <- min(L[[r]][paste0("DCBM_", loss[lq]),])
      ifelse(l.dcbm < l.sbm,
             paste0("DCSBM-",
                    K.CAND[which.min(L[[r]][paste0("DCBM_", loss[lq]),])]),
             paste0("SBM-",
                    K.CAND[which.min(L[[r]][paste0("SBM_", loss[lq]),])])
      )
    })
    
    
    obj2[[paste0(loss[lq], ".model")]] <- 
      modal(obj2[[paste0("Mod.K.hat.each.rep (", loss[lq], ")")]])
  }
  
  return(c(list('loss' = obj), 
           obj2))
}

################################################################################
################################################################################
##ECV with l2, bin.dev, and AUC
ECV.for.blockmodel <- function (A, max.K, cv = NULL, B = 3, holdout.p = 0.1, tau = 0, 
                                dc.est = 2, kappa = NULL) 
{
  n <- nrow(A)
  edge.index <- which(upper.tri(A))
  edge.n <- length(edge.index)
  holdout.index.list <- list()
  if (is.null(cv)) {
    holdout.n <- floor(holdout.p * edge.n)
    for (j in 1:B) {
      holdout.index.list[[j]] <- sample(x = edge.n, size = holdout.n)
    }
  }
  else {
    sample.index <- sample.int(edge.n)
    max.fold.num <- ceiling(edge.n/cv)
    fold.index <- rep(1:cv, each = max.fold.num)[edge.n]
    cv.index <- fold.index[sample.index]
    B <- cv
    for (j in 1:B) {
      holdout.index.list[[j]] <- which(cv.index == j)
    }
  }
  result <- lapply(holdout.index.list, holdout.evaluation.fast.all, 
                   A = A, max.K = max.K, tau = tau, dc.est = dc.est, p.sample = 1 - 
                     holdout.p, kappa = kappa)
  dc.block.err.mat <- dc.loglike.mat <- bin.dev.mat <- roc.auc.mat <- impute.err.mat <- block.err.mat <- loglike.mat <- matrix(0, 
                                                                                                                               nrow = B, ncol = max.K)
  sbm.auc.mat <- dc.auc.mat <- matrix(0, nrow = B, ncol = max.K)
  
  no.edge.seq <- rep(0, B)
  Omega.list <- A.list <- Imputed.A.list <- list()
  for (b in 1:B) {
    impute.err.mat[b, ] <- result[[b]]$impute.sq.err
    block.err.mat[b, ] <- result[[b]]$block.sq.err
    loglike.mat[b, ] <- result[[b]]$loglike
    roc.auc.mat[b, ] <- result[[b]]$roc.auc
    bin.dev.mat[b, ] <- result[[b]]$bin.dev
    no.edge.seq[b] <- result[[b]]$no.edge
    dc.block.err.mat[b, ] <- result[[b]]$dc.block.sq.err
    dc.loglike.mat[b, ] <- result[[b]]$dc.loglike
    
    sbm.auc.mat[b, ] <- result[[b]]$sbm.auc
    dc.auc.mat[b, ] <- result[[b]]$dc.auc
  }
  output <- list(impute.err = colMeans(impute.err.mat), l2 = colMeans(block.err.mat), 
                 dev = colSums(loglike.mat), auc = colMeans(sbm.auc.mat), 
                 dc.l2 = colMeans(dc.block.err.mat), dc.dev = colSums(dc.loglike.mat),
                 dc.auc = colMeans(dc.auc.mat),
                 sse = colMeans(impute.err.mat), auc.mat = roc.auc.mat, 
                 dev.mat = loglike.mat, l2.mat = block.err.mat, SSE.mat = impute.err.mat, 
                 auc.mat = sbm.auc.mat,
                 dc.dev.mat = dc.loglike.mat, dc.l2.mat = dc.block.err.mat,
                 dc.auc.mat = dc.auc.mat)
  if (min(output$dev) > min(output$dc.dev)) {
    dev.model <- paste("DCSBM", which.min(output$dc.dev), 
                       sep = "-")
  }
  else {
    dev.model <- paste("SBM", which.min(output$dev), sep = "-")
  }
  if (min(output$l2) > min(output$dc.l2)) {
    l2.model <- paste("DCSBM", which.min(output$dc.l2), sep = "-")
  }
  else {
    l2.model <- paste("SBM", which.min(output$l2), sep = "-")
  }
  if (min(output$auc) > min(output$dc.auc)) {
    auc.model <- paste("DCSBM", which.min(output$dc.auc), sep = "-")
  }
  else {
    auc.model <- paste("SBM", which.min(output$auc), sep = "-")
  }
  output$l2.model <- l2.model
  output$dev.model <- dev.model
  output$auc.model <- auc.model
  return(output)
}

holdout.evaluation.fast.all <- function(holdout.index,A,max.K,soft=TRUE,tau=0,dc.est=1,fast=FALSE,p.sample=1,kappa=NULL){
  n <- nrow(A)
  edge.index <- which(upper.tri(A))
  edge.n <- length(edge.index)
  A.new <- matrix(0,n,n)
  A.new[upper.tri(A.new)] <- A[edge.index]
  A.new[edge.index[holdout.index]] <- NA
  A.new <- A.new + t(A.new)
  degrees <- colSums(A.new,na.rm=TRUE)
  no.edge <- 0
  no.edge <- sum(degrees==0)
  
  Omega <- which(is.na(A.new))
  non.miss <- which(!is.na(A.new))
  #A.new[non.miss] <- A.new[non.miss] + 0.5
  SVD.result <- iter.SVD.core.fast.all(A.new,max.K,fast=TRUE,p.sample=p.sample)
  dc.block.sq.err <-  dc.loglike <- roc.auc <- bin.dev <- block.sq.err <- impute.sq.err <- loglike <- rep(0,max.K)
  sbm.auc <- dc.auc <- rep(0,max.K)
  
  for(k in 1:max.K){
    #print(k)
    ##print(fast)
    tmp.est <- SVD.result[[k]]
    A.approx <- tmp.est$A.thr
    impute.sq.err[k] <- sum((A.approx[Omega]-A[Omega])^2)
    response <- A[edge.index[holdout.index]]#A[Omega]
    predictors <- A.approx[edge.index[holdout.index]]#A.approx[Omega]
    #print("AUC calculation")
    ##print(system.time(tmp.roc <- pROC::roc(response=response,predictor=predictors)))
    ##print(length(unique(predictors)))
    #aa <- AUC::roc(predictions=predictors,labels=factor(response))
    aa <- 0
    #tmp.roc.smooth <- smooth(tmp.roc,method="binormal")
    #roc.auc[k] <- auc(aa)#as.numeric(tmp.roc$auc)
    roc.auc[k] <- 0
    ##print(tmp.roc$auc)
    ##print(auc(aa))
    #roc.auc[k] <- as.numeric(tmp.roc.smooth$auc)
    trunc.predictors <- predictors
    trunc.predictors[predictors>(1-1e-6)] <- 1-1e-6
    trunc.predictors[predictors<1e-6] <- 1e-6
    bin.dev[k] <- sum((response-trunc.predictors)^2)#-sum(response*log(trunc.predictors)) - sum((1-response)*log(1-trunc.predictors))
    if(k==1){
      pb <- (sum(A.new,na.rm=TRUE)+1)/(sum(!is.na(A.new)) -sum(!is.na(diag(A.new)))+1)
      if(pb < 1e-6) pb <- 1e-6
      if(pb > 1-1e-6) pb <- 1-1e-6
      A.Omega <- A[Omega]
      block.sq.err[k] <- sum((pb-A[Omega])^2)
      loglike[k] <- -sum(A.Omega*log(pb)) - sum((1-A.Omega)*log(1-pb))
      
    }
    
    #U.approx <- eigen(A.approx)$vectors[,1:k]
    #print("SBM calculation")
    ##print(k)
    ##print(dim(tmp.est$SVD$v))
    ptm <- proc.time()
    if(k==1) {U.approx <- matrix(tmp.est$SVD$v,ncol=k)}else{
      U.approx <- tmp.est$SVD$v[,1:k]
      if(tau>0){
        A.approx <- A.approx + tau*mean(colSums(A.approx))/n
        d.approx <- colSums(A.approx)
        L.approx <- diag(1/sqrt(d.approx))%*%A.approx%*%diag(1/sqrt(d.approx))
        A.approx.svd <- irlba(L.approx,nu=k,nv=k)
        U.approx <- A.approx.svd$v[,1:k]     
      }       
    }
    
    km <- kmeans(U.approx,centers=k,nstart=30,iter.max=30)
    B <- matrix(0,k,k)
    Theta <- matrix(0,n,k)
    for(i in 1:k){
      for(j in i:k){
        N.i <- which(km$cluster==i)
        N.j <- which(km$cluster==j)
        if(i!=j){
          B[i,j] <- B[j,i] <- (sum(A.new[N.i,N.j],na.rm=TRUE)+1)/(sum(!is.na(A.new[N.i,N.j]))+1)
        } else{
          ##print(max(N.i))
          ##print(max(N.j))
          ##print(dim(A.new))
          B[i,j] <- B[j,i] <- (sum(A.new[N.i,N.j],na.rm=TRUE)+1)/(sum(!is.na(A.new[N.i,N.j])) -sum(!is.na(diag(A.new[N.i,N.j])))+1)
        }
        
      }
      Theta[N.i,i] <- 1
    }
    P.hat <- Theta%*%B%*%t(Theta)
    diag(P.hat) <- 0
    block.sq.err[k] <- sum((P.hat[Omega]-A[Omega])^2)
    P.hat.Omega <- P.hat[Omega]
    A.Omega <- A[Omega]
    P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
    P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6
    loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
    sbm.auc[k] <- AUC(A.Omega, P.hat.Omega) ##SC addition
    ##print(oc.time() - ptm)
    #### Degree correct model
    V <- U.approx
    #print("DCSBM calculation")
    ptm <- proc.time()
    #V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
    if(k==1) {V.norms <- as.numeric(abs(V))}else{
      V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
    }
    
    iso.index <- which(V.norms==0)
    Psi <- V.norms
    Psi <- Psi / max(V.norms)
    inv.V.norms <- 1/V.norms
    inv.V.norms[iso.index] <- 1
    
    V.normalized <- diag(as.numeric(inv.V.norms))%*%V
    
    #V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
    #Psi <- V.norms
    #Psi <- Psi / max(V.norms)
    #V.normalized <- diag(1/V.norms)%*%V
    #Psi.outer <- outer(Psi,Psi)
    if(k==1){
      if(dc.est>1){
        B <- sum(A.new,na.rm=TRUE)+0.01
        
        partial.d <- colSums(A.new,na.rm=TRUE)
        partial.gd <- B
        phi <- rep(0,n)
        B.g <- partial.gd
        phi <- as.numeric(partial.d/B.g)
        B <- B/p.sample
        P.hat <- t(t(matrix(B,n,n)*phi)*phi)
        #P.hat <- diag(phi)%*%matrix(B,n,n)%*%diag(phi)
        diag(P.hat) <- 0
      }
      dc.block.sq.err[k] <- sum((pb-A[Omega])^2)
      P.hat.Omega <- P.hat[Omega]
      A.Omega <- A[Omega]
      P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
      P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6
      
      dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
      dc.auc[k] <- AUC(A.Omega, P.hat.Omega)
      
      
    }else{
      km <- kmeans(V.normalized,centers=k,nstart=30,iter.max=30)
      if(dc.est>1){
        B <- matrix(0,k,k)
        Theta <- matrix(0,n,k)
        for(i in 1:k){
          for(j in 1:k){
            N.i <- which(km$cluster==i)
            N.j <- which(km$cluster==j)
            B[i,j] <- sum(A.new[N.i,N.j],na.rm=TRUE)+0.01
          }
          Theta[N.i,i] <- 1
        }
        Theta <- Matrix(Theta,sparse=TRUE)
        partial.d <- colSums(A.new,na.rm=TRUE)
        partial.gd <- colSums(B)
        phi <- rep(0,n)
        B.g <- Theta%*%partial.gd
        phi <- as.numeric(partial.d/B.g)
        B <- B/p.sample
        tmp.int.mat <- Theta*phi
        P.hat <-as.matrix(tmp.int.mat%*%B%*%t(tmp.int.mat))
        #P.hat <- diag(phi)%*%Theta%*%B%*%t(Theta)%*%diag(phi)
        diag(P.hat) <- 0
      }
      dc.block.sq.err[k] <- sum((P.hat[Omega]-A[Omega])^2)
      P.hat.Omega <- P.hat[Omega]
      A.Omega <- A[Omega]
      P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
      P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6
      dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
      dc.auc[k] <- AUC(A.Omega, P.hat.Omega)
    }
    ###print(oc.time() - ptm)
    
    
    
  }
  return(list(impute.sq.err=impute.sq.err,block.sq.err=block.sq.err,loglike=loglike,roc.auc=roc.auc,no.edge=no.edge,dc.block.sq.err=dc.block.sq.err,dc.loglike=dc.loglike,bin.dev=bin.dev,
              sbm.auc = sbm.auc, dc.auc = dc.auc))
}

iter.SVD.core.fast.all <- function(A,Kmax,tol=1e-5,max.iter=100,sparse=TRUE,init=NULL,verbose=FALSE,tau=0,fast=FALSE,p.sample=1){
  if(sparse) A <- Matrix(A,sparse=TRUE)
  avg.p <- mean(as.numeric(A),na.rm=TRUE)
  cap <- 1#kappa*avg.p
  A[which(is.na(A))] <- 0
  A <- A/p.sample
  #svd.new <- svd(A,nu=K,nv=K)
  ##print("begin SVD")
  svd.new <- irlba(A,nu=Kmax,nv=Kmax)
  ##print("end SVD")
  result <- list()
  for(K in 1:Kmax){
    #print(K)
    if(K==1){
      A.new <- svd.new$d[1]*matrix(svd.new$u[,1],ncol=1)%*%t(matrix(svd.new$v[,1],ncol=1))
    }else{
      A.new <- A.new + svd.new$d[K]*matrix(svd.new$u[,K],ncol=1)%*%t(matrix(svd.new$v[,K],ncol=1))
    }
    A.new.thr <- A.new
    A.new.thr[A.new < 0+tau] <- 0+tau
    A.new.thr[A.new >cap] <- cap
    
    tmp.SVD <- list(u=svd.new$u[,1:K],v=svd.new$v[,1:K],d=svd.new$d[1:K])
    result[[K]] <- list(iter=NA,SVD=tmp.SVD,A=A.new,err.seq=NA,A.thr=A.new.thr)
  }
  return(result)
}

################################################################################
################################################################################
##NCV with l2, bin.dev, and AUC
NCV.for.blockmodel <- function (A, max.K, cv = 3, dc.est = 1) 
{
  dc.avg.se <- dc.avg.log <- dc.avg.auc <- 
    avg.se <- avg.log <- avg.auc <- rep(0, max.K)
  dc.avg.se[1] <- dc.avg.log[1] <- dc.avg.auc[1] <- 
    avg.se[1] <- avg.log[1] <- avg.auc[1] <- Inf
  dc.dev.mat <- dc.l2.mat <- dc.auc.mat <- 
    sbm.dev.mat <- sbm.l2.mat <- sbm.auc.mat <- matrix(0, cv, max.K)
  #####start from here - cv.evaluate/dc updated
  n <- nrow(A)
  sample.index <- sample.int(n)
  max.fold.num <- ceiling(n/cv)
  fold.index <- rep(1:cv, each = max.fold.num)[1:n]
  cv.index <- fold.index[sample.index]
  for (KK in 1:max.K) {
    dc.l2 <- l2 <- dc.log.like <- log.like <- sbm.auc <- dc.auc <- rep(0, cv)
    for (k in 1:cv) {
      holdout.index <- which(cv.index == k)
      train.index <- which(cv.index != k)
      #tmp.eval <- cv.evaluate(A, train.index, holdout.index, 
      #                        KK)
      #tmp.eval.dc <- cv.evaluate.DC(A, train.index, holdout.index, 
      #                              KK)
      tmp.all <- cv.evaluate.all(A, train.index, holdout.index, KK, dc.est)
      
      #log.like[k] <- tmp.eval$loglike
      sbm.l2.mat[k, KK] <- l2[k] <- tmp.all$l2
      sbm.dev.mat[k, KK] <- log.like[k] <- tmp.all$loglike
      sbm.auc.mat[k, KK] <- sbm.auc[k] <- tmp.all$auc
      dc.l2.mat[k, KK] <- dc.l2[k] <- tmp.all$dc.l2
      dc.dev.mat[k, KK] <- dc.log.like[k] <- tmp.all$dc.loglike
      dc.auc.mat[k, KK] <- dc.auc[k] <- tmp.all$dc.auc
    }
    avg.se[KK] <- mean(l2)
    avg.log[KK] <- mean(log.like)
    avg.auc[KK] <- mean(sbm.auc)
    dc.avg.se[KK] <- mean(dc.l2)
    dc.avg.log[KK] <- mean(dc.log.like)
    dc.avg.auc[KK] <- mean(dc.auc)
  }
  if (min(avg.log) > min(dc.avg.log)) {
    dev.model <- paste("DCSBM", which.min(dc.avg.log), sep = "-")
  }else {
    dev.model <- paste("SBM", which.min(avg.log), sep = "-")
  }
  if (min(avg.se) > min(dc.avg.se)) {
    l2.model <- paste("DCSBM", which.min(dc.avg.se), sep = "-")
  }else {
    l2.model <- paste("SBM", which.min(avg.se), sep = "-")
  }
  if (min(avg.auc) > min(dc.avg.auc)) {
    auc.model <- paste("DCSBM", which.min(dc.avg.auc), sep = "-")
  }else {
    auc.model <- paste("SBM", which.min(avg.auc), sep = "-")
  }
  return(list(dev = avg.log, l2 = avg.se, auc = avg.auc,
              dc.dev = dc.avg.log, dc.l2 = dc.avg.se, dc.auc = dc.avg.auc,
              sbm.l2.mat = sbm.l2.mat, sbm.dev.mat = sbm.dev.mat, sbm.auc.mat = sbm.auc.mat,
              dc.l2.mat = dc.l2.mat, dc.dev.mat = dc.dev.mat, dc.auc.mat = dc.auc.mat,
              l2.model = l2.model, 
              dev.model = dev.model,
              auc.model = auc.model))
}

cv.evaluate.all <- function(A,train.index,holdout.index,K,dc.est=1){
  n <- nrow(A)
  A.new <- A[c(train.index,holdout.index),c(train.index,holdout.index)]
  n.holdout <- length(holdout.index)
  n.train <- n-n.holdout
  A1 <- A.new[1:n.train,]
  A1.svd <- irlba(A1+0.001,nu=K,nv=K)
  dc.A1.svd <- irlba(A1,nu=K,nv=K)
  
  V <- dc.A1.svd$v[,1:K]
  if(K==1) {V.norms <- abs(V)}else{
    #V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
    V.norms <- sqrt(rowSums(V^2))
  }
  iso.index <- which(V.norms==0)
  Psi <- V.norms
  Psi <- Psi / max(V.norms)
  #Psi.outer <- outer(Psi,Psi)
  Psi.outer <- tcrossprod(Psi)
  inv.V.norms <- 1/V.norms
  inv.V.norms[iso.index] <- 1
  V.normalized <- crossprod(diag(inv.V.norms), V)
  
  if(K==1){
    A0 <- A1[1:n.train,1:n.train]
    pb <- sum(A0)/n.train^2
    if(pb < 1e-6) pb <- 1e-6
    if(pb > 1- 1e-6) pb <- 1-1e-6
    A.2 <- A.new[(n.train+1):n,(n.train+1):n]
    sum.index <- lower.tri(A.2)
    auc <- AUC(A.2[sum.index], rep(pb, length(A.2[sum.index])))
    loglike <- -sum(A.2[sum.index]*log(pb)) - sum((1-A.2[sum.index])*log(1-pb))
    l2 <- sum((A.2[sum.index]-pb)^2)
    
    N.1i <- 1:n.train
    N.2i <- (n.train+1):n
    dc.pb <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(sum(Psi.outer[N.1i,N.1i])/2 + sum(Psi.outer[N.1i,N.2i]) - sum(diag(Psi.outer))+1)
    #dc.P.hat.holdout <-  diag(Psi[(n.train+1):n])%*%matrix(1,ncol=(n-n.train),nrow=(n-n.train))%*%diag(Psi[(n.train+1):n])*dc.pb
    dc.P.hat.holdout <-  tcrossprod(crossprod(diag(Psi[(n.train+1):n]),matrix(1,ncol=(n-n.train),nrow=(n-n.train))),diag(Psi[(n.train+1):n]))*dc.pb
    dc.P.hat.holdout[dc.P.hat.holdout<1e-6] <- 1e-6
    dc.P.hat.holdout[dc.P.hat.holdout>(1-1e-6)] <- 1-1e-6
    #A.2 <- A.new[(n.train+1):n,(n.train+1):n]
    #sum.index <- lower.tri(A.2)
    dc.loglike <- -sum(A.2[sum.index]*log(dc.P.hat.holdout[sum.index])) - sum((1-A.2[sum.index])*log(1-dc.P.hat.holdout[sum.index]))
    dc.auc <- AUC(A.2[sum.index], dc.P.hat.holdout[sum.index])
    dc.l2 <- sum((A.2[sum.index]-dc.P.hat.holdout[sum.index])^2)
    return(list(loglike=loglike,l2=l2, auc = auc, 
                dc.loglike=dc.loglike,dc.l2=dc.l2,dc.auc=dc.auc,
                no.edge=NA,impute.err=NA))
  } 
  
  V <- A1.svd$v
  km <- kmeans(V,centers=K,nstart=30,iter.max=30)
  
  dc.km <- kmeans(V.normalized,centers=K,nstart=30,iter.max=30)
  
  degrees <- colSums(A1)
  no.edge <- sum(degrees==0)
  
  
  B <- dc.B <- matrix(0,K,K)
  
  tmp <- lapply(1:K, function(ii){
    N1 <- intersect(1:n.train,which(km$cluster==ii))
    N2 <- intersect((n.train+1):n,which(km$cluster==ii))
    
    dc.N1 <- intersect(1:n.train,which(dc.km$cluster==ii))
    dc.N2 <- intersect((n.train+1):n,which(dc.km$cluster==ii))
    
    return(list(N1 = N1, N2 = N2, 
                dc.N1 = dc.N1, dc.N2 = dc.N2))
  })
  
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      # N.1i <- intersect(1:n.train,which(km$cluster==i))
      # N.2i <- intersect((n.train+1):n,which(km$cluster==i))
      # N.1j <- intersect(1:n.train,which(km$cluster==j))
      # N.2j <- intersect((n.train+1):n,which(km$cluster==j))
      B[i,j] <- B[j,i] <- (
        sum(A.new[tmp[[i]]$N1,tmp[[j]]$N1]) + 
          sum(A.new[tmp[[i]]$N1,tmp[[j]]$N2]) + 
          sum(A.new[tmp[[j]]$N1,tmp[[i]]$N2])+1
      )/(
        length(tmp[[i]]$N1)*length(tmp[[j]]$N1)+
          length(tmp[[j]]$N1)*length(tmp[[i]]$N2)+
          length(tmp[[i]]$N1)*length(tmp[[j]]$N2)+1
      )
      #B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j])+1)/(length(N.1i)*length(N.1j)+length(N.1i)*length(N.2j)+1)
      
      # dc.N.1i <- intersect(1:n.train,which(dc.km$cluster==i))
      # dc.N.2i <- intersect((n.train+1):n,which(dc.km$cluster==i))
      # dc.N.1j <- intersect(1:n.train,which(dc.km$cluster==j))
      # dc.N.2j <- intersect((n.train+1):n,which(dc.km$cluster==j))
      dc.B[i,j] <- dc.B[j,i] <- (
        sum(A.new[tmp[[i]]$dc.N1,tmp[[j]]$dc.N1]) + 
          sum(A.new[tmp[[i]]$dc.N1,tmp[[j]]$dc.N2]) + 
          sum(A.new[tmp[[j]]$dc.N1,tmp[[i]]$dc.N2])+1)/(
            sum(Psi.outer[tmp[[i]]$dc.N1,tmp[[j]]$dc.N1]) + 
              sum(Psi.outer[tmp[[i]]$dc.N1,tmp[[j]]$dc.N2]) + 
              sum(Psi.outer[tmp[[j]]$dc.N1,tmp[[i]]$dc.N2])+1)
    }
  }
  #B <- B+t(B)
  Theta <- matrix(0,n,K)
  dc.Theta <- matrix(0,n,K)
  for(i in 1:K){
    # N.1i <- intersect(1:n.train,which(km$cluster==i))
    # N.2i <- intersect((n.train+1):n,which(km$cluster==i))
    B[i,i] <- (
      sum(A.new[tmp[[i]]$N1,tmp[[i]]$N1])/2 + 
        sum(A.new[tmp[[i]]$N1,tmp[[i]]$N2])+1)/(
          length(tmp[[i]]$N1)*(length(tmp[[i]]$N1)-1)/2+
            length(tmp[[i]]$N1)*length(tmp[[i]]$N2)+1)
    Theta[which(km$cluster==i),i] <- 1
    
    dc.B[i,i] <- (
      sum(A.new[tmp[[i]]$dc.N1,tmp[[i]]$dc.N1])/2 + 
        sum(A.new[tmp[[i]]$dc.N1,tmp[[i]]$dc.N2])+1)/(
          sum(Psi.outer[tmp[[i]]$dc.N1,tmp[[i]]$dc.N1])/2 + 
            sum(Psi.outer[tmp[[i]]$dc.N1,tmp[[i]]$dc.N2]) - 
            sum(diag(Psi.outer))+1)
    dc.Theta[which(dc.km$cluster==i),i] <- 1
    
  }
  
  #P.hat.holdout <- Theta[(n.train+1):n,]%*%B%*%t(Theta[(n.train+1):n,])
  P.hat.holdout <- tcrossprod(tcrossprod(Theta[(n.train+1):n,],B),Theta[(n.train+1):n,])
  P.hat.holdout[P.hat.holdout<1e-6] <- 1e-6
  P.hat.holdout[P.hat.holdout>(1-1e-6)] <- 1-1e-6
  A.2 <- A.new[(n.train+1):n,(n.train+1):n]
  sum.index <- lower.tri(A.2)
  loglike <- -sum(A.2[sum.index]*log(P.hat.holdout[sum.index])) - sum((1-A.2[sum.index])*log(1-P.hat.holdout[sum.index]))
  
  auc <- AUC(A.2[sum.index], P.hat.holdout[sum.index])
  
  l2 <- sum((A.2[sum.index]-P.hat.holdout[sum.index])^2)
  
  ## dc
  tmp.imt.mat <- dc.Theta[(n.train+1):n,]*Psi[(n.train+1):n]
  #dc.P.hat.holdout <-  tmp.imt.mat%*%dc.B%*%t(tmp.imt.mat)
  dc.P.hat.holdout <-  tcrossprod(tcrossprod(tmp.imt.mat,dc.B),tmp.imt.mat)
  #P.hat.holdout <-  diag(Psi[(n.train+1):n])%*%Theta[(n.train+1):n,]%*%B%*%t(Theta[(n.train+1):n,])%*%diag(Psi[(n.train+1):n])
  if(dc.est==2){
    dc.B <- matrix(0,K,K)
    dc.Theta <- matrix(0,n,K)
    A.new.na <- A.new
    A.new.na[(n.train+1):n,(n.train+1):n] <- NA
    for(i in 1:K){
      for(j in 1:K){
        N.i <- which(km$cluster==i)
        N.j <- which(km$cluster==j)
        dc.B[i,j] <- sum(A.new.na[N.i,N.j],na.rm=TRUE)+0.01
      }
      dc.Theta[N.i,i] <- 1
    }
    partial.d <- colSums(A.new.na,na.rm=TRUE)
    partial.gd <- colSums(dc.B)
    phi <- rep(0,n)
    B.g <- dc.Theta%*%partial.gd
    phi <- as.numeric(partial.d/B.g)
    #P.hat <- diag(phi)%*%dc.Theta%*%dc.B%*%t(dc.Theta)%*%diag(phi)
    P.hat <- tcrossprod(crossprod(diag(phi),tcrossprod(tcrossprod(dc.Theta,dc.B),dc.Theta)),diag(phi))
    diag(P.hat) <- 0
    dc.P.hat.holdout <- P.hat[(n.train+1):n,(n.train+1):n]
  }
  dc.P.hat.holdout[dc.P.hat.holdout<1e-6] <- 1e-6
  dc.P.hat.holdout[dc.P.hat.holdout>(1-1e-6)] <- 1-1e-6
  #A.2 <- A.new[(n.train+1):n,(n.train+1):n]
  #sum.index <- lower.tri(A.2)
  dc.loglike <- -sum(A.2[sum.index]*log(dc.P.hat.holdout[sum.index])) - sum((1-A.2[sum.index])*log(1-dc.P.hat.holdout[sum.index]))
  dc.auc <- AUC(A.2[sum.index], dc.P.hat.holdout[sum.index])
  dc.l2 <- sum((A.2[sum.index]-dc.P.hat.holdout[sum.index])^2)
  
  return(list(loglike=loglike,l2=l2, auc = auc,
              dc.loglike=dc.loglike,dc.l2=dc.l2,dc.auc=dc.auc,
              no.edge=no.edge,impute.err=NA))
}


################################################################################
################################################################################
## RDPG generation
sparse.RDPG.gen <- function(n, d, sparsity.multiplier = 1, ncore = 1){
  X <- matrix(runif(n*d), nrow = n, ncol = d)
  
  # X.norm <- X/sqrt(rowSums(X^2))
  
  P <- tcrossprod(X)
  P <- P*sparsity.multiplier/max(P)
  diag(P) <- 0
  
  stor <- do.call('rbind',
                  mclapply(1:(n-1), function(i) {
                    tmp <- which(rbinom(n-i, 1, P[i,(i+1):n]) == 1)
                    
                    if(length(tmp) == 0)
                      return(NULL)
                    else
                      return(cbind(rep(i, length(tmp)), i + tmp))
                  }, mc.cores = ncore))
  
  A <- sparseMatrix(i = stor[,1], j = stor[,2], dims = c(n,n), symmetric = T)
  
  return(list('A' = A, 'P' = P))
}


################################################################################
################################################################################
## Croissant for RDPG
croissant.rdpg <- function(A, d.cand, s, o, R,
                           laplace = F,
                           loss = c("l2", "bin.dev", "AUC"),
                           ncore = 1){
  n <- nrow(A)
  m <- (n-o)/s
  
  if(length(d.cand) == 1) d.cand <- 1:d.cand
  
  dmax <- max(d.cand)
  
  over <- lapply(1:R, function(ii) sample.int(n, o, F))
  non.over <- lapply(1:R, function(ii) sample((1:n)[-over[[ii]]], n-o, 
                                              replace = F))
  
  raw.ind <- cbind(rep(1:R, each = s), rep(1:s, R))
  colnames(raw.ind) <- c('r', 's')
  
  raw.out <- mclapply(1:nrow(raw.ind), function(ii){
    
    q <- raw.ind[ii, 's']
    r <- raw.ind[ii, 'r']
    
    sonn <- c(over[[r]], non.over[[r]][((q-1)*m+1):(q*m)])
    A.sonn <- A[sonn, sonn]
    
    if(!laplace){
      L.sonn <- A.sonn
    }else{
      degree <- rowSums(A.sonn)
      D <- sparseMatrix(i = 1:(o+m), j = 1:(o+m), x = 1/sqrt(degree))
      
      L.sonn <- tcrossprod(crossprod(D, A.sonn), D)
    }
    
    eig <- irlba::irlba(L.sonn, nv = dmax)
    U <- eig$v
    sigma.half <- diag(sqrt(abs(eig$d)))
    
    X.sonn <- tcrossprod(U, sigma.half)
    
    return(X.sonn)
    
  },mc.cores = ncore)
  
  match.out <- mclapply(1:nrow(raw.ind), function(ii){
    
    r <- raw.ind[ii, 'r']
    q <- raw.ind[ii, 's']
    #dd <- match.ind[ii, 'dd']
    
    #iip <- which(raw.ind[,'r'] == r & raw.ind[,'s'] == q)
    
    if(q == 1) return(raw.out[[ii]][-(1:o),])
    
    stand <- which(raw.ind[,'r'] == r & raw.ind[,'s'] == 1)
    
    proc.mat <- Procrustes(raw.out[[ii]][(1:o), ],
                           raw.out[[stand]][(1:o), ])$R
    
    X.rot <- raw.out[[ii]][-(1:o),] %*% proc.mat
    
    return(X.rot)
    
  }, mc.cores = ncore)
  
  non.size <- s*(s-1)/2
  ld <- length(d.cand)
  non.mat <- matrix(nrow = R*non.size*ld, ncol = 4)
  cc <- 1
  for(r in 1:R)
    for(dd in seq_along(d.cand))
      for(p in 1:(s-1))
        for(q in (p+1):s){
          non.mat[cc, ] <- c(r, dd, p, q)
          cc <- cc + 1
        }
  colnames(non.mat) <- c('r', 'dd', 'p', 'q')
  
  L.all <- mclapply(1:nrow(non.mat), function(ii){
    r <- non.mat[ii, 'r']
    dd <- non.mat[ii, 'dd']
    p <- non.mat[ii, 'p']
    q <- non.mat[ii, 'q']
    
    p.non <- non.over[[r]][((p-1)*m+1):(p*m)]
    q.non <- non.over[[r]][((q-1)*m+1):(q*m)]
    
    A.non <- A[p.non, q.non]
    
    L.temp <- matrix(0, nrow = length(loss), ncol = 1)
    row.names(L.temp) <- loss
    colnames(L.temp) <- as.character(d.cand[dd])
    
    ind1 <- which(raw.ind[, 'r'] == r & raw.ind[, 's'] == p)
    ind2 <- which(raw.ind[, 'r'] == r & raw.ind[, 's'] == q)
    
    P.hat <- tcrossprod(match.out[[ind1]][,1:d.cand[dd]],
                        match.out[[ind2]][,1:d.cand[dd]])
    
    message(sum(P.hat < 0 | P.hat > 1))
    
    P.hat[P.hat < 1e-6] <- 1e-6
    P.hat[P.hat > 1- 1e-6] <- 1 - 1e-6
    
    for(lq in seq_along(loss)){
      tmp.nm <- loss[lq]
      L.temp[tmp.nm, 1] <-  L.temp[tmp.nm, 1] +
        (do.call(loss[lq], list(A.non, P.hat)))/(s*(s-1)*0.5)
    }
    
    return(L.temp)},
    mc.cores = ncore)
  
  L <- list()
  
  for(r in 1:R){
    L[[r]] <- do.call('cbind', lapply(1:ld, function(dd){
      Reduce('+', L.all[which(non.mat[,'dd'] == dd & non.mat[,'r'] == r)])
    }))
    
    row.names(L[[r]]) <- loss
    colnames(L[[r]]) <- as.character(d.cand)
  }
  
  obj <- data.table(`Candidate Rank` = d.cand)
  
  for(lq in seq_along(loss))
    for(r in 1:R){
      obj[[paste0(loss[lq], "-Rep=", r)]] <- 
        as(rbind(L[[r]][loss[lq], ]), "vector")
    }
  
  obj2 <- list()
  
  obj2[["Candidate Rank"]] <- d.cand
  
  for(lq in seq_along(loss)){
    
    obj2[[paste0("d.hat.each.rep (", loss[lq], ")")]] <- sapply(1:R, function(r){
      l.rdpg <- d.cand[which.min(L[[r]][loss[lq],])]
    })
    
    
    obj2[[paste0(loss[lq], ".model")]] <- 
      modal(obj2[[paste0("d.hat.each.rep (", loss[lq], ")")]])
  }
  
  return(c(list('loss' = obj), 
           obj2))
  
}

##ECV for RDPG
missing.undirected.Rank.weighted.fast.all <- function(holdout.index,A,max.K,soft=FALSE,fast=fast,p.sample=1){
  n <- nrow(A)
  #A.new <- A
  #A.new[holdout.index] <- NA
  edge.index <- which(upper.tri(A))
  edge.n <- length(edge.index)
  A.new <- matrix(0,n,n)
  A.new[upper.tri(A.new)] <- A[edge.index]
  A.new[edge.index[holdout.index]] <- NA
  A.new <- A.new + t(A.new)
  diag(A.new) <- diag(A)
  degrees <- colSums(A.new,na.rm=TRUE)
  no.edge <- 0
  no.edge <- sum(degrees==0)
  
  Omega <- which(is.na(A.new))
  imputed.A <- list()
  sse <- roc.auc <- dev <- rep(0,max.K)
  SVD.result <- iter.SVD.core.fast.all(A.new,max.K,fast=TRUE,p.sample=p.sample)
  for(k in 1:max.K){
    print(k)
    tmp.est <- SVD.result[[k]]
    #if(k==1){
    #A.approx <- matrix(tmp.est$SVD$u,ncol=1)%*%t(matrix(tmp.est$SVD$v,ncol=1))*tmp.est$SVD$d[1]
    #}else{
    #   A.approx <- tmp.est$SVD$u%*%t(tmp.est$SVD$v*tmp.est$SVD$d)
    #}
    A.approx <- tmp.est$A
    response <- A[Omega]
    predictors <- A.approx[Omega]
    #aa <- AUC::roc(predictions=predictors,labels=factor(response))
    roc.auc[k] <- AUC(response, predictors)
    sse[k] <- mean((response-predictors)^2)
    predictors[predictors < 1e-6] <- 1e-6
    predictors[predictors > 1-1e-6] <- 1-1e-6
    dev[k] <- bin.dev(matrix(response, ncol = k), matrix(predictors, ncol = k))
    imputed.A[[k]] <- A.approx
  }
  return(list(imputed.A=imputed.A,Omega=Omega, roc.auc = roc.auc, sse=sse, dev = dev))
}

ECV.undirected.Rank <- function(A,max.K,B=3,holdout.p=0.1,soft=FALSE,fast=fast){
  n <- nrow(A)
  #edge.index <- 1:n^2
  #edge.n <- length(edge.index)
  edge.index <- which(upper.tri(A))
  edge.n <- length(edge.index)
  
  holdout.index.list <- list()
  
  holdout.n <- floor(holdout.p*edge.n)
  
  for(j in 1:B){
    holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
  }
  result <- lapply(holdout.index.list,
                   missing.undirected.Rank.weighted.fast.all,
                   A=A,max.K=max.K,soft=soft,fast=fast,p.sample=1-holdout.p)
  
  sse.mat <- roc.auc.mat <- dev.mat <- matrix(0,nrow=B,ncol=max.K)
  
  for(b in 1:B){
    roc.auc.mat[b,] <- result[[b]]$roc.auc
    sse.mat[b,] <- result[[b]]$sse
    dev.mat[b,] <- result[[b]]$dev
  }
  
  auc.seq <- colMeans(roc.auc.mat)
  #auc.sd <- apply(roc.auc.mat,2,sd)/sqrt(B)
  sse.seq <- colMeans(sse.mat)
  #sse.sd <- apply(sse.mat,2,sd)/sqrt(B)
  dev.seq <- colMeans(dev.mat)
  # return(list(sse=sse.seq,sse.sd=sse.sd))
  return(list(rank.sse=which.min(sse.seq), sse=sse.seq,
              rank.dev = which.min(dev.seq), dev = dev.seq,
              rank.auc=which.min(auc.seq),auc=auc.seq
  ))
}


## Latent space
latent.gen <- function(n, d, alpha = 1, sparsity = 1, ncore = 1){
  Z <- matrix(runif(n*d), nrow = n, ncol = d)
  
  stor <- do.call('rbind',
                  mclapply(1:(n-1), function(i) {
                    
                    logodds <- alpha - sapply((i+1):n, 
                                              \(jj) sqrt(sum((Z[i, ] -Z[jj, ])^2))
                    )
                    
                    pp <- sparsity*exp(logodds)/(1+exp(logodds))
                    
                    tmp <- which(rbinom(n-i, 1, pp) == 1)
                    
                    if(length(tmp) == 0)
                      return(NULL)
                    else
                      return(cbind(rep(i, length(tmp)), i + tmp))
                  }, mc.cores = ncore))
  
  A <- as(sparseMatrix(i = stor[,1], j = stor[,2], dims = c(n,n), symmetric = T),
          'dMatrix')
  
  return(list('A' = A, 'Z' = Z))
}


## Croissant for latent space
croissant.latent <- function(A, d.cand, s, o, R,
                             loss = c("l2", "bin.dev", "AUC"),
                             ncore = 1){
  n <- nrow(A)
  m <- (n-o)/s
  
  if(length(d.cand) == 1) d.cand <- 1:d.cand
  
  dmax <- max(d.cand)
  
  over <- lapply(1:R, function(ii) sample.int(n, o, F))
  non.over <- lapply(1:R, function(ii) sample((1:n)[-over[[ii]]], n-o, 
                                              replace = F))
  
  # raw.ind <- cbind(rep(1:R, each = s), rep(1:s, R))
  # colnames(raw.ind) <- c('r', 's')
  
  ld <- length(d.cand)
  raw.ind <- matrix(nrow = R*s*ld, ncol = 3)
  cc <- 1
  for(r in 1:R)
    for(dd in seq_along(d.cand))
      for(q in 1:s){
        raw.ind[cc, ] <- c(r, dd, q)
        cc <- cc + 1
      }
  colnames(raw.ind) <- c('r', 'dd', 'q')
  
  system.time(raw.out <- mclapply(1:nrow(raw.ind), function(ii){
    
    q <- raw.ind[ii, 'q']
    r <- raw.ind[ii, 'r']
    dd <- raw.ind[ii, 'dd']
    
    sonn <- c(over[[r]], non.over[[r]][((q-1)*m+1):(q*m)])
    A.sonn <- A[sonn, sonn]
    net.sonn <- as.network(A.sonn, matrix.type = "adjacency")
    
    out.lat <- ergmm(net.sonn ~ euclidean(d = d.cand[dd]), tofit = "mle")
    Z.sonn <- out.lat$mle$Z
    beta.sonn <- out.lat$mle$beta
    
    return(list(Z.hat = Z.sonn, beta.hat = beta.sonn))
    
  },mc.cores = ncore))
  
  match.out <- mclapply(1:nrow(raw.ind), function(ii){
    
    r <- raw.ind[ii, 'r']
    q <- raw.ind[ii, 'q']
    dd <- raw.ind[ii, 'dd']
    #dd <- match.ind[ii, 'dd']
    
    #iip <- which(raw.ind[,'r'] == r & raw.ind[,'s'] == q)
    
    if(q == 1) return(list('Z.rot' = raw.out[[ii]]$Z.hat[-(1:o),],
                           'beta.hat' = raw.out[[ii]]$beta.hat))
    
    stand <- which(raw.ind[,'r'] == r & raw.ind[,'q'] == 1 & 
                     raw.ind[,'dd'] == dd)
    
    proc.par <- Procrustes(cbind(raw.out[[ii]]$Z.hat[(1:o), ]),
                           cbind(raw.out[[stand]]$Z.hat[(1:o), ]),
                           translate = T,
                           dilate = F)
    
    Z.rot <- cbind(raw.out[[ii]]$Z.hat[-(1:o),]) %*% proc.par$R +
      matrix(proc.par$t, nrow = m, ncol = d.cand[dd])
    
    return(list('Z.rot' = Z.rot, 'beta.hat' = raw.out[[ii]]$beta.hat))
    
  }, mc.cores = ncore)
  
  non.size <- s*(s-1)/2
  ld <- length(d.cand)
  non.mat <- matrix(nrow = R*non.size*ld, ncol = 4)
  cc <- 1
  for(r in 1:R)
    for(dd in seq_along(d.cand))
      for(p in 1:(s-1))
        for(q in (p+1):s){
          non.mat[cc, ] <- c(r, dd, p, q)
          cc <- cc + 1
        }
  colnames(non.mat) <- c('r', 'dd', 'p', 'q')
  
  L.all <- mclapply(1:nrow(non.mat), function(ii){
    r <- non.mat[ii, 'r']
    dd <- non.mat[ii, 'dd']
    p <- non.mat[ii, 'p']
    q <- non.mat[ii, 'q']
    
    p.non <- non.over[[r]][((p-1)*m+1):(p*m)]
    q.non <- non.over[[r]][((q-1)*m+1):(q*m)]
    
    A.non <- A[p.non, q.non]
    
    L.temp <- matrix(0, nrow = length(loss), ncol = 1)
    row.names(L.temp) <- loss
    colnames(L.temp) <- as.character(d.cand[dd])
    
    ind1 <- which(raw.ind[, 'r'] == r & raw.ind[, 'q'] == p &
                    raw.ind[,'dd'] == dd)
    ind2 <- which(raw.ind[, 'r'] == r & raw.ind[, 'q'] == q &
                    raw.ind[,'dd'] == d.cand[dd])
    
    # P.hat <- tcrossprod(match.out[[ind1]][,1:d.cand[dd]],
    #                     match.out[[ind2]][,1:d.cand[dd]])
    
    Z1.hat <- match.out[[ind1]]$Z.rot
    Z2.hat <- match.out[[ind2]]$Z.rot
    beta.hat <- (match.out[[ind1]]$beta.hat + match.out[[ind2]]$beta.hat)/2
    
    log.hat <- beta.hat - cdist(Z1.hat, Z2.hat)
    
    P.hat <- exp(log.hat)/(1+exp(log.hat))
    
    message(sum(P.hat < 0 | P.hat > 1))
    
    # P.hat[P.hat < 1e-6] <- 1e-6
    # P.hat[P.hat > 1- 1e-6] <- 1 - 1e-6
    
    for(lq in seq_along(loss)){
      tmp.nm <- loss[lq]
      L.temp[tmp.nm, 1] <-  L.temp[tmp.nm, 1] +
        (do.call(loss[lq], list(A.non, P.hat)))/(s*(s-1)*0.5)
    }
    
    return(L.temp)},
    mc.cores = ncore)
  
  L <- list()
  
  for(r in 1:R){
    L[[r]] <- do.call('cbind', lapply(1:ld, function(dd){
      Reduce('+', L.all[which(non.mat[,'dd'] == dd & non.mat[,'r'] == r)])
    }))
    
    row.names(L[[r]]) <- loss
    colnames(L[[r]]) <- as.character(d.cand)
  }
  
  obj <- data.table(`Candidate Rank` = d.cand)
  
  for(lq in seq_along(loss))
    for(r in 1:R){
      obj[[paste0(loss[lq], "-Rep=", r)]] <- 
        as(rbind(L[[r]][loss[lq], ]), "vector")
    }
  
  obj2 <- list()
  
  obj2[["Candidate Rank"]] <- d.cand
  
  for(lq in seq_along(loss)){
    
    obj2[[paste0("d.hat.each.rep (", loss[lq], ")")]] <- sapply(1:R, function(r){
      l.rdpg <- d.cand[which.min(L[[r]][loss[lq],])]
    })
    
    
    obj2[[paste0(loss[lq], ".model")]] <- 
      modal(obj2[[paste0("d.hat.each.rep (", loss[lq], ")")]])
  }
  
  return(c(list('loss' = obj), 
           obj2))
  
}

## Param tune reg SP
croissant.tune.regsp <- function(A, K, tau.cand,
                                 DCBM = F,
                                 s, o, R,
                                 laplace = F,
                                 dc.est = 2,
                                 loss = c("l2", "bin.dev", "AUC"),
                                 ncore = 1){
  n <- nrow(A)
  m <- (n-o)/s
  
  L <- list()
  
  over <- lapply(1:R, function(ii) sample.int(n, o, F))
  non.over <- lapply(1:R, function(ii) sample((1:n)[-over[[ii]]], n-o, replace = F))
  
  raw.ind <- cbind(rep(1:R, each = s), rep(1:s, R))
  
  raw.out <- mclapply(1:nrow(raw.ind), function(ii){
    
    q <- raw.ind[ii, 2]
    r <- raw.ind[ii, 1]
    
    sonn <- c(over[[r]], non.over[[r]][((q-1)*m+1):(q*m)])
    A.sonn <- A[sonn, sonn]
    
    deg <- rowSums(A.sonn)
    avg.deg <- mean(deg)
    
    out.BM <- list()
    
    for(tt in seq_along(tau.cand)){
      A.sonn.tau <- A.sonn + tau.cand[tt]*avg.deg/n
      
      d.sonn.tau <- sparseMatrix( i = 1:(o+m), j = 1:(o+m),
                                  x = 1/sqrt(deg + tau.cand[tt]*avg.deg))
      
      L.sonn <- A.sonn.tau
      
      if(laplace){
        L.sonn <- tcrossprod(crossprod(d.sonn.tau, A.sonn.tau), d.sonn.tau)
        L.sonn[is.na(L.sonn)] <- 0
      }
      
      eig.max <- irlba::partial_eigen(x = L.sonn, n = K, 
                                      symmetric = T)$vectors
      
      if(K == 1){
        out.BM[[tt]] <- rep(1, o+m)
        next
      }
      
      rn.eig <- eig.max
      
      if(DCBM){
        rownorm <- sqrt(rowSums(eig.max^2))
        rownorm[rownorm == 0] <- 1
        
        rn.eig <- eig.max/rownorm
      }
      
      out.BM[[tt]] <- as.integer(pam(rn.eig, K,
                                     metric = "euclidean",
                                     do.swap = F, cluster.only = T,
                                     pamonce = 6))
    }
    
    return(list('BM' = out.BM))
    # 'psi' = psi.hat))
  },
  mc.cores = ncore)
  
  tau.size <- length(tau.cand)
  
  est.out <- mclapply(1:(tau.size*nrow(raw.ind)), function(ii){
    
    tt <- ii %% tau.size
    tt <- ifelse(tt == 0, tau.size, tt)
    
    rot <- ceiling(ii / tau.size)
    
    q <- raw.ind[rot, 2]
    r <- raw.ind[rot, 1]
    
    #message(paste0("Est. at s=",q, " started"))
    sonn <- c(over[[r]], non.over[[r]][((q-1)*m+1):(q*m)])
    A.sonn <- A[sonn, sonn]
    
    out.BM.std <- raw.out[raw.ind[,1] == r][[1]]$BM[[tt]]
    
    out.BM <- raw.out[raw.ind[,1] == r][[q]]$BM[[tt]]
    
    work.tau <- tau.cand[tt]
    
    if(K == 1){
      mat.BM <- rep(1, m)
      
      if(!DCBM){
        B.BM <- fast.SBM.est(A.sonn, rep(1,o+m), o+m, 1)
        mat.BM <- rep(1, m)
        psi.BM <- rep(1, m)
        
        return(list('gBM' = mat.BM, 'BBM' = B.BM,
                    'psiBM' = psi.BM))
      }
      
      if(DCBM){
        if(dc.est == 2){
          tmp <- fast.DCBM.est(A.sonn, rep(1,o+m), o+m, 1, o,
                               p.sample = 1)
        }else{
          #psi.hat <- raw.out[raw.ind[,1] == r][[q]]$psi[[k.cand]]
          tmp <- eigen.DCBM.est(A.sonn, rep(1,o+m), o+m, 1, o,
                                p.sample = 1)
        }
      }
      
      B.BM <- tmp$Bsum
      psi.BM <- tmp$psi
      mat.BM <- rep(1, m)
      
      return(list('gBM' = mat.BM,
                  'BBM' = B.BM,
                  'psiBM' = psi.BM))
    }
    
    E.BM.kc <- best.perm.label.match(out.BM[1:o], 
                                     out.BM.std[1:o],
                                     o, K)
    
    tmp.BM <- sparseMatrix(i = 1:(o+m), j = out.BM, 
                           dims = c((o+m), K))
    
    mat.BM <- as.vector(tcrossprod(tcrossprod(tmp.BM, E.BM.kc),
                                   rbind(1:K)))
    if(!DCBM){
      B.BM <- fast.SBM.est(A.sonn, mat.BM, o+m, K)
      mat.BM <- mat.BM[-(1:o)]
      psi.BM <- rep(1, m)
      
      return(list('gBM' = mat.BM, 'BBM' = B.BM,
                  'psiBM' = psi.BM))
    }
    
    if(dc.est == 2){
      tmp <- fast.DCBM.est(A.sonn, mat.BM, o+m, K, o, 
                           p.sample = 1)
    }else{
      #psi.hat <- raw.out[raw.ind[,1] == r][[q]]$psi[[k.cand]]
      tmp <- eigen.DCBM.est(A.sonn, mat.BM, o+m, K, o, 
                            p.sample = 1)
    }
    
    B.BM <- tmp$Bsum
    psi.BM <- tmp$psi
    mat.BM <- mat.BM[-(1:o)]
    
    message(paste0("Est. at s=",q, " finished"))
    
    return(list('gBM' = mat.BM,
                'BBM' = B.BM,
                'psiBM' = psi.BM))
  },
  mc.cores = ncore)
  
  g.BM <- list()
  B.BM <- list()
  psi.BM <- list()
  
  raw.mat <- cbind(raw.ind[rep(1:nrow(raw.ind), each = tau.size), ], 
                   rep(1:tau.size, nrow(raw.ind)))
  
  for(r in 1:R){
    g.BM[[r]] <- list()
    B.BM[[r]] <- list()
    psi.BM[[r]] <- list()
    
    for(tt in seq_along(tau.cand)){
      tmp.est <- est.out[which(raw.mat[,3] == tt & raw.mat[,1] == r)]
      B.BM[[r]][[tt]] <- 0
      g.BM[[r]][[tt]] <- list()
      psi.BM[[r]][[tt]] <- list()
      
      for(q in 1:s){
        B.BM[[r]][[tt]] <- B.BM[[r]][[tt]] + 
          tmp.est[[q]]$BBM/s
        
        g.BM[[r]][[tt]][[q]] <- tmp.est[[q]]$gBM
        
        psi.BM[[r]][[tt]][[q]] <- tmp.est[[q]]$psiBM
      }
    }
  }
  
  non.size <- s*(s-1)/2
  non.mat <- matrix(nrow = R*non.size*tau.size, ncol = 4)
  cc <- 1
  for(r in 1:R)
    for(tt in seq_along(tau.cand))
      for(p in 1:(s-1))
        for(q in (p+1):s){
          non.mat[cc, ] <- c(r, tt, p, q)
          cc <- cc + 1
        }
  
  L.all <- mclapply(1:nrow(non.mat), function(ii){
    r <- non.mat[ii, 1]
    tt <- non.mat[ii, 2]
    p <- non.mat[ii, 3]
    q <- non.mat[ii, 4]
    
    p.non <- non.over[[r]][((p-1)*m+1):(p*m)]
    q.non <- non.over[[r]][((q-1)*m+1):(q*m)]
    
    A.non <- A[p.non, q.non]
    
    L.temp <- matrix(0, nrow = length(loss), ncol = 1)
    row.names(L.temp) <- loss
    colnames(L.temp) <- as.character(tau.cand[tt])
    
    P.BM <- B.BM[[r]][[tt]][g.BM[[r]][[tt]][[p]],
                            g.BM[[r]][[tt]][[q]] ] *
      tcrossprod(psi.BM[[r]][[tt]][[p]],
                 psi.BM[[r]][[tt]][[q]])
    
    P.BM[P.BM < 1e-6] <- 1e-6
    P.BM[P.BM > 1- 1e-6] <- 1 - 1e-6
    
    for(lq in seq_along(loss)){
      tmp.nm <- loss[lq]
      L.temp[tmp.nm, 1] <-  L.temp[tmp.nm, 1] +
        do.call(loss[lq], list(as.numeric(A.non), P.BM))/(s*(s-1)*0.5)
    }
    
    return(L.temp)},
    mc.cores = ncore)
  
  L <- list()
  
  for(r in 1:R){
    L[[r]] <- do.call('cbind', lapply(1:tau.size, function(tt){
      Reduce('+', L.all[which(non.mat[,2] == tt & non.mat[,1] == r)])
    }))
    
    row.names(L[[r]]) <- loss
    colnames(L[[r]]) <- as.character(tau.cand)
  }
  
  obj <- data.table(`Candidate Tau` = tau.cand)
  
  for(lq in seq_along(loss))
    for(r in 1:R){
      obj[[paste0(loss[lq], "-Rep=", r)]] <- 
        as(rbind(L[[r]][loss[lq], ]), "vector")
    }
  
  obj2 <- list()
  
  obj2[["Candidate Tau"]] <- tau.cand
  
  for(lq in seq_along(loss)){
    
    obj2[[paste0("tau.hat.each.rep (", loss[lq], ")")]] <- sapply(1:R, function(r){
      tau.BM <- tau.cand[which.min(L[[r]][loss[lq],])]
    })
    
    
    obj2[[paste0(loss[lq], ".model")]] <- 
      mean(obj2[[paste0("tau.hat.each.rep (", loss[lq], ")")]])
  }
  
  return(c(list('loss' = obj), 
           obj2))
}































