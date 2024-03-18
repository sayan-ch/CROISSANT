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

net <- blockmodel.gen.fast(1000, 5, diag(0.2, 5, 5) + matrix(0.1, 5, 5),
                           ncore = 6)
A <- net$A; K <- 5; tau.cand <- 0.1*(0:10); DCBM = T
s = 3; o = 250; R = 3
laplace = F; dc.est = 2; loss = c("l2", "bin.dev", "AUC")
ncore = 6

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

library(poweRlaw)
net <- blockmodel.gen.fast(3000, 5, diag(.2, 5, 5) + matrix(.04, 5, 5),
                           psi = rplcon(3000, 1, 5),
                           ncore = 6)

mean(rowSums(net$A))

system.time(out <- croissant.tune.regsp(A = net$A, K = 5,
                                        tau.cand = seq(0, 2, by = 0.1), DCBM = T,
                                        laplace = T,
                                        s = 2, o = 2400, R = 6, ncore = 6))

library(randnet)

cls <- err <- list()

ii <- 1
laplace <- T

n = 3000

for(tau in c(seq(0, 2, by = 0.1), out$l2.model, out$bin.dev.model, out$AUC.model)){
  print(tau)
  
  deg <- rowSums(net$A)
  avg.deg <- mean(deg)
  
  A.tau <- net$A + tau * avg.deg / n
  
  L.tau <- A.tau
  
  if(laplace){
  d.tau <- sparseMatrix( i = 1:n, j = 1:n,
                              x = 1/sqrt(deg + tau*avg.deg))
  
    L.tau <- tcrossprod(crossprod(d.tau, A.tau), d.tau)
    L.tau[is.na(L.tau)] <- 0
  }
  
  eig <- irlba::partial_eigen(L.tau, n = 5, symmetric = T)$vectors
  rowsum <- sqrt(rowSums(eig ^ 2))
  
  rowsum[rowsum == 0] <- 1
  
  eig.rn <- eig / rowsum
  
  cls[[ii]] <- as.integer(pam(x = eig.rn, k = 5,
                              pamonce = 6,
                              do.swap = F,
                              cluster.only = T))
  
  err[[ii]] <-
    100 * mean(net$member != matched.lab(cls[[ii]], net$member))
  
  ii <- ii + 1
}

btt <- c(seq(0, 2, by = 0.1), out$l2.model, out$bin.dev.model, out$AUC.model)

plot(btt, do.call('c', err))

error <- do.call('c', err); error
btt[which.min(error)]





