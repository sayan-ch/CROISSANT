library(parallel)
library(Matrix)
library(dplyr)
#library(permute)
library(RSpectra)
library(cluster)
library(Rfast)
library(data.table)
#library(raster)
library(readr)
library(irlba)

sumFast <- function(X){
  if(is.vector(X))  return(sum(X))
  return(sum(rowSums(X)))
}

modal <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

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
  
  E <- best.perm.label.match(lab, fixed, K)
  
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
    B.sum[K,K] <- sumFast(A) + 1e-6
    
    if(psi.omit > 0){
      psi <- as.numeric(rowSums(A[-(1:psi.omit),])/B.sum[K,K])
      return(list(Bsum = B.sum/p.sample, psi = psi))
    }
    
    psi <- as.numeric(rowSums(A)/B.sum[K,K]) + 1e-6
    
    return(list(Bsum = B.sum/p.sample, psi = psi))
  }
  
  G <- lapply(1:K, function(k) which(g == k))
  
  for(k in 1:K){
    for(l in k:K){
      B.sum[k,l] <- B.sum[l,k] <- sumFast(A[G[[k]], G[[l]]]) + 1e-6
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

# auroc <- function(score, bool) {
#   n1 <- sum(!bool)
#   #n2 <- sum(bool)
#   n2 <- length(score) - n1
#   U  <- sum(rank(score)[!bool]) - n1 * (n1 + 1) / 2
#   return(1 - U / n1 / n2)
# }

AUC <- function(A, P){
  -Rfast::auc(as(A, "vector"), as(P, "vector"))
  #-auroc(as(P, 'vector'), as(A, 'vector'))
}



################################################################################
################################################################################
################################################################################
################################################################################

SBM.DCBM.fast <- function(A, K.CAND,
                          s, o, R, tau = 1, laplace = F,
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
    
    if(!laplace){
      L.sonn <- A.sonn.tau
    }else{
      L.sonn <- tcrossprod(crossprod(d.sonn.tau, A.sonn.tau), d.sonn.tau)
    }
    
    #eig.max <- eigs_sym(L.sonn, K.max, "LM")$vectors
    eig.max <- irlba(L.sonn, nu = K.max, nv = K.max)$v
    
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
                                          cluster.only = T,
                                          pamonce = 6))
      
      #psi.hat[[k.cand]] <- sqrt(rowSums(eig.max[,1:work.K]^2))
      rownorm <- sqrt(rowSums(eig.max[, 1:work.K]^2))
      
      rn.eig <- eig.max[,1:work.K]/rownorm
      
      #out.DCBM[[k.cand]] <- as.integer(pam(rn.eig, work.K,
      #                                     metric = "manhattan",
      #                                     cluster.only = T))
      
      out.DCBM[[k.cand]] <- as.integer(kmeans(rn.eig, work.K,
                                               nstart = 100,
                                               iter.max = 10^7)$cluster)
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
    
    out.SBM.std <- raw.out[raw.ind[,1] == r][[1]]$SBM[[1]]
    out.DCBM.std <- raw.out[raw.ind[,1] == r][[1]]$DCBM[[1]]
    
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
            do.call(loss[lq], list(A.non, P.SBM))/(s*(s-1)*0.5) 
        }
        next}
      else{
        for(lq in seq_along(loss)){
          tmp.nm <- paste(mod[mq], loss[lq], sep = "_")
          L.temp[tmp.nm, 1] <-  L.temp[tmp.nm, 1] +
            do.call(loss[lq], list(A.non, P.DCBM))/(s*(s-1)*0.5)
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

###############################################################################
nn <- 10000
KK <- 5
p.intra <- 0.05
p.inter <- 0.01
psps <- rep(1, nn)

KKCAND <- 10

BB <- diag(p.intra-p.inter, KK, KK) + matrix(p.inter, KK, KK)


setwd('~/new_model_select/SBM')
dir.create('tmp_sbm')
dir.create('final_sbm')

write_csv(data.table(`simulation` = 0, `model` = "", n = 0, K = 0,
                     `cand.K` = 0, B = "",
                     psi = 0, `sbm.fit.method` = "",
                     `dc.fit.method` = "",
                     `avg.deg` = 0,
                     s = 0, o = 0, R = 0,
                     `time` = 0,
                     `best.l2` = "",
                     `best.bd` = "",
                     `best.auc` = "",
                     `buffer` = "",
                     `l2.only.sbm` = "",
                     `l2.only.dcbm` = "",
                     `l2.both` = "",
                     `bd.only.sbm` = "",
                     `bd.only.dcbm` = "",
                     `bd.both` = "",
                     `auc.only.sbm` = "",
                     `auc.only.dcbm` = "",
                     `auc.both` = ""), 
          "final_sbm/sonnet_sbm_1-5.csv", append = F)


sonpar <- list(c(5, 2000, 1),
               c(5, 1000, 1),
               c(10, 1000, 1),
               c(5, 2000, 3),
               c(5, 1000, 3),
               c(10, 1000, 3),
               c(5, 2000, 5),
               c(5, 1000, 5),
               c(10, 1000, 5))

for(i in 1:100){
  system.time(net <- 
                    blockmodel.gen.fast(n = nn, K = KK, B = BB,
                                        psi = psps, ncore = 20))
  
  for(par in seq_along(sonpar)){
    ss <- sonpar[[par]][1]
    oo <- sonpar[[par]][2]
    RR <- sonpar[[par]][3]
    
    message("Simulation - ", i, ":\t SONNET running for ", ss, "-", oo, "-", RR)
    
    time.sonn <- system.time(out.sonn <- 
                    SBM.DCBM.fast(A = net$A, K.CAND = KKCAND, 
                                  s = ss, o = oo, R = RR,
                                  tau = 0,
                                  laplace = F,
                                  dc.est = 1,
                                  loss = c("l2", "bin.dev", "AUC"),
                                  ncore = 20)
                                )[3]
    
    l2.only <- as.matrix(out.sonn$loss %>% select(starts_with("l2")))
    
    l2.out <- apply(apply(l2.only, 2, function(xx){
      ll <- length(xx)
      
      c(`sbm` = which.min(xx[seq(from = 1, to = ll, by = 2)]),
      `dcbm` = which.min(xx[seq(from = 2, to = ll, by = 2)]))
    }), 1, modal)
    
    bd.only <- as.matrix(out.sonn$loss %>% select(starts_with("bin.dev")))
    
    bd.out <- apply(apply(bd.only, 2, function(xx){
      ll <- length(xx)
      
      c(`sbm` = which.min(xx[seq(from = 1, to = ll, by = 2)]),
        `dcbm` = which.min(xx[seq(from = 2, to = ll, by = 2)]))
    }), 1, modal)
    
    auc.only <- as.matrix(out.sonn$loss %>% select(starts_with("AUC")))
    
    auc.out <- apply(apply(auc.only, 2, function(xx){
      ll <- length(xx)
      
      c(`sbm` = which.min(xx[seq(from = 1, to = ll, by = 2)]),
        `dcbm` = which.min(xx[seq(from = 2, to = ll, by = 2)]))
    }), 1, modal)
    
    
    
    write_csv(data.table(`simulation` = i, `model` = "SBM", n = nn, K = KK,
               `cand.K` = KKCAND, B = paste0(p.intra, "-", p.inter),
               psi = 1, `sbm.fit.method` = "pam+nc",
               `dc.fit.method` = "km+ns100+rn+CL",
               `avg.deg` = mean(rowSums(net$A)),
               s = ss, o = oo, R = RR,
               `time` = time.sonn,
               `best.l2` = out.sonn$l2.model,
               `best.bd` = out.sonn$bin.dev.model,
               `best.auc` = out.sonn$AUC.model,
               `buffer` = "",
               `l2.only.sbm` = l2.out['sbm'],
               `l2.only.dcbm` = l2.out['dcbm'],
               `l2.both` = paste0(out.sonn$`Mod.K.hat.each.rep (l2)`, collapse = " "),
               `bd.only.sbm` = bd.out['sbm'],
               `bd.only.dcbm` = bd.out['dcbm'],
               `bd.both` = paste0(out.sonn$`Mod.K.hat.each.rep (bin.dev)`, collapse = " "),
               `auc.only.sbm` = auc.out['sbm'],
               `auc.only.dcbm` = auc.out['dcbm'],
               `auc.both` = paste0(out.sonn$`Mod.K.hat.each.rep (AUC)`, collapse = " ")
    ), "final_sbm/sonnet_sbm_1-5.csv", append = T)
    
    write_csv(data.table(simulation = i, model = "SBM",n = nn, K = KK,
               `cand.K` = KKCAND, B = paste0(p.intra, "-", p.inter),
               psi = 1, `sbm.fit.method` = "pam+nc",
               `dc.fit.method` = "km+ns100+rn+CL",
               `avg.deg` = mean(rowSums(net$A)),
               s = ss, o = oo, R = RR,
               out.sonn$loss),
              "tmp_sbm/sonnet_sbm_loss_1-5.csv", append = T)
    
    
   gc() 
  }
}








