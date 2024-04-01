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

##
library(latentnet)

data(sampson)
samp.fit <- ergmm(samplike ~ euclidean(d=2))

net <- blockmodel.gen.fast(100, 3, diag(0.4, 3,3) + matrix(0.1, 3,3))

adj <- as(net$A, 'dMatrix')

netnet <- as.network(x = adj, 
           matrix.type = "adjacency")

try <- ergmm(netnet ~ euclidean(d=3), tofit = "mkl")

Z.hat <- try$mcmc.mle$Z

plot(try)


## Generating latent space network
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

#####
net <- latent.gen(100, 4, 0, 1, 6)

netnet <- as.network(net$A, matrix.type = "adjacency")

system.time(try <- ergmm(netnet ~ euclidean(d=1)))

Z.hat <- try$mcmc.mle$Z
beta <- try$mcmc.mle$beta


plot(try)

time1 <- system.time(out1 <- croissant.latent(A = net$A, d.cand = 10, 
                                     s = 17, o = 15, R = 5,
                                     loss = c("l2", "bin.dev", "AUC"),
                                     ncore = 6))









