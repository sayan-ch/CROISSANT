devtools::source_url("https://github.com/sayan-ch/CROISSANT/blob/main/all_base.R?raw=TRUE")


###############################################################################
nn <- 10000
KK <- 20
p.intra <- 0.3
p.inter <- 0.1
psps <- rep(1, nn)

KKCAND <- 30

BB <- diag(p.intra-p.inter, KK, KK) + matrix(p.inter, KK, KK)


setwd('~/ecv/sbm')
if(!dir.exists("out_sbm"))
  dir.create("out_sbm")

write_csv(data.table(`simulation` = 0, `model` = "", n = 0, K = 0,
                     `cand.K` = 0, B = "",
                     psi = 0,
                     `avg.deg` = 0,
                     `time` = 0,
                     `time.st` = 0,
                     `best.l2` = "",
                     `best.bd` = "",
                     `best.auc` = "",
                     `best.l2.st` = "",
                     `best.bd.st` = "",
                     `best.auc.st` = "",
                     `buffer` = "",
                     `l2.only.sbm` = "",
                     `l2.only.dcbm` = "",
                     `l2.only.sbm.st` = "",
                     `l2.only.dcbm.st` = "",
                     `bd.only.sbm` = "",
                     `bd.only.dcbm` = "",
                     `bd.only.sbm.st` = "",
                     `bd.only.dcbm.st` = "",
                     `auc.only.sbm` = "",
                     `auc.only.dcbm` = "",
                     `auc.only.sbm.st` = "",
                     `auc.only.dcbm.st` = ""),
          "out_sbm/ecv_sbm_10-30.csv", append = F)


for(i in 1:100){
  system.time(net <- 
                blockmodel.gen.fast(n = nn, K = KK, B = BB,
                                    psi = psps, ncore = 20))
  
  
  tmp.out <- mclapply(1:20, function(sta){
    time.ecv <- system.time(out.ecv <- 
                    ECV.for.blockmodel(net$A, KKCAND, B = 3)
    )[3]
    
    l2.only.sbm <- paste0("SBM-",which.min(out.ecv$l2))
    l2.only.dcbm <- paste0("DCSBM-",which.min(out.ecv$dc.l2))
    l2.both <- out.ecv$l2.model
    
    bd.only.sbm <- paste0("SBM-",which.min(out.ecv$dev))
    bd.only.dcbm <- paste0("DCSBM-",which.min(out.ecv$dev))
    bd.both <- out.ecv$dev.model
    
    auc.only.sbm <- paste0("SBM-",which.max(out.ecv$auc))
    auc.only.dcbm <- paste0("DCSBM-",which.max(out.ecv$dc.auc))
    auc.both <- out.ecv$auc.model
    
    return(list(time.ecv = time.ecv, 
                l2.only.sbm = l2.only.sbm, l2.only.dcbm = l2.only.dcbm, 
                l2.both = l2.both,
                bd.only.sbm = bd.only.sbm, bd.only.dcbm = bd.only.dcbm, 
                bd.both = bd.both,
                auc.only.sbm = auc.only.sbm, auc.only.dcbm = auc.only.dcbm, 
                auc.both = auc.both))
    
  }, mc.cores = 20)
  
  l2.only.sbm <- l2.only.dcbm <- l2.both <- vector()
  bd.only.sbm <- bd.only.dcbm <- bd.both <- vector()
  auc.only.sbm <- auc.only.dcbm <- auc.both <- vector()
  time.ecv <- vector()
  
  for(sta in 1:20){
    time.ecv[sta] <- tmp.out[[sta]]$time.ecv
    
    l2.only.sbm[sta] <- tmp.out[[sta]]$l2.only.sbm
    l2.only.dcbm[sta] <- tmp.out[[sta]]$l2.only.dcbm
    l2.both[sta] <- tmp.out[[sta]]$l2.both
    
    bd.only.sbm[sta] <- tmp.out[[sta]]$bd.only.sbm
    bd.only.dcbm[sta] <- tmp.out[[sta]]$bd.only.dcbm
    bd.both[sta] <- tmp.out[[sta]]$bd.both
    
    auc.only.sbm[sta] <- tmp.out[[sta]]$auc.only.sbm
    auc.only.dcbm[sta] <- tmp.out[[sta]]$auc.only.dcbm
    auc.both[sta] <- tmp.out[[sta]]$auc.both
  }  
  
  
  write_csv(data.table(`simulation` = i, `model` = "sbm", n = nn, K = KK,
                       `cand.K` = KKCAND, B = paste0(p.intra, "-", p.inter),
                       psi = 1,
                       `avg.deg` = mean(rowSums(net$A)),
                       `time` = time.ecv[1],
                       `time.st` = sum(time.ecv),
                       `best.l2` = l2.both[1],
                       `best.bd` = bd.both[1],
                       `best.auc` = auc.both[1],
                       `best.l2.st` = modal(l2.both),
                       `best.bd.st` = modal(bd.both),
                       `best.auc.st` = modal(auc.both),
                       `buffer` = "",
                       `l2.only.sbm` = l2.only.sbm[1],
                       `l2.only.dcbm` = l2.only.dcbm[1],
                       `l2.only.sbm.st` = modal(l2.only.sbm),
                       `l2.only.dcbm.st` = modal(l2.only.dcbm),
                       `bd.only.sbm` = bd.only.sbm[1],
                       `bd.only.dcbm` = bd.only.dcbm[1],
                       `bd.only.sbm.st` = modal(bd.only.sbm),
                       `bd.only.dcbm.st` = modal(bd.only.dcbm),
                       `auc.only.sbm` = auc.only.sbm[1],
                       `auc.only.dcbm` = auc.only.dcbm[1],
                       `auc.only.sbm.st` = modal(auc.only.sbm),
                       `auc.only.dcbm.st` = modal(auc.only.dcbm)
  ), "out_sbm/ecv_sbm_10-30.csv", append = T)
   gc()
}








