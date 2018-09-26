


GRM <- R6Class("GRM",
                     inherit = UR,
                     public = list(
                       initialize = function(URobj){
                         private$ref       <- URobj$.__enclos_env__$private$ref
                         private$alt       <- URobj$.__enclos_env__$private$alt
                         private$ratio     <- URobj$.__enclos_env__$private$ref/(URobj$.__enclos_env__$private$ref + URobj$.__enclos_env__$private$alt)
                         private$chrom     <- URobj$.__enclos_env__$private$chrom
                         private$pos       <- URobj$.__enclos_env__$private$pos
                         private$SNP_Names <- URobj$.__enclos_env__$private$SNP_Names
                         private$indID     <- URobj$.__enclos_env__$private$indID
                         private$nSnps     <- URobj$.__enclos_env__$private$nSnps
                         private$nInd      <- URobj$.__enclos_env__$private$nInd
                         private$gform     <- URobj$.__enclos_env__$private$gform
                         private$AFrq      <- URobj$.__enclos_env__$private$AFrq
                         private$infilename<- URobj$.__enclos_env__$private$infilename
                         private$ploid     <- URobj$.__enclos_env__$private$ploid
                         private$pfreq     <- URobj$.__enclos_env__$pfreq
                         private$gfreq     <- URobj$.__enclos_env__$gfreq
                         private$ep        <- URobj$.__enclos_env__$ep
                       },
                       computeGRM = function(p = "pest", thres = 0.01){
                         
                         
                         if(p == "pest")
                           phat = private$pfreq
                         else if(p == "gest")
                           phat = private$gfreq
                         
                         snpsubset <- which(phat < 1-thres & phat > thres)
                         phat <- phat[snpsubset]
                         nSnps <- length(snpsubset)
                         nInd <- private$nInd
                         ep <- matrix(private$ep[snpsubset], nrow=nInd, ncol=nSnps, byrow=T)
                         ratio <- private$ratio[,snpsubset]
                         depth <- private$ref[,snpsubset] + private$alt[,snpsubset]
                         ## Compute the adjusted GRM
                         genon <- ratio - (2*ploid)*rep.int(phat, rep(nInd, nSnps))
                         genon[is.na(ratio)] <- 0
                         genon[depth < 2] <- 0
                         P0 <- matrix(phat,nrow=nInd,ncol=nSnps,byrow=T)
                         P1 <- 1-P0
                         P0[depth < 2] <- 0
                         P1[depth < 2] <- 0
                         div0 <- (2*ploid)*tcrossprod(P0,P1)
                         GRM <- (tcrossprod(genon/sqrt(1-4*ep*(1-ep))) - tcrossprod(sqrt((2*ploid*ep)^2*(1-4*P0*P1)/(1-4*ep*(1-ep)))))/div0
                         depth.temp <- depth
                         depth.temp[which(depth < 2)] <- 0
                         depth.temp <- 1/depth.temp
                         depth.temp2 <- depth.temp
                         depth.temp[is.infinite(depth.temp)] <- 1
                         depth.temp2[is.infinite(depth.temp2)] <- 0
                         depth.temp3 <- depth
                         depth.temp3[which(depth.temp3 < 2)] <- 0
                         depth.temp3[which(depth.temp3 > 1)] <- 1
                         adj <- depth.temp3*(2*ploid)^2*(P0*P1*(depth.temp + 4*ep*(1-ep)*(1-depth.temp)) +
                                                           ep*(ep+(1-ep)*depth.temp - 4*P0*P1))
                         diag(GRM) <- rowSums(  (genon^2 - adj)/((1-depth.temp2)*(1-4*ep*(1-ep))))/diag(div0)
                         return(GRM)
                       }
                     ),
                     private = list(
                       ratio  = NULL
                     )
)