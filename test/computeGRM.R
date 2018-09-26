

computeGRM = function(URobj, snpsubset=NULL, indsubset=NULL, phat=NULL, ephat=NULL, thres = 0.01){
  
  if(is.null(indsubset))
    indsubset <- 1:URobj$.__enclos_env__$private$nInd
  if(is.null(snpsubset))
    snpsubset <- 1:URobj$.__enclos_env__$private$nSnps
  if(is.null(phat))
    phat <- URobj$.__enclos_env__$private$pfreq[snpsubset]
  if(is.null(ephat))
    ephat = rep(0, length(snpsubset))

  temp <- which(phat < 1-thres & phat > thres)
  snpsubset <- snpsubset[temp]
  phat <- phat[temp]
  nSnps <- length(snpsubset)
  nInd <- length(indsubset)
  ep <- matrix(ephat[temp], nrow=nInd, ncol=nSnps, byrow=T)
  depth <- URobj$.__enclos_env__$private$ref[indsubset,snpsubset] + URobj$.__enclos_env__$private$alt[indsubset,snpsubset]
  ratio <- URobj$.__enclos_env__$private$ref[indsubset,snpsubset]/depth
  ploid = URobj$.__enclos_env__$private$ploid
  badSnps <- which(depth > 50)
  ## Compute the adjusted GRM
  genon <- ploid*(ratio - rep.int(phat, rep(nInd, nSnps)))
  genon[is.na(ratio)] <- 0
  genon[depth < 2] <- 0
  genon[badSnps] <- 0
  P0 <- matrix(phat,nrow=nInd,ncol=nSnps,byrow=T)
  P1 <- 1-P0
  P0[depth < 2] <- 0
  P1[depth < 2] <- 0
  P0[badSnps] <- 0
  P1[badSnps] <- 0
  div0 <- (ploid)*tcrossprod(P0,P1)
  GRM <- (tcrossprod(genon/sqrt(1-4*ep*(1-ep))) - tcrossprod(sqrt((ploid*ep)^2*(1-4*P0*P1)/(1-4*ep*(1-ep)))))/div0
  depth.temp <- depth
  depth.temp3 <- depth
  depth.temp[which(depth < 2)] <- 0
  depth.temp <- 1/depth.temp
  depth.temp2 <- depth.temp
  depth.temp[is.infinite(depth.temp)] <- 1
  depth.temp2[is.infinite(depth.temp2)] <- 0
  depth.temp3[which(depth.temp3 < 2)] <- 0
  depth.temp3[which(depth.temp3 > 1)] <- 1
  adj <- depth.temp3*(ploid)^2*(P0*P1*(depth.temp + 4*ep*(1-ep)*(1-depth.temp)) +
                                    ep*(ep+(1-ep)*depth.temp - 4*P0*P1))
  diag(GRM) <- rowSums(  (genon^2 - adj)/((1-depth.temp2)*(1-4*ep*(1-ep))))/diag(div0)
  return(GRM)
}
