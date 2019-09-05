
auto_sc <- function(data = NULL) {
  ns <- dim(data)[1]
  nv <- dim(data)[2]
  data_sd <- apply(data, 2, sd)
  null_var <- NULL
  if (any(data_sd == 0))
    null_var <- which(data_sd == 0)
  data_sd <- t(matrix(as.matrix(data_sd), nv, ns))
  data_mean <- t(matrix(as.matrix(colMeans(data)), nv, ns))
  results <- (data - data_mean) / data_sd
  if (!is.null(null_var)) 
    results[, null_var] <- 0
  return(results)
}

blk_sc <- function(data = NULL) {
  nv <- dim(data)[2]
  scaling_factor <- 1 / sqrt(nv - 1)
  data_scale <- auto_sc(data)
  data_scale <- data_scale * scaling_factor
  return(data_scale)
}

c.factor <- function(..., recursive = TRUE) unlist(list(...), recursive = recursive)

mbpca_host <- function(X, meta, instrument = FALSE, model = "consensus", no_pcs = 2){
  ns <- dim(X)[1]
  nvar <- dim(X)[2]
  no_factors <- NCOL(meta)
  if (instrument){
    lvls <- unlist(meta[[1]])
    unique_lvls <- unique(lvls)
    block_label <- unlist(meta[[2]])
    no_lvls <- length(unique_lvls)
    cnt <- 0
    block_id <- list()
    sub_matrix <- NA
    X_in <- matrix(0, no_lvls, 2)
    nPCs <- NA
    for (i in 1:no_lvls){
      block_id <- c(block_id, paste("Block: ", unique_lvls[i], sep = ""))
      nvars <- length(which(lvls == unique_lvls[i]))
      X_in[i, ] <- c(cnt + 1, cnt + nvars)
      cnt <- cnt + nvars
      if (is.na(nPCs))
        nPCs <- no_pcs
      else
        nPCs <- c(nPCs, no_pcs)
    }
    matrix_now <- X
    if (model == "consensus"){
      result <- cpca(matrix_now, X_in, nPCs)
      } else{ 
          if (model == "hierarchical"){
            result <- hpca(matrix_now, X_in, nPCs)
          } 
          else
            stop("Unsupported MB-PCA model! Available options are either consensus or hierarchical")
      }
    results <- list(model_names = "Data source blocking", models = result, block_label = block_label, block_id = block_id)
    return(results)
  }else{
    mb_matrices <- list()
    factor_name = names(meta)
    no_lvls <- rep(0, no_factors)
    for (i in 1:length(factor_name)){
      no_lvls[i] <- length(unique(meta[, factor_name[i]]))
    }
    model_names <- list()
    model_collection <- list()
    for (i in 1:no_factors){
      model_names[i] <- paste("Blocking for ", factor_name[i], sep = "")
      matrix_now <- NA
      unique_lvls <- unique(meta[, i])
      block_label <- NA
      nPCs <- NA
      for (ii in 1:length(unique_lvls)){
        cnt <- 1
        block_id <- list()
        for (iii in 1:no_factors){
          if (i == iii) next
          no_lvls_sub <- no_lvls[iii]
          X_in <- matrix(0, no_lvls_sub, 2)
          unique_lvls_sub <- unique(meta[, iii])
          sub_matrix <- NA
  
          for (iiii in 1:no_lvls_sub){
            block_id <- c(block_id, paste("Block: ", factor_name[iii], "=", 
                                          unique_lvls_sub[iiii], sep = ""))
            if (all(is.na(sub_matrix)))
              sub_matrix <- X[which((meta[, i] == unique_lvls[ii]) & meta[, iii] == unique_lvls_sub[iiii]) , ]
            else{
              ns1 <- dim(sub_matrix)[1]
              X_append <- X[which((meta[, i] == unique_lvls[ii]) & meta[, iii] == unique_lvls_sub[iiii]), ]
              ns2 <- dim(X_append)[1]
              if (ns1 == ns2)
                sub_matrix <- cbind(sub_matrix, X_append)
              else {
                if (ns1 > ns2){
                  no_missing <- ns1 - ns2
                  mat_add <- t(matrix(colMeans(X_append), nvar, no_missing))
                  X_append <- rbind(X_append, mat_add)
                } else{
                  no_remove <- ns2 - ns1
                  X_append <- X_append[-(1:no_remove), ]
                }
                sub_matrix <- cbind(sub_matrix, X_append)
              }
            }
            X_in[cnt, 1] <- (cnt - 1)*nvar + 1
            X_in[cnt, 2] <- cnt*nvar
            nPCs[cnt] <- no_pcs
            cnt <- cnt + 1
          }
        }
        if (all(is.na(matrix_now))){
          matrix_now <- sub_matrix
          block_label <- rep(unique_lvls[ii], dim(sub_matrix)[1])
  #        sub_matrix <- NA
        }
        else{
          matrix_now <- rbind(matrix_now, sub_matrix)
          block_label <- c.factor(block_label, rep(unique_lvls[ii], dim(sub_matrix)[1]))
   #       sub_matrix <- NA
        }
  
      }
      
      if (model == "consensus"){
        result <- cpca(matrix_now, X_in, nPCs)
      }
      else { if (model == "hierarchical") {
        result <- hpca(matrix_now, X_in, nPCs)
      }else{
        stop("Unsupported MB-PCA model! Available options are either consensus or hierarchical")
      }
      }
      model_collection[[i]] <- list(model_names = model_names[[i]], models = result, block_label = block_label, block_id = block_id)
    }
    return(model_collection)
  }
}

cpca <- function(data = NULL, Xin = NULL, nPC = NULL, xscale = TRUE, tol = 1e-8){
  maxiter <- 500
  data <- as.matrix(data)
  ns <- dim(data)[1]
  nv <- dim(data)[2]
  nb <- dim(Xin)[1]
  print(paste("number of samples = ", ns, "; number of variables = ", nv, "; number of blocks = ", nb))
  data_mean <- colMeans((data))
  data <- data - t(matrix(as.matrix(data_mean), nv, ns))

  if (xscale) {
    for (i in 1:nb){
      rowi <- Xin[i,1]:Xin[i,2]
      data_scale <- blk_sc(data[, rowi])
      data[, rowi] <- data_scale
    }
  }
  maxPC <- max(nPC)
  Tb <- matrix(0, ns, maxPC * nb)
  Pb <- matrix(0, nv, maxPC)
  Wt <- matrix(0, nb, maxPC)
  Tt <- matrix(0, ns, maxPC)
  ssq <- matrix(0, maxPC, 1 + nb)
  Rbo <- matrix(0, ns, maxPC * nb)
  Rbv <- matrix(0, nv, maxPC)
  Lbo <- matrix(0, ns, maxPC * nb)
  Lbv <- matrix(0, nv, maxPC)
  Lto <- matrix(0, ns, maxPC)
  ssqX <- matrix(0, nb + 1, 1)
  ssqX[1] <- sum(colSums(data^2))
  
  for (i in 1:nb){
    rowi <- Xin[i,1]:Xin[i,2]
    ssqX[i+1] <- sum(colSums(data[,rowi]^2))
  }
  
  eigval <- eigs(data %*% t(data), maxPC)
  v <- eigval$vectors
  
  for (i in 1:maxPC){
    iter <- 0
    Tt[, i] <- v[, i]
    t_old <- Tt[,i] * 100
    while ( sum((t_old - Tt[,i])^2) > tol & iter < maxiter){
      iter <- iter + 1
      t_old <- Tt[, i]
      for (ii in 1:nb) {
        if (nPC[ii] >= i) {
          rowi <- Xin[ii,1]:Xin[ii,2]
          coli <- (i-1) * nb + ii
          Pb[rowi, i] <- t(data[, rowi]) %*% Tt[, i] / sum(Tt[,i]^2)
          Pb[rowi, i] <- Pb[rowi, i] / sqrt(sum(Pb[rowi,i]^2))
          Tb[, coli] <- data[, rowi] %*% Pb[rowi, i] / sum(Pb[rowi, i]^2)
        }
      }
      index <- ((i - 1) * nb + 1) : (i*nb)
      Wt[, i] <- t(Tb[, index]) %*% Tt[, i] / sum(Tt[, i]^2)
      Wt[, i] <- Wt[, i] / sqrt(sum(Wt[, i]^2))
      Tt[, i] <- Tb[, index] %*% Wt[, i] / sum(Wt[, i]^2)
    }
  
    if (iter == maxiter) 
      print("Warning: maximum number of iterations reached before convergence")
    for (ii in 1:nb){
      if (nPC[ii] >= i){
        rowi <- Xin[ii,1] : Xin[ii, 2]
        Pb[rowi, i] <- t(data[, rowi]) %*% Tt[, i] / sum(Tt[, i]^2)
        data[, rowi] <- data[, rowi] - Tt[, i] %*% t(Pb[rowi, i])
      }
    }
    
    ssq[i, 1] <- (ssqX[1] - sum(colSums(data^2))) / ssqX[1]
    
    for (ii in 1:nb) {
      rowi <- Xin[ii,1] : Xin[ii, 2]
      coli <- (i-1)*nb + ii
      ssq[i, ii + 1] <- (ssqX[ii + 1] - sum(colSums(data[, rowi]^2))) / ssqX[ii + 1]
      Rbo[, coli] <- sqrt(rowSums(data[, rowi]^2))
      index <- seq(from = ii, to = (i-1)*nb + ii, by = nb)
      Lbo[, coli] <- diag(Tb[, index] %*% ginv(t(Tb[, index]) %*% Tb[, index]) %*% t(Tb[, index]))
    }
    Rbv[, i] <- as.matrix(sqrt(colSums(data^2)))
    Lbv[, i] <- diag(Pb[, 1:i] %*% t(Pb[, 1:i]))
    Lto[, i] <- diag(Tt[, 1:i] %*% ginv(t(Tt[, 1:i]) %*% Tt[, 1:i]) %*% t(Tt[, 1:i]))
  }
  results <- list(block_scores = Tb, block_loadings = Pb, super_scores = Tt, block_weight = Wt, 
                  ssq = ssq, Obj_residue = Rbo, blk_var_residue = Rbv, Obj_leverage = Lbo, 
                  var_leverage = Lbv, super_obj_leverages = Lto)
  return(results)
}

hpca <- function(data = NULL, Xin = NULL, nPC = NULL, xscale = TRUE, tol = 1e-8){
  maxiter <- 500
  ns <- dim(data)[1]
  nv <- dim(data)[2]
  nb <- dim(Xin)[1]
  maxPC <- max(nPC)
  
  data <- as.matrix(data)
  data_mean <- colMeans((data))
  data <- data - t(matrix(as.matrix(data_mean), nv, ns))
  
  if (xscale) {
    for (i in 1:nb){
      rowi <- Xin[i,1]:Xin[i,2]
      data_scale <- blk_sc(data[, rowi])
      data[, rowi] <- data_scale
    }
  }

  
  Tb <- matrix(0, ns, maxPC * nb)
  Pb <- matrix(0, nv, maxPC)
  Wt <- matrix(0, nb, maxPC)
  Tt <- matrix(0, ns, maxPC)
  ssq <- matrix(0, maxPC, 1 + nb)
  Rbo <- matrix(0, ns, maxPC * nb)
  Rbv <- matrix(0, nv, maxPC)
  Lbo <- matrix(0, ns, maxPC * nb)
  Lbv <- matrix(0, nv, maxPC)
  Lto <- matrix(0, ns, maxPC)
  ssqX <- matrix(0, nb + 1, 1)
  ssqX[1] <- sum(colSums(data^2))
  
  for (i in 1:nb) {
    rowi <- Xin[i,1] : Xin[i,2]
    ssqX[i + 1] <- sum(colSums(data[, rowi]^2))
  }
  
  for (i in 1:maxPC){
    iter <- 0
    eigval <- eigs(data %*% t(data), 1)
    v <- eigval$vectors
    Tt[, i] <- v
    t_old <- Tt[, i] * 100
    
    while ((sum((t_old - Tt[, i])^2) > tol) & (iter < maxiter)) {
      iter = iter + 1
      t_old <- Tt[, i]
      for (ii in 1:nb){
        if (nPC[ii] >= i){
          rowi <- Xin[ii, 1] : Xin[ii, 2]
          coli <- (i-1) * nb + ii
          Pb[rowi, i] <- t(data[, rowi]) %*% Tt[, i] / sum(Tt[, i]^2)
          Tb[, coli] <- data[, rowi] %*% Pb[rowi, i] / sum(Pb[,i]^2)
          Tb[, coli] <- Tb[, coli] / sqrt(sum(Tb[, coli]^2))
        }
      }
      
      index <- ((i - 1) * nb + 1) : (i * nb)
      Wt[, i] <- t(Tb[, index]) %*% Tt[, i] / sum(Tt[, i]^2)
      Tt[, i] <- Tb[, index] %*% Wt[, i] / sum(Wt[, i]^2)
      Tt[, i] <- Tt[, i] / sqrt(sum(Tt[, i]^2))
    }
    
    if (iter == maxiter)
      print("WARNING: maximum number of iterations reached before convergence")
    
    data <- data - Tt[, i] %*% t(Pb[, i])
    
    ssq[i ,1] <- (ssqX[1] - sum(colSums(data^2))) / ssqX[1]
    
    for (ii in 1:nb){
      rowi <- Xin[ii, 1] : Xin[ii, 2]
      coli <- (i-1) * nb + ii
      ssq[i, ii + 1] <- (ssqX[ii + 1] - sum(colSums(data[, rowi]^2))) / ssqX[ii + 1]
      Rbo[, coli] <- sqrt(sum(rowMeans(data[, rowi]^2)))
      index <- seq(from = ii, to = (i - 1) * nb + ii, by = nb)
      Lbo[, coli] <- diag(Tb[, index] %*% ginv(t(Tb[, index]) %*% Tb[, index]) %*% t(Tb[, index]))
    }
    
    Rbv[, i] <- as.matrix(sqrt(colSums((data^2))))
    Lbv[, i] <- diag(Pb[, 1:i] %*% t(Pb[, 1:i]))
    Lto[, i] <- diag(Tt[, 1:i] %*% ginv(t(Tt[, 1:i]) %*% Tt[, 1:i]) %*% t(Tt[, 1:i]))
  }
  
  results <- list(block_scores = Tb, block_loadings = Pb, super_scores = Tt, block_weight = Wt, 
                  ssq = ssq, Obj_residue = Rbo, blk_var_residue = Rbv, Obj_leverage = Lbo, 
                  var_leverage = Lbv, super_obj_leverages = Lto)
  return(results)
}

mbpls <- function(data = NULL, label = NULL, nPC = NULL,  Xin = NULL, Yin = NULL, xcentre = TRUE, Xscale = FALSE, tol = 1e-8) {
  maxiter <- 2000
  ns <- dim(data)[1]
  nv <- dim(data)[2]
  nvY <- dim(label)[2]
  maxLV <- max(nPC)
  maxb <- which(nPC == maxLV)
  
  if (length(maxb) > 1) maxb <- maxb[1]
  if (is.null(Yin)) Yin <- t(as.matrix(c(1, nvY)))
  if (xcentre) 
    data <- data - t(matrix(as.matrix(colMeans(data)), nv, ns))
  
  nbX <- dim(Xin)[1]
  nbY <- dim(Yin)[1]
  Tb <- matrix(0, ns, maxLV * nbX)
  Pb <- matrix(0, nv, maxLV)
  Wb <- matrix(0, nv, maxLV)
  Wt <- matrix(0, nbX, maxLV)
  Tt <- matrix(0, ns, maxLV)
  Ub <- matrix(0, ns, maxLV * nbY)
  Qb <- matrix(0, nvY, maxLV)
  Wu <- matrix(0, nbY, maxLV)
  Tu <- matrix(0, ns, maxLV)
  B <- matrix(0, ns, maxLV*nvY)
  ssq <- matrix(0, maxLV, nbX + nbY + 2)
  ssqXY <- matrix(0, nbX + nbY + 2, 1)
  ssqXY[1] <- sum(colMeans(data)^2)
  ssqXY[nbx + 2] <- sum(colSums(label^2))
  
  for (i in 1:nbX) {
    rowi <- Xin[i,1] : Xin[i,2]
    ssqXY[i+1] <- sum(colSums(data[, rowi]^2))
  }
  
  for (i in 1:nbY) {
    rowi <- Yin[i,1]:Yin[i,2]
    ssqXY[i + nbX + 2] <- sum(colSums(label[,rowi]^2))
  }
  
  for (i in 1:maxLV){
    iter <- 0
    Tu[, i] <- label[,1]
    Tt[, i] <- data[, Xin[maxb, 1]]
    t_old <- Tt[, i] * 100
    while((sum((t_old - Tt[, i])^2)) > tol & (iter < maxiter)) {
      iter <- iter + 1
      for (ii in 1:nbX) {
        if (nLV[ii] >= i) {
          rowi <- Xin[ii,1] : Xin[ii,2]
          coli <- (i-1) * nbX + ii
          Wb[rowi, i] <- t(data[, rowi]) %*% Tu[, i] / sum(Tu[, i]^2)
          Wb[rowi, i] <- Wb[rowi, i] / sqrt(sum(Wb[rowi, i]^2))
          Tb[, coli] <- data[, rowi] %*% Wb[rowi, i] / sum(Wb[rowi, i]^2)
        }
      }
      index <- ((i-1) * nbX + 1) : (i * nbX)
      Wt[, i] <- t(Tb[, index]) %*% Tu[, i] / sum(Tu[, i]^2)
      Wt[, i] <- Wt[, i] / sqrt(sum(Wt[, i]^2))
      Tt[, i] <- Tb[, index] %*% Wt[, i] / sum(Wt[, i]^2)
      for (ii in 1:nbY) {
        rowi <- Yin[aa, 1]:Yin[aa,2]
        coli <- (i-1) * nbY + ii
        Qb[rowi, i] <- t(label[, rowi]) %*% Tt[, i] / sum(Tt[, i]^2)
        Ub[, coli] <- label[, rowi] %*% Qb[rowi, i] / sum(Qb[rowi, i]^2)
      }
      index <- ((i-1)*nbY + 1) : (a*nbY)
      Wu[, i] <- t(Ub[, index]) %*% Tt[, i] / sum(Tt[, i]^2)
      Wu[, i] <- Wu[, i] / sqrt(sum(Wu[, i]^2))
      Tu[, i] <- Ub[, index] %*% Wu[, i] / sum(Wu[, i]^2)
    }
    if (iter == maxiter)
      print("Warning: maximum number of iterations reached before convergence")
    rowi <- Xin[i, 1]:Xin[i, 2]
    for (ii in i+1:nbX){
      rowi <- cbind(rowi, Xin[ii, 1]:Xin[ii, 2])
    }
    Pb[rowi, i] <- t(data[, rowi]) %*% Tt[, i] / sum(Tt[, i]^2)
    data <- data - Tt[, i] %*% t(Pb[, i])
    label <- label - Tt[, i] %*% t(Qb[, i])
    if (i > 1) {
      index <- ((i - 1) * nvY + 1) : (i * nvY)
      B[, index] <- B[, index - nvY]
    }
    index <- ((i - 1) * nvY + 1) : (i * nvY)
    B[, index] <- B[, index] + Wb[, i] %*% solve(t(Pb[, i]) %*% Wb[, i]) %*% t(Qb[, i])
    ssq[i, 1] <- (ssqXY[1] - sum(colSums(data^2))) / ssqXY[1]
    ssq[i, nbX + 2] <- (ssqXY(nbX + 2)- sum(colSums(label^2))) / ssqXY[nbX + 2]
    for (ii in 1:nbX) {
      rowi <- Xin[ii, 1]:Xin[ii, 2]
      coli <- (i-1)*nbX + aa
      ssq[i, ii + 1] <- (ssqXY(ii + 1) - sum(colSums(data[, rowi]^2))) / ssqXY[ii + 1]
    }
    for (ii in 1:nbY) {
      rowi <- Yin[ii,1]:Yin[ii,2]
      coli <- (i-1) * nbY + ii
      ssq[i, ii + 2 + nbX] <- (ssqXY[ii + 2 + nbX] - sum(colSums(Y[, rowi]^2))) / ssqXY[ii + 2 + nbX]
    }
  }
  YRegCoeff <- Wb %*% solve(t(Pb) %*% Wb) %*% t(Qb)
  results <- list(XBlockScores = Tb, XBlockLoadings = Pb, XWeights = Wb, XSupWeights = Wt, 
                  XSupScores = Tt, YBlockScores = Qb, YBlockWeights = Wu, YSupScores = Tu,
                  XRegCoeff = B, YRegCoeff = YRegCoeff, ssq = ssq)
}

mbpls_pred <- function(data = NULL, amean = NULL, model = NULL) {
  if (!is.null(amean))
    data <- data - t(matrix(as.matrix(amean), nv, ns))
  results <- data %*% model$YRegCoeff
  return(results)
}

plot.mbpca <- function(mb_model, model_names = NULL, block_label = NULL, block_id = NULL){
  super_scores <- mb_model$super_scores
  block_scores <- mb_model$block_scores
  block_loadings <- mb_model$block_loadings
  ssq <- mb_model$ssq
  no_blocks <- dim(ssq)[2] - 1
  unique_label <- unique(block_label)
  no_lvls <- length(unique_label)
  
  colors <- c("red", "blue", "green", "cyan", "yellow", "magenta", "black", "aquamarine", "orange", 
                "brown",  "coral", "chocolate", "plum", "darkgray", "gold", "salmon", "wheat", "gray")
  tev1 <- ssq[1, ]
  tev2 <- ssq[2, ] - ssq[1, ]
  
  # plot.new()
  par(mfrow = c(1, 1))
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(super_scores[, 1], super_scores[, 2], type = "n", main = paste(model_names, ": Super scores plot", sep = ""), 
       xlab = paste("PC 1, TEV = ", round(tev1[1]*10000)/100, "%", sep = ""),
       ylab = paste("PC 2, TEV = ", round(tev2[1]*10000)/100, "%", sep = "")
       )
  for (i in 1:no_lvls){
    points(super_scores[which(block_label == unique_label[i]), 1], 
           super_scores[which(block_label == unique_label[i]), 2],
           col = colors[i])
  }
  legend(x = 'topright', inset = c(-.3, 0), legend = unique_label, col = colors[1:no_lvls], pch = 1, bty = "n")
  if (no_blocks <= 3) {
    no_row <- 1
    no_col <- no_blocks
  } else{
    no_row <- round(sqrt(no_blocks))
    no_col <- ceiling(no_blocks / no_row)
  }
 # plot.new()
  par(mfrow = c(no_row, no_col))
  par(mar=c(5, 4, 4, 2) + 0.1, xpd=FALSE)
  for (i in 1:no_blocks){
    plot(block_scores[, i], block_scores[, i + no_blocks], type = "n", main = block_id[i], 
         xlab = paste("PC 1, TEV = ", round(tev1[i + 1]*10000)/100, "%", sep = ""), 
         ylab = paste("PC 2, TEV = ", round(tev2[i + 1]*10000)/100, "%", sep = ""))
    for (ii in 1:no_lvls){
      points(block_scores[which(block_label == unique_label[ii]), i], 
             block_scores[which(block_label == unique_label[ii]), i + no_blocks],
             col = colors[ii])
    }
  }

}
  

