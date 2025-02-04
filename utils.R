## blockNNGP and NNGP functions

meancov_nn <- function(Sigma, ind_obs, ind_neigblocks, indnum) {
  invC_nbi <- chol2inv(chol(Sigma[ind_neigblocks, ind_neigblocks]))
  B_bi <- (Sigma[ind_obs, ind_neigblocks]) %*% invC_nbi
  F_bi <- Sigma[ind_obs, ind_obs] - (B_bi %*% Sigma[ind_neigblocks, ind_obs])
  invFbi <- chol2inv(chol(F_bi))
  Bstar_bi <- matrix(0, length(ind_obs), dim(Sigma)[1])

  Bstar_bi[, ind_neigblocks] <- -B_bi
  Bstar_bi[indnum] <- 1

  val <- list(invFbi = invFbi, Bstar_bi = Bstar_bi)
  return(val)
}


PrecblockNNGP <- function(nloc, n.blocks, Sigma, nb, ind_obs1, num1, indb) {
  Fs_1 <- matrix(0, nloc, nloc)
  Bb <- matrix(0, nloc, nloc)
  cat("Sigma: ", length(Sigma), " //  ", nrow(Sigma), "x", ncol(Sigma), "\n")
  cat("ind_obs1: ", ind_obs1, "\n")
  Fs_1[1:nb[1], 1:nb[1]] <- chol2inv(chol(Sigma[ind_obs1, ind_obs1]))

  Bstar_bi <- matrix(0, length(ind_obs1), nloc)

  print(nloc == length(ind_obs1))
  diag(Bstar_bi) <- 1
  # diag(Bstar_bi) <- 1


  Bb[1:nb[1], ] <- Bstar_bi
  #    system.time(
  for (j in 2:n.blocks) {
    ress <- meancov_nn(Sigma, indb[[j - 1]][[1]], indb[[j - 1]][[2]], indb[[j - 1]][[4]])
    Fs_1[(nb[j - 1] + 1):(nb[j]), (nb[j - 1] + 1):(nb[j])] <- ress$invFbi
    Bb[(nb[j - 1] + 1):(nb[j]), ] <- ress$Bstar_bi
  }
  #    )

  Bbb <- as(Bb, "dgCMatrix")
  Fs_11 <- as(Fs_1, "dgCMatrix")

  invCs <- crossprod(Bbb, Fs_11) %*% Bbb

  return(invCs)
}



util.index <- function(i, blocks, AdjMatrix, newindex) {
  ind_obs <- which(blocks == i)
  ind_neig <- which(AdjMatrix[, i] == 1)
  ind_neigblocks <- NULL
  for (l in 1:length(ind_neig)) {
    ind_neigblocks1 <- which(blocks == ind_neig[l])
    ind_neigblocks <- c(ind_neigblocks, ind_neigblocks1)
  }
  index <- match(ind_neigblocks, newindex)
  num <- seq(1:length(ind_obs))
  index2 <- match(ind_obs, newindex)
  indnum <- cbind(num, index2)
  return(list(ind_obs, ind_neigblocks, index, indnum))
}



## NNGP functions

meancov_nn1 <- function(i, loc, AdjMatrix, Sigma) {
  ind_neig <- which(AdjMatrix[, i] == 1)
  invC_nbi <- chol2inv(chol(Sigma[ind_neig, ind_neig])) # inverse of Fbi
  B_bi <- (Sigma[i, ind_neig]) %*% invC_nbi
  F_bi <- Sigma[i, i] - (B_bi %*% Sigma[ind_neig, i])

  Bstar_bi <- matrix(0, 1, dim(loc)[1])
  Bstar_bi[, ind_neig] <- -B_bi
  Bstar_bi[i] <- 1


  val <- list(B_bi = B_bi, F_bi = F_bi, Bstar_bi = Bstar_bi)
  return(val)
}



Prec_NNGP <- function(loc, AdjMatrix, Sigma) {
  nloc <- dim(loc)[1]
  var <- matrix(NA, nloc, 1)
  Fs_1 <- matrix(0, nloc, nloc)
  Bb <- matrix(0, nloc, nloc)

  for (j in 1:nloc) {
    # print(j)
    if (j == 1) {
      var[j] <- Sigma[j, j]
      Fs_1[j, j] <- 1 / var[j]
      Bb[j, j] <- 1
    }
    if (j > 1) {
      res <- meancov_nn1(j, loc, AdjMatrix, Sigma)
      var[j] <- res$F_bi
      Fs_1[j, j] <- 1 / var[j]
      Bb[j, ] <- res$Bstar_bi
    }
  }

  m1 <- as(Bb, "dgCMatrix")
  m2 <- as(Fs_1, "dgCMatrix")
  invCs <- crossprod(m1, m2) %*% m1

  return(invCs)
}



hdist <- function(locMatrix) {
  nrow <- nrow(locMatrix)
  if (ncol(locMatrix) < 2) {
    stop("locMatrix must have at least two columns")
  }

  # Vetor de pontos
  ptsVec <- lapply(1:nrow, function(i) st_point(c(locMatrix[i, 1], locMatrix[i, 2])))

  # Matriz de distância (inicializa com zeros)
  distMatrix <- matrix(0, nrow = nrow, ncol = nrow)

  cat("Calculando matriz triangular de distâncias de Hausdorff, n =", nrow, "\n")
  for (i in 1:(nrow - 1)) {
    if (i %% 10 == 0) {
      print(paste("Linha", i))
    }
    for (j in (i + 1):nrow) {
      dist <- st_distance(ptsVec[[i]], ptsVec[[j]], which = "Hausdorff")
      distMatrix[i, j] <- dist
      distMatrix[j, i] <- dist
    }
  }

  print("Matriz de distâncias calculada com sucesso!")
  return(distMatrix)
}


hdist_sf <- function(locVec) {
  nrow <- length(locVec)
  distMatrix <- matrix(0, nrow = nrow, ncol = nrow)

  # Calcular matriz triangular de distâncias de Hausdorff
  cat("Calculando matriz triangular de distâncias de Hausdorff, n =", nrow, "\n")
  for (i in 1:(nrow - 1)) {
    if (i %% 10 == 0) {
      print(paste("Linha", i))
    }
    for (j in (i + 1):nrow) {
      geom_i <- locVec[[i]]
      geom_j <- locVec[[j]]

      dist <- st_distance(geom_i, geom_j, which = "Hausdorff")
      distMatrix[i, j] <- dist
      distMatrix[j, i] <- dist
    }
  }

  return(distMatrix)
}
