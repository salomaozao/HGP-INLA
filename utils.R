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
  Fs_1[1:nb[1], 1:nb[1]] <- chol2inv(chol(Sigma[ind_obs1, ind_obs1]))

  Bstar_bi <- matrix(0, length(ind_obs1), nloc)

  diag(Bstar_bi) <- 1

  Bb[1:nb[1], ] <- Bstar_bi
  #    system.time(
  for (j in 2:n.blocks) {
    ress <- meancov_nn(
      Sigma,
      indb[[j - 1]][[1]],
      indb[[j - 1]][[2]],
      indb[[j - 1]][[4]]
    )
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

hdist_sf <- function(locVec) {
  nrow <- nrow(locVec)
  distMatrix <- matrix(0, nrow = nrow, ncol = nrow)

  cat(
    "Calculando matriz triangular de distâncias de Hausdorff, n =",
    nrow,
    "\n"
  )
  for (i in 1:(nrow - 1)) {
    if (i %% 10 == 0) {
      print(paste("Linha", i))
    }
    for (j in (i + 1):nrow) {
      # Acessar a geometria das linhas i e j
      geom_i <- st_geometry(locVec[i, ])
      geom_j <- st_geometry(locVec[j, ])

      # Calcular distância de Hausdorff
      dist <- st_distance(geom_i, geom_j, which = "Hausdorff")
      distMatrix[i, j] <- as.numeric(dist)
      distMatrix[j, i] <- as.numeric(dist)
    }
  }
  return(distMatrix)
}

get_blocksdata = function(loc, sf, n.blocks, num.nb) {
  nloc <- dim(loc)[1]
  blocks <- NULL
  loc.blocks <- matrix(NA, n.blocks, 2)
  nb <- NULL # Pontos por bloco

  centroids = st_centroid(sf) # get centroids to sort the blocks
  centroids_coords = st_coordinates(centroids)
  points <- data.frame(x = centroids_coords[, 1], y = centroids_coords[, 2]) # Cria dataframe com coordenadas.
  tree <- kdtree(points) # Cria kd-tree
  treenew <- tree[1:(n.blocks - 1), ] # cria subsets de kd-tree para dividir
  blocks <- kdtree_blocks(treenew, n.blocks, loc) # atribui pontos aos blocos

  for (k in 1:n.blocks) {
    indblock <- which(blocks == k)
    loc.blocks[k, ] <- c(
      mean(centroids_coords[indblock, 1]),
      mean(centroids_coords[indblock, 2])
    )
  }

  ind <- sort.int(loc.blocks[, 2], index.return = TRUE)
  new.locblocks <- cbind(loc.blocks[ind$ix, 1], ind$x)
  blocks0 <- blocks

  for (i in (1:n.blocks)) {
    indi <- which(blocks0 == ind$ix[i])
    blocks[indi] <- i
  }

  if (n.blocks %in% c(8, 16)) indr <- 4
  if (n.blocks %in% c(32, 64)) indr <- 8
  if (n.blocks == 128) indr <- 16

  blocksr <- 1:n.blocks
  indexsort1 <- NULL
  for (j in 1:(n.blocks / indr)) {
    h1 <- new.locblocks[(((j - 1) * indr) + 1):(j * indr), ]
    indh1 <- sort.int(h1[, 1], index.return = TRUE)
    indexsort1 <- c(indexsort1, indh1$ix + ((j - 1) * indr))
  }

  blocks01 <- blocks
  for (i in (1:n.blocks)) {
    indi <- which(blocks01 == indexsort1[i])
    blocks[indi] <- i
  }
  #  build the adjacency matrix,
  sortloc <- new.locblocks[indexsort1, ]
  dist.mat <- as.matrix(dist(sortloc))

  AdjMatrix <- matrix(0, n.blocks, n.blocks)

  for (j in 2:n.blocks) {
    if (j <= num.nb + 1) {
      AdjMatrix[1:(j - 1), j] = 1
    } else {
      ind1 <- (sort(dist.mat[, j], index.return = TRUE))$ix
      ind <- (ind1[which(ind1 < j)])[1:num.nb]
      AdjMatrix[ind, j] = 1
    }
  }

  nloc <- dim(loc)[1]
  ind1 <- sort.int(blocks, index.return = TRUE)

  return(
    list(
      ind1 = ind1,
      AdjMatrix = AdjMatrix,
      blocks = blocks
    )
  )
}

get_precMatrixData = function(n.blocks, blocks, AdjMatrix, sf) {
  ### creating new indexes
  newindex <- NULL
  nb <- matrix(NA, n.blocks, 1)
  nb[1] <- length(which(blocks == 1))

  for (j in 1:n.blocks) {
    ind_obs <- which(blocks == j)
    newindex <- c(newindex, ind_obs)
    if (j > 1) {
      nbj <- length(ind_obs)
      nb[j] <- nb[j - 1] + nbj
    }
  }

  nloc <- dim(sf)[1]
  n <- nloc # Definindo n para usar no PrecblockNNGP
  ind_obs1 <- which(blocks == 1)
  num1 <- seq(1, length(ind_obs1))

  indb <- NULL
  for (k in 1:(n.blocks - 1)) {
    indb[[k]] <- util.index(k + 1, blocks, AdjMatrix, newindex)
  }

  ## mask for precision-blockNNGP
  coords.D <- hdist_sf(sf)
  C1 <- exp(-0.04 * coords.D) # for sparseMatrix

  invC <- PrecblockNNGP(n, n.blocks, C1, nb, ind_obs1, num1, indb) # create precision matrix
  invCsp <- as.matrix(invC) # 1 if points connect, 0 otherwise
  invCsp[which(invC > 0)] <- 1
  invCsp[which(invC < 0)] <- 1
  invCsp[which(invC == 0)] <- 0

  W = invCsp
  W <- as(W, "sparseMatrix")

  return(
    list(
      W = W,
      nb = nb,
      ind_obs1 = ind_obs1,
      num1 = num1,
      indb = indb,
      coords.D = coords.D
    )
  )
}

get_HGPdata = function(loc, sf, n.blocks, num.nb ) {
  block_struc = get_blocksdata(loc, sf, n.blocks, num.nb)
  ind1 = block_struc$ind1
  AdjMatrix = block_struc$AdjMatrix
  blocks <- block_struc$blocks

  blocks <- blocks[(ind1$ix)]
  sf <- sf[(ind1$ix), ]

  precMatrixData = get_precMatrixData(n.blocks, blocks, AdjMatrix, sf)

  return(
    list(
      model = inla.rgeneric.define(
        inla.rgeneric.blockNNGP.model,
        W = precMatrixData$W,
        n = n,
        n.blocks = n.blocks,
        nb = precMatrixData$nb,
        ind_obs1 = precMatrixData$ind_obs1,
        num1 = precMatrixData$num1,
        indb = precMatrixData$indb,
        coords.D = precMatrixData$coords.D
      ),
      W = precMatrixData$W,
      nb = precMatrixData$nb,
      ind_obs1 = precMatrixData$ind_obs1,
      indb = precMatrixData$indb,
      coords.D = precMatrixData$coords.D
    )
  )
}
