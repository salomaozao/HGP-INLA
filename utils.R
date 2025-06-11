source("Irregblock.R")


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

#' Calculate Precision Matrix for Block NNGP
#'
#' @param nloc Number of locations
#' @param n.blocks Number of blocks
#' @param Sigma Covariance matrix
#' @param nb Vector containing cumulative number of points per block
#' @param ind_obs1 Indices of observations in first block
#' @param num1 Sequence for first block
#' @param indb List of block information
#'
#' @return Precision matrix (inverse covariance matrix)

PrecblockNNGP <- function(nloc, n.blocks, Sigma, nb, ind_obs1, num1, indb) {
  #  Validation
  if (!is.matrix(Sigma)) {
    stop("Sigma must be a matrix")
  }
  if (nrow(Sigma) != ncol(Sigma)) {
    stop("Sigma must be a square matrix")
  }
  if (length(nb) < n.blocks) {
    stop("nb must have at least n.blocks elements")
  }

  precision_blocks <- matrix(0, nloc, nloc) # (?) use sparseMatrix?
  block_matrix <- matrix(0, nloc, nloc)

  precision_blocks[1:nb[1], 1:nb[1]] <- chol2inv(chol(Sigma[
    ind_obs1,
    ind_obs1
  ]))

  Bstar_bi <- matrix(0, length(ind_obs1), nloc)

  diag(Bstar_bi) <- 1

  block_matrix[1:nb[1], ] <- Bstar_bi
  #    system.time(
  for (j in 2:n.blocks) {
    r <- meancov_nn(
      Sigma,
      indb[[j - 1]][[1]],
      indb[[j - 1]][[2]],
      indb[[j - 1]][[4]]
    )

    start <- nb[j - 1] + 1
    end <- nb[j]
    precision_blocks[start:end, start:end] <- r$invFbi
    block_matrix[start:end, ] <- r$Bstar_bi
  }
  #    )

  Bmat <- as(block_matrix, "dgCMatrix")
  Fmat <- as(precision_blocks, "dgCMatrix")
  invCs = crossprod(Bmat, Fmat) %*% Bmat

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

meancov_nn <- function(Sigma, ind_obs, ind_neigblocks, indnum) {
  invC_nbi <- chol2inv(chol(Sigma[ind_neigblocks, ind_neigblocks]))
  B_bi <- Sigma[ind_obs, ind_neigblocks] %*% invC_nbi
  F_bi <- Sigma[ind_obs, ind_obs] - (B_bi %*% Sigma[ind_neigblocks, ind_obs])
  invFbi <- chol2inv(chol(F_bi))
  Bstar_bi <- matrix(0, length(ind_obs), ncol(Sigma))
  Bstar_bi[, ind_neigblocks] <- -B_bi
  Bstar_bi[indnum] <- 1
  list(invFbi = invFbi, Bstar_bi = Bstar_bi)
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


get_HGPdata = function(
  #loc,
  sf,
  n.blocks,
  num.nb,
  priors = list(a = 0, b = 10), # (?) boa definição?
  alpha = 1
) {
  get_blocksdata = function(sf, n.blocks, num.nb) {
    #nloc <- dim(loc)[1]
    blocks <- NULL
    loc.blocks <- matrix(NA, n.blocks, 2)
    nb <- NULL # Pontos por bloco

    centroids = st_centroid(sf) # get centroids to sort the blocks
    centroids_coords = st_coordinates(centroids)
    nloc <- nrow(centroids_coords)
    points <- data.frame(x = centroids_coords[, 1], y = centroids_coords[, 2]) # Cria dataframe com coordenadas.
    tree <- kdtree(points) # Cria kd-tree
    treenew <- tree[1:(n.blocks - 1), ] # cria subsets de kd-tree para dividir
    blocks <- kdtree_blocks(treenew, n.blocks, centroids_coords) # atribui pontos aos blocos

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

    if (n.blocks %in% c(8, 16)) {
      indr <- 4
    }
    if (n.blocks %in% c(32, 64)) {
      indr <- 8
    }
    if (n.blocks == 128) {
      indr <- 16
    }

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

    #nloc <- dim(loc)[1]
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

    n <- nrow(sf)
    ind_obs1 <- which(blocks == 1)
    num1 <- seq(1, length(ind_obs1))

    indb <- NULL
    for (k in 1:(n.blocks - 1)) {
      indb[[k]] <- util.index(k + 1, blocks, AdjMatrix, newindex)
    }

    ## mask for precision-blockNNGP
    sf <- st_transform(sf, 32611)
    coords.D <- st_distance(sf, which = "Hausdorff") |>
      units::set_units(value = "km") |>
      units::set_units(value = NULL) |>
      as.matrix()

    C1 <- exp(-0.04 * coords.D) # for sparseMatrix

    cat(
      "is C1 positive semi-definite:",
      matrixcalc::is.positive.semi.definite(C1, tol = 1e-8)
    ) # !

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

  if (!"sf" %in% class(sf)) {
    stop("The object is not a spatial data frame.")
  }

  #block_struc = get_blocksdata(loc, sf, n.blocks, num.nb)
  block_struc = get_blocksdata(sf, n.blocks, num.nb)
  ind1 = block_struc$ind1
  AdjMatrix = block_struc$AdjMatrix
  blocks <- block_struc$blocks

  blocks <- blocks[(ind1$ix)]
  sf <- sf[(ind1$ix), ]

  precMatrixData = get_precMatrixData(
    n.blocks,
    blocks,
    AdjMatrix,
    sf
  )

  return(
    list(
      model = inla.rgeneric.define(
        inla.rgeneric.blockNNGP.model,
        W = precMatrixData$W,
        n = nrow(sf),
        n.blocks = n.blocks,
        nb = precMatrixData$nb,
        ind_obs1 = precMatrixData$ind_obs1,
        num1 = precMatrixData$num1,
        indb = precMatrixData$indb,
        coords.D = precMatrixData$coords.D,
        alpha = alpha,
        a = priors$a,
        b = priors$b
      ),
      nb = precMatrixData$nb,
      blocks = blocks,
      ind_obs1 = precMatrixData$ind_obs1,
      indb = precMatrixData$indb,
      coords.D = precMatrixData$coords.D,
      order = ind1$ix
    )
  )
}
