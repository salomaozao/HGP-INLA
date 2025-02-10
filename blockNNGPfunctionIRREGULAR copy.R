## case: name to identify regular blockNNGP results
## loc: locations
## y: observed data
## X: covariates
## w: spatial random effects (used to compare possterior estimations with the true values.
## n.blocks: number of  blocks  (irregular)
## num.nb: number of neighbor blocks.
## dir.save: directory to save the results
## formula:
##  ‘y ~ 1 + x + f(idx, model = blockNNGP.model)’
##
## For different family distribution chage the inla() function on line 162.

blockNNGP = function(
  case = "irregular",
  loc,
  sf,
  y,
  X,
  w.obs,
  dir.save,
  n.partition,
  num.nb,
  coords.D
) {
  # createblocks = function (loc, n.blocks){
  #  create blocks via kdtree
  nloc <- dim(loc)[1]
  blocks <- NULL
  loc.blocks <- matrix(NA, n.blocks, 2)
  nb <- NULL # Pontos por bloco

  points <- data.frame(x = loc[, 1], y = loc[, 2]) # Cria dataframe com coordenadas.
  tree <- kdtree(points) # Cria kd-tree
  treenew <- tree[1:(n.blocks - 1), ] # cria subsets de kd-tree para dividir
  blocks <- kdtree_blocks(treenew, n.blocks, loc) # atribui pontos aos blocos

  

  for (k in 1:n.blocks) {
    indblock <- which(blocks == k)
    loc.blocks[k, ] <- c(mean(loc[indblock, 1]), mean(loc[indblock, 2]))
  }

  #  sort blocks,
  # blocks 	  	<- NULL
  # blocks 	  	<- kdtree_blocks(treenew, n.blocks, loc)

  #  reorder the data points
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
  dist.mat <- hdist(sortloc)

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
  #  reorder the data points according to the block order

  # g1 = graph.adjacency(AdjMatrix, mode = 'directed')

  ind1 <- sort.int(blocks, index.return = TRUE)
  loc <- loc[(ind1$ix), ]

  # return(list(
  #   blocks = blocks,          # Vetor de atribuição de blocos
  #   loc.blocks = sortloc,     # Centróides ordenados dos blocos
  #   AdjMatrix = AdjMatrix,    # Matriz de adjacência
  #   loc.sorted = loc,         # Localizações reordenadas
  #   idx.sort = ind1$ix        # Índices de ordenação original
  # ))
  # }

  blocks <- blocks[(ind1$ix)]

  y <- y[(ind1$ix)]
  X <- X[(ind1$ix), ]
  w <- w[(ind1$ix)]

  #   # Passo 1: Criar estrutura de blocos
  #   block_structure <- create_blocks(loc, n.partition, num.nb)

  #   # Passo 2: Reordenar dados conforme blocos
  #   y <- y[block_structure$idx.sort]
  #   X <- X[block_structure$idx.sort, ]
  #   w.obs <- w.obs[block_structure$idx.sort]
  #   loc <- block_structure$loc.sorted

  rnorm_n.obs <- rnorm_n.obs[(ind1$ix)]

  ### needed indexes to built the precision matrix of block_NNGP
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
  nloc <- dim(loc)[1]
  ind_obs1 <- which(blocks == 1)
  num1 <- seq(1:length(ind_obs1))

  indb <- NULL
  for (k in 1:(n.blocks - 1)) {
    indb[[k]] <- util.index(k + 1, blocks, AdjMatrix, newindex)
  }

  ## mask for precision-blockNNGP
  coords.D <- hdist(loc)
  C1 <- exp(-0.04 * coords.D)
  invC <- PrecblockNNGP(n, n.blocks, C1, nb, ind_obs1, num1, indb)
  invCsp <- as.matrix(invC)
  invCsp[which(invC > 0)] <- 1
  invCsp[which(invC < 0)] <- 1
  invCsp[which(invC == 0)] <- 0

  W = invCsp
  W <- as(W, "sparseMatrix")

  ###%%%%%%%%%%%%%%%%%%%%%%% END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## run model with INLA

  ## set the blockNNGP as latent effect
  blockNNGP.model <- inla.rgeneric.define(
    inla.rgeneric.blockNNGP.model,
    W = W,
    n = n,
    n.blocks = n.blocks,
    nb = nb,
    ind_obs1 = ind_obs1,
    num1 = num1,
    indb = indb,
    coords.D = coords.D
  )

  # response variable and covariates
  data1 <- data.frame(y = y, x = X[, 2])
  data1$idx <- 1:nrow(data1)

  # formula to fit the blockNNGP model
  # f() ensure that the latent effect is the blockNNGP
  f.blockNNGP <- y ~ 1 + x + f(idx, model = blockNNGP.model)

  # inla function to fit the model
  # The user can change the family and priors here.
  resf <- inla(f.blockNNGP, data = as.data.frame(data1), family = "gaussian")

  # Recovering the posterior mean estimation of nugget, marginal variance and phi parameters.
  # It depends on the internal representation of the hyperparameters.
  tau.est <- inla.emarginal(
    foo <- function(x) {
      1 / x
    },
    resf$marginals.hyperpar[[1]]
  )
  sigmasq.est <- inla.emarginal(
    foo <- function(x) {
      exp(-x)
    },
    resf$marginals.hyperpar[[2]]
  )
  phi.est = inla.emarginal(
    foo <- function(x) {
      1 - 1 / (1 + exp(x))
    },
    resf$marginals.hyperpar[[3]]
  )
  summary.theta <- c(tau.est, sigmasq.est, phi.est)
  print(c("tau.est", "sigmasq.est", "phi.est"))
  print(summary.theta)

  ## posterior mean of spatial random effects
  est.w = resf$summary.random$idx$mean

  return(resf)
}
