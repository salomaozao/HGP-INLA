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
  n.partitions,
  n.blocks,
  num.nb,
  coords.D
) {
  if (!'sf' %in% class(sf)) stop("sf is not a spatial data frame.")


  block_struc = createblocks(loc, sf, n.blocks, num.nb)
  ind1 = block_struc$ind1
  AdjMatrix = block_struc$AdjMatrix
  blocks <- block_struc$blocks  

  y <- y[(ind1$ix)]
  X <- X[(ind1$ix), ]
  w <- w[(ind1$ix)]
  blocks <- blocks[(ind1$ix)]
  orderedLoc <- sf[(ind1$ix), ] 

  # print(ind1$ix)
  # print(y)
  # print(X)
  # print(w)
  print(blocks)
  #   # Passo 1: Criar estrutura de blocos
  #   block_structure <- create_blocks(loc, n.partition, num.nb)

  #   # Passo 2: Reordenar dados conforme blocos
  # blocks <- blocks[(block_structure$idx.sort)]
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
  nloc <- dim(orderedLoc)[1]
  ind_obs1 <- which(blocks == 1)
  num1 <- seq(1:length(ind_obs1))

  indb <- NULL
  for (k in 1:(n.blocks - 1)) {
    indb[[k]] <- util.index(k + 1, blocks, AdjMatrix, newindex)
  }

  ## mask for precision-blockNNGP
  coords.D <- hdist_sf(orderedLoc)
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
