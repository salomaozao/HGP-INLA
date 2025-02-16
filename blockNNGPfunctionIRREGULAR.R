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
  num.nb
) {
  if (!'sf' %in% class(sf)) stop("sf is not a spatial data frame.")

  get_HGPdata = function(n.blocks, blocks, sf) {
    block_struc = get_blocksdata(loc, sf, n.blocks, num.nb)
    ind1 = block_struc$ind1
    AdjMatrix = block_struc$AdjMatrix
    blocks <- block_struc$blocks

    y <<- y[(ind1$ix)]
    X <<- X[(ind1$ix), ]
    w <<- w.obs[(ind1$ix)]
    blocks <- blocks[(ind1$ix)]
    sf <- sf[(ind1$ix), ]
    rnorm_n.obs <<- rnorm_n.obs[(ind1$ix)]

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

  HGPdata = get_HGPdata(n.blocks, NULL, sf) # Note que 'blocks' será definido dentro de get_HGPdata
  W = HGPdata$W
  nb = HGPdata$nb
  ind_obs1 = HGPdata$ind_obs1
  indb = HGPdata$indb
  coords.D = HGPdata$coords.D

  ###%%%%%%%%%%%%%%%%%%%%%%% END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## run model with INLA

  ## set the blockNNGP as latent effect
  blockNNGP.model <- HGPdata$model

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
  tau.est <- inla.emarginal(function(x) 1 / x, resf$marginals.hyperpar[[1]])
  sigmasq.est <- inla.emarginal(
    function(x) exp(-x),
    resf$marginals.hyperpar[[2]]
  )
  phi.est <- inla.emarginal(
    function(x) 1 - 1 / (1 + exp(x)),
    resf$marginals.hyperpar[[3]]
  )
  summary.theta <- c(tau.est, sigmasq.est, phi.est)
  print(c("tau.est", "sigmasq.est", "phi.est"))
  print(summary.theta)

  ## posterior mean of spatial random effects
  est.w = resf$summary.random$idx$mean

  return(resf)
}
