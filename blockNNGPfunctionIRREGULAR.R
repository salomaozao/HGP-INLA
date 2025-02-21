## Description:
##
## blockNNGP_reg() function  fits  block-NNGP models through IRREGULAR BLOCKS using INLA
# Arguments:
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
## For a list of possible alternatives and use ‘inla.doc’ for detailed
##  docs for individual families.
##
##  Value:
##      ‘inla’ returns an object of class ‘"inla"’.
##

blockNNGP = function(
  case = "irregular",
  loc,
  sf,
  y,
  X,
  E = NULL,
  w.obs,
  dir.save,
  n.partitions,
  n.blocks,
  num.nb
) {
  inla.setOption(num.threads = "2:1") # 4 threads principais + 1 para hyper
  if (!'sf' %in% class(sf)) stop("sf is not a spatial data frame.")

  # if this function is defined in utils, the parameters must be function(n.blocks, sf, y, X, w.obs, rnorm_n.obs, num.nb)
  get_HGPdata = function(n.blocks, sf) {
    block_struc = get_blocksdata(loc, sf, n.blocks, num.nb)
    ind1 = block_struc$ind1
    AdjMatrix = block_struc$AdjMatrix
    blocks <- block_struc$blocks

    y <- y[(ind1$ix)]
    X <- X[(ind1$ix), ]
    w <- w.obs[(ind1$ix)]
    blocks <- blocks[(ind1$ix)]
    sf <- sf[(ind1$ix), ]
    rnorm_n.obs <- rnorm_n.obs[(ind1$ix)]

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

  HGPdata = get_HGPdata(n.blocks, sf) # blocks é definido dentro de get_HGPdata
  W = HGPdata$W
  nb = HGPdata$nb
  ind_obs1 = HGPdata$ind_obs1
  indb = HGPdata$indb
  coords.D = HGPdata$coords.D

  ###%%%%%%%%%%%%%%%%%%%%%%% END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## run model with INLA

  ## set the blockNNGP as latent effect
  blockNNGP.model <- HGPdata$model

  # response variable and covariates FOR POINT DATA
  data1 <- data.frame(y = y, x = X[, 2])
  data1$idx <- 1:nrow(data1)

  # formula to fit the blockNNGP model
  # f() ensure that the latent effect is the blockNNGP
  f.blockNNGP <- y ~ 1 + x + f(idx, model = blockNNGP.model)

  # inla function to fit the model
  # The user can change the family and priors here.
  resf <- inla(
    f.blockNNGP,
    data = as.data.frame(data1),
    family = "poisson",
    E = E
  )

  print(names(resf$marginals.hyperpar))

  # Recovering the posterior mean estimation of nugget, marginal variance and phi parameters.
  tau.est <- inla.emarginal(function(x) 1 / x, resf$marginals.hyperpar[[1]])
  sigmasq.est <- inla.emarginal(
    function(x) exp(-x),
    resf$marginals.hyperpar[[2]]
  )
  # phi.est <- inla.emarginal(
  #   function(x) 1 - 1 / (1 + exp(x)),
  #   resf$marginals.hyperpar[[3]]
  # )
  # summary.theta <- c(tau.est, sigmasq.est, phi.est)
  # print(c("tau.est", "sigmasq.est", "phi.est"))

  summary.theta <- c(tau.est, sigmasq.est)
  print(c("tau.est", "sigmasq.est"))  # supposed to be sigmasq.est and rho.est?

  print(summary.theta)

  ## posterior mean of spatial random effects
  est.w = resf$summary.random$idx$mean

  return(resf)
}
