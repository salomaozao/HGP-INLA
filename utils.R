source("irregblock.R")


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

  # diag(Bstar_bi) <- 1 # (?) erro?
  # Bb[1:nb[1], ] <- Bstar_bi
  indices_to_set <- cbind(1:length(ind_obs1), ind_obs1)
  Bstar_bi[indices_to_set] <- 1

  Bb[1:nb[1], 1:nloc] <- Bstar_bi

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

show_results = function(res) {
  rho.est <- phi.est * (log(10))^(1 / alpha)
  cat("\n--- Resultados da Análise Espacial ---\n")
  cat("\nEfeitos Fixos:\n")
  cat(round(beta.est, 4))
  cat("\nHiperparâmetros Espaciais Estimados:\n")
  cat(sprintf("  - Variância Espacial (σ²): %.4f\n", sigmasq.est))
  cat(sprintf("  - Desvio Padrão Espacial (σ): %.4f\n", sqrt(sigmasq.est)))
  cat(sprintf("  - Parâmetro de Alcance (ρ): %.4f\n", rho.est))
  cat("\nCritérios de Ajuste do Modelo:\n")
  cat(sprintf("  - DIC: %.2f\n", resf$dic$dic))
  cat(sprintf("  - WAIC: %.2f\n", resf$waic$waic))
  cat(sprintf("  - Tempo de execução (segundos): %.2f\n", exec_time["elapsed"]))
  cat("-------------------------------------\n")
}

get_HGPdata = function(
  sf,
  n.blocks,
  num.nb,
  alpha = 1,
  priors = list(a = 0, b = 10)
) {
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
    ###############

    contagem_blocos <- as.data.frame(table(blocks)) %>%
      arrange(as.numeric(as.character(blocks)))

    colnames(contagem_blocos) <- c("Bloco", "Pontos")
    print(contagem_blocos) #
    # print(1:n.blocks)

    # print(which(blocks == 1))
    # Poisson: > 13143, 13144, 13145, ...
    # Gaussian: > 42, 163, 172, ...

    # print(length(blocks))
    # Poisson: > 26465
    # Gaussian: > 200

    # print(length(centroids_coords))
    # Poisson: > 268
    # Gaussian: > 400

    ###############

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

    # if (n.blocks %in% c(8, 16)) {
    #   indr <- 4
    # }
    # if (n.blocks %in% c(32, 64)) {
    #   indr <- 8
    # }
    # if (n.blocks == 128) {
    #   indr <- 16
    # }
    indr = sqrt(n.blocks)
    indexsort1 <- NULL

    for (j in 1:(n.blocks / indr)) {
      start = (((j - 1) * indr) + 1)
      end = (j * indr)

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

    coords.D <- st_distance(st_centroid(sf), which = "Hausdorff") |>
      # units::set_units(value = "km") |>
      units::set_units(value = NULL) |>
      as.matrix()

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
        coords.D = coords.D,
        Q = invCsp
      )
    )
  }

  if (!"sf" %in% class(sf)) {
    stop("sf is not a spatial data frame.")
  }

  loc <- st_coordinates(st_centroid(sf))

  # print(n.blocks)
  block_struc = get_blocksdata(loc, sf, n.blocks, num.nb)
  ind1 = block_struc$ind1
  AdjMatrix = block_struc$AdjMatrix
  blocks <- block_struc$blocks

  # bad code?
  # y <- y[(ind1$ix)]
  # X <- X[(ind1$ix), ]
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
        n = n,
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
      ind_obs1 = precMatrixData$ind_obs1,
      indb = precMatrixData$indb,
      coords.D = precMatrixData$coords.D,
      order = ind1$ix,
      blocks = blocks,
      Q = precMatrixData$Q
    )
  )
}


#' @title Summarize a custom HGP model from an INLA object
#'
#' @description Extracts and formats key results, including parameter estimates,
#'              95% credible intervals, and model fit statistics.
#'
#' @param resf A fitted model object from INLA.
#' @param alpha The smoothness parameter (nu) used in the PEC function.
#'              This is required to calculate the practical range rho.
#' @param ci_level The desired level for the credible interval (e.g., 0.95 for 95%).
#'
#' @return A list containing the extracted parameter summaries and fit statistics.
#'         The function also prints a formatted summary to the console.

summarize_inla_hgp_CA <- function(resf, alpha, ci_level = 0.95) {
  # --- Validação de Entradas ---
  if (!inherits(resf, "inla")) {
    stop("O objeto 'resf' não é um resultado válido da função INLA.")
  }
  if (missing(alpha)) {
    stop(
      "O argumento 'alpha' (parâmetro de suavidade nu) é obrigatório para calcular rho."
    )
  }

  # --- Preparação ---
  # Define os quantis para o intervalo de credibilidade
  lower_q <- (1 - ci_level) / 2
  upper_q <- 1 - lower_q
  quants <- c(lower_q, 0.5, upper_q) # Inferior, Mediana, Superior

  # --- 1. Efeitos Fixos (Intercepto) ---
  beta_summary <- resf$summary.fixed[
    "(Intercept)",
    c("0.025quant", "mean", "0.975quant")
  ]
  # Usamos a média como estimativa pontual para consistência com summary()
  beta_vals <- c(
    beta_summary[["0.025quant"]],
    beta_summary[["mean"]],
    beta_summary[["0.975quant"]]
  )

  # --- 2. Hiperparâmetros ---
  hyper_marginals <- resf$marginals.hyperpar

  # Parâmetro 1: Tau (desvio padrão do erro de medição)
  # O hiperparâmetro é a precisão (1/tau^2). A transformação é tau = 1/sqrt(precisão).
  prec_obs_summary <- inla.qmarginal(quants, hyper_marginals[[1]])
  tau_vals <- 1 / sqrt(prec_obs_summary)

  # Parâmetro 2: Sigma^2 (variância do efeito espacial)
  # O hiperparâmetro do rgeneric é log(1/sigma^2). A transformação é sigma^2 = exp(-log(1/sigma^2)).
  log_inv_sigmasq_summary <- inla.qmarginal(quants, hyper_marginals[[2]])
  sigmasq_vals <- exp(-log_inv_sigmasq_summary)

  # Parâmetro 3: Phi (parâmetro base do alcance)
  # O hiperparâmetro do rgeneric é logit(phi). A transformação é a logística (plogis).
  logit_phi_summary <- inla.qmarginal(quants, hyper_marginals[[3]])
  phi_vals <- plogis(logit_phi_summary)

  # --- 3. Parâmetro Derivado (Rho - Alcance Prático) ---
  rho_vals <- phi_vals * (log(10))^(1 / alpha)

  # --- 4. Critérios de Ajuste ---
  dic <- resf$dic$dic
  waic <- resf$waic$waic
  resf$waic
  # LOOIC só está disponível se control.compute = list(loo = TRUE) foi usado
  looic <- if (!is.null(resf$loo)) resf$loo$looic else NA

  # --- 5. Tempo de Execução ---
  cpu_time <- resf$cpu.used["Total"]

  # --- Formatação da Saída ---
  cat("\n--- Resumo do Modelo HGP via INLA ---\n")

  cat("\nEFEITOS FIXOS\n")
  cat(sprintf(
    "  - Intercepto (β0):      %.4f (%.4f, %.4f)\n",
    beta_vals[2],
    beta_vals[1],
    beta_vals[3]
  ))

  cat("\nHIPERPARÂMETROS\n")
  cat(sprintf(
    "  - D.P. do Erro (τ):       %.4f (%.4f, %.4f)\n",
    tau_vals[2],
    tau_vals[1],
    tau_vals[3]
  ))
  cat(sprintf(
    "  - Variância Espacial (σ²): %.4f (%.4f, %.4f)\n",
    sigmasq_vals[2],
    sigmasq_vals[1],
    sigmasq_vals[3]
  ))

  cat("\nPARÂMETRO DE ALCANCE\n")
  cat(sprintf(
    "  - Alcance Prático (ρ):    %.4f (%.4f, %.4f)\n",
    rho_vals[2],
    rho_vals[1],
    rho_vals[3]
  ))

  cat("\nCRITÉRIOS DE AJUSTE DO MODELO\n")
  cat(sprintf("  - DIC:   %.2f\n", dic))
  cat(sprintf("  - WAIC:  %.2f\n", waic))
  cat(sprintf("  - LOOIC: %.2f\n", looic))

  cat("\nTEMPO DE EXECUÇÃO\n")
  cat(sprintf("  - Tempo total (CPU): %.2f segundos\n", cpu_time))
  cat("-------------------------------------\n")

  # --- Retorno dos Resultados (Opcional) ---
  # Retorna uma lista para que os valores possam ser usados posteriormente
  invisible(
    list(
      fixed_effects = data.frame(
        median = beta_vals[2],
        lower = beta_vals[1],
        upper = beta_vals[3],
        row.names = "beta0"
      ),
      hyperparameters = data.frame(
        median = c(tau_vals[2], sigmasq_vals[2]),
        lower = c(tau_vals[1], sigmasq_vals[1]),
        upper = c(tau_vals[3], sigmasq_vals[3]),
        row.names = c("tau", "sigma_sq")
      ),
      range_parameter = data.frame(
        median = rho_vals[2],
        lower = rho_vals[1],
        upper = rho_vals[3],
        row.names = "rho"
      ),
      fit_criteria = c(DIC = dic, WAIC = waic, LOOIC = looic)
    )
  )
}
