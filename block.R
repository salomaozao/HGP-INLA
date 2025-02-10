createblocks = function(loc, n.blocks) {
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

  for (i in (1:n.blocks)) {
    indi <- which(blocks == indexsort1[i])
    blocks[indi] <- i
  }

  #  build the adjacency matrix,
  createAdjMatrix <- function(blocks, new.locblocks, indexsort1, num.nb) {
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
  }

  #  reorder the data points according to the block order
  ind1 <- sort.int(blocks, index.return = TRUE)
  loc <- loc[(ind1$ix), ]

  return(
    list(
      blocks = blocks, # Vetor de atribuição de blocos
      loc.blocks = sortloc, # Centróides ordenados dos blocos
      AdjMatrix = AdjMatrix, # Matriz de adjacência
      loc.sorted = loc, # Localizações reordenadas
      idx.sort = ind1$ix # Índices de ordenação original
    )
  )
}

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
