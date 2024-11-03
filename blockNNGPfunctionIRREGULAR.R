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

blockNNGP = function(case="irregular", loc,  y, X, w.obs , dir.save, n.partition, num.nb){

###%%%%%%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## The next codes implement the Adjacency matrix of blockNNGP and objects that the
## 'inla.rgeneric.blockNNGP.model' function needs, for specific number of
## blocks ( n.blocks) and num.nb neighbor blocks.
## The user does not need to change it.


nloc 	<- dim(loc)[1]
points 	<- data.frame(x=loc[ ,1], y=loc[ ,2])
tree 	<- kdtree(points)
treenew	<- tree[1:(n.blocks-1), ] 
blocks 	<- NULL
blocks 	<- kdtree_blocks(treenew, n.blocks, loc)
file1 	<- paste('/case',case,'_',nloc,'_',n.blocks,'_',num.nb,sep="")

png(paste(dir.save,'/case',case,'_',nloc,'_',n.blocks,'_',num.nb,"_NNGPblocks1_new.png",sep=""))
loc.blocks<- matrix(NA, n.blocks, 2)
nb 	<- NULL
plot(loc,  main= paste('M = ',n.blocks), pch='.')
for(k in 1:n.blocks){
	indblock 	<- which(blocks==k)
	nb	 	<- c(nb, length(indblock))
	loc.blocks[k,1] <- mean(loc[indblock, 1])
	loc.blocks[k,2] <-mean(loc[indblock, 2])
	points(loc[indblock,1:2],col=k,pch=19)
	text(loc[indblock,1]+0.01, loc[indblock,2], k, col='red', cex=0.8)
	}
dev.off()


#########################################################

blocks 	  	<- NULL
blocks 	  	<- kdtree_blocks(treenew, n.blocks, loc)
ind 	  	<- sort.int(loc.blocks[ ,2], index.return=TRUE)
x2 	  	<- ind$x
indexsort 	<- ind$ix
x1 	  	<- loc.blocks[indexsort, 1]
new.locblocks 	<- cbind(x1, x2)
blocks0 	<- blocks

for(i in(1:n.blocks)){
	indi 	     <- which(blocks0==ind$ix[i])
	blocks[indi] <- i
}



if(n.blocks == 8  | n.blocks == 16) indr = 4 
if(n.blocks == 32 | n.blocks == 64) indr = 8
if(n.blocks == 128 ) indr = 16 


blocksr 	<- 1:n.blocks
indexsort1	<- NULL
for(j in 1: (n.blocks/indr)){
	h1 		<- new.locblocks[(((j-1)*indr)+1):(j*indr), ]
	indh1 		<- sort.int(h1[,1], index.return = TRUE)
	indexsort1 	<- c(indexsort1, indh1$ix+ ((j-1)*indr))
}


blocks01 	<- blocks
for(i in(1:n.blocks)){
	indi 		<- which(blocks01 == indexsort1[i])
	blocks[indi] 	<- i
}


newf.locblocks 	<- new.locblocks[indexsort1, ]
x1		<- newf.locblocks[ , 1]
x2		<- newf.locblocks[ , 2]



sortloc 	<- cbind(x1, x2)
dist.mat 	<- rdist(sortloc)
AdjMatrix 	<-  matrix(0, n.blocks, n.blocks)
AdjMatrix[1,1] 	<-0


for (j in 2:n.blocks){
if (j <= num.nb+1) AdjMatrix[1:(j-1), j] = 1
if (j > num.nb+1){
 ind1 <-  (sort(dist.mat[ , j], index.return=TRUE))$ix
 ind <- (ind1[which(ind1 < j)])[1:num.nb]
 AdjMatrix[ind, j] = 1
}
}

AdjMatrix


g1 = graph.adjacency(AdjMatrix, mode = 'directed')
is.dag(g1)


ind1 	<- sort.int(blocks, index.return = TRUE)
loc 	<- loc[(ind1$ix), ]
blocks 	<- blocks[(ind1$ix)]
y 	<- y[(ind1$ix)]
X 	<- X[(ind1$ix), ]
w 	<- w[(ind1$ix)]
rnorm_n.obs <- rnorm_n.obs[(ind1$ix)]

### needed indexes to built the precision matrix of block_NNGP
  newindex     	<- NULL
  nb           	<- matrix(NA, n.blocks, 1)
  nb[1]        	<- length(which(blocks == 1))
  for (j in 1:n.blocks){
    ind_obs    	<- which(blocks == j)
    newindex   	<- c(newindex, ind_obs)
    if(j > 1){
    nbj 	<- length(ind_obs)
    nb[j] 	<- nb[j-1] + nbj
    }
  }
  nloc         	<- dim(loc)[1]
  ind_obs1      <- which(blocks == 1)
  num1          <- seq(1:length(ind_obs1))


indb <- NULL
for (k in 1:(n.blocks-1)){
indb[[k]] <- util.index(k+1, blocks, AdjMatrix, newindex)
}


## mask for precision-blockNNGP
coords.D 	<- rdist(loc)
C1 <-  exp(-0.04*coords.D)
invC   <-  PrecblockNNGP(n, n.blocks,C1,nb,ind_obs1,num1,indb)
invCsp <- as.matrix(invC)
invCsp[which(invC>0)] <- 1
invCsp[which(invC<0)] <- 1
invCsp[which(invC==0)] <- 0


W = invCsp
W <- as(W, "sparseMatrix")

###%%%%%%%%%%%%%%%%%%%%%%% END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## run model with INLA

## set the blockNNGP as latent effect
blockNNGP.model <- inla.rgeneric.define(inla.rgeneric.blockNNGP.model, W = W, n = n, 
					n.blocks = n.blocks, nb = nb, ind_obs1 = ind_obs1, 
					num1 = num1, indb = indb, coords.D = coords.D)

# response variable and covariates
data1		<- data.frame(y = y, x = X[,2])
data1$idx 	<- 1:nrow(data1)

# formula to fit the blockNNGP model
# f() ensure that the latent effect is the blockNNGP
f.blockNNGP 	<- y ~ 1 + x + f(idx, model = blockNNGP.model)


# inla function to fit the model
# The user can change the family and priors here.
resf <- inla(f.blockNNGP, data = as.data.frame(data1), family = "gaussian")

# Recovering the posterior mean estimation of nugget, marginal variance and phi parameters.
# It depends on the internal representation of the hyperparameters.
tau.est 	<- 1/resf$summary.hyperpar$mean[1]
sigmasq.est 	<- exp(resf$summary.hyperpar$mean[2])
phi.est = 30 - (29)/ (1 + exp(resf$summary.hyperpar$mean[3]))
summary.theta 	<- c(tau.est, sigmasq.est, phi.est)
print( c("tau.est","sigmasq.est","phi.est"))
print(summary.theta)


## posterior mean of spatial random  effects
est.w = resf$summary.random$idx$mean

## plot to compare true spatial random  effects and their posterior mean estimation.
png(paste(dir.save,file1,"_NNGPblocks12.png",sep=""))
par(mfrow=c(1,2))
# plot original spatial field
int.elev <- mba.surf(cbind(loc,w), 100, 100, extend=TRUE)$xyz.est
int.elev2 <- mba.surf(cbind(loc,est.w), 100, 100, extend=TRUE)$xyz.est
image.plot(int.elev , main='True w',
	   xaxs = 'r', yaxs = 'r',
           xlim = range(loc[,1]),
           ylim = range(loc[,2]),           zlim=c(min(int.elev$z,na.rm=TRUE)-1,max(int.elev$z,na.rm=TRUE)+1),
           xlab='Longitude', ylab='Latitude')

## plot estimated spatial field
image.plot(int.elev2  , main= paste('M = ',n.blocks,', nb=',num.nb),
	   xaxs = 'r', yaxs = 'r',
           xlim = range(loc[,1]),
           ylim = range(loc[,2]),
           zlim=c(min(int.elev$z,na.rm=TRUE)-1,max(int.elev$z,na.rm=TRUE)+1),
           xlab='Longitude', ylab='Latitude')
dev.off()

return(resf)

}



