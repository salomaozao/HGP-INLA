## Description:
##
## blockNNGP_reg() function  fits  NNGP models  using INLA
# Arguments:
## case: name to identify NNGP results
## loc: locations
## y: observed data
## X: covariates
## w: spatial random effects (used to compare possterior estimations with the true values.
## num.nb: number of neighbor observations.
## dir.save: directory to save the results.
## formula:
##  ‘y ~ 1 + x + f(idx, model = NNGP.model)’
##
## For different family distribution chage the inla() function on line 162.
## For a list of possible alternatives and use ‘inla.doc’ for detailed
##  docs for individual families.
##
##  Value:
##      ‘inla’ returns an object of class ‘"inla"’.
##


NNGP = function(case="NNGP", loc,  y, X, w,  dir.save, num.nb ){

###%%%%%%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## The next codes implement the Adjacency matrix of NNGP and objects that the
## 'inla.rgeneric.NNGP.model' function needs, for specific number of
## neigbors (num.nb).
## The user does not need to change it.

n.blocks <- 1
nloc <- dim(loc)[1]

case= 'simNNGP'

file1 <- paste('/case',case,'_',nloc,'_',num.nb,sep="")

coords <- loc 


nloc <- dim(loc)[1]     

ind <- sort.int(loc[,2], index.return=TRUE)

x2 <- ind$x
indexsort <- ind$ix
x1 <- loc[indexsort,1]


sortloc <- cbind(x1, x2)
dist.mat <- rdist(sortloc)

AdjMatrix <-  matrix(0, nloc, nloc)
AdjMatrix[1,1] <-0


for (j in 2:n){
if (j <= num.nb+1) AdjMatrix[1:(j-1),j] = 1
if (j > num.nb+1){
 ind1 <-  (sort(dist.mat[,j], index.return=TRUE))$ix
 ind <- (ind1[which(ind1<j)])[1:num.nb] 
 AdjMatrix[ind,j] = 1
}
}



g1 = graph.adjacency(AdjMatrix, mode='directed')
is.dag(g1)



## to make things simple let reorder loc and other stuff
rnorm_n.obs.orig <- rnorm_n.obs
w.orig 		<- w

ind1	 	<- sort.int(loc[,2], index.return=TRUE)
loc 		<- loc[(ind1$ix),]
rnorm_n.obs 	<- rnorm_n.obs[(ind1$ix)]
y 		<- y[(ind1$ix)]
w 		<- w[(ind1$ix)]
X 		<- X[(ind1$ix),]

Y=y
coords=loc

## mask for precision-NNGP
coords.D 	<- rdist(loc)
C1 <- 0.04 *exp(-phi*coords.D)
invC   <-  Prec_NNGP(coords,AdjMatrix,C1) 
invCsp <- as.matrix(invC)
invCsp[which(invC>0)] <- 1
invCsp[which(invC<0)] <- 1
invCsp[which(invC==0)] <- 0

W = invCsp
W <- as(W, "sparseMatrix")

###%%%%%%%%%%%%%%%%%%%%%%% END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## run model with INLA
library("INLA")
#INLA:::inla.dynload.workaround() 

## set the NNGP as latent effect
NNGP.model <- inla.rgeneric.define(inla.rgeneric.NNGP.model, W = W,coords.D=coords.D, AdjMatrix=AdjMatrix)

# response variable and covariates
data11<- data.frame(y=y,x =X[,2])
data11$idx <- 1:nrow(data11)

# formula to fit the NNGP model
# f() ensure that the latent effect is the NNGP
f.NNGP <- y ~ 1 + x + f(idx, model = NNGP.model)

# inla function to fit the model
# The user can change the family and priors here.
resf <- inla(f.NNGP, data = as.data.frame(data11), family = "gaussian")

# Recovering the posterior mean estimation of nugget, marginal variance and phi parameters.
# It depends on the internal representation of the hyperparameters.
tau.est = 1/resf$summary.hyperpar$mean[1]
sigmasq.est = exp(-resf$summary.hyperpar$mean[2])
phi.est = 30 - (29)/ (1 + exp(resf$summary.hyperpar$mean[3]))

summary.theta <- c(tau.est,sigmasq.est,phi.est)
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
           ylim = range(loc[,2]),
 zlim=c(min(int.elev$z,na.rm=TRUE)-1,max(int.elev$z,na.rm=TRUE)+1),
           xlab='Longitude', ylab='Latitude')

## plot estimated spatial field
image.plot(int.elev2  , main= paste('nb=',num.nb),
	   xaxs = 'r', yaxs = 'r',
           xlim = range(loc[,1]),
           ylim = range(loc[,2]),
           zlim=c(min(int.elev$z,na.rm=TRUE)-1,max(int.elev$z,na.rm=TRUE)+1),
           xlab='Longitude', ylab='Latitude')
dev.off()



return(resf)

}

