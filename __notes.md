TODO: 
- Entender de onde vem o theta    
- Entender parametrização

***
    
(?):
- RMSP is the root mean squared error of prediction
- CPP is the frequentist coverage of the 95% prediction interval
- Width is the width of the prediction interval
-  IS is the interval score.

*** 

Several R packages contain procedures for simulation (and conditional simulation) of Gaussian spatial processes, including Spatial and geoRglm. You need to supply a variogram model for the spatial autocorrelation (and both of these contain procedures to help estimate such a model from data)
also check out the RandomFields package, which has a very complete set of Gaussian simulation methods 