# Main function implementing the two stage algortihm 

fbioSelect = function(y, # response vector (only continuous responses acceptable)
                  id, # a vector with the unique ids corresponding to the response vector supplied
                  t.var, # the varying time or dosage under which the responses were recorded 
                  Z, # a matrix of low-dimensional covariates (object obtained from the model.matrix function)
                  X, # a data.frame of high-dimensional covariates design matrix subject to variable screening 
                  low.dim.cov.number = 1, # number of low-dimensional covariates, minimum is 1
                  threshold.size, # threshold used for the first step of the algorithm (variable screening step)
                  type.of.penalty = c("grLasso", "grMCP", "grSCAD", "gel", "cMCP"), # type of penalty used for the second step (group penalised regression)
                  degree = 3, # degree of the piecewise polynomial basis used—default is 3 for cubic splines.
                  intercept = FALSE, # if TRUE, an intercept is included in the basis; default is FALSE.
                  degrFreedom = NULL, # the degrees of freedom supplied for the basis expansion (splines::bs function)
                  knots = NULL, # the internal breakpoints that define the spline. The default is NULL, which results in a basis for ordinary polynomial regression.
                  weights = NULL, # regression weights, if null then independence covariance structure will be assumed
                  nfolds = 5, # number of folds for the cross-validation selecting the penalty size for the second step of the algorithm
                  seed = 1000){
 
  require(tidyverse) 
  require(foreach)
  require(grpreg)
  require(nlme)
  require(doParallel)
  
  ## add error functions
  registerDoParallel(cores = 3)
  set.seed(seed)
  
  Lgrid = 101 # number of grid points for the coefficient functions
  uniq_id = unique(id) # identify experimental units
  luniq_id = length(uniq_id) # number of experimental units in the data 
  ut.var = sort(unique(t.var)) # identify time/dosage levels
  lut.var = length(ut.var) # number of time/dosage levels in the experiment
  t.varGrid = seq(from = min(ut.var), to = max(ut.var), length = Lgrid) # create time/dosage grid
  high_dim_cov_names = colnames(X) # names of the high dimensional covariates in the data
  high_dim_cov_size = length(high_dim_cov_names) # number of high dimensional covariates subject to screening
  p = ncol(Z)-1 # number of the low-dimensional parameters without accounting for the intercept
  Interval = seq(min(ut.var), max(ut.var), by = 0.01) 
  
  ## Construct the basis function matrix for the low dimensional covariates
  BasisMtrx = BasisFunc.Approximation(designMtrx = Z, grid = t.var, degrFreedom = degrFreedom, 
                                      degree = degree, intercept = intercept, knots = knots)
  
  ## Make sure weights information is not NULL
  temp_datafrm = data.frame(y, t.var, id)
  if (is.null(weights)){
    ni = niFind(data = temp_datafrm, idLabel = "id")
    weights = 1/ni
  }  
  
  #################################################
  ## Step 1: Apply the variable screenign algorithm 
  ## Screening function ( iteratively fit the varying coefficient model for each high dimensional covariate separately and 
  ## calculate each high dimensional covariate's utility (mse) ) based on which the high dimensional covariates are ranked later on
  Ggroup = rep(1:(p+2), each = degrFreedom)
  ## Screening step
  mse_screening = foreach(i = 1:high_dim_cov_size, .combine = c, .inorder = T) %dopar%
    screening(i, BasisMtrx, X, y, Ggroup, luniq_id, uniq_id, lut.var, ut.var,p,
              Lgrid, t.varGrid, degrFreedom, intercept, degree, knots, t.var, weights)
  
  ## Rank the high dimensional covariates based on their utility
  geneRanking = data.frame(hdcov = high_dim_cov_names, score = mse_screening) %>% arrange(score) # rank the high dimensional covariates according to their utility
  ## Select the top ranked high dimensional covariates according to the threshold applied and the ranking based on their utility
  outranked_hdcov = geneRanking[1:threshold.size,] 
  
  ## Update the design matrix based on the variable screening step
  Z2 = cbind(Z, X[, as.character(outranked_hdcov$hdcov)]) 
  ## Update design matrix by adding the top ranked high dimensional covariates
  BasisMtrx2 = BasisFunc.Approximation(designMtrx =  Z2, grid = t.var, degrFreedom = degrFreedom, degree = degree, intercept = intercept, knots = knots)
  
  #################################################
  ## Step 2: Apply a group version of regularisation penalty to the selected high dimensional covariate set to 
  ##further decrease the number of covariates detected as being trully associated with the drug response
  group = c(rep(attributes(Z)$assign, each = degrFreedom)+ifelse(attributes(Z)$assign[1] == 0,1,0), 
            rep((max(attributes(Z)$assign+ifelse(attributes(Z)$assign[1] == 0,1,0))+1):ncol(Z2), each = degrFreedom))
  groupMultiplier = c(rep(0,ncol(Z)), rep(sqrt(degrFreedom), threshold.size)) 
  
  ## Cross-validation for penalty size selection
  lambda = cv.grpreg(BasisMtrx2[,-1], y, group[-1], nfolds = nfolds, penalty = type.of.penalty, weights = weights)$lambda.min
  ## Penalised regression model implementation
  fit = grpreg(BasisMtrx2[,-1], y, group[-1], penalty = type.of.penalty, lambda = lambda, weights = weights, group.multiplier = groupMultiplier)
  
  ## Estimate the coefficient functions
  b = data.frame(fit$beta)
  b[2:degrFreedom,] = b[2:degrFreedom,]+b[1,]
  names(b) = "fit.coefficients"
  b$group = group
  
  Beta = coeff.Functions.Estimation(p = nrow(outranked_hdcov), LgridA = Lgrid, gridA = Interval, degrFreedom = degrFreedom, degree = degree, intercept = intercept, b = b, knots = knots)
  coefFtest = data.frame(beta = as.vector(b$fit.coefficients), group = b$group)
  
  ## Detect the non-zero coefficient functions
  coefFnonZero = c()
  for (ir in 1:length(unique(coefFtest$group))){
    tmp = coefFtest %>% filter(group == ir)
    if (sum(tmp$beta) != 0){
      coefFnonZero = rbind(coefFnonZero, tmp)
    }
  }
  
  Ahdcov_est = outranked_hdcov$hdcov[unique(coefFnonZero$group)[-c(1:(low.dim.cov.number+1))]-(low.dim.cov.number+1)] # active high dimensional covariate set
 
  rslt = list(gene.ranking.marginal.screening = geneRanking, 
              active.hdcov.pen.regression = Ahdcov_est,
              coefFunctions.pen.regression = Beta, 
              regression.weights = weights,
              pen.regression.fit = fit)
  return(rslt)
}

BasisFunc.Approximation = function(designMtrx, grid, degrFreedom, degree, intercept, knots){
  require(splines)
  require(doParallel)
 
  B0 = bs(grid, df =  degrFreedom ,intercept = intercept, degree = degree, knots = knots)
  modsplinesA = foreach(a = iter(as.matrix(designMtrx[,-1]), by = "col"),  .combine = "cbind") %do% apply(B0, 2, function(x){x*a})
  modsplines = cbind(B0, modsplinesA)
  
  return(modsplines)
}

niFind = function(data, idLabel){
  
  freq = data.frame(table(data[, idLabel]))
  names(freq) = c(idLabel, "ni")
  newData = merge(data.frame(id = data[,idLabel]), freq, by = c(idLabel))
  return(newData$ni)
}

coeff.Functions.Estimation = function(p, LgridA, gridA, LgridB = NULL, gridB = NULL, degrFreedom, degree, intercept, b, knots){
  
  Beta = matrix(0, nrow = p+1, ncol = LgridA)
  if(!is.null(LgridB)){Betap = matrix(0, nrow = p+1, ncol = LgridB)}
  
  BB = bs(gridA, df = degrFreedom, intercept = intercept, degree = degree, knots = knots)
  if(!is.null(LgridB)){BBp = bs(gridB, df = degrFreedom, intercept = intercept, degree = degree, knots = knots)}
  for (ib in 1:(p+1)){
    ind = (ib*degrFreedom-(degrFreedom-1)):(ib*degrFreedom)
    Beta[ib,] = BB %*%b[ind,"fit.coefficients"]
    if(!is.null(LgridB)){Betap[ib,] = BBp %*%b[ind,"fit.coefficients"]}
  }
  rownames(Beta) = 1:(p+1)
  colnames(Beta) = gridA
  if(!is.null(LgridB)){
    rownames(Betap) = 1:(p+1)
    colnames(Betap) = gridB
  }
  
  if(!is.null(LgridB)){return(list(Beta = Beta, Betap = Betap))}else{return(Beta)}
}

screening = function(i,BasisMtrx, X, y, Ggroup, luniq_id, uniq_id, lut.var, ut.var,p,
                     Lgrid, t.varGrid, degrFreedom, intercept, degree, knots, t.var, weights){
  
  require(splines)
  
  B0 = bs(t.var, df =  degrFreedom ,intercept = intercept, degree = degree, knots = knots)
  
  BasisMtrx2 = cbind(BasisMtrx, apply(B0, 2, function(x){x*X[,i]}))
  
  fit = lm(y~BasisMtrx2[,-1], weights = weights)
  
  b = as.matrix(fit$coefficients)
  b[2:degrFreedom,] = b[2:degrFreedom,]+b[1,] 
  rownames(b) = Ggroup
  colnames(b) = "fit.coefficients"
  
  Beta = coeff.Functions.Estimation(p = p+1, LgridA = Lgrid, gridA = t.varGrid, LgridB = NULL, gridB = NULL, degrFreedom = degrFreedom,
                                     degree = degree, intercept = intercept, b = b, knots = knots)
  
  
  fitted_valTmp = t(apply(BasisMtrx2, 1, function(x){ x*t(b)}))
  fitted_values = rowSums(fitted_valTmp) 
  residualsTmp = y-fitted_values
  
  #print(paste(i))
  
  return((t(residualsTmp)*weights)%*%residualsTmp/luniq_id)
}

