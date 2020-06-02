#' Two-stage variable selection algorithm
#'
#' From a high-dimensional set of time-invariant covariates, selects those which are associated to the functional (repeated measures) response given some how-dimensional covariates.
#' @param y a response vector where only continuous responses are acceptable. Responses are expected to be time or dose varying.
#' @param id a vector including the identification codes for the response vector supplied.
#' @param t.var the varying time or dose under which the responses were recorded.
#' @param Z a matrix of low-dimensional covariates (should be an object obtained from the \code{\link{model.matrix}} function).
#' @param X a dataframe or matrix of high-dimensional covariates subject to variable selection.
#' @param low.dim.cov.number number of low-dimensional covariates, minimum is 1.
#' @param threshold.size threshold used for the variable screening step of the algorithm.
#' @param type.of.penalty type of penalty used for the second step of the algorithm (group penalised regression). Options of penalties can be found in the documentation of the \code{\link{grpreg}} documentation.
#' @param degree the degree of the polynomial basis used, default is 3 which corresponds to cubic splines.
#' @param intercept either TRUE or FALSE. See bs documentation. See \code{\link{bs}} documentation.
#' @param degrFreedom the degrees of freedom for the polynomial basis used. See \code{\link{bs}} documentation.
#' @param knots the position of interior knots used for the basis expansion, if NULL then the median is used. See \code{\link{bs}} documentation.
#' @param weights regression weights, if NULL then independence covariance structure is used for the repeated measures data.
#' @param nfolds number of folds for the cross-validation selecting the penalty size for the second step of the algorithm. See \code{\link{grpreg}} documentation.
#' @param seed the seed of R's random number generator. Default is 1000.
#' @param ncores number of cores used for parallel processing. Default is 3.
#' @return A list including the ranking of the high dimensional covariates based on the marginal fit (variable screening step), the high-dimensional covariates that have been selected from the second step (group penalised regression), the estimated coefficient functions, the regression weights used and the group penalised regression path (\code{\link{grpreg}} output).
#' @export
fbioSelect = function(y,id,t.var,Z,X,low.dim.cov.number=1,threshold.size,type.of.penalty,degree=3,intercept=FALSE,degrFreedom=NULL,knots = NULL,weights=NULL,nfolds=10,seed=1000, ncores = 3){

  if (requireNamespace("dplyr", quietly = TRUE)) {
    dplyr::arrange()
    dplyr::`%>%`
    dplyr::filter()
  }
  if (requireNamespace("foreach", quietly = TRUE)) {
    foreach::foreach()
    foreach::`%do%`
    foreach::`%dopar%`
  }
  if (requireNamespace("grpreg", quietly = TRUE)) {
    grpreg::grpreg()
    grpreg::cv.grpreg()
  }
  if (requireNamespace("doParallel", quietly = TRUE)) {
    doParallel::registerDoParallel()
  }

  registerDoParallel(cores = ncores)
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
  BasisMtrx = BasisFuncApproximation(designMtrx = Z, grid = t.var, degrFreedom = degrFreedom,
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
  geneRanking = data.frame(hdcov = high_dim_cov_names, score = mse_screening) %>% arrange("score") # rank the high dimensional covariates according to their utility
  ## Select the top ranked high dimensional covariates according to the threshold applied and the ranking based on their utility
  outranked_hdcov = geneRanking[1:threshold.size,]

  ## Update the design matrix based on the variable screening step
  Z2 = cbind(Z, X[, as.character(outranked_hdcov$hdcov)])
  ## Update design matrix by adding the top ranked high dimensional covariates
  BasisMtrx2 = BasisFuncApproximation(designMtrx =  Z2, grid = t.var, degrFreedom = degrFreedom, degree = degree, intercept = intercept, knots = knots)

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

  Beta = coeffFunctionEstimation(p = nrow(outranked_hdcov), LgridA = Lgrid, gridA = Interval, degrFreedom = degrFreedom, degree = degree, intercept = intercept, b = b, knots = knots)
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
