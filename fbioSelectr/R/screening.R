#' Variable screening
#'
#' Marginally fits a model including one of the high-dimensional covariates and the low-dimensional covariate matrix.
#' @param i iteration index.
#' @param BasisMtrx  the basis approximation matrix of the low-dimensional covariate design matrix.
#' @param X a dataframe or matrix of high-dimensional covariates subject to variable selection.
#' @param y a response vector where only continuous responses are acceptable. Responses are expected to be time or dose varying.
#' @param Ggroup the groups created from the basis expansion in the BasisMtrx matrix.
#' @param luniq_id length of unique ids in the data.
#' @param uniq_id unique ids in the data.
#' @param lut.var length of unique time or dose points.
#' @param ut.var the unique time or dose points.
#' @param p number of parameters in the data.
#' @param Lgrid number of grid points on which the coefficient functions are estimated.
#' @param t.varGrid time of dose grid point on which the coefficient functions are estimated.
#' @param degrFreedom the degrees of freedom for the polynomial basis used. See \code{\link{bs}} documentation.
#' @param degree the degree of the polynomial basis used, default is 3 which corresponds to cubic splines.See \code{\link{bs}} documentation.
#' @param intercept either TRUE or FALSE. See \code{\link{bs}} documentation.
#' @param knots the position of interior knots used for the basis expansion, if NULL then the median is used. See \code{\link{bs}} documentation.
#' @param t.var  the varying time or dose under which the responses were recorded.
#' @param weights regression weights.
#' @return Te weighted mean square error of model fit.
#' @export
screening = function(i,BasisMtrx, X, y, Ggroup, luniq_id, uniq_id, lut.var, ut.var,p,Lgrid, t.varGrid, degrFreedom, intercept, degree, knots, t.var, weights){

  if (requireNamespace("splines", quietly = TRUE)) {
    splines::bs()
  }

  B0 = bs(t.var, df =  degrFreedom ,intercept = intercept, degree = degree, knots = knots)
  BasisMtrx2 = cbind(BasisMtrx, apply(B0, 2, function(x){x*X[,i]}))

  fit = lm(y~BasisMtrx2[,-1], weights = weights)

  b = as.matrix(fit$coefficients)
  b[2:degrFreedom,] = b[2:degrFreedom,]+b[1,]
  rownames(b) = Ggroup
  colnames(b) = "fit.coefficients"

  Beta = coeffFunctionEstimation(p = p+1, LgridA = Lgrid, gridA = t.varGrid, degrFreedom = degrFreedom,degree = degree, intercept = intercept, b = b, knots = knots)

  fitted_valTmp = t(apply(BasisMtrx2, 1, function(x){ x*t(b)}))
  fitted_values = rowSums(fitted_valTmp)
  residualsTmp = y-fitted_values

  return((t(residualsTmp)*weights)%*%residualsTmp/luniq_id)
}

