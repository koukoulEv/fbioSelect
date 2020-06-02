#' Coefficinet function estimation
#'
#' Takes a design matrix and converts it to the corresponding basis function approximation.
#' @param p number of parameters in the data.
#' @param LgridA number of grid points on which the coefficient functions will be estimated.
#' @param gridA the grid points on which the coefficient functions will be estimated.
#' @param degrFreedom the degrees of freedom for the polynomial basis used. See \code{\link{bs}} documentation.
#' @param degree the degree of the polynomial basis used, default is 3 which corresponds to cubic splines. See \code{\link{bs}} documentation.
#' @param intercept either TRUE or FALSE. See \code{\link{bs}} documentation.
#' @param b the estimated coefficients' vector.
#' @param knots the position of interior knots used for the basis expansion, if NULL then the median is used. See \code{\link{bs}} documentation.
#' @return The estimated coefficient functions.
#' @export
coeffFunctionEstimation = function(p, LgridA, gridA, degrFreedom, degree, intercept, b, knots){
  if (requireNamespace("splines", quietly = TRUE)) {
    splines::bs()
  }
  Beta = matrix(0, nrow = p+1, ncol = LgridA)

  BB = bs(gridA, df = degrFreedom, intercept = intercept, degree = degree, knots = knots)
   for (ib in 1:(p+1)){
    ind = (ib*degrFreedom-(degrFreedom-1)):(ib*degrFreedom)
    Beta[ib,] = BB %*%b[ind,"fit.coefficients"]
    }
  rownames(Beta) = 1:(p+1)
  colnames(Beta) = gridA

  return(Beta)
}
