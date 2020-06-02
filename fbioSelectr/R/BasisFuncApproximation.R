#' Basis Function Approximation
#'
#' Takes a design matrix and converts it to the corresponding basis function approximation.
#' @param designMtrx the design matrix of time-invariant covariates.
#' @param grid a vector of the grid points for the basis function approximation (length equal to the number of rows of designMtrx).
#' @param degrFreedom the degrees of freedom for the polynomial basis used. See \code{\link{bs}} documentation.
#' @param degree the degree of the polynomial basis used, default is 3 which corresponds to cubic splines.
#' @param intercept either TRUE or FALSE. See bs documentation.
#' @param knots the position of interior knots used for the basis expansion, if NULL then the median is used. See bs documentation.
#' @return The basis approximation matrix of the imported design matrix.
#' @export
BasisFuncApproximation = function(designMtrx, grid, degrFreedom, degree, intercept, knots){

  if (requireNamespace("foreach", quietly = TRUE)) {
    foreach::foreach()
  }
  if (requireNamespace("iterators", quietly = TRUE)) {
    iterators::iter()
  }
  if (requireNamespace("splines", quietly = TRUE)) {
    splines::bs()
  }

  B0 = bs(grid, df =  degrFreedom ,intercept = intercept, degree = degree, knots = knots)
  modsplinesA = foreach(a = iter(as.matrix(designMtrx[,-1]), by = "col"),  .combine = "cbind") %do% apply(B0, 2, function(x){x*a})
  modsplines = cbind(B0, modsplinesA)

  return(modsplines)
}
