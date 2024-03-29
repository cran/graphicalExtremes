

#' Fit value(s) in interval
#' 
#' Fit value(s) in interval, all arguments are recycled where necessary.
#' 
#' @param x Numeric vector
#' @param xMin Numeric vector
#' @param xMax Numeric vector
#' @return Numeric vector
#' @keywords internal
fitInInterval <- function(x, xMin=-Inf, xMax=Inf){
  if(any(xMax<xMin)){
    stop('Make sure that xMax>=xMin!')
  }
  x <- pmax(x, xMin)
  x <- pmin(x, xMax)
  return(x)
}


# Replaces a principal submatrix of a Gamma matrix, preserving definiteness.
# Other entries are kept "heuristically similar".
# Zeros (instead of NAs) on the diagonal of `Gamma.fix` indicate which principal submatrix to replace.
# If not the entire submatrix in `Gamma.fix` is specified, [complete_Gamma_decomposable()] is used
replaceGammaSubMatrix <- function(Gamma.est, Gamma.fix){
  # check which are fixed
  ind <- which(!is.na(diag(Gamma.fix)))
  if(length(ind)==0){
    return(Gamma.est)
  }
  if(any(is.na(Gamma.fix[ind, ind]))){
    Gamma.fix[ind, ind] <- complete_Gamma(Gamma.fix[ind, ind], allowed_graph_type = 'decomposable')
  }
  if(length(ind) == ncol(Gamma.est)){
    return(Gamma.fix)
  }
  
  # naive attempt:
  G1 <- Gamma.est
  G1[ind, ind] <- Gamma.fix[ind, ind]
  if(is_valid_Gamma(G1)){
    return(G1)
  }

  # heuristic, but safe solution:
  k <- ind[1]
  M.est <- Gamma2Sigma(Gamma.est, k, check = FALSE)
  M.fix <- Gamma2Sigma(Gamma.fix, k, check = FALSE)
  M <- replaceSpdSubmatrix(M.est, M.fix)
  G <- Sigma2Gamma(M, k=k, check = FALSE)
  return(G)
}

# replaces part of a positive definite matrix, keeping definiteness
replaceSpdSubmatrix <- function(M.est, M.fix){
  indFix <- !is.na(diag(M.fix))
  indMod <- !indFix

  indFix <- which(indFix)
  indMod <- which(indMod)

  if(length(indFix) == 0){
    return(M.est)
  } else if(length(indMod) == 0){
    return(M.fix)
  }

  ind <- c(indFix, indMod)

  M2 <- M.est[ind, ind]

  indFix2 <- seq_along(indFix)
  indMod2 <- seq_along(indMod) + length(indFix2)

  C <- M.fix[indFix, indFix]

  L <- chol(M2)
  LC <- chol(C)

  L[indFix2, indFix2] <- LC

  M2 <- t(L) %*% L

  M <- matrix(0, NROW(M2), NROW(M2))

  M[ind, ind] <- M2

  return(M)
}

rdunif <- function(n, a, b){
  a + floor((b - a + 1) * stats::runif(n))
}

pdet <- function(M, tol=get_small_tol()){
  ev <- eigen(M, only.values = TRUE)$values
  prod(ev[abs(ev) > tol])
}


upper.tri.val <- function(M, diag=FALSE){
  M[upper.tri(M, diag)]
}


is_eq <- function(a, b, tol=get_small_tol()) {
  abs(a - b) < tol
}
is_greater <- function(a, b, tol=get_small_tol()) {
  a - b > tol
}
is_less <- function(a, b, tol=get_small_tol()) {
  is_greater(b, a, tol)
}
is_leq <- function(a, b, tol=get_small_tol()) {
  !is_greater(a, b, tol)
}
is_geq <- function(a, b, tol=get_small_tol()) {
  !is_less(a, b, tol)
}


makeUnitVector <- function(d, k){
    ek <- numeric(d)
    ek[k] <- 1
    return(ek)
}
makeOneVec <- function(d){
    v <- rep(1, d)
}
makeProjK <- function(d, k){
    ek <- makeUnitVector(d, k)
    oneVec <- rep(1, d)
    diag(d) - ek %*% t(oneVec)
}



#' Convert indices to numerical indices
#' 
#' Converts (possibly) logical indices to numerical ones.
#' Also ensures unique indices and sorts them if specified.
#' 
#' @param ind The numerical or logical index vector
#' @param n Max numerical index (used if `ind` is logical and might be recycled)
#' @param unique Whether to keep every (numerical) index at most once
#' @param sort Whether to sort the numerical indices
#' 
#' @return A numerical index vector
#' 
#' @keywords internal
make_numeric_indices <- function(ind, n=NULL, unique=TRUE, sort=TRUE){
  if(is.logical(ind)){
    if(is.null(n)){
      n <- length(ind)
    }
    ind0 <- seq_len(n)
    ind <- ind0[ind]
  }
  ind <- as.integer(ind)
  if(unique){
    ind <- unique(ind)
  }
  if(sort){
    ind <- sort(ind)
  }
  return(ind)
}

