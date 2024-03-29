

#' Parameter fitting for Huesler-Reiss graphical models
#' 
#' Fits the parameter matrix (variogram) of a multivariate Huesler-Reiss Pareto distribution
#' with a given graphical structure, using maximum-likelihood estimation
#' or the empirical variogram.
#' 
#' @param data Numeric \nxd matrix, where `n` is the
#' number of observations and `d` is the number of dimensions.
#'
#' @param graph Undirected, connected \[`igraph::graph`\] object with `d` vertices,
#' representing the graphical structure of the fitted Huesler-Reiss model.
#'
#' @param p Numeric between 0 and 1 or `NULL`. If `NULL` (default),
#' it is assumed that the `data` is already on a multivariate Pareto scale.
#' Else, `p` is used as the probability in the function [data2mpareto()]
#' to standardize the `data`.
#'
#' @param method One of `c('vario', 'ML')`, with `'vario'` as default, indicating
#' the method to be used for parameter estimation. See details.
#' 
#' @param handleCliques How to handle cliques and separators in the graph.
#' See details.
#' 
#' @param ... Arguments passed to [fmpareto_HR_MLE()]. Currently `cens`, `maxit`,
#' `optMethod`, and `useTheta` are supported.
#' 
#' @return The estimated parameter matrix.
#' 
#' @details
#' If `handleCliques='average'`, the marginal parameter matrix is estimated for
#' each maximal clique of the `graph` and then combined into a partial parameter
#' matrix by taking the average of entries from overlapping cliques. Lastly,
#' the full parameter matrix is computed using [complete_Gamma()].
#' 
#' If `handleCliques='full'`, first the full parameter matrix is estimated using the
#' specified `method` and then the non-edge entries are adjusted such that the
#' final parameter matrix has the graphical structure indicated by `graph`.
#' 
#' If `handleCliques='sequential'`, `graph` must be decomposable, and
#' `method='ML'` must be specified. The parameter matrix is first estimated on
#' the (recursive) separators and then on the rest of the cliques, keeping
#' previously estimated entries fixed.
#' 
#' If `method='ML'`, the computational cost is mostly influenced by the total size
#' of the graph (if `handleCliques='full'`) or the size of the cliques,
#' and can already take a significant amount of time for modest dimensions (e.g. `d=3`).
#' 
#' @family parameterEstimation
#' @export
fmpareto_graph_HR <- function(
  data,
  graph,
  p = NULL,
  method = c('vario', 'ML'),
  handleCliques = c( 'average','full', 'sequential'),
  ...
){
  # match args
  method <- match.arg(method)
  handleCliques <- match.arg(handleCliques)

  # standardize data
  if(!is.null(p)){
    data <- data2mpareto(data, p)
  }
  
  # call other functions, depending on handleCliques and method
  if(handleCliques == 'full'){
    # Estimate the entire parameter matrix at once
    if(method == 'ML'){
      Gamma <- fmpareto_HR_MLE(
        data,
        graph = graph,
        ...
      )
    } else{ # method == 'vario'
      variogram <- emp_vario(data)
      Gamma <- complete_Gamma(variogram, graph)
    }
  } else if(handleCliques == 'sequential'){
    # Estimate one clique after the other, fixing previously estimated entries.
    # Works only with MLE on decomposable graphs
    if(method != 'ML' || !is_decomposable_graph(graph)){
      stop('Arguments handleCliques="sequential" only works with decomposable graphs and method="ML"!')
    }
    Gamma <- fmpareto_graph_HR_clique_sequential(
      data,
      graph,
      ...
    )
  } else { # handleCliques == 'average'
    # Estimate each clique separately, taking the average on separators
    Gamma <- fmpareto_graph_HR_clique_average(
      data,
      graph,
      method,
      ...
    )
    if(!is_valid_Gamma(Gamma)){
      stop('Averaging on the separators did not yield a valid variogram matrix!')
    }
  }
  return(Gamma)
}

#' HR Parameter fitting - Helper functions
#' 
#' Helper functions called by [`fmpareto_HR_MLE`].
#' 
#' @rdname fmpareto_HR_helpers
#' @keywords internal
fmpareto_graph_HR_clique_average <- function(
  data,
  graph,
  method = c('ML', 'vario'),
  ...
){
  method <- match.arg(method)
  
  graph <- check_graph(graph)
  cliques <- igraph::max_cliques(graph)
  if(method == 'vario'){
    subGammas <- parallel::mclapply(
      mc.cores = get_mc_cores(),
      cliques,
      function(cli){
        data.cli <- mparetomargins(data, cli)
        emp_vario(data.cli)
      }
    )
  } else {
    subGammas <- parallel::mclapply(
      mc.cores = get_mc_cores(),
      cliques,
      function(cli){
        data.cli <- mparetomargins(data, cli)
        tmp <- fmpareto_HR_MLE(
          data.cli,
          ...
        )
        if(is.null(tmp$Gamma)){
          stop('MLE did not find a valid Gamma for clique: ', paste0(cli, collapse = ','))
        }
        return(tmp$Gamma)
      }
    )
  }
  Gamma_partial <- combine_clique_estimates_by_averaging(cliques, subGammas)
  Gamma <- complete_Gamma(Gamma_partial, graph)
  return(Gamma)
}

#' @rdname fmpareto_HR_helpers
#' @keywords internal
fmpareto_graph_HR_clique_sequential <- function(
  data,
  graph,
  ...
){
  # Check that not `useTheta=TRUE`
  argsMLE <- list(...)
  if(!is.null(argsMLE$useTheta) && argsMLE$useTheta){
    stop('Sequential handling of cliques only works when optimizing on Gamma-level!')
  }
  
  # check inputs
  d <- ncol(data)
  graph <- check_graph(graph, graph_type = 'decomposable', nVertices = d)
  
  # Compute cliques and (recursive) separators
  layers <- get_cliques_and_separators(graph, sortIntoLayers = TRUE)
  
  # Estimate entries corresponding to separators/cliques
  Ghat <- matrix(NA, d, d)
  for(cliques in layers){
    # Compute estimate for each new clique/separator
    subGammas <- parallel::mclapply(
      mc.cores = get_mc_cores(),
      cliques,
      function(cli){
        # get margins data
        data.cli <- mparetomargins(data, cli)
        
        # find (already) fixed entries
        G.cli <- Ghat[cli, cli]
        par.cli <- matrix2par(G.cli)
        fixParams.cli <- !is.na(par.cli)
        
        # get initial parameters that agree with the fixed ones (heuristic, close to empirical variogram):
        G0 <- emp_vario(data.cli)
        G1 <- replaceGammaSubMatrix(G0, G.cli)
        init.cli <- matrix2par(G1)
        
        # estimate parameters
        opt <- fmpareto_HR_MLE(
          data = data.cli,
          init = init.cli,
          fixParams = fixParams.cli,
          useTheta = FALSE,
          ...
        )
        if(is.null(opt$Gamma)){
          stop('MLE did not converge for clique: ', cli)
        }
        return(opt$Gamma)
      }
    )

    # Fill newly computed entries in Ghat
    for(i in seq_along(cliques)){
      cli <- cliques[[i]]
      Ghat[cli, cli] <- subGammas[[i]]
    }
  }
  
  # Complete non-edge entries
  G_comp <- complete_Gamma(Ghat, graph)
  return(G_comp)
}

#' @rdname fmpareto_HR_helpers
#' @keywords internal
combine_clique_estimates_by_averaging <- function(cliques, subGammas){
  d <- do.call(max, cliques)
  Gamma <- matrix(0, d, d)
  overlaps <- matrix(0, d, d)
  for(i in seq_along(cliques)){
    cli <- cliques[[i]]
    Gamma[cli, cli] <- Gamma[cli, cli] + subGammas[[i]]
    overlaps[cli, cli] <- overlaps[cli, cli] + 1
  }
  Gamma[overlaps == 0] <- NA
  Gamma <- Gamma / overlaps
  return(Gamma)
}



#' Estimation of the variogram matrix \eGamma of a Huesler-Reiss distribution
#'
#' Estimates the variogram of the Huesler-Reiss distribution empirically.
#'
#' @param data Numeric \nxd matrix, where `n` is the
#' number of observations and `d` is the dimension.
#' @param k Integer between 1 and `d`. Component of the multivariate
#' observations that is conditioned to be larger than the threshold `p`.
#' If `NULL` (default), then an average over all `k` is returned.
#' @param p Numeric between 0 and 1 or `NULL`. If `NULL` (default),
#' it is assumed that the `data` are already on multivariate Pareto scale. Else,
#' `p` is used as the probability in the function [data2mpareto()]
#' to standardize the `data`.
#' 
#' @details
#' `emp_vario_pairwise` calls `emp_vario` for each pair of observations.
#' This is more robust if the data contains many `NA`s, but can take rather long.
#'
#' @return Numeric \dxd matrix. The estimated variogram of the Huesler-Reiss distribution.
#' 
#' @examples
#' G <- generate_random_Gamma(d=5)
#' y <- rmpareto(n=100, par=G)
#' Ghat <- emp_vario(y)
#' 
#' @rdname emp_vario
#' @family parameterEstimation
#' @export
emp_vario <- function(data, k = NULL, p = NULL) {

  # helper ####
  G.fun <- function(i, data) {
    idx <- which(data[, i] > 1)
    if (length(idx) > 1) {
      xx <- Sigma2Gamma(stats::cov(log(data[idx, ])), full = TRUE, check = FALSE)
    } else {
      xx <- matrix(NA, d, d)
    }
    return(xx)
  }

  # body ####
  if (!is.matrix(data)) {
    stop("The data should be a matrix")
  }
  if (ncol(data) <= 1) {
    stop("The data should be a matrix with at least two columns.")
  }

  d <- ncol(data)
  if (!is.null(p)) {
    data <- data2mpareto(data, p)
  }

  if (!is.null(k)) {
    G <- G.fun(k, data)

    if (any(is.na(G))) {
      warning(paste(
        "Produced NA matrix since there are no exceedances in the component k =",
        k
      ))
    }
  } else {

    # take the average
    row_averages <- rowMeans(sapply(1:d, FUN = function(i) {
      G.fun(i, data)
    }), na.rm = TRUE)
    G <- matrix(row_averages, nrow = d, ncol = d)
  }

  return(G)
}


#' @param verbose Print verbose progress information
#' @rdname emp_vario
#' @export
emp_vario_pairwise <- function(data, k = NULL, p = NULL, verbose = FALSE){
  vcat <- function(...){
    if(verbose){
      cat(...)
    }
  }
  d <- ncol(data)
  vario <- matrix(0, d, d)
  vario[] <- NA
  pct <- -Inf
  vcat('Computing emp_vario_pairwise, reporting progress:\n')
  for(i in seq_len(d)){
    vario[i,i] <- 0
    if(i + 1 > d){
      next
    }
    for(j in (i+1):d){
      data_ij <- data2mpareto(data[,c(i,j)], p, na.rm=TRUE)
      vario_ij <- emp_vario(data_ij, p=NULL)
      vario[i,j] <- vario_ij[1,2]
      vario[j,i] <- vario_ij[1,2]
      pct1 <- floor(100 * sum(!is.na(vario)) / length(vario))
      if(pct1 > pct){
        pct <- pct1
        vcat(pct, '% ', sep='')
      }
    }
  }
  cat('\nDone.\n')
  return(vario)
}


#' Empirical estimation of extremal correlation matrix \eChi
#'
#' Estimates empirically the matrix of bivariate extremal correlation coefficients \eChi.
#'
#' @inheritParams emp_chi_multdim
#'
#' @return Numeric matrix \dxd. The matrix contains the
#' bivariate extremal coefficients \eqn{\chi_{ij}}, for \eqn{i, j = 1, ..., d}.
#' 
#' 
#' @details
#' `emp_chi_pairwise` calls `emp_chi` for each pair of observations.
#' This is more robust if the data contains many `NA`s, but can take rather long.
#'
#' @examples
#' n <- 100
#' d <- 4
#' p <- .8
#' Gamma <- cbind(
#'   c(0, 1.5, 1.5, 2),
#'   c(1.5, 0, 2, 1.5),
#'   c(1.5, 2, 0, 1.5),
#'   c(2, 1.5, 1.5, 0)
#' )
#'
#' set.seed(123)
#' my_data <- rmstable(n, "HR", d = d, par = Gamma)
#' emp_chi(my_data, p)
#' 
#' @rdname emp_chi
#' @family parameterEstimation
#' @export
emp_chi <- function(data, p = NULL) {
  if (!is.matrix(data)) {
    stop("The data should be a matrix")
  }
  if (ncol(data) <= 1) {
    stop("The data should be a matrix with at least two columns.")
  }

  if (!is.null(p)) {
    data <- data2mpareto(data, p)
  }

  n <- nrow(data)
  d <- ncol(data)

  ind <- data > 1
  ind_mat <- matrix(colSums(ind), byrow = TRUE, ncol = d, nrow = d)
  crossprod(ind, ind) / (1 / 2 * (ind_mat + t(ind_mat)))
}

#' @param verbose Print verbose progress information
#' @rdname emp_chi
#' @export
emp_chi_pairwise <- function(data, p = NULL, verbose=FALSE){
  vcat <- function(...){
    if(verbose){
      cat(...)
    }
  }
  d <- ncol(data)
  chi <- matrix(0, d, d)
  chi[] <- NA
  pct <- -Inf
  vcat('Computing emp_chi_pairwise, reporting progress:\n')
  for(i in seq_len(d)){
    chi[i,i] <- 1
    if(i + 1 > d){
      next
    }
    for(j in (i+1):d){
      data_ij <- data2mpareto(data[,c(i,j)], p, na.rm=TRUE)
      chi_ij <- emp_chi(data_ij, NULL)
      chi[i,j] <- chi_ij[1,2]
      chi[j,i] <- chi_ij[1,2]
      pct1 <- floor(100 * sum(!is.na(chi)) / length(chi))
      if(pct1 > pct){
        pct <- pct1
        vcat(pct, '% ', sep='')
      }
    }
  }
  cat('\nDone.\n')
  return(chi)
}



#' Empirical estimation of extremal correlation \eChi
#'
#' Estimates the `d`-dimensional extremal correlation coefficient \eChi empirically.
#'
#' @param data Numeric \nxd matrix, where `n` is the
#' number of observations and `d` is the dimension.
#' @param p Numeric scalar between 0 and 1 or `NULL`. If `NULL` (default),
#' it is assumed that the `data` are already on multivariate Pareto scale. Else,
#' `p` is used as the probability in [data2mpareto()] to standardize the `data`.
#'
#' @return Numeric scalar. The empirical `d`-dimensional extremal correlation coefficient \eChi
#' for the `data`.
#' @examples
#' n <- 100
#' d <- 2
#' p <- .8
#' G <- cbind(
#'   c(0, 1.5),
#'   c(1.5, 0)
#' )
#'
#' set.seed(123)
#' my_data <- rmstable(n, "HR", d = d, par = G)
#' emp_chi_multdim(my_data, p)
#' 
#' @family parameterEstimation
#' @export
emp_chi_multdim <- function(data, p = NULL) {
  if (!is.matrix(data)) {
    stop("The data should be a matrix")
  }
  if (ncol(data) <= 1) {
    stop("The data should be a matrix with at least two columns.")
  }

  d <- ncol(data)
  data <- stats::na.omit(data)

  if (!is.null(p)) {
    data <- data2mpareto(data, p)
  }



  rowmin <- apply(data, 1, min)
  chi <- mean(sapply(1:ncol(data), FUN = function(i) {
    mean(rowmin[which(data[, i] > 1)] > 1)
  }), na.rm = TRUE)
  return(chi)
}



#' Compute Huesler-Reiss log-likelihood, AIC, and BIC
#'
#' Computes (censored) Huesler-Reiss log-likelihood, AIC, and BIC values.
#'
#' @param data Numeric \nxd matrix. It contains
#' observations following a multivariate HR Pareto distribution.
#'
#' @param graph An \[`igraph::graph`\] object or `NULL`. The `graph` must be undirected and
#' connected. If no graph is specified, the complete graph is used.
#'
#' @param Gamma Numeric \nxd matrix.
#' It represents a variogram matrix \eGamma.
#'
#' @param cens Boolean. If true, then censored log-likelihood is computed.
#' By default, `cens = FALSE`.
#'
#' @param p Numeric between 0 and 1 or `NULL`. If `NULL` (default),
#' it is assumed that the `data` are already on multivariate Pareto scale.
#'  Else, `p` is used as the probability in the function [data2mpareto()]
#' to standardize the `data`.
#'
#' @return Numeric vector `c("loglik"=..., "aic"=..., "bic"=...)` with the evaluated
#' log-likelihood, AIC, and BIC values.
#'
#' @family parameterEstimation
#' @export
loglik_HR <- function(data, p = NULL, graph = NULL, Gamma, cens = FALSE){
  if (!is.null(p)) {
    data <- data2mpareto(data, p)
  }

  loglik <- logLH_HR(
    data = data,
    Gamma = Gamma,
    cens = cens
  )

  if(is.null(graph)){
    d <- ncol(data)
    n_edges <- d*(d-1)/2
  } else{
    n_edges <-  igraph::ecount(graph)
  }
  n <- nrow(data)

  aic <- 2 * n_edges - 2 * loglik

  bic <- log(n) * n_edges - 2 * loglik

  return(c("loglik" = loglik, "aic" = aic, "bic" = bic))
}

