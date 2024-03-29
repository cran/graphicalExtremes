

#' Number of cores to be used in parallel computations
#' 
#' Helper function that returns the number of cores to be used in parallel computations.
#' Will always be 1 on Windows. On other systems, this value can be set using
#' `setOption('graphicalExtremes.mc.cores', ...)`.
#' 
#' @param overwrite Use this value (if it is valid and not on Windows)
#' @return An integer to be used as number of cores
#' 
#' @seealso [`graphicalExtremes-package`]
#' @family getOptions
#' @export
get_mc_cores <- function(overwrite = NULL){
  # Always 1 on Windows
  if(.Platform$OS.typ == 'windows'){
    return(1L)
  }
  # Use overwrite if specified
  if(is.numeric(overwrite) && length(overwrite) >= 1){
    nc_overwrite <- fitInInterval(overwrite[1], 1, Inf)
    return(as.integer(nc_overwrite))
  }
  # Try to use package option
  nc_option <- getOption('graphicalExtremes.mc.cores')
  if(is.numeric(nc_option) && length(nc_option) >= 1){
    nc_option <- fitInInterval(nc_option[1], 1, Inf)
    return(as.integer(nc_option))
  }
  # Try to detect
  nc_detected <- parallel::detectCores()
  if(!is.na(nc_detected)){
    return(nc_detected)
  }
  # Fall back to 1
  return(1L)
}


#' Tolerances to be used in computations
#' 
#' Helper function that returns the tolerance to be used in internal computations.
#' 
#' There are two different tolerances used in the package, for details see
#' [`graphicalExtremes-package`]. The default values for these tolerances can be
#' set using the options `"graphicalExtremes.tol.small"` and
#' `"graphicalExtremes.tol.large"`.
#' 
#' @param overwrite `NULL` or numeric scalar. If specified, use this value
#' instead of the option value.
#' @return A non-negative numerical scalar
#' 
#' @rdname get_tol
#' @seealso [`graphicalExtremes-package`]
#' @family getOptions
#' @export
get_small_tol <- function(overwrite = NULL){
  # Use overwrite if specified
  if(is.numeric(overwrite) && length(overwrite) >= 1){
    tol_overwrite <- fitInInterval(overwrite[1], 0, Inf)
    return(tol_overwrite)
  }
  # Try package option
  tol_option <- getOption('graphicalExtremes.tol.small')
  if(is.numeric(tol_option) && length(tol_option) >= 1){
    tol_option <- fitInInterval(tol_option[1], 0, Inf)
    return(tol_option)
  }
  # Fall back to some default value
  return(1e-12)
}


#' @rdname get_tol
#' @export
get_large_tol <- function(overwrite = NULL){
  # Use overwrite if specified
  if(is.numeric(overwrite) && length(overwrite) >= 1){
    tol_overwrite <- fitInInterval(overwrite[1], 0, Inf)
    return(tol_overwrite)
  }
  # Try package option
  tol_option <- getOption('graphicalExtremes.tol.large')
  if(is.numeric(tol_option) && length(tol_option) >= 1){
    tol_option <- fitInInterval(tol_option[1], 0, Inf)
    return(tol_option)
  }
  # Fall back to some default value
  return(1e-9)
}



#' Get alert function
#' 
#' Get a function that can be used to alert the user of invalid inputs.
#' Returns the value implied by the `overwrite` argument,
#' or the option `"graphicalExtremes.default.alert"`,
#' falling back to `warning()` if neither is specified.
#' 
#' @param overwrite `NULL` or `TRUE` to read the option value,
#' `FALSE` to return a dummy function,
#' or a function that takes an arbitrary number of strings as arguments (e.g. `stop()`).
#' 
#' @return A function that takes an arbitrary number of strings as arguments.
#' 
#' @seealso [`graphicalExtremes-package`]
#' @family getOptions
#' @export
get_alert_function <- function(overwrite = NULL){
  OPTION_NAME_DEFAULT_ALERT <- 'graphicalExtremes.default.alert'

  # Handle overwrite
  if(is.null(overwrite) || isTRUE(overwrite)){
    # read from options below
  }else if(is.function(overwrite)){
    return(overwrite)
  }else if(isFALSE(overwrite)){
    return(ignore)
  }else {
    stop('`overwrite` must be a function, boolean, or NULL')
  }

  # Get from options
  alert_option <- getOption(OPTION_NAME_DEFAULT_ALERT)
  if(is.function(alert_option)){
    # return as is below
  }else if(is.null(alert_option) || isTRUE(alert_option)){
    alert_option <- warning
  }else if(isFALSE(overwrite)){
    alert_option <- ignore
  }else{
    stop('option "', OPTION_NAME_DEFAULT_ALERT, '" must a function, boolean, or NULL')
  }

  # Return
  return(alert_option)
}

# Helper function that ignores all input arguments
ignore <- function(...){invisible()}
