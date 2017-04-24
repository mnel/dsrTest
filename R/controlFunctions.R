#' @title Control Function for Dobson Method Confidence Intervals
#' 
#' @description Provides a list of arguments to pass to 
#' \code{\link[exactci]{poisson.exact}}.
#' @seealso \code{\link[exactci]{poisson.exact}}
#' @return a list with values 
#' - `midp`
#' - `tsmethod`
#' @param midp logical, use mid-p values? Currently only permitted
#' where `tsmethod = "central"`.
#' @param tsmethod `character` giving two-sided method
#' @export
dobsonControl <- function(midp = FALSE, 
                          tsmethod = c("central","minlike","blaker")){
  tsmethod <- match.arg(tsmethod)
  list(midp = midp, tsmethod = tsmethod)
}

#' @title Control Function for Asymptotic Method Confidence Intervals
#' 
#' @description Specify the transformation to apply to the distribution
#' of the MLE.
#'
#' @param trans Transformation apply to the MLE distribution. 
#'
#' @return
#' A list with values
#' - `trans`
#' @export
asymptoticControl <- function(trans = c("none","log")){
  trans <- match.arg(trans)
  list(trans = trans)
}


#' @title Control Function for Gamma Method Confidence Intervals
#' @description Provides a list of arguments to pass to 
#' \code{\link[asht]{wspoissonTest}}.
#' @seealso \code{\link[asht]{wspoissonTest}}
#' @param midp logical. Use mid-p confidence distribution method? Currently
#' only implemented where `wmtype = "max"``
#' @param nmc Calculation method when `midp=TRUE`.
#' @param wmtype type of modification for the Gamma confidence interval.
#' @param unirootTolFactor tolerance factor for uniroot where `midp = TRUE`
#' and `nmc = 0`.
#' @return
#' A list of arguments to pass to \code{\link[asht]{wspoissonTest}}.
#' If `midp = TRUE`, with values
#' - `midp`
#' - `nmc`
#' - `unirootTolFactor`
#' 
#' If `midp = FALSE`, with values
#' - `wmtype`
#' @export
gammaControl <- function(midp = FALSE, nmc = 0,  
                         wmtype = c("max", "mean", "minmaxavg", "tcz"), 
                         unirootTolFactor=10^(-6)){
  wmtype <- match.arg(wmtype)
  if(midp) {
    RVAL <- list(midp = midp, nmc = nmc, 
                 unirootTolFactor = unirootTolFactor)
  } else {
      RVAL = list(wmtype = wmtype)
    }
  RVAL
}

#' Control Function for Beta Method for Confidence Intervals
#'
#' @description In the future, this may provide arguments
#' that will permit modifications to the Beta method.
#'
#' @param ... Currently ignored
#'
#' @return an empty list
#'
#' @export
betaControl <- function(...) list()
  