#' @title Control Function for Dobson Method Confidence Intervals
#' 
#' @description Provides a list of arguments to pass to 
#' \code{\link[exactci]{poisson.exact}}.
#' @seealso \code{\link[exactci]{poisson.exact}}
#' @return a list with values 
#' \item{\code{midp}}{}
#' \item{\code{tsmethod}}{}
#' @param midp logical, use mid-p values? Currently only permitted
#' where `tsmethod = "central"`.
#' @param tsmethod `character` giving two-sided method
#' @param ... Currently ignored..
#' @export
dobsonControl <- function(midp = FALSE,
                          tsmethod = c("central", "minlike", "blaker"),
                          ...){
  tsmethod <- match.arg(tsmethod)
  list(midp = midp, tsmethod = tsmethod)
}

#' @title Control Function for Asymptotic Method Confidence Intervals
#' 
#' @description Specify the transformation to apply to the distribution
#' of the MLE.
#'
#' @param trans Transformation apply to the MLE distribution. 
#' @param ... Currently ignored.
#' @return
#' A list with values
#' \item{\code{trans}}{}
#' @export
asymptoticControl <- function(trans = c("none", "log", "loglog", "logit"), ...){
  trans <- match.arg(trans)
  list(trans = trans)
}

#' @title Control Function for Gamma Method Confidence Intervals
#' @description Provides a list of arguments to pass to 
#' \code{\link[asht]{wspoissonTest}}.
#' @seealso \code{\link[asht]{wspoissonTest}}
#' @param midp logical. Use mid-p confidence distribution method? Currently
#' only implemented where `wmtype = "max"`
#' @param nmc Calculation method when `midp = TRUE`.
#' @param wmtype type of modification for the Gamma confidence interval.
#' @param unirootTolFactor tolerance factor for uniroot where `midp = TRUE`
#' and `nmc = 0`.
#' @param ... Currently ignored.
#' @return
#' A list of arguments to pass to \code{\link[asht]{wspoissonTest}}.
#' If `midp = TRUE`, with values
#' \item{\code{midp}}{}
#' \item{\code{nmc}}{}
#' \item{\code{unirootTolFactor}}{}
#' 
#' If `midp = FALSE`, with values
#' \item{\code{wmtype}}{}
#' @export
gammaControl <- function(midp = FALSE, nmc = 0,
                         wmtype = c("max", "mean", "minmaxavg", "tcz"),
                         unirootTolFactor = 1e-06, ...){
  wmtype <- match.arg(wmtype)
  if (midp) {
    RVAL <- list(midp = midp, nmc = nmc,
                 unirootTolFactor = unirootTolFactor)
  } else {
      RVAL <- list(wmtype = wmtype)
    }
  RVAL
}

#' Control Function for Beta Method for Confidence Intervals
#'
#' @description Modification to the Beta method. The options are `"none"`
#' or the same modifications as applied to the Gamma Method 
#' (see \code{\link{gammaControl}}) are implemented. `wmtype="none"` and
#' `wmtype="tcz"` have been investigated by Tiwari et al (2006) and
#' Ng et al (2008).
#' @references 
#' Tiwari RC, Clegg LX, & Zou Z (2006). 'Efficient interval estimation 
#' for age-adjusted cancer rates.' 
#' *Statistical Methods in Medical Research* **15**: 547--569. 
#' \doi{doi:10.1177/0962280206070621}
#' 
#' Ng HKT, Filardo, G & Zheng G (2008). 'Confidence interval estimating 
#' procedures for standardized incidence rates.' 
#' *Computational Statistics and Data Analysis* **52** 3501--3516.
#' \doi{doi:10.1016/j.csda.2007.11.004}
#' @param wmtype `character` type of modification to the Beta Confidence
#' Interval
#' @param ... Currently ignored.
#'
#' @return a list with values
#' \item{wmtype}{modificator to Beta Confidence Interval to implement}
#'
#' @export
betaControl <- function(wmtype = c("none", "tcz", "max", "mean",
                                   "minmaxavg"), ...) {
  wmtype <- match.arg(wmtype)
  list(wmtype = wmtype)
}
