#' @title Tests and Confidence Intervals on Directly Standardized Rates
#' @seealso \code{\link[asht]{wspoissonTest}}, 
#' \code{\link[exactci]{poisson.exact}}, 
#' \code{\link{gammaControl}}, 
#' \code{\link{dobsonControl}}, 
#' \code{\link{asymptoticControl}}, 
#' \code{\link{betaControl}}
#' @description A number of methods have been proposed for calculating
#' confidence intervals for directly standardized rates. Ng et al (2008),
#' compare a number of methods, some of which are implemented here. 
#' The default uses the Gamma method by Fay and Feuer (1997) and 
#' implemented in \code{\link[asht]{wspoissonTest}}. 
#' 
#' @details
#' Four classes of method have been implemented here:
#' 
#' \describe{
#' \item{\code{"gamma"}}{Calls \code{\link[asht]{wspoissonTest}}. By default
#' uses the Gamma Method proposed by Fay and Feuer (1997). 
#' Modifications proposed by Tiwari et al (2006) and Fay and Kim (2017) 
#' also implemented - see \code{\link{gammaControl}}.} 
#' \item{\code{"asymptotic"}}{Using the normal approximation of the 
#' MLE or transformed MLE distribition - see \code{\link{asymptoticControl}}}
#' \item{\code{"dobson"}}{Uses the method proposed by Dobson et al (1991). 
#' Estimating the confidence interval on the unweighted sum is done by calling
#' \code{\link[exactci]{poisson.exact}} - both the exact method and
#' a mid-p method are possible - see \link{dobsonControl}.}
#' \item{\code{"beta"}}{Methods based on the beta distribution by Tiwari et
#' al (2006) - see \link{betaControl}.}
#' \item{\code{"bootstrap"}}{Approximate Bootstrap method by Swift (1995). 
#' P-values are estimated by solving for p.}
#' }
#' 
#' For each method there is a `control` function that will return a list of
#' parameters that can be used to define sub-types of each of the broad groups
#' @references 
#' Dobson, AJ, Kuulasmaa, K, Eberle, E and Scherer, J (1991) 
#' 'Confidence intervals for weighted sums of Poisson parameters', 
#' *Statistics in Medicine*, **10**: 457--462.
#' \doi{doi:10.1002/sim.4780100317}
#' 
#' Swift, MB (1995). 'Simple confidence intervals for 
#' standardized rates based on the approximate bootstrap method', 
#' *Statistics in Medicine*, **14**, 1875--1888.
#' \doi{doi:10.1002/sim.4780141704}. 
#' 
#' Fay MP & Feuer EJ (1997). 'Confidence intervals for directly 
#' standardized rates: a method based on the gamma distribution.
#' Statistics in Medicine*. **16**: 791--801.
#' \url{https://doi.org/10.1002/(SICI)1097-0258(19970415)16:7<791::AID-SIM500>3.0.CO;2-\%23}
#' 
#' Tiwari RC, Clegg LX, & Zou Z (2006). 'Efficient interval estimation 
#' for age-adjusted cancer rates.' 
#' *Statistical Methods in Medical Research* **15**: 547--569. 
#' \doi{doi:10.1177/0962280206070621}
#' 
#' Ng HKT, Filardo, G & Zheng G (2008). 'Confidence interval estimating 
#' procedures for standardized incidence rates.' 
#' *Computational Statistics and Data Analysis* **52** 3501--3516.
#' \doi{doi:10.1016/j.csda.2007.11.004}
#' @param x a vector of strata-specific counts. 
#' @param n a vector of strata-specific time bases for counts.
#' @param w a vector of strata-specific weights (or standard populations).
#' @param null.value a null hypothesis value of the directly 
#'  rate, if NULL no test is done. If not NULL, provide in rate per mult.
#' @param alternative type of alternative hypothesis.
#' @param conf.level confidence level for the returned 
#' confidence interval.
#' @param mult a factor to multiply the estimate and 
#' confidence intervals by, to give rates per mult.
#' @param method Method used to perform the test and construct the 
#' confidence interval. See details.
#' @param control list of arguments / type of modification used for
#' each method. See details and relevant `"xxxxControl"` documentation
#' @return a list with class `"htest"` containing the following 
#' components: 
#' \item{statistic}{number of strata or summands: 
#' \code{k = length(x)}}
#' \item{parameter}{mult}
#' \item{p.value}{p-value, set to \code{NA} if \code{null.value = NULL}}
#' \item{conf.int}{confidence interval on the true directly
#' standardized rate}
#' \item{estimate}{directly standardized rate}
#' \item{null.value}{null hypothesis value for the DSR}
#' \item{alternative}{alternative hypothesis type}
#' \item{method}{description of the method}
#' \item{data.name}{description of the data}
#' 
#' @importFrom stats qnorm qlnorm qbeta pbeta pnorm uniroot qlogis plogis
#' @importFrom exactci poisson.exact
#' @importFrom asht wspoissonTest
#' @importFrom loglognorm qloglognorm ploglognorm
#' @importFrom utils packageVersion
#' @export
dsrTest <- function (x, n, w, null.value = NULL,
  alternative = c("two.sided", "less", "greater"),
  conf.level = 0.95, mult = 1,
  method = c("gamma", "asymptotic", "dobson", "beta", "bootstrap"),
  control = list()){
  # Functions for use later
  ciBound <- function(alternative, ci){
    # return the appropriate CI from c(alpha, 1-alpha)
    switch(alternative,
        less = c(0, ci[2]), greater = c(ci[1], Inf), two.sided = ci)
  }
  pz <- function(alternative, null.value, mean, sd){
  # return p.value from zscore (or NA where null.value is NULL)
      p.value <- NA
    if (!is.null(null.value)){
      zstat <- (null.value - mean) / sd
      p.value <- switch(alternative,
        less = stats::pnorm(zstat, lower.tail = FALSE),
        greater = stats::pnorm(zstat, lower.tail = TRUE),
        two.sided = 2 * (stats::pnorm(-abs(zstat))))
    }
  p.value
  }
  scaleNull <- function(null.value, mult, nullv = NULL, fun = identity){
  # scale null.value if not null
      r <- if (is.null(null.value)) nullv else fun(null.value / mult)
      r
  }
  methodName <- function(using, ...){
  # create a method name  (save a bit of typing)
    paste("Directly standardized rate: ", using, ..., sep = " ")
  }
  # - Argument checking
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  if (conf.level < 0.5)
    stop("conf.level must be greater than 0.5")
  alpha <- 1 - conf.level
  if (alternative == "two.sided")
    alpha <- alpha / 2
  if (length(x) != length(w) || length(x) != length(n))
    stop("'x', 'n' and 'w' must be of the same length")
  if (anyNA(x))
    stop("missing values in 'x'")
  if (anyNA(n))
    stop("missing values in 'n'")
  if (anyNA(w))
    stop("missing values in 'w'")
  # controls
  control <- switch(method,
    dobson = do.call(dobsonControl, control),
    asymptotic = do.call(asymptoticControl, control),
    gamma = do.call(gammaControl, control),
    beta = do.call(betaControl, control),
    control)
  # estimated values
  # rescale weights to include timebases
  W <- (w / sum(w)) / n
  # dsr
  y <- sum(W * x)
  # variance of the dsr
  v <- sum(x * W ^ 2)
  # sum of x
  X <- sum(x)
  # CI and p values by method
  # -- Dobson - uses exactci::poisson.exact
  if (method == "dobson"){
    # Inverse (scale and shift) the null.value
    # exactci::poisson.exact requires non-NULL r
    nv <-
      if (is.null(null.value))
        1
    else
      X + ( (null.value / mult) - y) * sqrt(X / v)
    htest <- do.call(exactci::poisson.exact,
        c(list(x = X, r = nv, alternative = alternative,
          conf.level = conf.level), control))
    CINT <- ciBound(alternative,
                    y + sqrt(v / X ) * (htest[["conf.int"]] - X))
    p.value <- if (is.null(null.value)) NA else htest[["p.value"]]
    dname <- htest[["data.name"]]
    method <- methodName(
      "Dobson's method for Weighted Sum of Poissons with",
      htest[["method"]])
  }
  # Asymptotic
   if (method == "asymptotic"){
     control <- do.call(asymptoticControl, control)
     if (control[["trans"]] == "none"){
       CINT <- ciBound(alternative,
                stats::qnorm(c(alpha, 1 - alpha), y, sqrt(v)))
       p.value <- pz(alternative,
                     scaleNull(null.value, mult), y, sqrt(v))
       method <- methodName("Asymptotic method for Weighted Sum of Poissons",
                            "normal approximation of MLE")
     }
     if (control[["trans"]] == "log"){
       # variance of ln[y]
       vstar <- v / (y ^ 2)
       CINT <- ciBound(alternative,
        stats::qlnorm(c(alpha, 1 - alpha), log(y), sqrt(vstar)))
       zlstat <- scaleNull(null.value, mult, fun = log)
       p.value <- pz(alternative, zlstat, log(y), sqrt(vstar))
       method <- methodName("Asymptotic method for Weighted Sum of Poissons",
        "normal approximation of the log-transformed MLE")
     }
     if (control[["trans"]] == "loglog"){
       # variance of ln[ ln[-y]]
       vstar <- v / ( (y * log(y)) ^ 2)
       llog <- function(x) log(-log(x))
       CINT <- ciBound(alternative,
         loglognorm::qloglognorm(
           c(alpha, 1 - alpha), llog(y), sqrt(vstar)))
       nv <- scaleNull(null.value, mult)
       p.value <-
         if (is.null(null.value))
           NA
       else {
         pAL <- loglognorm::ploglognorm(nv, llog(y), sqrt(vstar))
         switch(alternative,
                less = 1 - pAL, greater = pAL,
                two.sided = min(1, 2 * pAL, 2 * pAG))
       }
       method <- methodName(
         "Asymptotic method for Weighted Sum of Poissons",
         "normal approximation of the log-log-transformed MLE")
     }
     if (control[["trans"]] == "logit"){
       # variance of logit[y]]
       vstar <- v / ( (y * (1 - y)) ^ 2)
       ql <- stats::qlogis
       CINT <- ciBound(alternative,
                 stats::plogis(stats::qnorm(
                 c(alpha, 1 - alpha), ql(y), sqrt(vstar))))
       zlstat <- scaleNull(null.value, mult, fun = ql)
       p.value <- pz(alternative, zlstat, ql(y), sqrt(vstar))
       method <- methodName(
         "Asymptotic method for Weighted Sum of Poissons",
         "normal approximation of the logit-transformed MLE")
     }
  }
  # gamma
  # uses asht::wsPoissonTest
  if (method == "gamma"){
    nv <- scaleNull(null.value, mult)
    htest <- do.call(asht::wspoissonTest, c(
      list(x = x, w = W, nullValue = nv, alternative = alternative,
           conf.level = conf.level), control))
    CINT <- htest[["conf.int"]]
    pv <- htest[["p.value"]]
    # approximately deal with asht issue that should be fixed soon
    # (only for midp = TRUE)
    p.value <-
      if (control[["midp"]] &&
          alternative  != "two.sided" &&
          utils::packageVersion("asht") <= "0.9")
        1 - pv
    else
      pv
    method <- methodName(htest[["method"]])
  }
  # beta
  if (method == "beta"){
    # different wmtype
    if (control[["wmtype"]] == "none"){
      yb <- y
      vb <- v
      modifierText <- ""
    }
    if (control[["wmtype"]] == "tcz"){
      yb <- y + mean(W)
      vb <- v + mean(W ^ 2)
      modifierText <- "[with Tiwari-Clegg-Zou modification]"
    }
    if (control[["wmtype"]] == "mean"){
      .wm <- mean(W)
      yb <- y + .wm
      vb <- v + .wm ^ 2
      modifierText <- "[with wm=mean(w) modification]"
    }
    if (control[["wmtype"]] == "minmaxavg"){
      .wm <- mean(c(min(W), max(W)))
      yb <- y + .wm
      vb <- v +  .wm ^ 2
      modifierText <- "[with wm=mean(min(w),max(w)) modification]"
    }
    if (control[["wmtype"]] == "max"){
      .wm <- max(W)
      yb <- y + .wm
      vb <- v + .wm ^ 2
      modifierText <- "[with wm=max(w) modification]"
    }
    a <- yb * (yb * (1 - yb) / vb - 1)
    b <- (1 - yb) * (yb * (1 - yb) / vb - 1)
    nv <- scaleNull(null.value, mult)
    CINT <- ciBound(alternative,
                    stats::qbeta(c(alpha, 1 - alpha), a, b))
    p.value <-
      if (is.null(null.value))
        NA
     else {
       pAL <- stats::pbeta(nv, a, b, lower.tail = FALSE)
       pAG <- stats::pbeta(nv, a, b, lower.tail = TRUE)
       switch(alternative,
         less = pAL, greater = pAG,
         two.sided = min(1, 2 * pAL, 2 * pAG))
     }
    method <- methodName("Beta method for Weighted Sum of Poissons",
                         modifierText)
  }
  # - Approximate bootstrap
  if (method == "bootstrap"){
    z0 <- a <- sum(x * (W ^ 3)) / (6 * v ^ (3 / 2))
    pFunction <- function(p){
      num <- z0 + stats::qnorm(p)
      y  + sqrt(v) * num / ( (1 - a * num) ^ 2)
    }
    CINT <- ciBound(alternative, pFunction(c(alpha, 1 - alpha)))
    if (is.null(null.value)){
      p.value <- NA
    } else {
      nv <- scaleNull(null.value, mult)
      f.AL <- function(p) pFunction(1 - p) - nv
      f.AG <- function(p)  (pFunction(p) - nv)
      pAL <- stats::uniroot(f.AL, interval = c(1e-09, 1 - 1e-09))$root
      pAG <- stats::uniroot(f.AG, interval = c(1e-09, 1 - 1e-09))$root
      p.value <- switch(alternative, less = pAL, greater = pAG,
        two.sided = min(1, 2 * pAG, 2 * pAL))
    }
    method <- methodName(
      "Approximate bootstrap method for Weighted Sum of Poissons")
  }
  # Set up htest attributes
  attr(CINT, "conf.level") <- conf.level
  names(y) <- "Directly standardarized rate"
  if (mult != 1)
    names(y) <- paste0("Directly standardized rate per ", mult)
  k <- length(x)
  names(k) <- "number of summands"
  parms <- mult
  names(parms)  <- "Multiplier"
  if (!is.null(null.value))
    names(null.value) <- "Directly standardized rate"
  xname <- paste0("x = ", deparse(substitute(x)))
  timebasename <- paste0("time bases: n = ", deparse(substitute(n)))
  wname <- paste0("weights: w = ", deparse(substitute(w)))
  dname <- paste(xname, timebasename, wname, sep = ", ")
  # returned object
  structure(
    list(
      statistic = k,
      parameter = parms,
      p.value = p.value,
      conf.int = CINT * mult,
      estimate = y * mult,
      null.value = null.value,
      alternative = alternative,
      method = method,
      data.name = dname
    ),
    class = "htest")
}
