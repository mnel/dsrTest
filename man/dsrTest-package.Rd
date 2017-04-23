\name{dsrTest-package}
\docType{package}
\title{\packageTitle{dsrTest}}
\description{\packageDescription{dsrTest}}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{dsrTest}
  \packageIndices{dsrTest}

  The function \code{dsrTest()} implements a number of different 
  parameterizations and methods for computing confidence intervals on 
  directly standardized rates. These methods are described in detail in
  Fay and Feuer (1997), Tiwari et al (2006), Ng et al (2008) and Fay
  and Kim (2017).
  
  The various \code{<method>Control()} functions allow modifications to
  the general approaches. 
}
\author{
  \packageAuthor{dsrTest}
  
  Maintainer: \packageMaintainer{dsrTest}
}
\references{
Dobson, AJ, Kuulasmaa, K, Eberle, E and Scherer, J (1991)
'Confidence intervals for weighted sums of Poisson parameters',
\emph{Statistics in Medicine}, \strong{10}: 457--462. \doi{doi:10.1002/sim.4780100317}


Swift, MB (1995)  'Simple confidence intervals for standardized rates based on the approximate bootstrap method', \emph{Statistics in Medicine}, \strong{14}, 1875--1888. \doi{doi:10.1002/sim.4780141704}.

Fay MP & Feuer EJ (1997) 'Confidence intervals for directly standardized rates: a method based on the gamma distribution.' \emph{Statistics in Medicine}. \strong{16}: 791--801. \url{https://doi.org/10.1002/(SICI)1097-0258(19970415)16:7<791::AID-SIM500>3.0.CO;2-\%23}

Tiwari RC, Clegg LX, & Zou Z (2006)  'Efficient interval estimation for age-adjusted cancer rates.' \emph{Statistical Methods in Medical Research} \strong{15}: 547--569. \doi{doi:10.1177/0962280206070621}

Ng HKT, Filardo, G & Zheng G (2008) 'Confidence interval estimating procedures for standardized incidence rates.' \emph{Computational Statistics and Data Analysis} \strong{52}: 3501--3516. \doi{doi:10.1016/j.csda.2007.11.004}

Fay, MP & Kim S (2017) 'Confidence intervals for directly standardized rates using mid-p gamma intervals.' \emph{Biometrical Journal}  \strong{59}(2): 377--387. \doi{doi:10.1002/bimj.201600111}

}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{package}
\seealso{
  \code{\link[asht]{wspoissonTest}}
}
\examples{
## using the example from asht::wspoissonTest
## birth data on Down's syndrome from Michigan, 1950-1964
## see Table II  of Fay and Feuer (1997)
## xfive = counts for mothers who have had 5 or more children
## nfive and ntotal are number of live births 
xfive<-c(0, 8, 63, 112, 262, 295)
nfive<-c(327, 30666, 123419, 149919, 104088, 34392)
ntotal<-c(319933, 931318, 786511, 488235, 237863, 61313)
## use mult =10^5 to give rates per 100,000
## gamma method of Fay and Feuer (1997) is default
dsrTest(xfive, nfive, ntotal, method = "gamma", mult = 1e5)
## Dobson et al (1991)
dsrTest(xfive, nfive, ntotal, method = "dobson", mult = 1e5)
## Asymptotic with log transformation
dsrTest(xfive, nfive, ntotal, method = "asymptotic", mult = 1e5,
control = list(trans = "log"))
}
