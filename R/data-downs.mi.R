#' @title Downs' syndrome cases and of total live births 
#' by maternal age and birth order, Michigan, 1950-1964.
#' 
#' @description
#' This data reproduces table 14.4 in Fleiss (1981) and which is
#' drawn from a large-scale study in Michigan 1950-1964 of the effect
#' of Maternal age and Birth order on Down Syndrome and Leukemia 
#' (Stark and Mantell, 1966).
#' @format
#' This `data.frame` contains the following columns:
#' \describe{
#' \item{\code{Age}}{The age category of the mothers.}
#' \item{\code{BirthOrder}}{The birth order.}
#' \item{\code{Cases}}{The number of cases of Down's Sydrome.}
#' \item{\code{Births}}{The number of live births.}
#' \item{\code{Standard}}{A "standard" population, derived as the total
#' number of births in each age category}
#' }
#' @source 
#' The data were obtained from table 14.4 (p 249) in
#' 
#' Fleiss, JL (1981) *Statistical Methods for Rates and Proportions*, 
#' Wiley, New York.
#' 
#' The original study is
#' 
#' Stark CR and Mantel N (1966) 'Effects of maternal age and 
#' birth order on the risk of mongolism and leukemia' 
#' *J Natl Cancer Inst* **37** (5) 687--698. \doi{doi:10.1093/jnci/37.5.687}
"downs.mi"