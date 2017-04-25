context("Sense checking confidence intervals")

xfive<-c(0, 8, 63, 112, 262, 295)
nfive<-c(327, 30666, 123419, 149919, 104088, 34392)
ntotal<-c(319933, 931318, 786511, 488235, 237863, 61313)

methods_list  <- list(
  gamma = list(
    list(wmtype = 'max'),
    list(midp = TRUE),
    list(wmtype = "tcz"),
    list(wmtype = "mean"),
    list(wmtype = "minmaxavg")),
  asymptotic = list(
    list(trans = "none"),
    list(trans = "log")),
  dobson = list(
    list(midp = FALSE), 
    list(midp = TRUE)),
  beta = list(list()),
  bootstrap = list(list())
)

all_ts <- mapply(dsrTest::dsrTest,
                      method = rep(names(methods_list), times = lengths(methods_list)),
                      control = do.call(c, unname(methods_list)),
                      MoreArgs = list(mult = 1e5, x = xfive, n = nfive, w = ntotal),
                      SIMPLIFY = FALSE)

all_less <- mapply(dsrTest::dsrTest,
                 method = rep(names(methods_list), times = lengths(methods_list)),
                 control = do.call(c, unname(methods_list)),
                 MoreArgs = list(mult = 1e5, x = xfive, n = nfive, w = ntotal,
                                 alternative = "less"),
                 SIMPLIFY = FALSE)

all_greater <- mapply(dsrTest::dsrTest,
                   method = rep(names(methods_list), times = lengths(methods_list)),
                   control = do.call(c, unname(methods_list)),
                   MoreArgs = list(mult = 1e5, x = xfive, n = nfive, w = ntotal,
                                   alternative = "greater"),
                   SIMPLIFY = FALSE)

test_that("Estimate is within confidence intervals for all methods",{
  expect_true(all(
    sapply(all_less, function(x) findInterval(x$estimate, x$conf.int)))==1L)
  expect_true(all(
    sapply(all_ts, function(x) findInterval(x$estimate, x$conf.int)))==1L)
  expect_true(all(
    sapply(all_greater, function(x) findInterval(x$estimate, x$conf.int)))==1L)
})


test_that("Confidence Intervals are ordered appropriately",
          {
            all(!sapply(all_less, function(x) is.unsorted(x$confint)))
            all(!sapply(all_greater, function(x) is.unsorted(x$confint)))
            all(!sapply(all_ts, function(x) is.unsorted(x$confint)))
          }
)