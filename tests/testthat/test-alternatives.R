context("Sense checking confidence intervals")

xfive<-c(0, 8, 63, 112, 262, 295)
nfive<-c(327, 30666, 123419, 149919, 104088, 34392)
ntotal<-c(319933, 931318, 786511, 488235, 237863, 61313)



test_that("Confidence Intervals are ordered appropriately",
  {
    expect_false(is.unsorted(dsrTest(xfive, nfive, ntotal, 
      method = 'beta')[['conf.int']]))
    expect_false(is.unsorted(dsrTest(xfive, nfive, ntotal, 
      method = 'bootstrap')[['conf.int']]))
    expect_false(is.unsorted(dsrTest(xfive, nfive, ntotal, 
      method = 'dobson')[['conf.int']]))
    expect_false(is.unsorted(dsrTest(xfive, nfive, ntotal, 
      method = 'dobson', control = list(midp = TRUE))[['conf.int']]))
    expect_false(is.unsorted(dsrTest(xfive, nfive, ntotal, 
      method = 'asymptotic')[['conf.int']]))
    expect_false(is.unsorted(dsrTest(xfive, nfive, ntotal, 
      method = 'asymptotic', control = list(trans = 'log'))[['conf.int']]))
  }
)
