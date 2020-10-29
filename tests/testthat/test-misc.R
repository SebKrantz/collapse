context("Misc")

# rm(list = ls())

m <- na_insert(qM(mtcars))

test_that("descr, pwcor, cov, Nobs, misc", {
  expect_visible(descr(wlddev))
  expect_visible(descr(GGDC10S))
  expect_visible(pwcor(nv(wlddev)))
  expect_visible(pwcor(nv(GGDC10S)))
  expect_visible(pwcov(nv(wlddev)))
  expect_visible(pwcov(nv(GGDC10S)))
  expect_visible(pwNobs(wlddev))
  expect_visible(pwNobs(GGDC10S))

  expect_visible(descr(m))
  expect_visible(pwcor(m))
  expect_visible(pwcov(m))
  expect_visible(pwNobs(m))

})
