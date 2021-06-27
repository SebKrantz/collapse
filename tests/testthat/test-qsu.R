context("qsu")

# rm(list = ls())

# Make more tests!! See also collapse general TODO !
test_that("qsu gives errors for wrong input", {
  expect_error(qsu(wlddev$year, 2:4))
  expect_error(qsu(wlddev$year, pid = 2:4))
  expect_error(qsu(wlddev, 2:4))
  expect_error(qsu(wlddev, pid = 2:4))
  expect_error(qsu(wlddev$year, letters))
  expect_error(qsu(wlddev$year, pid = letters))
  expect_error(qsu(wlddev, letters))
  expect_error(qsu(wlddev, pid = letters))

  expect_error(qsu(wlddev, ~ iso3c + bla))
  expect_error(qsu(wlddev, pid = ~ iso3c + bla))

  expect_visible(qsu(wlddev, PCGDP ~ region + income))
  expect_visible(qsu(wlddev, pid = PCGDP ~ region + income))

  expect_equal(qsu(wlddev, PCGDP ~ region + income, ~ iso3c), qsu(wlddev, ~ region + income, pid = PCGDP ~ iso3c))

  expect_error(qsu(wlddev, cols = 9:14))
  expect_error(qsu(wlddev, cols = c("PCGDP","bla")))
})
