context("psmat and psacf")



# rm(list = ls())

options(warn = -1)

test_that("psmat works as intended", {
  expect_identical(psmat(wlddev$PCGDP, wlddev$iso3c, wlddev$year), psmat(wlddev, PCGDP ~ iso3c, ~ year))
  expect_identical(psmat(wlddev$PCGDP, wlddev$iso3c, wlddev$year), psmat(wlddev[9], wlddev$iso3c, wlddev$year))
  expect_identical(psmat(wlddev[9:12], wlddev$iso3c, wlddev$year), psmat(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c, ~ year))
  expect_identical(psmat(wlddev[9:12], wlddev$iso3c, wlddev$year), psmat(wlddev, ~ iso3c, ~ year, cols = 9:12))
  # without year
  expect_identical(psmat(wlddev$PCGDP, wlddev$iso3c), psmat(wlddev, PCGDP ~ iso3c))
  expect_identical(psmat(wlddev$PCGDP, wlddev$iso3c), psmat(wlddev[9], wlddev$iso3c))
  expect_identical(psmat(wlddev[9:12], wlddev$iso3c), psmat(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c))
  expect_identical(psmat(wlddev[9:12], wlddev$iso3c), psmat(wlddev, ~ iso3c, cols = 9:12))
  # only nid's
  expect_identical(psmat(wlddev$PCGDP, 216), psmat(wlddev[9], 216))
  expect_identical(psmat(wlddev[9:12], 216), psmat(wlddev, 216, cols = 9:12))

  # TRANSPOSE
  expect_identical(psmat(wlddev$PCGDP, wlddev$iso3c, wlddev$year, transpose = TRUE),  `attr<-`(t(psmat(wlddev$PCGDP, wlddev$iso3c, wlddev$year)), "transpose", TRUE))
  expect_identical(psmat(wlddev$PCGDP, wlddev$iso3c, wlddev$year, transpose = TRUE), psmat(wlddev, PCGDP ~ iso3c, ~ year, transpose = TRUE))
  expect_identical(psmat(wlddev$PCGDP, wlddev$iso3c, wlddev$year, transpose = TRUE), psmat(wlddev[9], wlddev$iso3c, wlddev$year, transpose = TRUE))
  expect_identical(psmat(wlddev[9:12], wlddev$iso3c, wlddev$year, transpose = TRUE), psmat(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c, ~ year, transpose = TRUE))
  expect_identical(psmat(wlddev[9:12], wlddev$iso3c, wlddev$year, transpose = TRUE), psmat(wlddev, ~ iso3c, ~ year, cols = 9:12, transpose = TRUE))
  # without year
  expect_identical(psmat(wlddev$PCGDP, wlddev$iso3c, transpose = TRUE), psmat(wlddev, PCGDP ~ iso3c, transpose = TRUE))
  expect_identical(psmat(wlddev$PCGDP, wlddev$iso3c, transpose = TRUE), psmat(wlddev[9], wlddev$iso3c, transpose = TRUE))
  expect_identical(psmat(wlddev[9:12], wlddev$iso3c, transpose = TRUE), psmat(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c, transpose = TRUE))
  expect_identical(psmat(wlddev[9:12], wlddev$iso3c, transpose = TRUE), psmat(wlddev, ~ iso3c, cols = 9:12, transpose = TRUE))
  # only nid's
  expect_identical(psmat(wlddev$PCGDP, 216, transpose = TRUE), psmat(wlddev[9], 216, transpose = TRUE))
  expect_identical(psmat(wlddev[9:12], 216, transpose = TRUE), psmat(wlddev, 216, cols = 9:12, transpose = TRUE))

  # LIST-OUTPUT
  expect_true(is.array(psmat(wlddev[9:12], wlddev$iso3c, wlddev$year)))
  expect_true(is.list(psmat(wlddev[9:12], wlddev$iso3c, wlddev$year, array = FALSE)))
  expect_identical(psmat(wlddev[9:12], wlddev$iso3c, wlddev$year, array = FALSE), psmat(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c, ~ year, array = FALSE))
  expect_identical(psmat(wlddev[9:12], wlddev$iso3c, wlddev$year, array = FALSE), psmat(wlddev, ~ iso3c, ~ year, cols = 9:12, array = FALSE))
  # without year
  expect_identical(psmat(wlddev[9:12], wlddev$iso3c, array = FALSE), psmat(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c, array = FALSE))
  expect_identical(psmat(wlddev[9:12], wlddev$iso3c, array = FALSE), psmat(wlddev, ~ iso3c, cols = 9:12, array = FALSE))
  # only nid's
  expect_identical(psmat(wlddev[9:12], 216, array = FALSE), psmat(wlddev, 216, cols = 9:12, array = FALSE))
})

test_that("psacf works as intended", {
  x <- na_rm(wlddev$PCGDP)
  expect_equal(unclass(psacf(x, rep(1,length(x)), seq_along(x), lag.max = 12, plot = FALSE))[1:4], unclass(acf(x, lag.max = 12, plot = FALSE))[1:4], tolerance = 1e-3)
  expect_equal(unclass(psacf(x, rep(1,length(x)), seq_along(x), type = "covariance", lag.max = 12, gscale = FALSE, plot = FALSE))[1:4], unclass(acf(x, type = "covariance", lag.max = 12, plot = FALSE))[1:4], tolerance = 1e-3)
  expect_equal(unclass(pspacf(x, rep(1,length(x)), seq_along(x), lag.max = 12, plot = FALSE))[1:4], unclass(pacf(x, lag.max = 12, plot = FALSE))[1:4], tolerance = 1e-3)
  dat <- na_omit(get_vars(wlddev, c(9:10,12)))
  expect_equal(unclass(psacf(dat, rep(1,nrow(dat)), seq_row(dat), lag.max = 12, plot = FALSE))[1:4], unclass(acf(dat, lag.max = 12, plot = FALSE))[1:4], tolerance = 1e-2)
  expect_equal(unclass(psacf(dat, rep(1,nrow(dat)), seq_row(dat), type = "covariance", lag.max = 12, gscale = FALSE, plot = FALSE))[1:4], unclass(acf(dat, type = "covariance", lag.max = 12, plot = FALSE))[1:4], tolerance = 1e-2)
  # expect_equal(unclass(pspacf(dat, rep(1,nrow(dat)), seq_row(dat), lag.max = 12, plot = FALSE))[1:4], unclass(pacf(dat, lag.max = 12, plot = FALSE))[1:4], tolerance = 1e-2) # This is strange !!!!

  expect_equal(unclass(psacf(wlddev$PCGDP, wlddev$iso3c, wlddev$year, plot = FALSE))[1:4], unclass(psacf(wlddev, PCGDP ~ iso3c, ~ year, plot = FALSE))[1:4])
  expect_equal(unclass(psacf(wlddev$PCGDP, wlddev$iso3c, wlddev$year, plot = FALSE))[1:4], unclass(psacf(wlddev[9], wlddev$iso3c, wlddev$year, plot = FALSE))[1:4])
  expect_equal(unclass(psacf(wlddev[9:12], wlddev$iso3c, wlddev$year, plot = FALSE))[1:4], unclass(psacf(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c, ~ year, plot = FALSE))[1:4])
  expect_equal(unclass(psacf(wlddev[9:12], wlddev$iso3c, wlddev$year, plot = FALSE))[1:4], unclass(psacf(wlddev, ~ iso3c, ~ year, cols = 9:12, plot = FALSE))[1:4])
  # equality
  expect_equal(unclass(psacf(wlddev$PCGDP, wlddev$iso3c, plot = FALSE))[1:4], unclass(psacf(wlddev$PCGDP, wlddev$iso3c, wlddev$year, plot = FALSE))[1:4])
  # without year
  expect_equal(unclass(psacf(wlddev$PCGDP, wlddev$iso3c, plot = FALSE))[1:4], unclass(psacf(wlddev, PCGDP ~ iso3c, plot = FALSE))[1:4])
  expect_equal(unclass(psacf(wlddev$PCGDP, wlddev$iso3c, plot = FALSE))[1:4], unclass(psacf(wlddev[9], wlddev$iso3c, plot = FALSE))[1:4])
  expect_equal(unclass(psacf(wlddev[9:12], wlddev$iso3c, plot = FALSE))[1:4], unclass(psacf(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c, plot = FALSE))[1:4])
  expect_equal(unclass(psacf(wlddev[9:12], wlddev$iso3c, plot = FALSE))[1:4], unclass(psacf(wlddev, ~ iso3c, cols = 9:12, plot = FALSE))[1:4])
})

test_that("pspacf works as intended", {
  expect_equal(unclass(pspacf(wlddev$PCGDP, wlddev$iso3c, wlddev$year, plot = FALSE))[1:4], unclass(pspacf(wlddev, PCGDP ~ iso3c, ~ year, plot = FALSE))[1:4])
  expect_equal(unclass(pspacf(wlddev$PCGDP, wlddev$iso3c, wlddev$year, plot = FALSE))[1:4], unclass(pspacf(wlddev[9], wlddev$iso3c, wlddev$year, plot = FALSE))[1:4])
  expect_equal(unclass(pspacf(wlddev[9:12], wlddev$iso3c, wlddev$year, plot = FALSE))[1:4], unclass(pspacf(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c, ~ year, plot = FALSE))[1:4])
  expect_equal(unclass(pspacf(wlddev[9:12], wlddev$iso3c, wlddev$year, plot = FALSE))[1:4], unclass(pspacf(wlddev, ~ iso3c, ~ year, cols = 9:12, plot = FALSE))[1:4])
  # equality
  expect_equal(unclass(pspacf(wlddev$PCGDP, wlddev$iso3c, plot = FALSE))[1:4], unclass(pspacf(wlddev$PCGDP, wlddev$iso3c, wlddev$year, plot = FALSE))[1:4])
  # without year
  expect_equal(unclass(pspacf(wlddev$PCGDP, wlddev$iso3c, plot = FALSE))[1:4], unclass(pspacf(wlddev, PCGDP ~ iso3c, plot = FALSE))[1:4])
  expect_equal(unclass(pspacf(wlddev$PCGDP, wlddev$iso3c, plot = FALSE))[1:4], unclass(pspacf(wlddev[9], wlddev$iso3c, plot = FALSE))[1:4])
  expect_equal(unclass(pspacf(wlddev[9:12], wlddev$iso3c, plot = FALSE))[1:4], unclass(pspacf(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c, plot = FALSE))[1:4])
  expect_equal(unclass(pspacf(wlddev[9:12], wlddev$iso3c, plot = FALSE))[1:4], unclass(pspacf(wlddev, ~ iso3c, cols = 9:12, plot = FALSE))[1:4])
})

test_that("psmat gives errors for wrong input", {
  # wrong lengths
  expect_error(psmat(wlddev$PCGDP, wlddev$iso3c[-1], wlddev$year))
  expect_error(psmat(wlddev$PCGDP, wlddev$iso3c, wlddev$year[-1]))
  expect_error(psmat(wlddev[9:12], wlddev$iso3c[-1], wlddev$year))
  expect_error(psmat(wlddev[9:12], wlddev$iso3c, wlddev$year[-1]))
  # without year
  expect_error(psmat(wlddev$PCGDP, wlddev$iso3c[-1]))
  expect_error(psmat(wlddev[9:12], wlddev$iso3c[-1]))
  # only nid's
  expect_error(psmat(wlddev$PCGDP, 218))
  expect_error(psmat(wlddev[9:12], 218))

  # wrong formula
  expect_error(psmat(wlddev, PCGDP2 ~ iso3c, ~ year))
  expect_error(psmat(wlddev, PCGDP ~ iso3c2, ~ year))
  expect_error(psmat(wlddev, PCGDP ~ iso3c, ~ year2))
  expect_error(psmat(wlddev, PCGDP + LIFEEX + GINI + ODA + bla ~ iso3c, ~ year))
  expect_error(psmat(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c + bla, ~ year))
  expect_error(psmat(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c, ~ year + bla))
  # without year
  expect_error(psmat(wlddev, PCGDP2 ~ iso3c))
  expect_error(psmat(wlddev, PCGDP ~ iso3c2))
  expect_error(psmat(wlddev, PCGDP + LIFEEX + GINI + ODA + bla ~ iso3c, ~ year))
  expect_error(psmat(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c + bla, ~ year))
  # cols
  expect_error(psmat(wlddev, ~ iso3c, ~ year, cols = 14))
  expect_error(psmat(wlddev, ~ iso3c, ~ year, cols = "bla"))
  expect_visible(psmat(wlddev, ~ iso3c, ~ year, cols = sapply(wlddev, is.numeric)))
  expect_error(psmat(wlddev, ~ iso3c, ~ year, cols = sapply(wlddev, is.numeric)[-1]))
})

test_that("psacf gives errors for wrong input", {
  # wrong lengths
  expect_error(psacf(wlddev$PCGDP, wlddev$iso3c[-1], wlddev$year, plot = FALSE))
  expect_error(psacf(wlddev$PCGDP, wlddev$iso3c, wlddev$year[-1], plot = FALSE))
  expect_error(psacf(wlddev[9:12], wlddev$iso3c[-1], wlddev$year, plot = FALSE))
  expect_error(psacf(wlddev[9:12], wlddev$iso3c, wlddev$year[-1], plot = FALSE))
  # without year
  expect_error(psacf(wlddev$PCGDP, wlddev$iso3c[-1], plot = FALSE))
  expect_error(psacf(wlddev[9:12], wlddev$iso3c[-1], plot = FALSE))
  # this should give error...
  expect_error(psacf(wlddev$PCGDP, 218, plot = FALSE))
  expect_error(psacf(wlddev[9:12], 218, plot = FALSE))

  # wrong formula
  expect_error(psacf(wlddev, PCGDP2 ~ iso3c, ~ year, plot = FALSE))
  expect_error(psacf(wlddev, PCGDP ~ iso3c2, ~ year, plot = FALSE))
  expect_error(psacf(wlddev, PCGDP ~ iso3c, ~ year2, plot = FALSE))
  expect_error(psacf(wlddev, PCGDP + LIFEEX + GINI + ODA + bla ~ iso3c, ~ year, plot = FALSE))
  expect_error(psacf(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c + bla, ~ year, plot = FALSE))
  expect_error(psacf(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c, ~ year + bla, plot = FALSE))
  # without year
  expect_error(psacf(wlddev, PCGDP2 ~ iso3c, plot = FALSE))
  expect_error(psacf(wlddev, PCGDP ~ iso3c2, plot = FALSE))
  expect_error(psacf(wlddev, PCGDP + LIFEEX + GINI + ODA + bla ~ iso3c, ~ year, plot = FALSE))
  expect_error(psacf(wlddev, PCGDP + LIFEEX + GINI + ODA ~ iso3c + bla, ~ year, plot = FALSE))
  # cols
  expect_error(psacf(wlddev, ~ iso3c, ~ year, cols = 14, plot = FALSE))
  expect_error(psacf(wlddev, ~ iso3c, ~ year, cols = "bla", plot = FALSE))
  expect_visible(psacf(wlddev, ~ iso3c, ~ year, cols = sapply(wlddev, is.numeric), plot = FALSE))
  expect_error(psacf(wlddev, ~ iso3c, ~ year, cols = sapply(wlddev, is.numeric)[-1], plot = FALSE))
})

options(warn = 1)
