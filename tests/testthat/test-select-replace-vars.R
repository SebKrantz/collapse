context("select_replace vars")

test_that("selecting vars works well", {
  expect_identical(get_vars(wlddev, 4:8), wlddev[4:8])
  expect_identical(get_vars(wlddev, -(4:8)), wlddev[-(4:8)])
  expect_identical(get_vars(wlddev, sapply(wlddev, is.numeric)), wlddev[sapply(wlddev, is.numeric)])
  expect_identical(get_vars(wlddev, c("iso3c","PCGDP","ODA")), wlddev[c("iso3c","PCGDP","ODA")])
  expect_identical(get_vars(wlddev, is.factor), wlddev[sapply(wlddev, is.factor)])

  expect_identical(num_vars(wlddev), wlddev[sapply(wlddev, is.numeric)])
  expect_identical(cat_vars(wlddev), wlddev[sapply(wlddev, is.categorical)])
  expect_identical(char_vars(wlddev), wlddev[sapply(wlddev, is.character)])
  expect_identical(fact_vars(wlddev), wlddev[sapply(wlddev, is.factor)])
  expect_identical(Date_vars(wlddev), wlddev[sapply(wlddev, is.Date)])
})

test_that("replacing vars works well", {
  wlddevold <- wlddev
  get_vars(wlddev, 4:8) <- get_vars(wlddev, 4:8)
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  get_vars(wlddev, -(4:8)) <- get_vars(wlddev, -(4:8))
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  get_vars(wlddev, -(4:8)) <- as.list(get_vars(wlddev, -(4:8)))
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  get_vars(wlddev, c("iso3c","PCGDP","ODA")) <- get_vars(wlddev, c("iso3c","PCGDP","ODA"))
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  get_vars(wlddev, sapply(wlddev, is.numeric)) <- get_vars(wlddev, sapply(wlddev, is.numeric))
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  get_vars(wlddev, is.factor) <- get_vars(wlddev, is.factor)
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  num_vars(wlddev) <- num_vars(wlddev)
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  cat_vars(wlddev) <- cat_vars(wlddev)
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  char_vars(wlddev) <- char_vars(wlddev)
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  fact_vars(wlddev) <- fact_vars(wlddev)
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  logi_vars(wlddev) <- logi_vars(wlddev)
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  Date_vars(wlddev) <- Date_vars(wlddev)
  expect_identical(wlddevold, wlddev)
})

test_that("select vars errors for wrong input", {
  expect_error(get_vars(wlddev, 13))
  expect_error(get_vars(wlddev, 1:13))
  expect_error(get_vars(wlddev, -13))
  expect_error(get_vars(wlddev, c("PCGDP","ODA3")))
  expect_error(get_vars(wlddev, c(sapply(wlddev, is.numeric), TRUE)))
  expect_error(get_vars(wlddev, sapply(wlddev, is.numeric)[-1]))
})

test_that("replace vars errors for wrong input", {
  expect_error(get_vars(wlddev, 13) <- wlddev[12])
  expect_error(get_vars(wlddev, "ODA3") <- wlddev[12])
  expect_error(get_vars(wlddev, -13) <- wlddev[12])
  expect_error(get_vars(wlddev, 11:12) <- wlddev[12])
  expect_error(get_vars(wlddev, 9:12) <- wlddev[8:12])
  expect_invisible(get_vars(wlddev, 12) <- wlddev$ODA)
  expect_error(get_vars(wlddev, c(sapply(wlddev, is.numeric), TRUE)) <- wlddev)
  expect_error(get_vars(wlddev, sapply(wlddev, is.numeric)[-1]) <- wlddev)
})
