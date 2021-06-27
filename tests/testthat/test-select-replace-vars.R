context("select, replace or add vars")

# rm(list = ls())

test_that("selecting vars works well", {
  expect_identical(get_vars(wlddev, 4:8), wlddev[4:8])
  expect_identical(get_vars(wlddev, -(4:8)), wlddev[-(4:8)])
  expect_identical(get_vars(wlddev, sapply(wlddev, is.numeric)), wlddev[sapply(wlddev, is.numeric)])
  expect_identical(get_vars(wlddev, c("iso3c","PCGDP","ODA")), wlddev[c("iso3c","PCGDP","ODA")])
  expect_identical(get_vars(wlddev, "D", regex = TRUE), wlddev[c("OECD","PCGDP","ODA")])
  expect_identical(get_vars(wlddev, c("D","L"), regex = TRUE), wlddev[c("OECD","PCGDP","LIFEEX","ODA")])
  expect_identical(get_vars(wlddev, is.factor), wlddev[sapply(wlddev, is.factor)])

  expect_identical(num_vars(wlddev), wlddev[sapply(wlddev, is.numeric)])
  expect_identical(cat_vars(wlddev), wlddev[sapply(wlddev, is_categorical)])
  expect_identical(char_vars(wlddev), wlddev[sapply(wlddev, is.character)])
  expect_identical(fact_vars(wlddev), wlddev[sapply(wlddev, is.factor)])
  expect_identical(date_vars(wlddev), wlddev[sapply(wlddev, is_date)])
})

test_that("replacing vars works well", {
  wlddevold <- wlddev
  get_vars(wlddev, 4:8) <- get_vars(wlddev, 4:8)
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  fselect(wlddev, PCGDP:GINI) <- fselect(wlddev, PCGDP:GINI)
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  get_vars(wlddev, -(4:8)) <- get_vars(wlddev, -(4:8))
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  fselect(wlddev, -(PCGDP:GINI)) <- fselect(wlddev, -(PCGDP:GINI))
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  get_vars(wlddev, -(4:8)) <- as.list(get_vars(wlddev, -(4:8)))
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  get_vars(wlddev, c("iso3c","PCGDP","ODA")) <- get_vars(wlddev, c("iso3c","PCGDP","ODA"))
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  get_vars(wlddev, "D", regex = TRUE) <- get_vars(wlddev, "D", regex = TRUE)
  expect_identical(wlddevold, wlddev)

  wlddevold <- wlddev
  get_vars(wlddev, c("D","L"), regex = TRUE) <- get_vars(wlddev, c("D","L"), regex = TRUE)
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
  date_vars(wlddev) <- date_vars(wlddev)
  expect_identical(wlddevold, wlddev)
})

test_that("adding vars works well", {
  wlddev1 <- wlddev2 <- wlddev
  temp <- STD(get_vars(wlddev, 9:12))
  add_vars(wlddev1) <- temp
  wlddev2[names(temp)] <- temp
  expect_identical(wlddev1, wlddev2)

  wlddev1 <- wlddev
  temp <- STD(get_vars(wlddev, 9:12))
  add_vars(wlddev1, "front") <- temp
  expect_identical(wlddev1, add_vars(temp, wlddev))

  wlddev1 <- wlddev
  temp <- STD(get_vars(wlddev, 9:13))
  add_vars(wlddev1, c(10,12,14,16,18)) <- temp
  expect_true(all_identical(wlddev1,  add_vars(wlddev, temp, pos = c(10,12,14,16,18)),
                            add_vars(gv(wlddev, 1:9), gv(temp, 1), gv(wlddev, 10), gv(temp, 2),
                                     gv(wlddev, 11), gv(temp, 3), gv(wlddev, 12), gv(temp, 4), gv(wlddev, 13), gv(temp, 5))))

})

test_that("replacing with or adding atomic elements works well", {
  wlddev1 <- wlddev2 <- wlddev
  get_vars(wlddev1, 9) <- wlddev$PCGDP
  expect_identical(wlddev1, wlddev)

  get_vars(wlddev1, 9) <- qM(wlddev[9:12])
  wlddev2[9] <- qM(wlddev[9:12])
  expect_identical(wlddev1, wlddev2)

  wlddev1 <- wlddev2 <- wlddev
  add_vars(wlddev1) <- wlddev$PCGDP
  expect_identical(wlddev1, cbind(wlddev2, wlddev["PCGDP"]))

  wlddev1 <- wlddev2 <- wlddev
  add_vars(wlddev1) <- qM(wlddev[9:12])
  wlddev2["wlddev[9:12]"] <- qM(wlddev[9:12]) # formerly wlddev2["qM(wlddev[9:12])"], but no longer using deparse..
  expect_identical(wlddev1, wlddev2)

  wlddev1 <- wlddev2 <- wlddev
  add_vars(wlddev1, "front") <- wlddev$PCGDP
  expect_identical(wlddev1, add_vars(wlddev, wlddev$PCGDP, pos = 1))

  wlddev1 <- wlddev2 <- wlddev
  add_vars(wlddev1, "front") <- qM(wlddev[9:12])
  expect_identical(wlddev1, add_vars(wlddev, qM(wlddev[9:12]), pos = 1))
})


test_that("empty selections work well", {
  expect_equal(cat_vars(mtcars), mtcars[0L])
  expect_equal(char_vars(mtcars), mtcars[0L])
  expect_equal(fact_vars(mtcars), mtcars[0L])
  expect_equal(logi_vars(mtcars), mtcars[0L])
  expect_equal(get_vars(mtcars, is.character), mtcars[0L])
  expect_equal(get_vars(mtcars, 0L), mtcars[0L])
  expect_error(get_vars(mtcars, NULL))
})


test_that("select vars errors for wrong input", {
  expect_error(get_vars(wlddev, 14))
  expect_error(get_vars(wlddev, 1:14))
  expect_error(get_vars(wlddev, -14))
  expect_error(get_vars(wlddev, c("PCGDP","ODA3")))
  # expect_warning(get_vars(wlddev, "bla", regex = TRUE)) # Better give error
  expect_error(get_vars(wlddev, c(sapply(wlddev, is.numeric), TRUE)))
  expect_error(get_vars(wlddev, sapply(wlddev, is.numeric)[-1]))
})

test_that("replace vars errors for wrong input", {
  expect_error(get_vars(wlddev, 14) <- wlddev[12])
  expect_error(get_vars(wlddev, "ODA3") <- wlddev[12])
  expect_error(get_vars(wlddev, "bla", regex = TRUE) <- wlddev[12])
  expect_error(get_vars(wlddev, -14) <- wlddev[12])
  expect_error(get_vars(wlddev, 11:12) <- wlddev[12])
  expect_error(get_vars(wlddev, 9:12) <- wlddev[8:12])
  expect_invisible(get_vars(wlddev, 12) <- wlddev$ODA)
  expect_error(get_vars(wlddev, 12) <- wlddev$ODA[-1])
  expect_error(get_vars(wlddev, 12) <- qM(wlddev[9:12])[-1, ])
  expect_error(get_vars(wlddev, c(sapply(wlddev, is.numeric), TRUE)) <- wlddev)
  expect_error(get_vars(wlddev, sapply(wlddev, is.numeric)[-1]) <- wlddev)
})

test_that("add_vars errors for wrong input", {
  expect_error(add_vars(wlddev, 15) <- wlddev[12])
  expect_error(add_vars(wlddev, "ODA3") <- wlddev[12])

  expect_error(add_vars(wlddev) <- qM(wlddev[9:12])[-1, ])
  expect_error(add_vars(wlddev, "front") <- qM(wlddev[9:12])[-1, ])
  expect_error(add_vars(wlddev, 8) <- qM(wlddev[9:12])[-1, ])

  expect_error(add_vars(wlddev) <- wlddev[-1, 9:12])
  expect_error(add_vars(wlddev, "front") <- wlddev[-1, 9:12])
  expect_error(add_vars(wlddev, 8) <- wlddev[-1, 9:12])

  expect_error(add_vars(wlddev, 12) <- wlddev[9:12])
  expect_error(add_vars(wlddev, 9:12) <- wlddev[9:10])
})

test_that("fselect errors for wrong input", {
    expect_visible(fselect(mtcars, 1))
    expect_error(fselect(mtcars, "bla"))
    expect_visible(fselect(mtcars, "mpg"))
    expect_error(fselect(mtcars, mpg:bla))
    expect_error(fselect(mtcars, mpg > cyl))
    expect_error(fselect(mtcars, ~mpg))
})


test_that("fselect works properly", {
  expect_equal(fselect(mtcars, mpg, 2), mtcars[1:2])
  expect_equal(fselect(mtcars, mpg:vs), mtcars[1:8])
  expect_equal(names(fselect(mtcars, bla = mpg, cyl:vs)), c("bla", names(mtcars)[2:8]))
  expect_invisible(fselect(wlddev, -PCGDP) <- fselect(wlddev, -PCGDP))
})


test_that("no problems with numeric values", {
  expect_equal(fselect(mtcars, 1), mtcars[1])
  expect_equal(get_vars(mtcars, 1), mtcars[1])
  expect_equal(gv(mtcars, 1), mtcars[1])

  expect_invisible(fselect(mtcars, 1) <- mtcars[1])
  expect_invisible(get_vars(mtcars, 1) <- mtcars[1])
  expect_invisible(gv(mtcars, 1) <- mtcars[1])
  expect_invisible(av(mtcars, pos = 1) <- mtcars[1])
})
