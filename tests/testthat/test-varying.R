context("varying")

if(!is.null(attributes(identical(FALSE, TRUE)))) stop("OECD label issue")

# rm(list = ls())

if(identical(Sys.getenv("NCRAN"), "TRUE")) pwlddev <- eval(parse(text = paste0("plm", ":", ":", "pdata.frame(wlddev, index = c('iso3c', 'year'))")))
gwlddev <- fgroup_by(wlddev, iso3c)
wdm <- qM(`cat_vars<-`(wlddev, dapply(cat_vars(wlddev), qG)))
g <- GRP(wlddev, ~ region + year)

test_that("vector, matrix and data.frame methods work as intended", {

  expect_true(all(dapply(wlddev, varying)))
  expect_true(all(varying(wlddev)))
  expect_true(all(varying(wdm)))
  expect_true(is.atomic(varying(wlddev, drop = TRUE)))
  expect_true(is.atomic(varying(wdm, drop = TRUE)))
  expect_true(is.data.frame(varying(wlddev, drop = FALSE)))
  expect_true(is.matrix(varying(wdm, drop = FALSE)))

  expect_true(all_identical(dapply(wlddev, varying), varying(wlddev), varying(wdm)))
  expect_true(all_identical(dapply(wlddev, varying, drop = FALSE), varying(wlddev, drop = FALSE), qDF(varying(wdm, drop = FALSE))))
  if(identical(Sys.getenv("NCRAN"), "TRUE")) {
  expect_equal(dapply(unattrib(wlddev), varying, wlddev$iso3c), c(FALSE,FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE))
  expect_true(all_identical(dapply(wlddev, varying, wlddev$iso3c), varying(wlddev, wlddev$iso3c),  varying(wdm, wlddev$iso3c)))
  expect_true(all_identical(dapply(wlddev, varying, wlddev$iso3c, drop = FALSE), varying(wlddev, wlddev$iso3c, drop = FALSE),  qDF(varying(wdm, wlddev$iso3c, drop = FALSE))))
  }
  expect_true(all_identical(qM(dapply(wlddev, varying, wlddev$iso3c, any_group = FALSE)),
                            qM(varying(wlddev,  wlddev$iso3c, any_group = FALSE)),
                            varying(wdm,  wlddev$iso3c, any_group = FALSE)))

  expect_true(all_identical(qM(dapply(wlddev, varying, wlddev$iso3c, any_group = FALSE, drop = FALSE)),
                            qM(varying(wlddev,  wlddev$iso3c, any_group = FALSE, drop = FALSE)),
                            varying(wdm,  wlddev$iso3c, any_group = FALSE, drop = FALSE)))

  expect_true(all_identical(qM(dapply(wlddev, varying, wlddev$iso3c, any_group = FALSE, use.g.names = FALSE)),
                            qM(varying(wlddev,  wlddev$iso3c, any_group = FALSE, use.g.names = FALSE)),
                            varying(wdm,  wlddev$iso3c, any_group = FALSE, use.g.names = FALSE)))

  expect_true(all_identical(qM(dapply(wlddev, varying, wlddev$iso3c, any_group = FALSE, use.g.names = FALSE, drop = FALSE)),
                            qM(varying(wlddev,  wlddev$iso3c, any_group = FALSE, use.g.names = FALSE, drop = FALSE)),
                            varying(wdm,  wlddev$iso3c, any_group = FALSE, use.g.names = FALSE, drop = FALSE)))

  # With grouping objects...
  if(identical(Sys.getenv("NCRAN"), "TRUE")) {
  expect_equal(dapply(unattrib(wlddev), varying, g), c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE))
  expect_true(all_identical(dapply(wlddev, varying, g), varying(wlddev, g),  varying(wdm, g)))
  expect_true(all_identical(dapply(wlddev, varying, g, drop = FALSE), varying(wlddev, g, drop = FALSE),  qDF(varying(wdm, g, drop = FALSE))))
  }
  expect_true(all_identical(qM(dapply(wlddev, varying, g, any_group = FALSE)),
                            qM(varying(wlddev,  g, any_group = FALSE)),
                            varying(wdm,  g, any_group = FALSE)))

  expect_true(all_identical(qM(dapply(wlddev, varying, g, any_group = FALSE, drop = FALSE)),
                            qM(varying(wlddev,  g, any_group = FALSE, drop = FALSE)),
                            varying(wdm,  g, any_group = FALSE, drop = FALSE)))

  expect_true(all_identical(qM(dapply(wlddev, varying, g, any_group = FALSE, use.g.names = FALSE)),
                            qM(varying(wlddev,  g, any_group = FALSE, use.g.names = FALSE)),
                            varying(wdm,  g, any_group = FALSE, use.g.names = FALSE)))

  expect_true(all_identical(qM(dapply(wlddev, varying, g, any_group = FALSE, use.g.names = FALSE, drop = FALSE)),
                            qM(varying(wlddev,  g, any_group = FALSE, use.g.names = FALSE, drop = FALSE)),
                            varying(wdm,  g, any_group = FALSE, use.g.names = FALSE, drop = FALSE)))


})


test_that("data.frame method formula and cols work as intended", {
  expect_equal(varying(wlddev, cols = 2:5), varying(get_vars(wlddev, 2:5)))
  expect_equal(varying(wlddev, cols = c("PCGDP","country")), varying(get_vars(wlddev, c("PCGDP","country"))))
  expect_equal(varying(wlddev, cols = is.numeric), varying(num_vars(wlddev)))

  expect_equal(varying(wlddev, ~iso3c), varying(fselect(wlddev, -iso3c), wlddev$iso3c))
  expect_equal(varying(wlddev, PCGDP + country ~ iso3c), varying(fselect(wlddev, PCGDP, country), wlddev$iso3c))
  expect_equal(varying(wlddev, PCGDP + country ~ iso3c), varying(wlddev, ~ iso3c, cols = c("PCGDP", "country")))
  expect_equal(varying(wlddev, ~iso3c, any_group = FALSE), varying(fselect(wlddev, -iso3c), wlddev$iso3c, any_group = FALSE))
  expect_equal(varying(wlddev, PCGDP + country ~ iso3c, any_group = FALSE), varying(fselect(wlddev, PCGDP, country), wlddev$iso3c, any_group = FALSE))
  expect_equal(varying(wlddev, PCGDP + country ~ iso3c, any_group = FALSE), varying(wlddev, ~ iso3c, cols = c("PCGDP", "country"), any_group = FALSE))

  expect_equal(varying(wlddev, ~region + year), varying(fselect(wlddev, -region, -year), g))
  expect_equal(varying(wlddev, PCGDP + country ~ region + year), varying(fselect(wlddev, PCGDP, country), g))
  expect_equal(varying(wlddev, PCGDP + country ~ region + year), varying(wlddev, ~ region + year, cols = c("PCGDP", "country")))
  expect_equal(varying(wlddev, ~region + year, any_group = FALSE), varying(fselect(wlddev, -region, -year),g, any_group = FALSE))
  expect_equal(varying(wlddev, PCGDP + country ~ region + year, any_group = FALSE), varying(fselect(wlddev, PCGDP, country), g, any_group = FALSE))
  expect_equal(varying(wlddev, PCGDP + country ~ region + year, any_group = FALSE), varying(wlddev, ~ region + year, cols = c("PCGDP", "country"), any_group = FALSE))

  expect_error(varying(wlddev, ~ iso3c2))
  expect_error(varying(wlddev, PCGDP + country ~ iso3c2))
  expect_error(varying(wlddev, PCGDP + country2 ~ iso3c))
  expect_error(varying(wlddev, ~ iso3c, cols = c("PCGDP", "country2")))

  expect_error(varying(wlddev, ~ region2 + year))
  expect_error(varying(wlddev, PCGDP + country ~ region2 + year))
  expect_error(varying(wlddev, PCGDP + country2 ~ region3 + year))
  expect_error(varying(wlddev, ~ region + year, cols = c("PCGDP", "country2")))
})

if(identical(Sys.getenv("NCRAN"), "TRUE")) {

test_that("pseries and pdata.frame methods work as intended", {
  # pdata.frame
  expect_equal(unattrib(varying(pwlddev)), c(FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE))
  expect_equal(unattrib(varying(pwlddev, effect = "iso3c")), c(FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE))
  expect_equal(unattrib(varying(pwlddev, effect = 2L)), c(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE))
  expect_equal(unattrib(varying(pwlddev, effect = "year")), c(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE))

  expect_true(is.atomic(varying(pwlddev, drop = TRUE)))
  expect_true(is.data.frame(varying(pwlddev, drop = FALSE)))
  expect_true(is.data.frame(varying(pwlddev, any_group = FALSE)))

  atrapply <- function(X, FUN, ...) {
    res <- vector("list", fncol(X))
    for(i in seq_col(X)) {
      res[[i]] <- FUN(X[[i]], ...)
    }
    res
  }

  # Making sure fselect and get_vars etc. work properly.
  expect_identical(attributes(fselect(pwlddev, country:POP)), attributes(pwlddev))
  expect_identical(attributes(get_vars(pwlddev, seq_col(pwlddev))), attributes(pwlddev))

  # pseries
  expect_equal(unlist(atrapply(fselect(pwlddev, -iso3c), varying)), c(FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE))
  expect_equal(unlist(atrapply(fselect(pwlddev, -iso3c), varying, effect = "iso3c")), c(FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE))
  expect_equal(unlist(atrapply(fselect(pwlddev, -year), varying, effect = 2L)), c(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE))
  expect_equal(unlist(atrapply(fselect(pwlddev, -year), varying, effect = "year")), c(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE))

  expect_equal(varying(pwlddev$PCGDP), varying(wlddev$PCGDP, wlddev$iso3c))
  expect_equal(varying(pwlddev$PCGDP, any_group = FALSE), varying(wlddev$PCGDP, wlddev$iso3c, any_group = FALSE))
  expect_equal(varying(pwlddev$PCGDP, any_group = FALSE, use.g.names = FALSE), varying(wlddev$PCGDP, wlddev$iso3c, any_group = FALSE, use.g.names = FALSE))
  expect_equal(lengths(varying(pwlddev, any_group = FALSE), FALSE), lengths(atrapply(fselect(pwlddev, -iso3c), varying, any_group = FALSE)))

  # pdata.frame works like data.frame
  expect_identical(unattrib(varying(fselect(wlddev, -iso3c), wlddev$iso3c, any_group = FALSE)),
                   unattrib(varying(pwlddev, any_group = FALSE)))

  expect_identical(unattrib(varying(fselect(wlddev, -iso3c), wlddev$iso3c, any_group = FALSE, drop = FALSE)),
                   unattrib(varying(pwlddev, any_group = FALSE, drop = FALSE)))

  expect_identical(unattrib(varying(fselect(wlddev, -iso3c), wlddev$iso3c, any_group = FALSE, use.g.names = FALSE)),
                   unattrib(varying(pwlddev, any_group = FALSE, use.g.names = FALSE)))

  expect_identical(unattrib(varying(fselect(wlddev, -iso3c), wlddev$iso3c, any_group = FALSE, use.g.names = FALSE, drop = FALSE)),
                   unattrib(varying(pwlddev, any_group = FALSE, use.g.names = FALSE, drop = FALSE)))

  expect_identical(unattrib(varying(fselect(wlddev, -year), wlddev$year, any_group = FALSE)),
                   unattrib(varying(pwlddev, any_group = FALSE, effect = "year")))

  expect_identical(unattrib(varying(fselect(wlddev, -year), wlddev$year, any_group = FALSE, drop = FALSE)),
                   unattrib(varying(pwlddev, any_group = FALSE, drop = FALSE, effect = "year")))

  expect_identical(unattrib(varying(fselect(wlddev, -year), wlddev$year, any_group = FALSE, use.g.names = FALSE)),
                   unattrib(varying(pwlddev, any_group = FALSE, use.g.names = FALSE, effect = "year")))

  expect_identical(unattrib(varying(fselect(wlddev, -year), wlddev$year, any_group = FALSE, use.g.names = FALSE, drop = FALSE)),
                   unattrib(varying(pwlddev, any_group = FALSE, use.g.names = FALSE, drop = FALSE, effect = "year")))

})

}

test_that("grouped_df method works as intended", {

  expect_equal(unattrib(varying(gwlddev)), c(FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE))
  expect_true(is.atomic(varying(gwlddev, drop = TRUE)))
  expect_true(is.data.frame(varying(gwlddev, drop = FALSE)))
  expect_true(is.data.frame(varying(gwlddev, any_group = FALSE)))

  expect_identical(names(varying(gwlddev)), names(wlddev)[-2L])
  expect_identical(names(varying(get_vars(gwlddev, 9:12))), names(wlddev)[9:12])

  expect_identical(names(varying(gwlddev, any_group = FALSE)), c("iso3c", names(wlddev)[-2L]))
  expect_identical(names(varying(gwlddev, any_group = FALSE, keep.group_vars = FALSE)), names(wlddev)[-2L])
  expect_identical(names(varying(get_vars(gwlddev, 9:12), any_group = FALSE)), c("iso3c", names(wlddev)[9:12]))
  expect_identical(names(varying(get_vars(gwlddev, 9:12), any_group = FALSE, keep.group_vars = FALSE)), names(wlddev)[9:12])

  # grouped_df works like data.frame
  expect_identical(unattrib(varying(fselect(wlddev, -iso3c), wlddev$iso3c, any_group = FALSE)),
                   unattrib(varying(gwlddev, any_group = FALSE, keep.group_vars = FALSE)))

  expect_identical(unattrib(varying(fselect(wlddev, -iso3c), wlddev$iso3c, any_group = FALSE, drop = FALSE)),
                   unattrib(varying(gwlddev, any_group = FALSE, drop = FALSE, keep.group_vars = FALSE)))

  expect_identical(unclass(varying(fselect(wlddev, -iso3c), wlddev$iso3c, any_group = FALSE, use.g.names = FALSE)),
                   unclass(fungroup(varying(gwlddev, any_group = FALSE, use.g.names = FALSE, keep.group_vars = FALSE))))

  expect_identical(unclass(varying(fselect(wlddev, -iso3c), wlddev$iso3c, any_group = FALSE, drop = FALSE)),
                   unclass(fungroup(varying(gwlddev, any_group = FALSE, use.g.names = TRUE, drop = FALSE, keep.group_vars = FALSE))))

})
