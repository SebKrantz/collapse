context("collapse and sf")

library(sf)
nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)

test_that("sf methods work properly", {
  expect_visible(nc %>% fgroup_by(AREA))
  expect_visible(nc %>% fgroup_by(AREA) %>% fgroup_vars)
  expect_visible(descr(nc))
  expect_visible(qsu(nc))
  expect_visible(varying(nc))
  expect_equal(names(`nv<-`(nc, NULL)), c("NAME", "FIPS", "geometry"))
  nv(nc) <- NULL
  expect_equal(tfmv(nc, is.numeric, log), tfmv(nc, is.numeric, log, apply = FALSE))
  expect_equal(length(nc %>% gby(NAME) %>% varying), length(nc) - 2L)
  expect_true(is.data.frame(nc %>% gby(NAME) %>% varying(any_group = FALSE)))
  expect_visible(funique(nc, cols = 1))
})
