context("collapse and sf")

if(Sys.getenv("NMAC") == "TRUE") {

library(sf)
nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)

test_that("sf methods work properly", {
  expect_visible(nc %>% fgroup_by(AREA))
  expect_visible(nc %>% fgroup_by(AREA) %>% fgroup_vars)
  expect_visible(descr(nc))
  expect_visible(qsu(nc))
  expect_visible(varying(nc))
  expect_true(any(names(num_vars(nc)) == "geometry"))
  expect_true(any(names(fselect(nc, AREA, NAME:FIPSNO)) == "geometry"))
  expect_true(any(names(gv(nc, c("AREA", "NAME", "FIPS", "FIPSNO"))) == "geometry"))
  expect_true(any(names(fsubset(nc, AREA > fmean(AREA), AREA, NAME:FIPSNO)) == "geometry"))
  expect_true(any(names(ss(nc, 1:10, c("AREA", "NAME", "FIPS", "FIPSNO"))) == "geometry"))
  expect_true(inherits(rsplit(nc, AREA ~ SID74)[[1L]], "sf"))
  expect_equal(names(`nv<-`(nc, NULL)), c("NAME", "FIPS", "geometry"))
  # nv(nc) <- NULL
  expect_equal(tfmv(nc, is.numeric, log), tfmv(nc, is.numeric, log, apply = FALSE))
  expect_equal(length(nc %>% gby(NAME) %>% varying), length(nc) - 2L)
  expect_true(is.data.frame(nc %>% gby(NAME) %>% varying(any_group = FALSE)))
  expect_visible(funique(nc, cols = 1))
  expect_true(length(fcompute(nc, log_AREA = log(AREA))) == 2L)
  expect_true(length(fcomputev(nc, "AREA", log)) == 2L)
  expect_true(length(fcomputev(nc, "AREA", log, keep = "PERIMETER")) == 3L)
  expect_true(length(fcomputev(nc, "AREA", fscale, apply = FALSE)) == 2L)
  expect_true(length(fcomputev(nc, "AREA", fscale, apply = FALSE, keep = "PERIMETER")) == 3L)
  expect_true(inherits(nc %>% fgroup_by(SID74) %>%
                         fsummarise(AREA_Ag = fsum(AREA),
                                    Perimeter_Ag = fmedian(PERIMETER),
                                    geometry = st_union(geometry)), "sf"))
})

}
