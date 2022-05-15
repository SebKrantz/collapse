context("indexing")

wldi <- iby(wlddev, country, year)
library(magrittr)

test_that("unindexing and reindexing work well", {
  expect_equal(wlddev, unindex(wldi))
  expect_equal(wlddev$PCGDP, unindex(wldi$PCGDP))
  expect_equal(wlddev$region, unindex(wldi$region))

  expect_equal(wldi, reindex(wldi))
  expect_equal(wldi$PCGDP, reindex(wldi$PCGDP))
  expect_equal(wldi$region, reindex(wldi$region))

  expect_equal(wldi, reindex(wlddev, ix(wldi)))
  expect_equal(wldi$PCGDP, reindex(wldi$PCGDP, ix(wldi$PCGDP)))
  expect_equal(wldi$region, reindex(wldi$region, ix(wldi$region)))

})

test_that("subsetting works well", {

  expect_equal(fsubset(wldi, iso3c %in% c("KEN", "USA", "CHN")),
               findex_by(fsubset(wlddev, iso3c %in% c("KEN", "USA", "CHN")), country, year))

  expect_equal(fsubset(wldi, iso3c %in% c("KEN", "USA", "CHN"), country, year, PCGDP, POP),
               findex_by(fsubset(wlddev, iso3c %in% c("KEN", "USA", "CHN"), country, year, PCGDP, POP), country, year))

  expect_equal(wldi[wldi$iso3c %in% c("KEN", "USA", "CHN"), ] %>% setRownames(),
               ss(wlddev, wlddev$iso3c %in% c("KEN", "USA", "CHN")) %>% findex_by(country, year) %>% dapply(`attr<-`, "label", NULL))

  expect_true(all_obj_equal(wldi[.c(country, year, PCGDP, POP)],
                            wldi[, .c(country, year, PCGDP, POP)],
               wlddev[.c(country, year, PCGDP, POP)] %>% findex_by(country, year)))

  expect_equal(wldi$PCGDP[5:1000], reindex(wlddev$PCGDP[5:1000], ix(wldi)[5:1000, ]))
  expect_equal(wldi$PCGDP[100], wlddev$PCGDP[100])
  expect_equal(wldi$PCGDP[[100]], wlddev$PCGDP[[100]])

})

library(data.table)
wlddt <- qDT(wlddev)
wldidt <- wlddt %>% findex_by(iso3c, year)

test_that("indexed data.table works well", {

  expect_equal(unindex(wldidt[1:1000]), wlddt[1:1000])
  expect_equal(unindex(wldidt[year > 2000]), wlddt[year > 2000])
  expect_equal(wldidt[, .(sum_PCGDP = sum(PCGDP, na.rm = TRUE)), by = country],
               wlddt[, .(sum_PCGDP = sum(PCGDP, na.rm = TRUE)), by = country])

  expect_equal(wldidt[, lapply(.SD, sum, na.rm = TRUE), by = country, .SDcols = .c(PCGDP, LIFEEX)],
               wlddt[, lapply(.SD, sum, na.rm = TRUE), by = country, .SDcols = .c(PCGDP, LIFEEX)])

  expect_equal(wldidt[year > 2000, .(sum_PCGDP = sum(PCGDP, na.rm = TRUE)), by = country],
               wlddt[year > 2000, .(sum_PCGDP = sum(PCGDP, na.rm = TRUE)), by = country])

  # 'unclass' because of 'invisible' class
  expect_equal(unclass(unindex(wldidt[, PCGDP_growth_5Y := G(PCGDP, 5, power = 1/5)])),
               unclass(wlddt[, PCGDP_growth_5Y := G(PCGDP, 5, 1, iso3c, year, power = 1/5)]))

  expect_equal(unindex(wldidt[, PCGDP_growth_5Y := G(PCGDP, 5, power = 1/5)][1:5]),
               wlddt[, PCGDP_growth_5Y := G(PCGDP, 5, 1, iso3c, year, power = 1/5)][1:5])

  expect_equal(unindex(wldidt[, PCGDP_growth_5Y := G(PCGDP, 5, power = 1/5)][, .(sum_PCGDP = sum(PCGDP, na.rm = TRUE)), by = country]),
               wlddt[, PCGDP_growth_5Y := G(PCGDP, 5, 1, iso3c, year, power = 1/5)][, .(sum_PCGDP = sum(PCGDP, na.rm = TRUE)), by = country])

  expect_equal(unclass(unindex(wldidt[, .c(PCGDP_growth_5Y, LIFEEX_growth_5Y) := lapply(slt(.SD, PCGDP, LIFEEX), G, 5, power = 1/5)])),
               unclass(wlddt[, .c(PCGDP_growth_5Y, LIFEEX_growth_5Y) := lapply(slt(.SD, PCGDP, LIFEEX), G, 5, 1, iso3c, year, power = 1/5)]))

})


test_that("data selection by type works well", {
  for (FUN in list(num_vars, cat_vars, char_vars, logi_vars, fact_vars, date_vars))
     expect_equal(names(FUN(wlddev)), names(FUN(wldi)))
})

test_that("descriptives work well", {
  expect_equal(descr(wlddev), `attr<-`(descr(wldi), "name", "wlddev"))
  expect_equal(qsu(wlddev, pid = wlddev$country), qsu(wldi))
  expect_equal(varying(wlddev, by = ~country), varying(wldi))
  expect_equal(qtable(r = wlddev$region, i = wlddev$income), qtable(r = wldi$region, i = wldi$income))
  expect_equal(pwcor(nv(wlddev)), pwcor(nv(wldi)))
})

test_that("Id variables are properly preserved in operator methods", {
  wld1i <- findex_by(fsubset(wlddev, iso3c %==% "DEU"), year)
  GGDCii <- findex_by(GGDC10S, Variable, Country, Year)
  GGDCi <- findex_by(GGDC10S, Variable, Country, Year, interact.ids = FALSE)
  for(FUN in list(L, F, D, Dlog, G, B, W, STD)) {
      expect_identical(names(FUN(wld1i, cols = "PCGDP", stub = FALSE)), c("year", "PCGDP"))
      expect_identical(names(FUN(wld1i, cols = "PCGDP", keep.ids = FALSE, stub = FALSE)), "PCGDP")
      expect_identical(names(FUN(wldi, cols = "PCGDP", stub = FALSE)), c("country", "year", "PCGDP"))
      expect_identical(names(FUN(wldi, cols = "PCGDP", keep.ids = FALSE, stub = FALSE)), "PCGDP")
      expect_identical(names(FUN(GGDCi, cols = "SUM", stub = FALSE)), c("Country", "Variable", "Year", "SUM"))
      expect_identical(names(FUN(GGDCi, cols = "SUM", keep.ids = FALSE, stub = FALSE)), "SUM")
      expect_identical(names(FUN(GGDCii, cols = "SUM", stub = FALSE)), c("Country", "Variable", "Year", "SUM"))
      expect_identical(names(FUN(GGDCii, cols = "SUM", keep.ids = FALSE, stub = FALSE)), "SUM")
  }

})


