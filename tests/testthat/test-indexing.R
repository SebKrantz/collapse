context("indexing")

wldi <- iby(wlddev, country, year)

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
