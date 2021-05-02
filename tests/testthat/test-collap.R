context("collap")

# rm(list = ls())

options(warn = -1)

g <- GRP(wlddev, ~ country + decade)

Mode <- function(x, na.rm = FALSE) {
  if(na.rm) x <- x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# TODO: What about other return options and weighted multi-function aggregation ? And what about grouped_df method..

test_that("collap performs as intended in simple uses", {
  expect_equal(collap(mtcars, mtcars$cyl, keep.by = FALSE), fmean(mtcars, mtcars$cyl, use.g.names = FALSE))
  expect_equal(collap(mtcars, mtcars[2], keep.by = FALSE), fmean(mtcars, mtcars$cyl, use.g.names = FALSE))
  expect_equal(collap(mtcars, ~cyl), fmean(mtcars, mtcars$cyl, use.g.names = FALSE))
  expect_equal(collap(mtcars, ~cyl, keep.by = FALSE), fmean(mtcars[-2], mtcars$cyl, use.g.names = FALSE))
  expect_equal(collap(iris, ~Species, keep.by = FALSE), fmean(iris[-5], iris$Species, use.g.names = FALSE))
  expect_equal(collap(airquality, ~Month, keep.by = FALSE), fmean(airquality[-5], airquality$Month, use.g.names = FALSE))

  expect_equal(collap(wlddev, ~ country + decade, keep.col.order = FALSE),
               cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))
  expect_equal(collap(wlddev, ~ country + decade),
               cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))[order(c(1,5,4,9:12,2:3,6:8))])
  expect_equal(collap(wlddev, ~ country + decade, keep.col.order = FALSE, keep.by = FALSE),
               cbind(fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))
  expect_equal(collap(wlddev, ~ country + decade, keep.by = FALSE),
               cbind(fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))[order(c(4,9:12,2:3,6:8))])

  expect_equal(collap(wlddev, g, keep.by = FALSE),
               cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))[order(c(1,5,4,9:12,2:3,6:8))])

  expect_equal(collap(wlddev, wlddev[c("country","decade")], keep.by = FALSE),
               cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))[order(c(1,5,4,9:12,2:3,6:8))])

})

storage.mode(wlddev$year) <- "numeric" # Makes sure that no mismatch as year is aggregated
test_that("collap preserves data attributes", {
  expect_identical(lapply(collap(wlddev, ~country), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collap(wlddev, ~country)), vclasses(wlddev))
  expect_identical(vtypes(collap(wlddev, ~country)), vtypes(wlddev))
  expect_identical(lapply(collap(wlddev, ~iso3c), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collap(wlddev, ~iso3c)), vclasses(wlddev))
  expect_identical(vtypes(collap(wlddev, ~iso3c)), vtypes(wlddev))
  expect_identical(lapply(collap(wlddev, ~date), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collap(wlddev, ~date)), vclasses(wlddev))
  expect_identical(vtypes(collap(wlddev, ~date)), vtypes(wlddev))
  expect_identical(lapply(collap(wlddev, ~country + decade), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collap(wlddev, ~country + decade)), vclasses(wlddev))
  expect_identical(vtypes(collap(wlddev, ~country + decade)), vtypes(wlddev))
})

test_that("collap performs as intended in simple uses with base/stats functions", {
  expect_equal(collap(mtcars, mtcars$cyl, sum, keep.by = FALSE), fsum(mtcars, mtcars$cyl, use.g.names = FALSE))
  expect_equal(collap(mtcars, ~cyl, mean.default), fmean(mtcars, mtcars$cyl, use.g.names = FALSE))
  expect_equal(collap(mtcars, ~cyl, mean), fmean(mtcars, mtcars$cyl, use.g.names = FALSE)) # error !!

  expect_equal(collap(mtcars, mtcars[2], sum, keep.by = FALSE), fsum(mtcars, mtcars$cyl, use.g.names = FALSE))
  expect_equal(collap(mtcars, ~cyl, sum, keep.by = FALSE), fsum(mtcars[-2], mtcars$cyl, use.g.names = FALSE))
  expect_equal(collap(iris, ~Species, sum, keep.by = FALSE), fsum(iris[-5], iris$Species, use.g.names = FALSE))
  expect_equal(collap(airquality, ~Month, sum, na.rm = TRUE, keep.by = FALSE), fsum(airquality[-5], airquality$Month, use.g.names = FALSE))

  expect_equal(collap(wlddev, ~ country + decade, sum, Mode, na.rm = TRUE, keep.col.order = FALSE),
               cbind(g$groups, BY(get_vars(wlddev, c(4,9:12)), g, sum, na.rm = TRUE, use.g.names = FALSE),
                     BY(get_vars(wlddev, c(2:3,6:8)), g, Mode, na.rm = TRUE, use.g.names = FALSE)))
  expect_equal(collap(wlddev, ~ country + decade, sum, Mode, na.rm = TRUE),
               cbind(g$groups, BY(get_vars(wlddev, c(4,9:12)), g, sum, na.rm = TRUE, use.g.names = FALSE),
                     BY(get_vars(wlddev, c(2:3,6:8)), g, Mode, na.rm = TRUE, use.g.names = FALSE))[order(c(1,5,4,9:12,2:3,6:8))])
})

test_that("collap using 2-sided formula or cols performs as intended", {
  expect_equal(collap(mtcars, mpg ~ cyl, keep.by = FALSE), fmean(mtcars["mpg"], mtcars$cyl, use.g.names = FALSE))
  expect_equal(collap(mtcars, mpg ~ cyl, keep.by = FALSE, cols = 300:1000), fmean(mtcars["mpg"], mtcars$cyl, use.g.names = FALSE)) # cols is ignored, as should be
  expect_equal(collap(mtcars, ~ cyl, keep.by = FALSE, cols = 1), fmean(mtcars["mpg"], mtcars$cyl, use.g.names = FALSE))
  expect_equal(collap(mtcars, wt + mpg ~ cyl + vs + am, keep.by = FALSE), fmean(mtcars[c("mpg","wt")], mtcars[c("cyl","vs","am")], use.g.names = FALSE))
  expect_equal(collap(mtcars, ~ cyl + vs + am, keep.by = FALSE, cols = c(6,1)), fmean(mtcars[c("mpg","wt")], mtcars[c("cyl","vs","am")], use.g.names = FALSE))
  expect_equal(collap(iris, Sepal.Length + Sepal.Width ~ Species, keep.by = FALSE), fmean(iris[1:2], iris$Species, use.g.names = FALSE))
  expect_equal(collap(airquality, ~ Month, keep.by = FALSE), fmean(airquality[-5], airquality$Month, use.g.names = FALSE))
  expect_equal(collap(airquality, ~ Month, keep.by = FALSE, cols = 1:3), fmean(airquality[1:3], airquality$Month, use.g.names = FALSE))

  expect_equal(collap(wlddev, ~ country + decade, cols = 9:12), collap(wlddev, PCGDP + LIFEEX + GINI + ODA ~ country + decade))
  expect_equal(collap(wlddev, ~ country + decade, cols = 9:12), collap(wlddev, ~ country + decade, cols = 9:12, keep.col.order = FALSE))
  expect_equal(collap(wlddev, ~ country + decade, cols = c(2:3,6:8)), collap(wlddev, iso3c + date + region + income + OECD ~ country + decade))
  expect_false(identical(collap(wlddev, ~ country + decade, cols = c(2:3,6:8)), collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.col.order = FALSE)))

  expect_equal(collap(wlddev, ~ country + decade, cols = 9:12, keep.by = FALSE), collap(wlddev, PCGDP + LIFEEX + GINI + ODA ~ country + decade, keep.by = FALSE))
  expect_equal(collap(wlddev, ~ country + decade, cols = 9:12, keep.by = FALSE), collap(wlddev, ~ country + decade, cols = 9:12, keep.col.order = FALSE, keep.by = FALSE))
  expect_equal(collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.by = FALSE), collap(wlddev, iso3c + date + region + income + OECD ~ country + decade, keep.by = FALSE))
  expect_equal(collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.by = FALSE), collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.col.order = FALSE, keep.by = FALSE))

  expect_equal(collap(wlddev, g, cols = 9:12, keep.by = FALSE), collap(wlddev, PCGDP + LIFEEX + GINI + ODA ~ country + decade, keep.by = FALSE))
  expect_equal(collap(wlddev, g, cols = 9:12, keep.by = FALSE), collap(wlddev, ~ country + decade, cols = 9:12, keep.col.order = FALSE, keep.by = FALSE))
  expect_equal(collap(wlddev, g, cols = c(2:3,6:8), keep.by = FALSE), collap(wlddev, iso3c + date + region + income + OECD ~ country + decade, keep.by = FALSE))
  expect_equal(collap(wlddev, g, cols = c(2:3,6:8), keep.by = FALSE), collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.col.order = FALSE, keep.by = FALSE))

  expect_equal(collap(wlddev, wlddev[c("country","decade")], cols = 9:12, keep.by = FALSE), collap(wlddev, PCGDP + LIFEEX + GINI + ODA ~ country + decade, keep.by = FALSE))
  expect_equal(collap(wlddev, wlddev[c("country","decade")], cols = 9:12, keep.by = FALSE), collap(wlddev, ~ country + decade, cols = 9:12, keep.col.order = FALSE, keep.by = FALSE))
  expect_equal(collap(wlddev, wlddev[c("country","decade")], cols = c(2:3,6:8), keep.by = FALSE), collap(wlddev, iso3c + date + region + income + OECD ~ country + decade, keep.by = FALSE))
  expect_equal(collap(wlddev, wlddev[c("country","decade")], cols = c(2:3,6:8), keep.by = FALSE), collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.col.order = FALSE, keep.by = FALSE))

})

test_that("collap multi-function aggreagtion performs as intended", {
  expect_equal(collap(wlddev, ~ country + decade, list(fmean, fmedian), keep.col.order = FALSE, give.names = FALSE),
               cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))
  expect_equal(collap(wlddev, ~ country + decade, list(fmean, fmedian), list(fmode, flast), keep.col.order = FALSE, give.names = FALSE),
               cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE), flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))
  # with column ordering:
  expect_equal(unname(collap(wlddev, ~ country + decade, list(fmean, fmedian))),
               unname(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))[order(c(1,5,4,9:12,4,9:12,2:3,6:8))]))
  expect_equal(unname(collap(wlddev, ~ country + decade, list(fmean, fmedian), list(fmode, flast))),
               unname(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                      fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE), flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))[order(c(1,5,4,9:12,4,9:12,2:3,6:8,2:3,6:8))]))
})

test_that("collap custom aggreagtion performs as intended", {
  expect_equal(unname(collap(wlddev, ~ country + decade,
                      custom = list(fmean = 9:12, fsd = 9:10, fmode = 7:8), keep.col.order = FALSE)),
               unname(cbind(g$groups, fmean(wlddev[9:12], g, use.g.names = FALSE),
                                      fsd(wlddev[9:10], g, use.g.names = FALSE),
                                      fmode(wlddev[7:8], g, use.g.names = FALSE))))
  expect_equal(unname(collap(wlddev, ~ country + decade,
                             custom = list(fmean = 9:12, fsd = 9:10, fmode = 7:8))),
               unname(cbind(g$groups, fmean(wlddev[9:12], g, use.g.names = FALSE),
                            fsd(wlddev[9:10], g, use.g.names = FALSE),
                            fmode(wlddev[7:8], g, use.g.names = FALSE))[order(c(1,5,9:12,9:10,7:8))]))

  expect_equal(collap(wlddev, ~ country + decade,
               custom = list(fmean = 9:12, fsd = 9:10, fmode = 7:8)),
               collap(wlddev, ~ country + decade,
               custom = list(fmean = 9:12, fsd = c("PCGDP","LIFEEX"), fmode = 7:8)))
  expect_equal(collap(wlddev, ~ country + decade,
               custom = list(fmean = "PCGDP", fsd = c("LIFEEX","GINI"), flast = "date")),
               collap(wlddev, ~ country + decade,
               custom = list(fmean = "PCGDP", fsd = 10:11, flast = "date")))

  expect_equal(collap(wlddev, g,
                      custom = list(fmean = 9:12, fsd = 9:10, fmode = 7:8)),
               collap(wlddev, ~ country + decade,
                      custom = list(fmean = 9:12, fsd = c("PCGDP","LIFEEX"), fmode = 7:8)))
  expect_equal(collap(wlddev, g,
                      custom = list(fmean = "PCGDP", fsd = c("LIFEEX","GINI"), flast = "date")),
               collap(wlddev, g,
                      custom = list(fmean = "PCGDP", fsd = 10:11, flast = "date")))
  expect_equal(names(collap(wlddev, g,
                      custom = list(fmean = c(GDP = "PCGDP"), fsd = c("LIFEEX", GN = "GINI"), flast = "date"),
                      keep.by = FALSE, keep.col.order = FALSE)),
               .c(GDP, LIFEEX, GN, date))

})

test_that("collap weighted aggregations work as intended", {
  # Not keeping order ...
  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, keep.col.order = FALSE),
                 add_vars(g$groups,
                  fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                  fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                  fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, keep.col.order = FALSE, keep.by = FALSE),
               add_vars(fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, keep.col.order = FALSE, keep.w = FALSE),
               add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, keep.col.order = FALSE, keep.by = FALSE, keep.w = FALSE),
               add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  # keeping order ...
  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA),
               add_vars(g$groups,
                        fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[names(wlddev)])

  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, keep.by = FALSE),
               add_vars(fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[setdiff(names(wlddev), g$group.vars)])

  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, keep.w = FALSE),
               add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[setdiff(names(wlddev), "ODA")])

  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, keep.by = FALSE, keep.w = FALSE),
               add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[setdiff(names(wlddev), c(g$group.vars, "ODA"))])


})

test_that("collap multi-function aggreagtion with weights performs as intended", {
  expect_equal(collap(wlddev, ~ country + decade, list(fmean, fsd), w = ~ ODA, keep.col.order = FALSE, give.names = FALSE),
               cbind(g$groups, fsum(get_vars(wlddev, 12), g, use.g.names = FALSE), fmean(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                     fsd(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$ODA, use.g.names = FALSE)))
  expect_equal(collap(wlddev, ~ country + decade, list(fmean, fsd), list(fmode, flast), w = ~ ODA, keep.col.order = FALSE, give.names = FALSE),
               cbind(g$groups, fsum(get_vars(wlddev, 12), g, use.g.names = FALSE), fmean(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                     fsd(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$ODA, use.g.names = FALSE), flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))
  # with column ordering:
  expect_equal(unname(collap(wlddev, ~ country + decade, list(fmean, fsd), w = ~ ODA, wFUN = list(fsum, fmax))),
               unname(cbind(g$groups,
                            fsum(get_vars(wlddev, 12), g, use.g.names = FALSE),
                            fmax(get_vars(wlddev, 12), g, use.g.names = FALSE),
                            fmean(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                            fsd(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                            fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$ODA, use.g.names = FALSE))[order(c(1,5,12,12,4,9:11,4,9:11,2:3,6:8))]))
  expect_equal(unname(collap(wlddev, ~ country + decade, list(fmean, fsd), list(fmode, flast), w = ~ ODA, wFUN = list(fsum, fmax))),
               unname(cbind(g$groups,
                            fsum(get_vars(wlddev, 12), g, use.g.names = FALSE),
                            fmax(get_vars(wlddev, 12), g, use.g.names = FALSE),
                            fmean(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                            fsd(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                            fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$ODA, use.g.names = FALSE),
                            flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))[order(c(1,5,12,12,4,9:11,4,9:11,2:3,6:8,2:3,6:8))]))
})

test_that("collap weighted customized aggregation works as intended", {
  v1 <- c("year","PCGDP","LIFEEX","GINI")
  v2 <- c("iso3c","date","region","income", "OECD")
  # Not keeping order ...
  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, give.names = FALSE),
               add_vars(g$groups,
                        fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, keep.by = FALSE, give.names = FALSE),
               add_vars(fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, keep.w = FALSE, give.names = FALSE),
               add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, keep.by = FALSE, keep.w = FALSE, give.names = FALSE),
               add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  # keeping order ...
  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, custom = list(fmean = v1, fmode = v2), give.names = FALSE),
               add_vars(g$groups,
                        fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[names(wlddev)])

  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, keep.by = FALSE, custom = list(fmean = v1, fmode = v2), give.names = FALSE),
               add_vars(fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[setdiff(names(wlddev), g$group.vars)])

  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, keep.w = FALSE, custom = list(fmean = v1, fmode = v2), give.names = FALSE),
               add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[setdiff(names(wlddev), "ODA")])

  expect_equal(collap(wlddev, ~ country + decade, w = ~ ODA, keep.by = FALSE, keep.w = FALSE, custom = list(fmean = v1, fmode = v2), give.names = FALSE),
               add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[setdiff(names(wlddev), c(g$group.vars, "ODA"))])


})

test_that("collap gives informative errors", {
  expect_error(collap(~cyl, ~cyl)) # nah, need to give error in qDF
  expect_error(collap(wlddev, 1:3)) # only gives error in fmean.. a bit late..
  expect_error(collap(wlddev, "country")) # same thing
  expect_error(collap(wlddev, ~ country1))
  expect_error(collap(wlddev, ~ country, w = ~bla))
  expect_error(collap(wlddev, ~ country, w = ~ODA, wFUN = sum))
  expect_error(collap(wlddev, ~ country + year + bla))
  expect_error(collap(wlddev, bla ~ country))
  expect_warning(collap(wlddev, ~ country, bla = 1)) # passes to fmean.data.frame which give the error.
  # expect_error(collap(wlddev, ~ country, sum, cols = 9:12, bla = 1)) # This is an issue, sum(1:3, bla = 1) does not give an error
  expect_error(collap(wlddev, mtcars$cyl)) # again fmean error..
  expect_error(collap(wlddev, ~iso3c, cols = 9:13))
  # expect_error(collap(wlddev, ~iso3c, cols = 0:1)) # no error..
  expect_error(collap(wlddev, ~iso3c, cols = c("PCGDP","bla")))
  expect_error(collap(wlddev, ~iso3c, cols = c("PCGDP","LIFEEX1")))
  expect_error(collap(wlddev, ~iso3c, custom = ~ PCGDP))
  expect_error(collap(wlddev, ~iso3c, custom = list(fmean, fmode)))
  expect_error(collap(wlddev, ~iso3c, custom = list(fmean = 9:13, fmode = 4:6)))
  expect_error(collap(wlddev, ~iso3c, custom = list(fmean = 9:12, 4:6)))
  expect_error(collap(wlddev, ~iso3c, custom = list(fmean = 9:12, fmode2 = 4:6)))
  expect_error(collap(wlddev, ~iso3c, custom = list(fmean = 9:12, fmode = c("GINI","bla"))))
  expect_error(collap(wlddev, ~iso3c, custom = list(fmean = 9:12, fmode = c("GINI","PCGDP2"))))
})

# Note: one more thing to test is performance with vector-valued functions...


# Testing collapv
v <- c(1, 5)

test_that("collapv performs as intended in simple uses", {
  expect_equal(collapv(mtcars, 2), fmean(mtcars, mtcars$cyl, use.g.names = FALSE))
  expect_equal(collapv(mtcars, 2, keep.by = FALSE), fmean(mtcars[-2], mtcars$cyl, use.g.names = FALSE))
  expect_equal(collapv(iris, "Species", keep.by = FALSE), fmean(iris[-5], iris$Species, use.g.names = FALSE))
  expect_equal(collapv(airquality, "Month", keep.by = FALSE), fmean(airquality[-5], airquality$Month, use.g.names = FALSE))

  expect_equal(collapv(wlddev, v, keep.col.order = FALSE),
               cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))
  expect_equal(collapv(wlddev, v),
               cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))[order(c(1,5,4,9:12,2:3,6:8))])
  expect_equal(collapv(wlddev, v, keep.col.order = FALSE, keep.by = FALSE),
               cbind(fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))
  expect_equal(collapv(wlddev, v, keep.by = FALSE),
               cbind(fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))[order(c(4,9:12,2:3,6:8))])

})

storage.mode(wlddev$year) <- "numeric" # Makes sure that no mismatch as year is aggregated
test_that("collapv preserves data attributes", {
  expect_identical(lapply(collapv(wlddev, 1), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collapv(wlddev, 1)), vclasses(wlddev))
  expect_identical(vtypes(collapv(wlddev, 1)), vtypes(wlddev))
  expect_identical(lapply(collapv(wlddev, "iso3c"), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collapv(wlddev, "iso3c")), vclasses(wlddev))
  expect_identical(vtypes(collapv(wlddev, "iso3c")), vtypes(wlddev))
  expect_identical(lapply(collapv(wlddev, "date"), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collapv(wlddev, "date")), vclasses(wlddev))
  expect_identical(vtypes(collapv(wlddev, "date")), vtypes(wlddev))
  expect_identical(lapply(collapv(wlddev, v), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collapv(wlddev, v)), vclasses(wlddev))
  expect_identical(vtypes(collapv(wlddev, v)), vtypes(wlddev))
})

test_that("collapv performs as intended in simple uses with base/stats functions", {
  expect_equal(collapv(mtcars, "cyl", mean.default), fmean(mtcars, mtcars$cyl, use.g.names = FALSE))
  expect_equal(collapv(mtcars, "cyl", mean), fmean(mtcars, mtcars$cyl, use.g.names = FALSE))

  expect_equal(collapv(mtcars, 2, sum, keep.by = FALSE), fsum(mtcars[-2], mtcars$cyl, use.g.names = FALSE))
  expect_equal(collapv(iris, 5, sum, keep.by = FALSE), fsum(iris[-5], iris$Species, use.g.names = FALSE))
  expect_equal(collapv(airquality, "Month", sum, na.rm = TRUE, keep.by = FALSE), fsum(airquality[-5], airquality$Month, use.g.names = FALSE))

  expect_equal(collapv(wlddev, v, sum, Mode, na.rm = TRUE, keep.col.order = FALSE),
               cbind(g$groups, BY(get_vars(wlddev, c(4,9:12)), g, sum, na.rm = TRUE, use.g.names = FALSE),
                     BY(get_vars(wlddev, c(2:3,6:8)), g, Mode, na.rm = TRUE, use.g.names = FALSE)))
  expect_equal(collapv(wlddev, v, sum, Mode, na.rm = TRUE),
               cbind(g$groups, BY(get_vars(wlddev, c(4,9:12)), g, sum, na.rm = TRUE, use.g.names = FALSE),
                     BY(get_vars(wlddev, c(2:3,6:8)), g, Mode, na.rm = TRUE, use.g.names = FALSE))[order(c(1,5,4,9:12,2:3,6:8))])
})

test_that("collapv using cols performs as intended", {
  expect_equal(collapv(mtcars, 2, keep.by = FALSE, cols = 1), fmean(mtcars["mpg"], mtcars$cyl, use.g.names = FALSE))
  expect_equal(collapv(mtcars, c("cyl", "vs", "am"), keep.by = FALSE, cols = c(6,1)), fmean(mtcars[c("mpg","wt")], mtcars[c("cyl","vs","am")], use.g.names = FALSE))
  expect_equal(collapv(airquality, "Month", keep.by = FALSE), fmean(airquality[-5], airquality$Month, use.g.names = FALSE))
  expect_equal(collapv(airquality, "Month", keep.by = FALSE, cols = 1:3), fmean(airquality[1:3], airquality$Month, use.g.names = FALSE))

  expect_equal(collapv(wlddev, v, cols = 9:12), collap(wlddev, PCGDP + LIFEEX + GINI + ODA ~ country + decade))
  expect_equal(collapv(wlddev, v, cols = 9:12), collap(wlddev, ~ country + decade, cols = 9:12, keep.col.order = FALSE))
  expect_equal(collapv(wlddev, v, cols = c(2:3,6:8)), collap(wlddev, iso3c + date + region + income + OECD ~ country + decade))
  expect_false(identical(collapv(wlddev, v, cols = c(2:3,6:8)), collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.col.order = FALSE)))

  expect_equal(collapv(wlddev, v, cols = 9:12, keep.by = FALSE), collap(wlddev, PCGDP + LIFEEX + GINI + ODA ~ country + decade, keep.by = FALSE))
  expect_equal(collapv(wlddev, v, cols = 9:12, keep.by = FALSE), collap(wlddev, ~ country + decade, cols = 9:12, keep.col.order = FALSE, keep.by = FALSE))
  expect_equal(collapv(wlddev, v, cols = c(2:3,6:8), keep.by = FALSE), collap(wlddev, iso3c + date + region + income + OECD ~ country + decade, keep.by = FALSE))
  expect_equal(collapv(wlddev, v, cols = c(2:3,6:8), keep.by = FALSE), collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.col.order = FALSE, keep.by = FALSE))

})

test_that("collapv multi-function aggreagtion performs as intended", {
  expect_equal(collapv(wlddev, v, list(fmean, fmedian), keep.col.order = FALSE, give.names = FALSE),
               cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))
  expect_equal(collapv(wlddev, v, list(fmean, fmedian), list(fmode, flast), keep.col.order = FALSE, give.names = FALSE),
               cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE), flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))
  # with column ordering:
  expect_equal(unname(collapv(wlddev, v, list(fmean, fmedian))),
               unname(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                            fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))[order(c(1,5,4,9:12,4,9:12,2:3,6:8))]))
  expect_equal(unname(collapv(wlddev, v, list(fmean, fmedian), list(fmode, flast))),
               unname(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:12)), g, use.g.names = FALSE),
                            fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE), flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))[order(c(1,5,4,9:12,4,9:12,2:3,6:8,2:3,6:8))]))
})

test_that("collapv custom aggreagtion performs as intended", {
  expect_equal(unname(collapv(wlddev, v,
                             custom = list(fmean = 9:12, fsd = 9:10, fmode = 7:8), keep.col.order = FALSE)),
               unname(cbind(g$groups, fmean(wlddev[9:12], g, use.g.names = FALSE),
                            fsd(wlddev[9:10], g, use.g.names = FALSE),
                            fmode(wlddev[7:8], g, use.g.names = FALSE))))
  expect_equal(unname(collapv(wlddev, v,
                             custom = list(fmean = 9:12, fsd = 9:10, fmode = 7:8))),
               unname(cbind(g$groups, fmean(wlddev[9:12], g, use.g.names = FALSE),
                            fsd(wlddev[9:10], g, use.g.names = FALSE),
                            fmode(wlddev[7:8], g, use.g.names = FALSE))[order(c(1,5,9:12,9:10,7:8))]))

  expect_equal(collapv(wlddev, v,
                      custom = list(fmean = 9:12, fsd = 9:10, fmode = 7:8)),
               collapv(wlddev, v,
                      custom = list(fmean = 9:12, fsd = c("PCGDP","LIFEEX"), fmode = 7:8)))
  expect_equal(collapv(wlddev, v,
                      custom = list(fmean = "PCGDP", fsd = c("LIFEEX","GINI"), flast = "date")),
               collapv(wlddev, v,
                      custom = list(fmean = "PCGDP", fsd = 10:11, flast = "date")))

  expect_equal(collapv(wlddev, v,
                      custom = list(fmean = 9:12, fsd = 9:10, fmode = 7:8)),
               collapv(wlddev, v,
                      custom = list(fmean = 9:12, fsd = c("PCGDP","LIFEEX"), fmode = 7:8)))
  expect_equal(collapv(wlddev, v,
                      custom = list(fmean = "PCGDP", fsd = c("LIFEEX","GINI"), flast = "date")),
               collapv(wlddev, v,
                      custom = list(fmean = "PCGDP", fsd = 10:11, flast = "date")))

})

test_that("collapv weighted aggregations work as intended", {
  # Not keeping order ...
  expect_equal(collapv(wlddev, v, w = "ODA", keep.col.order = FALSE),
               add_vars(g$groups,
                        fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  expect_equal(collapv(wlddev, v, w = "ODA", keep.col.order = FALSE, keep.by = FALSE),
               add_vars(fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  expect_equal(collapv(wlddev, v, w = "ODA", keep.col.order = FALSE, keep.w = FALSE),
               add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  expect_equal(collapv(wlddev, v, w = "ODA", keep.col.order = FALSE, keep.by = FALSE, keep.w = FALSE),
               add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  # keeping order ...
  expect_equal(collapv(wlddev, v, w = "ODA"),
               add_vars(g$groups,
                        fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[names(wlddev)])

  expect_equal(collapv(wlddev, v, w = "ODA", keep.by = FALSE),
               add_vars(fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[setdiff(names(wlddev), g$group.vars)])

  expect_equal(collapv(wlddev, v, w = "ODA", keep.w = FALSE),
               add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[setdiff(names(wlddev), "ODA")])

  expect_equal(collapv(wlddev, v, w = "ODA", keep.by = FALSE, keep.w = FALSE),
               add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[setdiff(names(wlddev), c(g$group.vars, "ODA"))])


})

test_that("collapv multi-function aggreagtion with weights performs as intended", {
  expect_equal(collapv(wlddev, v, list(fmean, fsd), w = "ODA", keep.col.order = FALSE, give.names = FALSE),
               cbind(g$groups, fsum(get_vars(wlddev, 12), g, use.g.names = FALSE), fmean(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                     fsd(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$ODA, use.g.names = FALSE)))
  expect_equal(collapv(wlddev, v, list(fmean, fsd), list(fmode, flast), w = "ODA", keep.col.order = FALSE, give.names = FALSE),
               cbind(g$groups, fsum(get_vars(wlddev, 12), g, use.g.names = FALSE), fmean(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                     fsd(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$ODA, use.g.names = FALSE), flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))
  # with column ordering:
  expect_equal(unname(collapv(wlddev, v, list(fmean, fsd), w = "ODA", wFUN = list(fsum, fmax))),
               unname(cbind(g$groups,
                            fsum(get_vars(wlddev, 12), g, use.g.names = FALSE),
                            fmax(get_vars(wlddev, 12), g, use.g.names = FALSE),
                            fmean(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                            fsd(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                            fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$ODA, use.g.names = FALSE))[order(c(1,5,12,12,4,9:11,4,9:11,2:3,6:8))]))
  expect_equal(unname(collapv(wlddev, v, list(fmean, fsd), list(fmode, flast), w = "ODA", wFUN = list(fsum, fmax))),
               unname(cbind(g$groups,
                            fsum(get_vars(wlddev, 12), g, use.g.names = FALSE),
                            fmax(get_vars(wlddev, 12), g, use.g.names = FALSE),
                            fmean(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                            fsd(get_vars(wlddev, c(4,9:11)), g, wlddev$ODA, use.g.names = FALSE),
                            fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$ODA, use.g.names = FALSE),
                            flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))[order(c(1,5,12,12,4,9:11,4,9:11,2:3,6:8,2:3,6:8))]))
})

test_that("collapv weighted customized aggregation works as intended", {
  v1 <- c("year","PCGDP","LIFEEX","GINI")
  v2 <- c("iso3c","date","region","income", "OECD")
  # Not keeping order ...
  expect_equal(collapv(wlddev, v, w = "ODA", custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, give.names = FALSE),
               add_vars(g$groups,
                        fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  expect_equal(collapv(wlddev, v, w = "ODA", custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, keep.by = FALSE, give.names = FALSE),
               add_vars(fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  expect_equal(collapv(wlddev, v, w = "ODA", custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, keep.w = FALSE, give.names = FALSE),
               add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  expect_equal(collapv(wlddev, v, w = "ODA", custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, keep.by = FALSE, keep.w = FALSE, give.names = FALSE),
               add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE)))

  # keeping order ...
  expect_equal(collapv(wlddev, v, w = "ODA", custom = list(fmean = v1, fmode = v2), give.names = FALSE),
               add_vars(g$groups,
                        fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[names(wlddev)])

  expect_equal(collapv(wlddev, v, w = "ODA", keep.by = FALSE, custom = list(fmean = v1, fmode = v2), give.names = FALSE),
               add_vars(fsum(get_vars(wlddev, "ODA"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[setdiff(names(wlddev), g$group.vars)])

  expect_equal(collapv(wlddev, v, w = "ODA", keep.w = FALSE, custom = list(fmean = v1, fmode = v2), give.names = FALSE),
               add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[setdiff(names(wlddev), "ODA")])

  expect_equal(collapv(wlddev, v, w = "ODA", keep.by = FALSE, keep.w = FALSE, custom = list(fmean = v1, fmode = v2), give.names = FALSE),
               add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI")), g, wlddev$ODA, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$ODA, use.g.names = FALSE))[setdiff(names(wlddev), c(g$group.vars, "ODA"))])


})

test_that("collapv gives informative errors", {
  expect_error(collapv(~cyl, ~cyl)) # nah, need to give error in qDF
  expect_error(collapv(wlddev, ~ country)) # same thing
  expect_error(collapv(wlddev, 13))
  expect_error(collapv(wlddev, 1, w = 13))
  expect_error(collapv(wlddev, 1, w = "bla"))
  expect_error(collapv(wlddev, 1, w = 12, wFUN = sum))
  expect_error(collapv(wlddev, c(1,0)))
  expect_error(collapv(wlddev, c(1,13)))
  expect_warning(collapv(wlddev, 1, bla = 1)) # passes to fmean.data.frame which give the error.
  expect_error(collapv(wlddev, 2, cols = 9:13))
  expect_error(collapv(wlddev, 2, cols = c("PCGDP","bla")))
  expect_error(collapv(wlddev, 2, cols = c("PCGDP","LIFEEX1")))
  expect_error(collapv(wlddev, 2, custom = ~ PCGDP))
  expect_error(collapv(wlddev, 2, custom = list(fmean, fmode)))
  expect_error(collapv(wlddev, 2, custom = list(fmean = 9:13, fmode = 4:6)))
  expect_error(collapv(wlddev, 2, custom = list(fmean = 9:12, 4:6)))
  expect_error(collapv(wlddev, 2, custom = list(fmean = 9:12, fmode2 = 4:6)))
  expect_error(collapv(wlddev, 2, custom = list(fmean = 9:12, fmode = c("GINI","bla"))))
  expect_error(collapv(wlddev, 2, custom = list(fmean = 9:12, fmode = c("GINI","PCGDP2"))))
})


options(warn = 1)
