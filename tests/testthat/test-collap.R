context("collap")

if(!is.null(attributes(identical(FALSE, TRUE)))) stop("OECD label issue")

bsum <- base::sum
bmean <- base::mean

# rm(list = ls())

options(warn = -1)

g <- GRP(wlddev, ~ country + decade)

oa <- function(x) setAttrib(unattrib(x), attributes(x)[c("names", "row.names", "class")])
# Should use above, but sometimes still gives errors
if(Sys.getenv("NCRAN") != "TRUE") oa <- function(x) setNames(unattrib(x), names(x))

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

  expect_equal(oa(collap(wlddev, ~ country + decade, keep.col.order = FALSE)),
               oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))))
  expect_equal(oa(collap(wlddev, ~ country + decade)),
               oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))[order(c(1,5,4,9:13,2:3,6:8))])
  expect_equal(oa(collap(wlddev, ~ country + decade, keep.col.order = FALSE, keep.by = FALSE)),
               oa(cbind(fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))))
  expect_equal(oa(collap(wlddev, ~ country + decade, keep.by = FALSE)),
               oa(cbind(fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))[order(c(4,9:13,2:3,6:8))])

  expect_equal(oa(collap(wlddev, g, keep.by = FALSE)),
               oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))[order(c(1,5,4,9:13,2:3,6:8))])

  expect_equal(oa(collap(wlddev, wlddev[c("country","decade")], keep.by = FALSE)),
               oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))[order(c(1,5,4,9:13,2:3,6:8))])

})

test_that("collap preserves data attributes", {
  expect_identical(lapply(collap(wlddev, ~country), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collap(wlddev, ~country, fmin)), vclasses(wlddev))
  expect_identical(vtypes(collap(wlddev, ~country, fmax)), vtypes(wlddev))
  expect_identical(lapply(collap(wlddev, ~iso3c), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collap(wlddev, ~iso3c, fmin)), vclasses(wlddev))
  expect_identical(vtypes(collap(wlddev, ~iso3c, fmax)), vtypes(wlddev))
  expect_identical(lapply(collap(wlddev, ~date), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collap(wlddev, ~date, fmin)), vclasses(wlddev))
  expect_identical(vtypes(collap(wlddev, ~date, fmax)), vtypes(wlddev))
  expect_identical(lapply(collap(wlddev, ~country + decade), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collap(wlddev, ~country + decade, fmin)), vclasses(wlddev))
  expect_identical(vtypes(collap(wlddev, ~country + decade, fmax)), vtypes(wlddev))
})

# if(Sys.getenv("NCRAN") == "TRUE")
test_that("collap performs as intended in simple uses with base/stats functions", {
  expect_equal(oa(collap(mtcars, mtcars$cyl, bsum, keep.by = FALSE)),
               oa(fsum(mtcars, mtcars$cyl, use.g.names = FALSE)))
  expect_equal(oa(collap(mtcars, ~cyl, mean.default)),
               oa(fmean(mtcars, mtcars$cyl, use.g.names = FALSE)))
  expect_equal(oa(collap(mtcars, ~cyl, bmean)),
               oa(fmean(mtcars, mtcars$cyl, use.g.names = FALSE)))

  expect_equal(oa(collap(mtcars, mtcars[2], bsum, keep.by = FALSE)),
               oa(fsum(mtcars, mtcars$cyl, use.g.names = FALSE)))
  expect_equal(oa(collap(mtcars, ~cyl, bsum, keep.by = FALSE)),
               oa(fsum(mtcars[-2], mtcars$cyl, use.g.names = FALSE)))
  expect_equal(unattrib(collap(iris, ~Species, bsum, keep.by = FALSE)),
               unattrib(fsum(iris[-5], iris$Species, use.g.names = FALSE)))
  expect_equal(oa(collap(airquality, ~Month, bsum, na.rm = TRUE, keep.by = FALSE)),
               oa(fsum(airquality[-5], airquality$Month, use.g.names = FALSE)))

  expect_equal(oa(collap(wlddev, ~ country + decade, bsum, Mode, na.rm = TRUE, keep.col.order = FALSE)),
               oa(cbind(g$groups, BY(get_vars(wlddev, c(4,9:13)), g, bsum, na.rm = TRUE, use.g.names = FALSE),
                     BY(get_vars(wlddev, c(2:3,6:8)), g, Mode, na.rm = TRUE, use.g.names = FALSE))))
  expect_equal(oa(collap(wlddev, ~ country + decade, bsum, Mode, na.rm = TRUE)),
               oa(cbind(g$groups, BY(get_vars(wlddev, c(4,9:13)), g, bsum, na.rm = TRUE, use.g.names = FALSE),
                     BY(get_vars(wlddev, c(2:3,6:8)), g, Mode, na.rm = TRUE, use.g.names = FALSE)))[order(c(1,5,4,9:13,2:3,6:8))])
})

test_that("collap using 2-sided formula or cols performs as intended", {
  expect_equal(oa(collap(mtcars, mpg ~ cyl, keep.by = FALSE)),
               oa(fmean(mtcars["mpg"], mtcars$cyl, use.g.names = FALSE)))
  expect_equal(oa(collap(mtcars, mpg ~ cyl, keep.by = FALSE, cols = 300:1000)),
               oa(fmean(mtcars["mpg"], mtcars$cyl, use.g.names = FALSE))) # cols is ignored, as should be
  expect_equal(oa(collap(mtcars, ~ cyl, keep.by = FALSE, cols = 1)),
               oa(fmean(mtcars["mpg"], mtcars$cyl, use.g.names = FALSE)))
  expect_equal(oa(collap(mtcars, wt + mpg ~ cyl + vs + am, keep.by = FALSE)),
               oa(fmean(mtcars[c("mpg","wt")], mtcars[c("cyl","vs","am")], use.g.names = FALSE)))
  expect_equal(oa(collap(mtcars, ~ cyl + vs + am, keep.by = FALSE, cols = c(6,1))),
               oa(fmean(mtcars[c("mpg","wt")], mtcars[c("cyl","vs","am")], use.g.names = FALSE)))
  expect_equal(oa(collap(iris, Sepal.Length + Sepal.Width ~ Species, keep.by = FALSE)),
               oa(fmean(iris[1:2], iris$Species, use.g.names = FALSE)))
  expect_equal(oa(collap(airquality, ~ Month, keep.by = FALSE)),
               oa(fmean(airquality[-5], airquality$Month, use.g.names = FALSE)))
  expect_equal(oa(collap(airquality, ~ Month, keep.by = FALSE, cols = 1:3)),
               oa(fmean(airquality[1:3], airquality$Month, use.g.names = FALSE)))

  expect_equal(oa(collap(wlddev, ~ country + decade, cols = 9:13)),
               oa(collap(wlddev, PCGDP + LIFEEX + GINI + ODA + POP ~ country + decade)))
  expect_equal(oa(collap(wlddev, ~ country + decade, cols = 9:13)),
               oa(collap(wlddev, ~ country + decade, cols = 9:13, keep.col.order = FALSE)))
  expect_equal(oa(collap(wlddev, ~ country + decade, cols = c(2:3,6:8))),
               oa(collap(wlddev, iso3c + date + region + income + OECD ~ country + decade)))
  expect_false(identical(oa(collap(wlddev, ~ country + decade, cols = c(2:3,6:8))),
                         oa(collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.col.order = FALSE))))

  expect_equal(oa(collap(wlddev, ~ country + decade, cols = 9:12, keep.by = FALSE)),
               oa(collap(wlddev, PCGDP + LIFEEX + GINI + ODA ~ country + decade, keep.by = FALSE)))
  expect_equal(oa(collap(wlddev, ~ country + decade, cols = 9:13, keep.by = FALSE)),
               oa(collap(wlddev, ~ country + decade, cols = 9:13, keep.col.order = FALSE, keep.by = FALSE)))
  expect_equal(oa(collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.by = FALSE)),
               oa(collap(wlddev, iso3c + date + region + income + OECD ~ country + decade, keep.by = FALSE)))
  expect_equal(oa(collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.by = FALSE)),
               oa(collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.col.order = FALSE, keep.by = FALSE)))

  expect_equal(oa(collap(wlddev, g, cols = 9:12, keep.by = FALSE)),
               oa(collap(wlddev, PCGDP + LIFEEX + GINI + ODA ~ country + decade, keep.by = FALSE)))
  expect_equal(oa(collap(wlddev, g, cols = 9:13, keep.by = FALSE)),
               oa(collap(wlddev, ~ country + decade, cols = 9:13, keep.col.order = FALSE, keep.by = FALSE)))
  expect_equal(oa(collap(wlddev, g, cols = c(2:3,6:8), keep.by = FALSE)),
               oa(collap(wlddev, iso3c + date + region + income + OECD ~ country + decade, keep.by = FALSE)))
  expect_equal(oa(collap(wlddev, g, cols = c(2:3,6:8), keep.by = FALSE)),
               oa(collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.col.order = FALSE, keep.by = FALSE)))

  expect_equal(oa(collap(wlddev, wlddev[c("country","decade")], cols = 9:12, keep.by = FALSE)),
               oa(collap(wlddev, PCGDP + LIFEEX + GINI + ODA ~ country + decade, keep.by = FALSE)))
  expect_equal(oa(collap(wlddev, wlddev[c("country","decade")], cols = 9:13, keep.by = FALSE)),
               oa(collap(wlddev, ~ country + decade, cols = 9:13, keep.col.order = FALSE, keep.by = FALSE)))
  expect_equal(oa(collap(wlddev, wlddev[c("country","decade")], cols = c(2:3,6:8), keep.by = FALSE)),
               oa(collap(wlddev, iso3c + date + region + income + OECD ~ country + decade, keep.by = FALSE)))
  expect_equal(oa(collap(wlddev, wlddev[c("country","decade")], cols = c(2:3,6:8), keep.by = FALSE)),
               oa(collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.col.order = FALSE, keep.by = FALSE)))

})

test_that("collap multi-function aggreagtion performs as intended", {
  expect_equal(oa(collap(wlddev, ~ country + decade, list(fmean, fmedian), keep.col.order = FALSE, give.names = FALSE)),
               oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))))
  if(Sys.getenv("NCRAN") == "TRUE")
  expect_equal(oa(collap(wlddev, ~ country + decade, list(fmean, fmedian), list(fmode, flast), keep.col.order = FALSE, give.names = FALSE)),
               oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE), flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))))
  # with column ordering:
  expect_equal(unname(oa(collap(wlddev, ~ country + decade, list(fmean, fmedian)))),
               unname(oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))[order(c(1,5,4,9:13,4,9:13,2:3,6:8))])))
  expect_equal(unname(oa(collap(wlddev, ~ country + decade, list(fmean, fmedian), list(fmode, flast)))),
               unname(oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                      fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE), flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))[order(c(1,5,4,9:13,4,9:13,2:3,6:8,2:3,6:8))]))
})

test_that("collap custom aggreagtion performs as intended", {
  expect_equal(unname(oa(collap(wlddev, ~ country + decade,
                      custom = list(fmean = 9:13, fsd = 9:10, fmode = 7:8), keep.col.order = FALSE))),
               unname(oa(cbind(g$groups, fmean(wlddev[9:13], g, use.g.names = FALSE),
                                      fsd(wlddev[9:10], g, use.g.names = FALSE),
                                      fmode(wlddev[7:8], g, use.g.names = FALSE)))))
  expect_equal(unname(oa(collap(wlddev, ~ country + decade,
                             custom = list(fmean = 9:13, fsd = 9:10, fmode = 7:8)))),
               unname(oa(cbind(g$groups, fmean(wlddev[9:13], g, use.g.names = FALSE),
                            fsd(wlddev[9:10], g, use.g.names = FALSE),
                            fmode(wlddev[7:8], g, use.g.names = FALSE)))[order(c(1,5,9:13,9:10,7:8))]))

  expect_equal(oa(collap(wlddev, ~ country + decade,
               custom = list(fmean = 9:13, fsd = 9:10, fmode = 7:8))),
               oa(collap(wlddev, ~ country + decade,
               custom = list(fmean = 9:13, fsd = c("PCGDP","LIFEEX"), fmode = 7:8))))
  expect_equal(oa(collap(wlddev, ~ country + decade,
               custom = list(fmean = "PCGDP", fsd = c("LIFEEX","GINI"), flast = "date"))),
               oa(collap(wlddev, ~ country + decade,
               custom = list(fmean = "PCGDP", fsd = 10:11, flast = "date"))))

  expect_equal(oa(collap(wlddev, g,
                      custom = list(fmean = 9:13, fsd = 9:10, fmode = 7:8))),
               oa(collap(wlddev, ~ country + decade,
                      custom = list(fmean = 9:13, fsd = c("PCGDP","LIFEEX"), fmode = 7:8))))
  expect_equal(oa(collap(wlddev, g,
                      custom = list(fmean = "PCGDP", fsd = c("LIFEEX","GINI"), flast = "date"))),
               oa(collap(wlddev, g,
                      custom = list(fmean = "PCGDP", fsd = 10:11, flast = "date"))))
  expect_equal(names(collap(wlddev, g,
                      custom = list(fmean = c(GDP = "PCGDP"), fsd = c("LIFEEX", GN = "GINI"), flast = "date"),
                      keep.by = FALSE, keep.col.order = FALSE)),
               .c(GDP, LIFEEX, GN, date))

})

test_that("collap weighted aggregations work as intended", {
  # Not keeping order ...
  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, keep.col.order = FALSE)),
               oa(add_vars(g$groups,
                  fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                  fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                  fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, keep.col.order = FALSE, keep.by = FALSE)),
               oa(add_vars(fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, keep.col.order = FALSE, keep.w = FALSE)),
               oa(add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, keep.col.order = FALSE, keep.by = FALSE, keep.w = FALSE)),
               oa(add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  # keeping order ...
  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP)),
               oa(add_vars(g$groups,
                        fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[names(wlddev)])

  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, keep.by = FALSE)),
               oa(add_vars(fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[setdiff(names(wlddev), g$group.vars)])

  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, keep.w = FALSE)),
               oa(add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[setdiff(names(wlddev), "POP")])

  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, keep.by = FALSE, keep.w = FALSE)),
               oa(add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[setdiff(names(wlddev), c(g$group.vars, "POP"))])


})

if(Sys.getenv("NCRAN") == "TRUE")
test_that("collap multi-function aggreagtion with weights performs as intended", {

  expect_equal(oa(collap(wlddev, ~ country + decade, list(fmean, fsd), w = ~ POP, keep.col.order = FALSE, give.names = FALSE)),
               oa(cbind(g$groups, fsum(get_vars(wlddev, 13), g, use.g.names = FALSE), fmean(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                     fsd(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(oa(collap(wlddev, ~ country + decade, list(fmean, fsd), list(fmode, flast), w = ~ POP, keep.col.order = FALSE, give.names = FALSE)),
               oa(cbind(g$groups, fsum(get_vars(wlddev, 13), g, use.g.names = FALSE), fmean(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                     fsd(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$POP, use.g.names = FALSE), flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))))

  # with column ordering:
  expect_equal(unname(oa(collap(wlddev, ~ country + decade, list(fmean, fsd), w = ~ POP, wFUN = list(fsum, fmax)))),
               unname(oa(cbind(g$groups,
                            fsum(get_vars(wlddev, 13), g, use.g.names = FALSE),
                            fmax(get_vars(wlddev, 13), g, use.g.names = FALSE),
                            fmean(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                            fsd(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                            fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$POP, use.g.names = FALSE)))[order(c(1,5,13,13,4,9:12,4,9:12,2:3,6:8))]))

  expect_equal(unname(oa(collap(wlddev, ~ country + decade, list(fmean, fsd), list(fmode, flast), w = ~ POP, wFUN = list(fsum, fmax)))),
               unname(oa(cbind(g$groups,
                            fsum(get_vars(wlddev, 13), g, use.g.names = FALSE),
                            fmax(get_vars(wlddev, 13), g, use.g.names = FALSE),
                            fmean(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                            fsd(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                            fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$POP, use.g.names = FALSE),
                            flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))[order(c(1,5,13,13,4,9:12,4,9:12,2:3,6:8,2:3,6:8))]))
})

v1 <- c("year","PCGDP","LIFEEX","GINI","ODA")
v2 <- c("iso3c","date","region","income", "OECD")
test_that("collap weighted customized aggregation works as intended", {
  # Not keeping order ...
  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, give.names = FALSE)),
               oa(add_vars(g$groups,
                        fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(unattrib(collap(wlddev, ~ country + decade, w = ~ POP, custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, keep.by = FALSE, give.names = FALSE)),
               unattrib(add_vars(fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, keep.w = FALSE, give.names = FALSE)),
               oa(add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, keep.by = FALSE, keep.w = FALSE, give.names = FALSE)),
               oa(add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  # keeping order ...
  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, custom = list(fmean = v1, fmode = v2), give.names = FALSE)),
               oa(add_vars(g$groups,
                        fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[names(wlddev)])

  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, keep.by = FALSE, custom = list(fmean = v1, fmode = v2), give.names = FALSE)),
               oa(add_vars(fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[setdiff(names(wlddev), g$group.vars)])

  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, keep.w = FALSE, custom = list(fmean = v1, fmode = v2), give.names = FALSE)),
               oa(add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[setdiff(names(wlddev), "POP")])

  expect_equal(oa(collap(wlddev, ~ country + decade, w = ~ POP, keep.by = FALSE, keep.w = FALSE, custom = list(fmean = v1, fmode = v2), give.names = FALSE)),
               oa(add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[setdiff(names(wlddev), c(g$group.vars, "POP"))])


})

test_that("collap gives informative errors", {
  expect_error(collap(~cyl, ~cyl)) # nah, need to give error in qDF
  expect_error(collap(wlddev, 1:3)) # only gives error in fmean.. a bit late..
  expect_error(collap(wlddev, "country")) # same thing
  expect_error(collap(wlddev, ~ country1))
  expect_error(collap(wlddev, ~ country, w = ~bla))
  expect_error(collap(wlddev, ~ country, w = ~POP, wFUN = bsum))
  expect_error(collap(wlddev, ~ country + year + bla))
  expect_error(collap(wlddev, bla ~ country))
  expect_warning(collap(wlddev, ~ country, bla = 1)) # passes to fmean.data.frame which give the error.
  # expect_error(collap(wlddev, ~ country, bsum, cols = 9:13, bla = 1)) # This is an issue, bsum(1:3, bla = 1) does not give an error
  expect_error(collap(wlddev, mtcars$cyl)) # again fmean error..
  expect_error(collap(wlddev, ~iso3c, cols = 9:14))
  # expect_error(collap(wlddev, ~iso3c, cols = 0:1)) # no error..
  expect_error(collap(wlddev, ~iso3c, cols = c("PCGDP","bla")))
  expect_error(collap(wlddev, ~iso3c, cols = c("PCGDP","LIFEEX1")))
  expect_error(collap(wlddev, ~iso3c, custom = ~ PCGDP))
  expect_error(collap(wlddev, ~iso3c, custom = list(fmean, fmode)))
  expect_error(collap(wlddev, ~iso3c, custom = list(fmean = 9:14, fmode = 4:6)))
  expect_error(collap(wlddev, ~iso3c, custom = list(fmean = 9:13, 4:6)))
  expect_error(collap(wlddev, ~iso3c, custom = list(fmean = 9:13, fmode2 = 4:6)))
  expect_error(collap(wlddev, ~iso3c, custom = list(fmean = 9:13, fmode = c("GINI","bla"))))
  expect_error(collap(wlddev, ~iso3c, custom = list(fmean = 9:13, fmode = c("GINI","PCGDP2"))))
})

# Note: one more thing to test is performance with vector-valued functions...


# Testing collapv
v <- c(1, 5)

test_that("collapv performs as intended in simple uses", {
  expect_equal(oa(collapv(mtcars, 2)),
               oa(fmean(mtcars, mtcars$cyl, use.g.names = FALSE)))
  expect_equal(oa(collapv(mtcars, 2, keep.by = FALSE)),
               oa(fmean(mtcars[-2], mtcars$cyl, use.g.names = FALSE)))
  expect_equal(oa(collapv(iris, "Species", keep.by = FALSE)),
               oa(fmean(iris[-5], iris$Species, use.g.names = FALSE)))
  expect_equal(oa(collapv(airquality, "Month", keep.by = FALSE)),
               oa(fmean(airquality[-5], airquality$Month, use.g.names = FALSE)))

  expect_equal(oa(collapv(wlddev, v, keep.col.order = FALSE)),
               oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))))

  expect_equal(oa(collapv(wlddev, v)),
               oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))[order(c(1,5,4,9:13,2:3,6:8))])
  expect_equal(oa(collapv(wlddev, v, keep.col.order = FALSE, keep.by = FALSE)),
               oa(cbind(fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))))
  expect_equal(oa(collapv(wlddev, v, keep.by = FALSE)),
               oa(cbind(fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))[order(c(4,9:13,2:3,6:8))])

  expect_equal(names(collapv(wlddev, v,
                            custom = list(fmean = c(GDP = "PCGDP"), fsd = c("LIFEEX", GN = "GINI"), flast = "date"),
                            keep.by = FALSE, keep.col.order = FALSE)),
               .c(GDP, LIFEEX, GN, date))

})

test_that("collapv preserves data attributes", {
  expect_identical(lapply(collapv(wlddev, 1), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collapv(wlddev, 1, fmax)), vclasses(wlddev))
  expect_identical(vtypes(collapv(wlddev, 1, fmin)), vtypes(wlddev))
  expect_identical(lapply(collapv(wlddev, "iso3c"), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collapv(wlddev, "iso3c", fmax)), vclasses(wlddev))
  expect_identical(vtypes(collapv(wlddev, "iso3c", fmin)), vtypes(wlddev))
  expect_identical(lapply(collapv(wlddev, "date"), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collapv(wlddev, "date", ffirst)), vclasses(wlddev))
  expect_identical(vtypes(collapv(wlddev, "date", flast)), vtypes(wlddev))
  expect_identical(lapply(collapv(wlddev, v), attributes), lapply(wlddev, attributes))
  expect_identical(vclasses(collapv(wlddev, v, flast)), vclasses(wlddev))
  expect_identical(vtypes(collapv(wlddev, v, ffirst)), vtypes(wlddev))
})

# if(Sys.getenv("NCRAN") == "TRUE")
test_that("collapv performs as intended in simple uses with base/stats functions", {
  expect_equal(oa(collapv(mtcars, "cyl", mean.default)),
               oa(fmean(mtcars, mtcars$cyl, use.g.names = FALSE)))
  expect_equal(oa(collapv(mtcars, "cyl", bmean)),
               oa(fmean(mtcars, mtcars$cyl, use.g.names = FALSE)))

  expect_equal(oa(collapv(mtcars, 2, bsum, keep.by = FALSE)),
               oa(fsum(mtcars[-2], mtcars$cyl, use.g.names = FALSE)))
  expect_equal(oa(collapv(iris, 5, bsum, keep.by = FALSE)),
               oa(fsum(iris[-5], iris$Species, use.g.names = FALSE)))
  expect_equal(oa(collapv(airquality, "Month", bsum, na.rm = TRUE, keep.by = FALSE)),
               oa(fsum(airquality[-5], airquality$Month, use.g.names = FALSE)))

  expect_equal(oa(collapv(wlddev, v, bsum, Mode, na.rm = TRUE, keep.col.order = FALSE)),
               oa(cbind(g$groups, BY(get_vars(wlddev, c(4,9:13)), g, bsum, na.rm = TRUE, use.g.names = FALSE),
                     BY(get_vars(wlddev, c(2:3,6:8)), g, Mode, na.rm = TRUE, use.g.names = FALSE))))
  expect_equal(oa(collapv(wlddev, v, bsum, Mode, na.rm = TRUE)),
               oa(cbind(g$groups, BY(get_vars(wlddev, c(4,9:13)), g, bsum, na.rm = TRUE, use.g.names = FALSE),
                     BY(get_vars(wlddev, c(2:3,6:8)), g, Mode, na.rm = TRUE, use.g.names = FALSE)))[order(c(1,5,4,9:13,2:3,6:8))])
})

test_that("collapv using cols performs as intended", {
  expect_equal(oa(collapv(mtcars, 2, keep.by = FALSE, cols = 1)),
               oa(fmean(mtcars["mpg"], mtcars$cyl, use.g.names = FALSE)))
  expect_equal(oa(collapv(mtcars, c("cyl", "vs", "am"), keep.by = FALSE, cols = c(6,1))),
               oa(fmean(mtcars[c("mpg","wt")], mtcars[c("cyl","vs","am")], use.g.names = FALSE)))
  expect_equal(oa(collapv(airquality, "Month", keep.by = FALSE)),
               oa(fmean(airquality[-5], airquality$Month, use.g.names = FALSE)))
  expect_equal(oa(collapv(airquality, "Month", keep.by = FALSE, cols = 1:3)),
               oa(fmean(airquality[1:3], airquality$Month, use.g.names = FALSE)))

  expect_equal(oa(collapv(wlddev, v, cols = 9:12)),
               oa(collap(wlddev, PCGDP + LIFEEX + GINI + ODA ~ country + decade)))
  expect_equal(oa(collapv(wlddev, v, cols = 9:13)),
               oa(collap(wlddev, ~ country + decade, cols = 9:13, keep.col.order = FALSE)))
  expect_equal(oa(collapv(wlddev, v, cols = c(2:3,6:8))),
               oa(collap(wlddev, iso3c + date + region + income + OECD ~ country + decade)))
  expect_false(identical(collapv(wlddev, v, cols = c(2:3,6:8)),
                         collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.col.order = FALSE)))

  expect_equal(oa(collapv(wlddev, v, cols = 9:12, keep.by = FALSE)),
               oa(collap(wlddev, PCGDP + LIFEEX + GINI + ODA ~ country + decade, keep.by = FALSE)))
  expect_equal(oa(collapv(wlddev, v, cols = 9:13, keep.by = FALSE)),
               oa(collap(wlddev, ~ country + decade, cols = 9:13, keep.col.order = FALSE, keep.by = FALSE)))
  expect_equal(oa(collapv(wlddev, v, cols = c(2:3,6:8), keep.by = FALSE)),
               oa(collap(wlddev, iso3c + date + region + income + OECD ~ country + decade, keep.by = FALSE)))
  expect_equal(oa(collapv(wlddev, v, cols = c(2:3,6:8), keep.by = FALSE)),
               oa(collap(wlddev, ~ country + decade, cols = c(2:3,6:8), keep.col.order = FALSE, keep.by = FALSE)))

})

test_that("collapv multi-function aggreagtion performs as intended", {

  expect_equal(oa(collapv(wlddev, v, list(fmean, fmedian), keep.col.order = FALSE, give.names = FALSE)),
               oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))))
  expect_equal(oa(collapv(wlddev, v, list(fmean, fmedian), list(fmode, flast), keep.col.order = FALSE, give.names = FALSE)),
               oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE), flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))))
  # with column ordering:
  expect_equal(unname(oa(collapv(wlddev, v, list(fmean, fmedian)))),
               unname(oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                            fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))[order(c(1,5,4,9:13,4,9:13,2:3,6:8))]))
  expect_equal(unname(oa(collapv(wlddev, v, list(fmean, fmedian), list(fmode, flast)))),
               unname(oa(cbind(g$groups, fmean(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE), fmedian(get_vars(wlddev, c(4,9:13)), g, use.g.names = FALSE),
                            fmode(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE), flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))[order(c(1,5,4,9:13,4,9:13,2:3,6:8,2:3,6:8))]))
})

test_that("collapv custom aggreagtion performs as intended", {
  expect_equal(unname(oa(collapv(wlddev, v,
                             custom = list(fmean = 9:13, fsd = 9:10, fmode = 7:8), keep.col.order = FALSE))),
               unname(oa(cbind(g$groups, fmean(wlddev[9:13], g, use.g.names = FALSE),
                            fsd(wlddev[9:10], g, use.g.names = FALSE),
                            fmode(wlddev[7:8], g, use.g.names = FALSE)))))
  expect_equal(unname(oa(collapv(wlddev, v,
                             custom = list(fmean = 9:13, fsd = 9:10, fmode = 7:8)))),
               unname(oa(cbind(g$groups, fmean(wlddev[9:13], g, use.g.names = FALSE),
                            fsd(wlddev[9:10], g, use.g.names = FALSE),
                            fmode(wlddev[7:8], g, use.g.names = FALSE)))[order(c(1,5,9:13,9:10,7:8))]))

  expect_equal(oa(collapv(wlddev, v,
                      custom = list(fmean = 9:13, fsd = 9:10, fmode = 7:8))),
               oa(collapv(wlddev, v,
                      custom = list(fmean = 9:13, fsd = c("PCGDP","LIFEEX"), fmode = 7:8))))
  expect_equal(oa(collapv(wlddev, v,
                      custom = list(fmean = "PCGDP", fsd = c("LIFEEX","GINI"), flast = "date"))),
               oa(collapv(wlddev, v,
                      custom = list(fmean = "PCGDP", fsd = 10:11, flast = "date"))))

  expect_equal(oa(collapv(wlddev, v,
                      custom = list(fmean = 9:13, fsd = 9:10, fmode = 7:8))),
               oa(collapv(wlddev, v,
                      custom = list(fmean = 9:13, fsd = c("PCGDP","LIFEEX"), fmode = 7:8))))
  expect_equal(oa(collapv(wlddev, v,
                      custom = list(fmean = "PCGDP", fsd = c("LIFEEX","GINI"), flast = "date"))),
               oa(collapv(wlddev, v,
                      custom = list(fmean = "PCGDP", fsd = 10:11, flast = "date"))))

})

test_that("collapv weighted aggregations work as intended", {
  # Not keeping order ...
  expect_equal(oa(collapv(wlddev, v, w = "POP", keep.col.order = FALSE)),
               oa(add_vars(g$groups,
                        fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(oa(collapv(wlddev, v, w = "POP", keep.col.order = FALSE, keep.by = FALSE)),
               oa(add_vars(fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(oa(collapv(wlddev, v, w = "POP", keep.col.order = FALSE, keep.w = FALSE)),
               oa(add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(oa(collapv(wlddev, v, w = "POP", keep.col.order = FALSE, keep.by = FALSE, keep.w = FALSE)),
               oa(add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  # keeping order ...
  expect_equal(oa(collapv(wlddev, v, w = "POP")),
               oa(add_vars(g$groups,
                        fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[names(wlddev)])

  expect_equal(unattrib(collapv(wlddev, v, w = "POP", keep.by = FALSE)),
               unattrib(add_vars(fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))[setdiff(names(wlddev), g$group.vars)]))

  expect_equal(oa(collapv(wlddev, v, w = "POP", keep.w = FALSE)),
               oa(add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[setdiff(names(wlddev), "POP")])

  expect_equal(oa(collapv(wlddev, v, w = "POP", keep.by = FALSE, keep.w = FALSE)),
               oa(add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[setdiff(names(wlddev), c(g$group.vars, "POP"))])


})

if(Sys.getenv("NCRAN") == "TRUE")
test_that("collapv multi-function aggreagtion with weights performs as intended", {

  expect_equal(oa(collapv(wlddev, v, list(fmean, fsd), w = "POP", keep.col.order = FALSE, give.names = FALSE)),
               oa(cbind(g$groups, fsum(get_vars(wlddev, 13), g, use.g.names = FALSE), fmean(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                     fsd(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(oa(collapv(wlddev, v, list(fmean, fsd), list(fmode, flast), w = "POP", keep.col.order = FALSE, give.names = FALSE)),
               oa(cbind(g$groups, fsum(get_vars(wlddev, 13), g, use.g.names = FALSE), fmean(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                     fsd(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                     fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$POP, use.g.names = FALSE), flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE))))
  # with column ordering:
  expect_equal(unname(oa(collapv(wlddev, v, list(fmean, fsd), w = "POP", wFUN = list(fsum, fmax)))),
               unname(oa(cbind(g$groups,
                            fsum(get_vars(wlddev, 13), g, use.g.names = FALSE),
                            fmax(get_vars(wlddev, 13), g, use.g.names = FALSE),
                            fmean(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                            fsd(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                            fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$POP, use.g.names = FALSE)))[order(c(1,5,13,13,4,9:12,4,9:12,2:3,6:8))]))

  expect_equal(unattrib(collapv(wlddev, v, list(fmean, fsd), list(fmode, flast), w = "POP", wFUN = list(fsum, fmax))),
               unattrib(cbind(g$groups,
                            fsum(get_vars(wlddev, 13), g, use.g.names = FALSE),
                            fmax(get_vars(wlddev, 13), g, use.g.names = FALSE),
                            fmean(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                            fsd(get_vars(wlddev, c(4,9:12)), g, wlddev$POP, use.g.names = FALSE),
                            fmode(get_vars(wlddev, c(2:3,6:8)), g, wlddev$POP, use.g.names = FALSE),
                            flast(get_vars(wlddev, c(2:3,6:8)), g, use.g.names = FALSE)))[order(c(1,5,13,13,4,9:12,4,9:12,2:3,6:8,2:3,6:8))])
})

v1 <- c("year","PCGDP","LIFEEX","GINI","ODA")
v2 <- c("iso3c","date","region","income", "OECD")
test_that("collapv weighted customized aggregation works as intended", {
  # Not keeping order ...
  expect_equal(oa(collapv(wlddev, v, w = "POP", custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, give.names = FALSE)),
               oa(add_vars(g$groups,
                        fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(unattrib(collapv(wlddev, v, w = "POP", custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, keep.by = FALSE, give.names = FALSE)),
               unattrib(add_vars(fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(oa(collapv(wlddev, v, w = "POP", custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, keep.w = FALSE, give.names = FALSE)),
               oa(add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  expect_equal(oa(collapv(wlddev, v, w = "POP", custom = list(fmean = v1, fmode = v2), keep.col.order = FALSE, keep.by = FALSE, keep.w = FALSE, give.names = FALSE)),
               oa(add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE))))

  # keeping order ...
  expect_equal(oa(collapv(wlddev, v, w = "POP", custom = list(fmean = v1, fmode = v2), give.names = FALSE)),
               oa(add_vars(g$groups,
                        fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[names(wlddev)])

  expect_equal(oa(collapv(wlddev, v, w = "POP", keep.by = FALSE, custom = list(fmean = v1, fmode = v2), give.names = FALSE)),
               oa(add_vars(fsum(get_vars(wlddev, "POP"), g, use.g.names = FALSE),
                        fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[setdiff(names(wlddev), g$group.vars)])

  expect_equal(oa(collapv(wlddev, v, w = "POP", keep.w = FALSE, custom = list(fmean = v1, fmode = v2), give.names = FALSE)),
               oa(add_vars(g$groups, fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[setdiff(names(wlddev), "POP")])

  expect_equal(oa(collapv(wlddev, v, w = "POP", keep.by = FALSE, keep.w = FALSE, custom = list(fmean = v1, fmode = v2), give.names = FALSE)),
               oa(add_vars(fmean(get_vars(wlddev, c("year","PCGDP","LIFEEX","GINI","ODA")), g, wlddev$POP, use.g.names = FALSE),
                        fmode(get_vars(wlddev, c("iso3c","date","region","income", "OECD")), g, wlddev$POP, use.g.names = FALSE)))[setdiff(names(wlddev), c(g$group.vars, "POP"))])


})

test_that("collapv gives informative errors", {
  expect_error(collapv(~cyl, ~cyl)) # nah, need to give error in qDF
  expect_error(collapv(wlddev, ~ country)) # same thing
  expect_error(collapv(wlddev, 14))
  expect_error(collapv(wlddev, 1, w = 14))
  expect_error(collapv(wlddev, 1, w = "bla"))
  expect_error(collapv(wlddev, 1, w = 13, wFUN = bsum))
  expect_error(collapv(wlddev, c(1,0)))
  expect_error(collapv(wlddev, c(1,14)))
  expect_warning(collapv(wlddev, 1, bla = 1)) # passes to fmean.data.frame which give the error.
  expect_error(collapv(wlddev, 2, cols = 9:14))
  expect_error(collapv(wlddev, 2, cols = c("PCGDP","bla")))
  expect_error(collapv(wlddev, 2, cols = c("PCGDP","LIFEEX1")))
  expect_error(collapv(wlddev, 2, custom = ~ PCGDP))
  expect_error(collapv(wlddev, 2, custom = list(fmean, fmode)))
  expect_error(collapv(wlddev, 2, custom = list(fmean = 9:14, fmode = 4:6)))
  expect_error(collapv(wlddev, 2, custom = list(fmean = 9:14, 4:6)))
  expect_error(collapv(wlddev, 2, custom = list(fmean = 9:14, fmode2 = 4:6)))
  expect_error(collapv(wlddev, 2, custom = list(fmean = 9:13, fmode = c("GINI","bla"))))
  expect_error(collapv(wlddev, 2, custom = list(fmean = 9:13, fmode = c("GINI","PCGDP2"))))
})


options(warn = 1)
