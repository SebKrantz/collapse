context("pivot")

library(data.table)
mtcDT <- qDT(mtcars)
mtcnaDT <- qDT(na_insert(mtcars))
irisDT <- qDT(iris)
wldDT <- qDT(wlddev)
GGDCDT <- qDT(GGDC10S)

rmnic <- function(x) {
  if(!length(fci <- fact_vars(x, "indices"))) return(x)
  for (i in fci) oldClass(x[[i]]) <- setdiff(oldClass(x[[i]]), "na.included")
  x
}

test_that("long pivots work properly", {
  # No id's
  expect_identical(rmnic(pivot(mtcDT)), melt(mtcDT, measure.vars = seq_along(mtcDT)))
  expect_identical(rmnic(pivot(mtcDT, values = 3:11)), melt(mtcDT, measure.vars = 3:11))
  expect_identical(rmnic(pivot(mtcnaDT, na.rm = TRUE)), melt(mtcnaDT, measure.vars = seq_along(mtcnaDT), na.rm = TRUE))
  expect_identical(rmnic(pivot(mtcnaDT, values = 3:11, na.rm = TRUE)), melt(mtcnaDT, measure.vars = 3:11, na.rm = TRUE))

  expect_identical(names(pivot(gv(wlddev, 9:10), labels = TRUE)), c("variable", "label", "value"))
  expect_identical(names(pivot(gv(wlddev, 9:10), labels = "bla")), c("variable", "bla", "value"))
  expect_identical(names(pivot(gv(wlddev, 9:10), labels = TRUE, na.rm = TRUE)), c("variable", "label", "value"))
  expect_identical(names(pivot(gv(wlddev, 9:10), labels = "bla", na.rm = TRUE)), c("variable", "bla", "value"))
  expect_warning(pivot(mtcnaDT, check.dups = TRUE))

  # with ids
  expect_identical(rmnic(pivot(irisDT, "Species")), melt(irisDT, "Species"))
  expect_identical(rmnic(setLabels(pivot(wldDT, 1:8), NULL)), setLabels(melt(wldDT, 1:8), NULL))
  expect_identical(rmnic(setLabels(pivot(wldDT, 1:8, na.rm = TRUE), NULL)), setLabels(melt(wldDT, 1:8, na.rm = TRUE), NULL))
  expect_warning(pivot(irisDT, "Species", check.dups = TRUE))
  # with labels
  expect_identical(names(pivot(wldDT, c("iso3c", "year"), values = 9:10, labels = TRUE)), c("iso3c", "year", "variable", "label", "value"))
  expect_identical(names(pivot(wldDT, c("iso3c", "year"), values = 9:10, names = list("var", "val"), labels = "lab")), c("iso3c", "year", "var", "lab", "val"))
  expect_identical(names(pivot(wldDT, c("iso3c", "year"), values = 9:10, names = list(value = "val"), labels = "lab")), c("iso3c", "year", "variable", "lab", "val"))
  expect_identical(names(pivot(wldDT, c("iso3c", "year"), values = 9:10, names = list(variable = "var"), labels = "lab")), c("iso3c", "year", "var", "lab", "value"))

})


# test_that("wide pivots work properly", {
#
#   pivot(wlddev, "iso3c", "PCGDP", "year", how = "wider")
#   pivot(wlddev, "iso3c", "PCGDP", "year", how = "wider", check.dups = TRUE, na.rm = TRUE, sort = c("ids", "names"))
#
#   pivot(wlddev, "iso3c", "PCGDP", "year", "decade", how = "wider", check.dups = TRUE, na.rm = TRUE, sort = c("ids", "names"))
#   pivot(wlddev, "iso3c", c("PCGDP", "LIFEEX"), "year", "decade", how = "wider", check.dups = TRUE, na.rm = TRUE, sort = c("ids", "names"))
#
#   pivot(wlddev, "iso3c", c("PCGDP", "LIFEEX"), "year", "decade", how = "wider", check.dups = TRUE, na.rm = TRUE, sort = c("ids", "names"), transpose = c("cols", "names"))
#
# })
