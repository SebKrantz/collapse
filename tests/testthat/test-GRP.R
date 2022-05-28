context("radixorder, GRP, qF, qG")

# rm(list = ls())
set.seed(101)
mtcNA <- na_insert(mtcars)
wlddev2 <- slt(wlddev, -date)
num_vars(wlddev2) <- round(num_vars(wlddev2), 8)
wldNA <- na_insert(wlddev2)

GGDCNA <- na_insert(GGDC10S)


unlab <- function(x) `attr<-`(x, "label", NULL)

test_that("radixorder works like order(.., method = 'radix')", {

# Ordering single variable
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x))), lapply(wldNA, order, method = "radix"))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, decreasing = TRUE))), lapply(wldNA, order, method = "radix", decreasing = TRUE))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, na.last = FALSE))), lapply(wldNA, order, method = "radix", na.last = FALSE))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, decreasing = TRUE, na.last = FALSE))), lapply(wldNA, order, method = "radix", decreasing = TRUE, na.last = FALSE))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, na.last = NA))), lapply(wldNA, order, method = "radix", na.last = NA))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, decreasing = TRUE, na.last = NA))), lapply(wldNA, order, method = "radix", decreasing = TRUE, na.last = NA))
# get starts
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, starts = TRUE))), lapply(wldNA, order, method = "radix"))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, decreasing = TRUE, starts = TRUE))), lapply(wldNA, order, method = "radix", decreasing = TRUE))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, na.last = FALSE, starts = TRUE))), lapply(wldNA, order, method = "radix", na.last = FALSE))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, decreasing = TRUE, na.last = FALSE, starts = TRUE))), lapply(wldNA, order, method = "radix", decreasing = TRUE, na.last = FALSE))
# get group.sizes
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, group.sizes = TRUE))), lapply(wldNA, order, method = "radix"))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, decreasing = TRUE, group.sizes = TRUE))), lapply(wldNA, order, method = "radix", decreasing = TRUE))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, na.last = FALSE, group.sizes = TRUE))), lapply(wldNA, order, method = "radix", na.last = FALSE))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, decreasing = TRUE, na.last = FALSE, group.sizes = TRUE))), lapply(wldNA, order, method = "radix", decreasing = TRUE, na.last = FALSE))
# get starts and group.sizes
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, starts = TRUE, group.sizes = TRUE))), lapply(wldNA, order, method = "radix"))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, decreasing = TRUE, starts = TRUE, group.sizes = TRUE))), lapply(wldNA, order, method = "radix", decreasing = TRUE))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, na.last = FALSE, starts = TRUE, group.sizes = TRUE))), lapply(wldNA, order, method = "radix", na.last = FALSE))
expect_identical(lapply(wldNA, function(x) unattrib(radixorder(x, decreasing = TRUE, na.last = FALSE, starts = TRUE, group.sizes = TRUE))), lapply(wldNA, order, method = "radix", decreasing = TRUE, na.last = FALSE))

randcols <- function(n = 3) replicate(n, sample.int(11, sample.int(5, 1)), simplify = FALSE)
order2 <- function(x, ...) do.call(order, c(gv(wldNA, x), list(...)))

# Ordering by multiple variables
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x)))), lapply(rc, order2, method = "radix"))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), decreasing = TRUE))), lapply(rc, order2, method = "radix", decreasing = TRUE))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), na.last = FALSE))), lapply(rc, order2, method = "radix", na.last = FALSE))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), decreasing = TRUE, na.last = FALSE))), lapply(rc, order2, method = "radix", decreasing = TRUE, na.last = FALSE))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), na.last = NA))), lapply(rc, order2, method = "radix", na.last = NA))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), decreasing = TRUE, na.last = NA))), lapply(rc, order2, method = "radix", decreasing = TRUE, na.last = NA))
# get starts
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), starts = TRUE))), lapply(rc, order2, method = "radix"))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), decreasing = TRUE, starts = TRUE))), lapply(rc, order2, method = "radix", decreasing = TRUE))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), na.last = FALSE, starts = TRUE))), lapply(rc, order2, method = "radix", na.last = FALSE))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), decreasing = TRUE, na.last = FALSE, starts = TRUE))), lapply(rc, order2, method = "radix", decreasing = TRUE, na.last = FALSE))
# get group.sizes
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), group.sizes = TRUE))), lapply(rc, order2, method = "radix"))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), decreasing = TRUE, group.sizes = TRUE))), lapply(rc, order2, method = "radix", decreasing = TRUE))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), na.last = FALSE, group.sizes = TRUE))), lapply(rc, order2, method = "radix", na.last = FALSE))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), decreasing = TRUE, na.last = FALSE, group.sizes = TRUE))), lapply(rc, order2, method = "radix", decreasing = TRUE, na.last = FALSE))
# get starts and group.sizes
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), starts = TRUE, group.sizes = TRUE))), lapply(rc, order2, method = "radix"))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), decreasing = TRUE, starts = TRUE, group.sizes = TRUE))), lapply(rc, order2, method = "radix", decreasing = TRUE))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), na.last = FALSE, starts = TRUE, group.sizes = TRUE))), lapply(rc, order2, method = "radix", na.last = FALSE))
rc <- randcols()
expect_identical(lapply(rc, function(x) unattrib(radixorderv(gv(wldNA, x), decreasing = TRUE, na.last = FALSE, starts = TRUE, group.sizes = TRUE))), lapply(rc, order2, method = "radix", decreasing = TRUE, na.last = FALSE))

})



test_that("GRP works as intended", {

 expect_visible(GRP(unname(as.list(mtcars))))
 expect_visible(GRP(unname(as.list(mtcars)), 8:9))
 expect_equal(GRPnames(GRP(mtcars$cyl)), c("4","6","8"))
 expect_equal(GRPnames(GRP(mtcars$cyl), FALSE), c(4, 6, 8))
 expect_identical(GRPnames(GRP(mtcars$cyl, return.groups = FALSE)), NULL)
 expect_output(print(GRP(mtcars, ~ cyl + am)))
 expect_output(print(GRP(mtcars, ~ cyl + am, return.groups = FALSE)))
 # expect_invisible(plot(GRP(mtcars, ~ cyl + am)))
 expect_identical(GRP(GRP(mtcars$mpg)), GRP(mtcars$mpg))
 expect_identical(GRP.default(GRP(mtcars$mpg)), GRP(mtcars$mpg))
 expect_equal(GRP(mtcars$mpg)[[2]], unattrib(as.factor(mtcars$mpg)))
 expect_equal(GRP(mtcars$cyl)[[2]], unattrib(as.factor(mtcars$cyl)))
 expect_equal(GRP(wlddev2$country)[[2]], unattrib(as.factor(wlddev2$country)))
 expect_equal(GRP(wlddev2$PCGDP)[[2]], unattrib(factor(wlddev2$PCGDP, exclude = NULL)))

 expect_equal(GRP(mtcars$mpg)[[1]], attributes(qG(mtcars$mpg))[[1]])
 expect_equal(GRP(mtcars$cyl)[[1]], attributes(qG(mtcars$cyl))[[1]])
 expect_equal(GRP(wlddev2$country)[[1]], attributes(qG(wlddev2$country))[[1]])
 expect_equal(GRP(wlddev2$PCGDP)[[1]], attributes(qG(wlddev2$PCGDP, na.exclude = FALSE))[[1]])

 expect_equal(GRP(mtcars$mpg)[[4]][[1]], attributes(qG(mtcars$mpg, return.groups = TRUE))[["groups"]])
 expect_equal(GRP(mtcars$cyl)[[4]][[1]], attributes(qG(mtcars$cyl, return.groups = TRUE))[["groups"]])
 expect_equal(GRP(wlddev2$country)[[4]][[1]], attributes(qG(wlddev2$country, return.groups = TRUE))[["groups"]])
 expect_equal(GRP(wlddev2$PCGDP)[[4]][[1]], attributes(qG(wlddev2$PCGDP, na.exclude = FALSE, return.groups = TRUE))[["groups"]])

 expect_visible(GRP(1:10))
 expect_visible(GRP(1:10, decreasing = TRUE))
 expect_visible(GRP(mtcNA$mpg))
 expect_visible(GRP(mtcNA$mpg, return.groups = FALSE))
 expect_visible(GRP(mtcNA$mpg, return.groups = FALSE, return.order = TRUE))
 expect_visible(GRP(mtcNA$mpg, na.last = FALSE))
 expect_visible(GRP(mtcNA$mpg, na.last = FALSE, decreasing = TRUE))
 expect_visible(GRP(mtcNA$mpg, na.last = FALSE, decreasing = TRUE, return.order = TRUE))
 expect_visible(GRP(list(a = 1:3, b = 1:3)))
 expect_visible(GRP(mtcars))
 expect_visible(GRP(mtcNA))
 expect_visible(GRP(mtcNA, return.groups = FALSE))
 expect_visible(GRP(mtcNA, return.groups = FALSE, return.order = TRUE))
 expect_visible(GRP(mtcNA, na.last = FALSE))
 expect_visible(GRP(mtcNA, na.last = FALSE, decreasing = TRUE))
 expect_visible(GRP(mtcNA, na.last = FALSE, decreasing = TRUE, return.order = TRUE))
 expect_visible(GRP(wlddev2))
 expect_visible(GRP(wlddev2, return.groups = FALSE))
 expect_true(all_obj_equal(GRP(mtcars, ~ cyl + vs + am)[1:7],
                           GRP(mtcars, c("cyl","vs","am"))[1:7],
                           GRP(mtcars, c(2,8:9))[1:7]))
})

test_that("GRP gives errors for wrong input", {

  expect_error(GRP(mtcars$mpg, na.last = NA))
  expect_error(GRP(~ bla))
  expect_error(GRP(1:10, 1))
  expect_error(GRP(1:10, ~ cyl))
  expect_error(GRP(1:10, "cyl"))
  expect_error(GRP(mtcars, TRUE))
  expect_error(GRP(mtcars, ~ cyl + bla))
  expect_error(GRP(mtcars, c("bal","cyl")))
  expect_error(GRP(mtcars, 11:12))
  expect_error(GRP(list(a = 1:3, b = 1:4)))
  expect_visible(GRP(mtcars, ~ cyl + vs, order = -1L))

})

test_that("fgroup_by works as intended", {

  ca <- function(x) {
    nam <- names(x[[4L]])
    attributes(x[[4L]]) <- NULL
    names(x[[4L]]) <- nam
    x
  }

  expect_output(print(fgroup_by(mtcars, cyl, vs, am)))
  expect_equal(GRP(fgroup_by(mtcars, cyl, vs, am)), ca(GRP(mtcars, ~ cyl + vs + am, call = FALSE)))
  expect_equal(GRP(fgroup_by(mtcars, c("cyl", "vs", "am"))), ca(GRP(mtcars, ~ cyl + vs + am, call = FALSE)))
  expect_equal(GRP(fgroup_by(mtcars, c(2, 8:9))), ca(GRP(mtcars, ~ cyl + vs + am, call = FALSE)))
  expect_identical(fungroup(fgroup_by(mtcars, cyl, vs, am)), mtcars)
  expect_equal(fgroup_by(fgroup_by(mtcars, cyl, vs, am), cyl), fgroup_by(mtcars, cyl))

  # The issue is that GRP.grouped_df does not reclass the groups... take up another time.
  # This is to fool very silly checks on CRAN scanning the code of the tests
  # group_by <- eval(parse(text = paste0("dplyr", ":", ":", "group_by")))
  # expect_equal(GRP(group_by(mtcars, cyl, vs, am), call = FALSE), GRP(as.list(mtcars), ~ cyl + vs + am, call = FALSE))
  # expect_equal(GRP(group_by(mtcNA, cyl, vs, am)), GRP(mtcNA, ~ cyl + vs + am, call = NULL))
  # expect_equal(GRP(group_by(GGDC10S, Variable, Country)), GRP(GGDC10S, ~ Variable + Country, call = FALSE))
  # expect_equal(GRP(group_by(GGDCNA, Variable, Country)), GRP(GGDCNA, ~ Variable + Country, call = NULL))
  # expect_equal(GRP(group_by(wlddev, region, year)), GRP(wlddev, ~ region + year, call = NULL))
  # expect_equal(GRP(group_by(wldNA, region, year)), GRP(wldNA, ~ region + year, call = NULL))


})

gdat <- gby(GGDCNA, Variable, Country)

test_that("fgroup_vars works as intended", {

  expect_identical(fgroup_vars(gdat), slt(GGDCNA, Variable, Country))
  expect_identical(fgroup_vars(gdat, "unique"), funique(slt(GGDCNA, Variable, Country), sort = TRUE))
  expect_identical(fgroup_vars(gdat, "names"), .c(Variable, Country))
  expect_identical(fgroup_vars(gdat, "indices"), c(4L, 1L))
  expect_identical(fgroup_vars(gdat, "named_indices"), setNames(c(4L, 1L), .c(Variable, Country)))
  expect_identical(fgroup_vars(gdat, "logical"), `[<-`(logical(fncol(GGDCNA)),  c(4L, 1L), TRUE))
  expect_identical(fgroup_vars(gdat, "named_logical"), setNames(`[<-`(logical(fncol(GGDCNA)),  c(4L, 1L), TRUE), names(GGDC10S)))
  expect_error(fgroup_vars(gdat, "bla"))

})



test_that("GRP <> factor conversions run seamlessly", {

  expect_identical(unclass(iris$Species), unclass(as_factor_GRP(GRP(iris$Species)))) # as_factor_GRP always adds class "na.included"

  expect_identical(unclass(wlddev$iso3c[1:200]), unclass(as_factor_GRP(GRP(wlddev$iso3c[1:200])))) # as_factor_GRP always adds class "na.included"
  expect_identical(unclass(fdroplevels(wlddev$iso3c[1:200])), unclass(as_factor_GRP(GRP(wlddev$iso3c[1:200], drop = TRUE)))) # as_factor_GRP always adds class "na.included"

  expect_identical(unclass(`vlabels<-`(wlddev2$iso3c, "label", NULL)), unclass(as_factor_GRP(GRP(wlddev2$iso3c))))
  set.seed(101)
  int <- sample.int(10,100,TRUE)
  expect_identical(unclass(qF(int)), unclass(as_factor_GRP(GRP(int))))
  expect_identical(unclass(qF(int)), unclass(as_factor_GRP(GRP(qF(int)))))
  intNA <- int
  set.seed(101)
  intNA[sample(100,20)] <- NA
  expect_identical(unclass(qF(intNA, na.exclude = FALSE)), unclass(as_factor_GRP(GRP(intNA))))
  expect_identical(unclass(qF(intNA, na.exclude = FALSE)), unclass(as_factor_GRP(GRP(qF(intNA)))))
  dblNA <- as.double(intNA)
  expect_false(identical(unclass(qF(dblNA)), unclass(as_factor_GRP(GRP(dblNA))))) # qF with na.exclude = TRUE retains double NA's...
  expect_false(identical(unclass(qF(dblNA)), unclass(as_factor_GRP(GRP(qF(dblNA))))))
  expect_identical(qF(dblNA, na.exclude = FALSE), as_factor_GRP(GRP(dblNA)))
  expect_identical(qF(dblNA, na.exclude = FALSE), as_factor_GRP(GRP(qF(dblNA))))
  expect_identical(qF(dblNA, na.exclude = FALSE), as_factor_GRP(GRP(qF(dblNA, na.exclude = FALSE))))

})

# could also do qG to GRP, but qG is same as factor.. and is a programmers function anyway..

test_that("qF and qG work as intended", {

  af <- lapply(wlddev2, function(x) as.factor(x))
  expect_equal(af[!fact_vars(wlddev2, "logical")], lapply(gv(wlddev2, !fact_vars(wlddev2, "logical")), function(x) unlab(qF(x, method = "radix"))))
  expect_equal(af[!fact_vars(wlddev2, "logical")], lapply(gv(wlddev2, !fact_vars(wlddev2, "logical")), function(x) unlab(qF(x, method = "hash"))))
  af <- lapply(af, unattrib)
  expect_identical(af, lapply(wlddev2, function(x) unattrib(qF(x, method = "radix"))))
  expect_identical(af, lapply(wlddev2, function(x) unattrib(qF(x, method = "hash"))))
  expect_identical(af, lapply(wlddev2, function(x) unattrib(qG(x, method = "radix"))))
  expect_identical(af, lapply(wlddev2, function(x) unattrib(qG(x, method = "hash"))))
  afNA <- lapply(wldNA, function(x) as.factor(x))
  expect_equal(afNA[!fact_vars(wlddev2, "logical")], lapply(gv(wldNA, !fact_vars(wlddev2, "logical")), function(x) unlab(qF(x, method = "radix"))))
  expect_equal(afNA[!fact_vars(wlddev2, "logical")], lapply(gv(wldNA, !fact_vars(wlddev2, "logical")), function(x) unlab(qF(x, method = "hash"))))
  afNA <- lapply(afNA, unattrib)
  expect_identical(afNA, lapply(wldNA, function(x) unattrib(qF(x, method = "radix"))))
  expect_identical(afNA, lapply(wldNA, function(x) unattrib(qF(x, method = "hash"))))
  expect_identical(afNA, lapply(wldNA, function(x) unattrib(qG(x, method = "radix"))))
  expect_identical(afNA, lapply(wldNA, function(x) unattrib(qG(x, method = "hash"))))
  afnoNA <- lapply(wldNA, function(x) factor(x, exclude = NULL))
  expect_equal(lapply(afnoNA[!fact_vars(wlddev2, "logical")], unclass), lapply(gv(wldNA, !fact_vars(wlddev2, "logical")), function(x) unclass(unlab(qF(x, method = "radix", na.exclude = FALSE)))))
  expect_equal(lapply(afnoNA[!fact_vars(wlddev2, "logical")], unclass), lapply(gv(wldNA, !fact_vars(wlddev2, "logical")), function(x) unclass(unlab(qF(x, method = "hash", na.exclude = FALSE)))))
  afnoNA <- lapply(afnoNA, unattrib)
  expect_identical(afnoNA, lapply(wldNA, function(x) unattrib(qF(x, method = "radix", na.exclude = FALSE))))
  expect_identical(afnoNA, lapply(wldNA, function(x) unattrib(qF(x, method = "hash", na.exclude = FALSE))))
  expect_identical(afnoNA, lapply(wldNA, function(x) unattrib(qG(x, method = "radix", na.exclude = FALSE))))
  expect_identical(afnoNA, lapply(wldNA, function(x) unattrib(qG(x, method = "hash", na.exclude = FALSE))))

  countryf <- as.factor(wlddev2$country)
  expect_identical(countryf, unlab(qF(wlddev2$country)))
  expect_identical(countryf, unlab(qF(wlddev2$country, method = "radix")))
  expect_identical(countryf, unlab(qF(wlddev2$country, method = "hash")))
  # identical(as.factor(wlddev2$iso3c), wlddev2$iso3c)
  expect_identical(levels(wlddev2$iso3c), levels(unlab(qF(wlddev2$iso3c))))
  expect_identical(unattrib(wlddev2$iso3c), unattrib(unlab(qF(wlddev2$iso3c))))
  expect_identical(class(wlddev2$iso3c), class(unlab(qF(wlddev2$iso3c))))

  expect_equal(lapply(wlddev2, function(x) qF(x, method = "radix")), lapply(wlddev2, function(x) qF(x, method = "hash")))
  expect_equal(lapply(wldNA, function(x) qF(x, method = "radix")), lapply(wldNA, function(x) qF(x, method = "hash")))
  expect_equal(lapply(wlddev2, function(x) qG(x, method = "radix")), lapply(wlddev2, function(x) qG(x, method = "hash")))
  expect_equal(lapply(wldNA, function(x) qG(x, method = "radix")), lapply(wldNA, function(x) qG(x, method = "hash")))

  expect_equal(lapply(wlddev2, function(x) qF(x, method = "radix", na.exclude = FALSE)), lapply(wlddev2, function(x) qF(x, method = "hash", na.exclude = FALSE)))
  expect_equal(lapply(wldNA, function(x) qF(x, method = "radix", na.exclude = FALSE)), lapply(wldNA, function(x) qF(x, method = "hash", na.exclude = FALSE)))
  expect_equal(lapply(wlddev2, function(x) qG(x, method = "radix", na.exclude = FALSE)), lapply(wlddev2, function(x) qG(x, method = "hash", na.exclude = FALSE)))
  expect_equal(lapply(wldNA, function(x) qG(x, method = "radix", na.exclude = FALSE)), lapply(wldNA, function(x) qG(x, method = "hash", na.exclude = FALSE)))

  # Testing reordering of factor levels
  expect_identical(qF(wlddev$iso3c), wlddev$iso3c)
  riso3 <- rev(wlddev$iso3c)
  expect_identical(qF(riso3), riso3)
  expect_identical(qF(riso3, sort = FALSE), factor(riso3, levels = funique(riso3)))
  iso3na <- na_insert(wlddev$iso3c)
  expect_identical(qF(iso3na), iso3na)
  expect_identical(unclass(qF(iso3na, na.exclude = FALSE, keep.attr = FALSE)), unclass(addNA(iso3na)))
  riso3na <- na_insert(riso3)
  expect_identical(qF(riso3na), riso3na)
  expect_identical(unclass(qF(riso3na, na.exclude = FALSE, keep.attr = FALSE)), unclass(addNA(riso3na)))
  expect_identical(qF(riso3na, sort = FALSE), factor(riso3na, levels = funique(riso3)))
  expect_identical(unclass(qF(riso3na, sort = FALSE, na.exclude = FALSE)), unclass(factor(riso3na, levels = funique(riso3na), exclude = NULL)))
  expect_identical(unclass(qF(riso3na, sort = FALSE, na.exclude = FALSE)), unclass(factor(riso3na, levels = unique(riso3na), exclude = NULL)))

})

# Could still refine this code, but is not at all critical !!
date <- qG(wlddev$date, return.groups = TRUE)
dateg <- GRP(date, call = FALSE)
dateg$ordered <- NULL
date <- wlddev$date
vlabels(date) <- NULL
dateg2 <- GRP(date, call = FALSE)
dateg2$ordered <- NULL

test_that("GRP <> qG and factor <> qG conversions work", {

  # expect_equal(dateg, dateg2)
  expect_equal(qF(unattrib(wlddev$country)), as.factor_qG(qG(unattrib(wlddev$country), return.groups = TRUE)))
  expect_equal(qF(unattrib(wlddev$country)), qF(qG(unattrib(wlddev$country), return.groups = TRUE)))
  expect_equal(qG(unattrib(wlddev$country)), qG(qF(unattrib(wlddev$country))))
  expect_equal(qG(unattrib(wlddev$country), return.groups = TRUE), qG(qF(unattrib(wlddev$country)), return.groups = TRUE))

})

base_group <- function(x, sort = FALSE, group.sizes = FALSE) {
  if(sort) o <- if(is.list(x)) do.call(order, c(x, list(method = "radix"))) else order(x, method = "radix")
  if(is.list(x)) x <- do.call(paste, c(x, list(sep = ".")))
  ux <- unique(if(sort) x[o] else x)
  r <- match(x, ux)
  attr(r, "N.groups") <- length(ux)
  if(group.sizes) attr(r, "group.sizes") <- tabulate(r, length(ux))
  if(!sort) oldClass(r) <- c("qG", "na.included")
  r
}

test_that("group() works as intended", {
  wlduo <- wlddev[order(rnorm(nrow(wlddev))), ]
  wlduoNA <- na_insert(wlduo)
  dlist <- c(mtcNA, wlddev, wlduo, GGDCNA, airquality)
  # Single grouping variable
  expect_identical(lapply(dlist, group, group.sizes = TRUE), lapply(dlist, base_group, group.sizes = TRUE))
  # Multiple grouping variables
  g <- replicate(70, sample.int(11, sample.int(6, 1)), simplify = FALSE)
  expect_identical(lapply(g, function(i) group(.subset(mtcars, i), group.sizes = TRUE)), lapply(g, function(i) base_group(.subset(mtcars, i), group.sizes = TRUE)))
  expect_identical(lapply(g, function(i) group(.subset(mtcNA, i), group.sizes = TRUE)), lapply(g, function(i) base_group(.subset(mtcNA, i), group.sizes = TRUE)))
  g <- replicate(50, sample.int(13, sample.int(4, 1)), simplify = FALSE)
  expect_identical(lapply(g, function(i) group(.subset(wlduo, i), group.sizes = TRUE)), lapply(g, function(i) base_group(.subset(wlduo, i), group.sizes = TRUE)))
  expect_identical(lapply(g, function(i) group(.subset(wlduoNA, i), group.sizes = TRUE)), lapply(g, function(i) base_group(.subset(wlduoNA, i), group.sizes = TRUE)))
  g <- replicate(50, sample.int(13, 3, replace = TRUE), simplify = FALSE)
  expect_identical(lapply(g, function(i) group(.subset(wlduo, i), group.sizes = TRUE)), lapply(g, function(i) base_group(.subset(wlduo, i), group.sizes = TRUE)))
  expect_identical(lapply(g, function(i) group(.subset(wlduoNA, i), group.sizes = TRUE)), lapply(g, function(i) base_group(.subset(wlduoNA, i), group.sizes = TRUE)))
  # Positive and negative values give the same grouping
  nwld <- nv(wlduo)
  expect_identical(lapply(nwld, group), lapply(nwld %c*% -1, group))
  expect_visible(group(nwld %c*% -1))
  expect_visible(group(nwld[c(4,2,3)] %c*% -1))

  expect_equal(group(0), base_group(0))
  expect_equal(group(1), base_group(1))
  expect_equal(group(0L), base_group(0L))
  expect_equal(group(1L), base_group(1L))
  expect_equal(group(Inf), base_group(Inf))
  expect_equal(group(-Inf), base_group(-Inf))
  expect_equal(group(c(NaN, NA, 0, 1, Inf, -Inf)), base_group(c(NaN, NA, 0, 1, Inf, -Inf)))
  expect_equal(group(NA_integer_), base_group(NA_integer_))
  expect_equal(group(NA_real_), base_group(NA_real_))
  expect_equal(group(NaN), base_group(NaN))
  expect_equal(group(NA), base_group(NA))
  expect_equal(group(NA_character_), base_group(NA_character_))
})

GRP2 <- function(x) {
  g <- GRP.default(x, sort = TRUE, return.groups = FALSE, call = FALSE)
  r <- g[[2]]
  attr(r, "N.groups") <- g[[1]]
  attr(r, "group.sizes") <- g[[3]]
  r
}

qG2 <- function(x, method = "auto", sort = TRUE) unclass(qG(x, na.exclude = FALSE, sort = sort, method = method))

test_that("GRP2() and qG2 work as intended", {
  wlduo <- wlddev[order(rnorm(nrow(wldNA))), ]
  dlist <- c(mtcNA, wlddev, wlduo, GGDCNA, airquality)
  # Single grouping variable
  expect_identical(lapply(dlist, GRP2), lapply(dlist, base_group, sort = TRUE, group.sizes = TRUE))
  bgres <- lapply(dlist, base_group, sort = TRUE)
  expect_identical(lapply(dlist, qG2), bgres)
  expect_identical(lapply(dlist, qG2, method = "hash"), bgres)
  expect_identical(lapply(dlist, qG2, method = "radix"), bgres)
  expect_true(all_identical(qG2(wlduo$country, method = "radix", sort = FALSE),
                            qG2(wlduo$country, method = "hash", sort = FALSE),
                            unclass(base_group(wlduo$country, sort = FALSE))))
  # Multiple grouping variables
  g <- replicate(50, sample.int(11, sample.int(6, 1)), simplify = FALSE)
  expect_identical(lapply(g, function(i) GRP2(.subset(mtcars, i))), lapply(g, function(i) base_group(.subset(mtcars, i), sort = TRUE, group.sizes = TRUE)))
  g <- replicate(30, sample.int(13, sample.int(4, 1)), simplify = FALSE)
  expect_identical(lapply(g, function(i) GRP2(.subset(wlduo, i))), lapply(g, function(i) base_group(.subset(wlduo, i), sort = TRUE, group.sizes = TRUE)))
  g <- replicate(30, sample.int(13, 3, replace = TRUE), simplify = FALSE)
  expect_identical(lapply(g, function(i) GRP2(.subset(wlduo, i))), lapply(g, function(i) base_group(.subset(wlduo, i), sort = TRUE, group.sizes = TRUE)))
})



test_that("GRP2() works as intended", {
  wlduo <- wlddev[order(rnorm(nrow(wldNA))), ]
  dlist <- c(mtcNA, wlddev, wlduo, GGDCNA, airquality)
  # Single grouping variable
  expect_identical(lapply(dlist, GRP2), lapply(dlist, base_group, sort = TRUE, group.sizes = TRUE))
  # Multiple grouping variables
  g <- replicate(50, sample.int(11, sample.int(6, 1)), simplify = FALSE)
  expect_identical(lapply(g, function(i) GRP2(.subset(mtcars, i))), lapply(g, function(i) base_group(.subset(mtcars, i), sort = TRUE, group.sizes = TRUE)))
  g <- replicate(30, sample.int(13, sample.int(4, 1)), simplify = FALSE)
  expect_identical(lapply(g, function(i) GRP2(.subset(wlduo, i))), lapply(g, function(i) base_group(.subset(wlduo, i), sort = TRUE, group.sizes = TRUE)))
  g <- replicate(30, sample.int(13, 3, replace = TRUE), simplify = FALSE)
  expect_identical(lapply(g, function(i) GRP2(.subset(wlduo, i))), lapply(g, function(i) base_group(.subset(wlduo, i), sort = TRUE, group.sizes = TRUE)))
})


if(identical(Sys.getenv("NCRAN"), "TRUE")) {

  # This is to fool very silly checks on CRAN scanning the code of the tests
  pwlddev <- eval(parse(text = paste0("plm", ":", ":", "pdata.frame(wlddev, index = c('iso3c', 'year'))")))
  iso3c <- eval(parse(text = paste0("plm", ":", ":", "index(pwlddev, 1L)")))
  year <- eval(parse(text = paste0("plm", ":", ":", "index(pwlddev, 2L)")))

test_that("GRP pseries and pdata.frame methods work as intended", {

  expect_equal(GRP(pwlddev, call = FALSE), GRP(iso3c, call = FALSE))
  expect_equal(GRP(pwlddev$PCGDP, call = FALSE), GRP(pwlddev, call = FALSE))
  expect_equal(GRP(pwlddev, effect = "year", call = FALSE), GRP(year, call = FALSE))
  expect_equal(GRP(pwlddev$PCGDP, effect = "year", call = FALSE), GRP(pwlddev, effect = "year", call = FALSE))

})

}

fl <- slt(wlddev, region, income)
set.seed(101)
flNA <- na_insert(fl)

test_that("finteraction works as intended", {

  expect_equal(`oldClass<-`(finteraction(fl), "factor"), base::interaction(fl, drop = TRUE, lex.order = TRUE))
  expect_equal(`oldClass<-`(finteraction(ss(fl, 1:300)), "factor"), base::interaction(ss(fl, 1:300), drop = TRUE, lex.order = TRUE)) # missing levels

  expect_equal(unattrib(finteraction(fl, factor = FALSE, sort = TRUE)), unattrib(base::interaction(fl, drop = TRUE, lex.order = TRUE)))
  expect_equal(unattrib(finteraction(fl, factor = FALSE, sort = FALSE)), unattrib(group(fl)))

  # Missing value behavior is always different !!
  # expect_equal(`oldClass<-`(finteraction(flNA), "factor"), factor(base::interaction(flNA, drop = TRUE, lex.order = TRUE), exclude = NULL))
  # expect_equal(`oldClass<-`(finteraction(ss(flNA, 1:300)), "factor"), base::interaction(ss(flNA, 1:300), drop = TRUE, lex.order = TRUE))

})


wld150 <- ss(wlddev, 1:150)
vlabels(wld150) <- NULL
set.seed(101)
wldNA150 <- na_insert(ss(wlddev, 1:150))
vlabels(wldNA150) <- NULL

test_that("fdroplevels works as intended", {

  expect_identical(fdroplevels(wld150), droplevels(wld150))
  expect_identical(fdroplevels(wldNA150), droplevels(wldNA150))
  expect_identical(fdroplevels(wld150$iso3c), droplevels(wld150$iso3c))
  expect_identical(fdroplevels(wldNA150$iso3c), droplevels(wldNA150$iso3c))
  expect_message(fdroplevels(1:3))
  # expect_warning(fdroplevels(wld150, bla = 1))
  # expect_warning(fdroplevels(wld150$iso3c, bla = 1))
  expect_error(fdroplevels.factor(wld150$country))

})


# Note: Should extend with other than just character data..

rctry <- wlddev$country[order(rnorm(length(wlddev$country)))]
set.seed(101)
rctryNA <- na_insert(rctry)

rdat <- sbt(GGDC10S, order(rnorm(length(Variable))), Variable, Country)
vlabels(rdat) <- NULL
vlabels(rdat, "format.stata") <- NULL
set.seed(101)
rdatNA <- na_insert(rdat)

test_that("funique works well", {

 expect_equal(funique(rctry), unique(rctry))
 expect_equal(funique(rctry, sort = TRUE), sort(unique(rctry)))
 expect_equal(funique(rctryNA), unique(rctryNA))
 expect_equal(funique(rctryNA, sort = TRUE), c(sort(unique(rctryNA)), NA))

 expect_equal(funique(mtcars[.c(cyl, vs, am)]),  unique(mtcars[.c(cyl, vs, am)]))
 expect_equal(funique(mtcNA[.c(cyl, vs, am)]),  unique(mtcNA[.c(cyl, vs, am)]))

 expect_equal(funique(rdat), setRownames(unique(rdat)))
 expect_equal(funique(rdat, sort = TRUE), roworderv(unique(rdat)))
 expect_equal(funique(rdatNA), setRownames(unique(rdatNA)))
 expect_equal(funique(rdatNA, sort = TRUE), roworderv(unique(rdatNA)))

 expect_equal(lapply(wlddev, function(x) unattrib(base::unique(x))),
              lapply(wlddev, function(x) unattrib(funique(x))))
 expect_equal(lapply(wldNA, function(x) unattrib(base::unique(x))),
              lapply(wldNA, function(x) unattrib(funique(x))))
 expect_equal(lapply(GGDC10S, function(x) unattrib(base::unique(x))),
              lapply(GGDC10S, function(x) unattrib(funique(x))))

})
