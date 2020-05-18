context("GRP, qF, qG")

mtcNA <- na_insert(mtcars)
wlddev2 <- slt(wlddev, -date)
num_vars(wlddev2) <- round(num_vars(wlddev2), 8)
wldNA <- na_insert(wlddev2)

unlab <- function(x) `attr<-`(x, "label", NULL)

test_that("GRP works as intended", {

 expect_visible(GRP(unname(as.list(mtcars))))
 expect_visible(GRP(unname(as.list(mtcars)), 8:9))
 expect_visible(GRP(mtcars$mpg))
 expect_equal(GRP(mtcars$mpg)[[2]], unattrib(as.factor(mtcars$mpg)))
 expect_equal(GRP(mtcars$cyl)[[2]], unattrib(as.factor(mtcars$cyl)))
 expect_equal(GRP(wlddev2$country)[[2]], unattrib(as.factor(wlddev2$country)))
 expect_equal(GRP(wlddev2$PCGDP)[[2]], unattrib(factor(wlddev2$PCGDP, exclude = NULL)))

 expect_equal(GRP(mtcars$mpg)[[1]], attributes(qG(mtcars$mpg))[[1]])
 expect_equal(GRP(mtcars$cyl)[[1]], attributes(qG(mtcars$cyl))[[1]])
 expect_equal(GRP(wlddev2$country)[[1]], attributes(qG(wlddev2$country))[[1]])
 expect_equal(GRP(wlddev2$PCGDP)[[1]], attributes(qG(wlddev2$PCGDP, na.exclude = FALSE))[[1]])

 expect_equal(GRP(mtcars$mpg)[[4]][[1]], attributes(qG(mtcars$mpg, return.groups = TRUE))[[2]])
 expect_equal(GRP(mtcars$cyl)[[4]][[1]], attributes(qG(mtcars$cyl, return.groups = TRUE))[[2]])
 expect_equal(GRP(wlddev2$country)[[4]][[1]], attributes(qG(wlddev2$country, return.groups = TRUE))[[2]])
 expect_equal(GRP(wlddev2$PCGDP)[[4]][[1]], attributes(qG(wlddev2$PCGDP, na.exclude = FALSE, return.groups = TRUE))[[2]])

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

  expect_error(GRP(~ bla))
  expect_error(GRP(1:10, 1))
  expect_error(GRP(1:10, ~ cyl))
  expect_error(GRP(1:10, "cyl"))
  expect_error(GRP(mtcars, TRUE))
  expect_error(GRP(mtcars, ~ cyl + bla))
  expect_error(GRP(mtcars, c("bal","cyl")))
  expect_error(GRP(mtcars, 11:12))
  expect_error(GRP(list(a = 1:3, b = 1:4)))

})

test_that("GRP <> factor conversions run seamlessly", {

  expect_identical(unclass(iris$Species), unclass(as.factor.GRP(GRP(iris$Species)))) # as.factor.GRP always adds class "na.included"
  expect_identical(unclass(`vlabels<-`(wlddev2$iso3c, "label", NULL)), unclass(as.factor.GRP(GRP(wlddev2$iso3c))))
  int <- sample.int(10,100,TRUE)
  expect_identical(unclass(qF(int)), unclass(as.factor.GRP(GRP(int))))
  expect_identical(unclass(qF(int)), unclass(as.factor.GRP(GRP(qF(int)))))
  intNA <- int
  intNA[sample(100,20)] <- NA
  expect_identical(unclass(qF(intNA, na.exclude = FALSE)), unclass(as.factor.GRP(GRP(intNA))))
  expect_identical(unclass(qF(intNA, na.exclude = FALSE)), unclass(as.factor.GRP(GRP(qF(intNA)))))
  dblNA <- as.double(intNA)
  expect_false(identical(unclass(qF(dblNA)), unclass(as.factor.GRP(GRP(dblNA))))) # qF with na.exclude = TRUE retains double NA's...
  expect_false(identical(unclass(qF(dblNA)), unclass(as.factor.GRP(GRP(qF(dblNA))))))
  expect_identical(qF(dblNA, na.exclude = FALSE), as.factor.GRP(GRP(dblNA)))
  expect_identical(qF(dblNA, na.exclude = FALSE), as.factor.GRP(GRP(qF(dblNA))))
  expect_identical(qF(dblNA, na.exclude = FALSE), as.factor.GRP(GRP(qF(dblNA, na.exclude = FALSE))))

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

  # countryf <- as.factor(wlddev2$country)
  # expect_identical(countryf, unlab(qF(wlddev2$country)))
  # expect_identical(countryf, unlab(qF(wlddev2$country, method = "radix")))
  # identical(as.factor(wlddev2$iso3c), wlddev2$iso3c)
  # expect_identical(levels(wlddev2$iso3c), levels(unlab(qF(wlddev2$iso3c))))
  # expect_identical(unattrib(wlddev2$iso3c), unattrib(unlab(qF(wlddev2$iso3c))))
  # expect_identical(class(wlddev2$iso3c), class(unlab(qF(wlddev2$iso3c))))

  expect_equal(lapply(wlddev2, function(x) qF(x, method = "radix")), lapply(wlddev2, function(x) qF(x, method = "hash")))
  expect_equal(lapply(wldNA, function(x) qF(x, method = "radix")), lapply(wldNA, function(x) qF(x, method = "hash")))
  expect_equal(lapply(wlddev2, function(x) qG(x, method = "radix")), lapply(wlddev2, function(x) qG(x, method = "hash")))
  expect_equal(lapply(wldNA, function(x) qG(x, method = "radix")), lapply(wldNA, function(x) qG(x, method = "hash")))

  expect_equal(lapply(wlddev2, function(x) qF(x, method = "radix", na.exclude = FALSE)), lapply(wlddev2, function(x) qF(x, method = "hash", na.exclude = FALSE)))
  expect_equal(lapply(wldNA, function(x) qF(x, method = "radix", na.exclude = FALSE)), lapply(wldNA, function(x) qF(x, method = "hash", na.exclude = FALSE)))
  expect_equal(lapply(wlddev2, function(x) qG(x, method = "radix", na.exclude = FALSE)), lapply(wlddev2, function(x) qG(x, method = "hash", na.exclude = FALSE)))
  expect_equal(lapply(wldNA, function(x) qG(x, method = "radix", na.exclude = FALSE)), lapply(wldNA, function(x) qG(x, method = "hash", na.exclude = FALSE)))

})
