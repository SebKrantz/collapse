context("flag / L / F")

if(!is.null(attributes(identical(FALSE, TRUE)))) stop("OECD label issue")

# rm(list = ls())
# TODO: test computations on irregular time series and panels
set.seed(101)
x <- abs(10*rnorm(100))
xNA <- x
xNA[sample.int(100, 20)] <- NA
f <- as.factor(rep(1:10, each = 10))
t <- as.factor(rep(1:10, 10))

data <- setRownames(wlddev[wlddev$iso3c %in% c("BLZ","IND","USA","SRB","GRL"), ])
g <- GRP(droplevels(data$iso3c))
td <- as.factor(data$year)
dataNA <- na_insert(data)
m <- qM(data)
suppressWarnings(storage.mode(m) <- "numeric")
mNAc <- qM(dataNA)
mNA <- mNAc
suppressWarnings(storage.mode(mNA) <- "numeric")

# Creatung unordered data:
o = order(rnorm(100))
xuo = x[o]
xNAuo = xNA[o]
fuo = f[o]
tuo = t[o]
t2uo = seq_len(100)[o]
o = order(o)

od = order(rnorm(length(td)))
muo = m[od, ]
datauo = data[od, ]
guo = as_factor_GRP(g)[od]
tduo = td[od]
t2duo = seq_along(od)[od]
od = order(od)

baselag <- function(x, n = 1) c(rep(NA_real_, n), x[1:(length(x)-n)])
baselead <- function(x, n = 1) c(rep(NA_real_, n), x[1:(length(x)-n)])


# flag

test_that("flag performs like baselag", {
  expect_equal(flag(1:10), baselag(1:10))
  expect_equal(flag(1:10, 2), baselag(1:10, 2))
  expect_equal(flag(-1:1), baselag(-1:1))
  expect_equal(flag(x), baselag(x))
  expect_equal(flag(x, 2), baselag(x, 2))
  expect_equal(flag(xNA), baselag(xNA))
  expect_equal(flag(xNA, 2), baselag(xNA, 2))
  expect_equal(flag(m, stubs = FALSE), dapply(m, baselag))
  expect_equal(flag(m, 2, stubs = FALSE), dapply(m, baselag, 2))
  expect_equal(flag(mNA, stubs = FALSE), dapply(mNA, baselag))
  expect_equal(flag(mNA, 2, stubs = FALSE), dapply(mNA, baselag, 2))
  expect_equal(flag(num_vars(data), stubs = FALSE), dapply(num_vars(data), baselag))
  expect_equal(flag(num_vars(data), 2, stubs = FALSE), dapply(num_vars(data), baselag, 2))
  expect_equal(flag(num_vars(dataNA), stubs = FALSE), dapply(num_vars(dataNA), baselag))
  expect_equal(flag(num_vars(dataNA), 2, stubs = FALSE), dapply(num_vars(dataNA), baselag, 2))
  expect_equal(flag(x, 1, f), BY(x, f, baselag, use.g.names = FALSE))
  expect_equal(flag(x, 2, f), BY(x, f, baselag, 2, use.g.names = FALSE))
  expect_equal(flag(xNA, 1, f), BY(xNA, f, baselag, use.g.names = FALSE))
  expect_equal(flag(xNA, 2, f), BY(xNA, f, baselag, 2, use.g.names = FALSE))
  expect_equal(flag(m, 1, g, stubs = FALSE), BY(m, g, baselag, use.g.names = FALSE))
  expect_equal(flag(m, 2, g, stubs = FALSE), BY(m, g, baselag, 2, use.g.names = FALSE))
  expect_equal(flag(mNA, 1, g, stubs = FALSE), BY(mNA, g, baselag, use.g.names = FALSE))
  expect_equal(flag(mNA, 2, g, stubs = FALSE), BY(mNA, g, baselag, 2, use.g.names = FALSE))
  expect_equal(flag(num_vars(data), 1, g, stubs = FALSE), BY(num_vars(data), g, baselag, use.g.names = FALSE))
  expect_equal(flag(num_vars(data), 2, g, stubs = FALSE), BY(num_vars(data), g, baselag, 2, use.g.names = FALSE))
  expect_equal(flag(num_vars(dataNA), 1, g, stubs = FALSE), BY(num_vars(dataNA), g, baselag, use.g.names = FALSE))
  expect_equal(flag(num_vars(dataNA), 2, g, stubs = FALSE), BY(num_vars(dataNA), g, baselag, 2, use.g.names = FALSE))
  # Adding time-variable: Computing fully identified panel-lags !!
  expect_equal(flag(x, 1, f, t), BY(x, f, baselag, use.g.names = FALSE))
  expect_equal(flag(x, 2, f, t), BY(x, f, baselag, 2, use.g.names = FALSE))
  expect_equal(flag(xNA, 1, f, t), BY(xNA, f, baselag, use.g.names = FALSE))
  expect_equal(flag(xNA, 2, f, t), BY(xNA, f, baselag, 2, use.g.names = FALSE))
  expect_equal(flag(m, 1, g, td, stubs = FALSE), BY(m, g, baselag, use.g.names = FALSE))
  expect_equal(flag(m, 2, g, td, stubs = FALSE), BY(m, g, baselag, 2, use.g.names = FALSE))
  expect_equal(flag(mNA, 1, g, td, stubs = FALSE), BY(mNA, g, baselag, use.g.names = FALSE))
  expect_equal(flag(mNA, 2, g, td, stubs = FALSE), BY(mNA, g, baselag, 2, use.g.names = FALSE))
  expect_equal(flag(num_vars(data), 1, g, td, stubs = FALSE), BY(num_vars(data), g, baselag, use.g.names = FALSE))
  expect_equal(flag(num_vars(data), 2, g, td, stubs = FALSE), BY(num_vars(data), g, baselag, 2, use.g.names = FALSE))
  expect_equal(flag(num_vars(dataNA), 1, g, td, stubs = FALSE), BY(num_vars(dataNA), g, baselag, use.g.names = FALSE))
  expect_equal(flag(num_vars(dataNA), 2, g, td, stubs = FALSE), BY(num_vars(dataNA), g, baselag, 2, use.g.names = FALSE))
})

test_that("flag performs (panel-) vector lags and leads without errors", {
  expect_visible(flag(1:10, -2:2))
  expect_visible(flag(1:10, 1:2))
  expect_visible(flag(1:10, -1:-2))
  expect_visible(flag(1:10, 0))

  expect_visible(flag(xNA, -2:2))
  expect_visible(flag(xNA, 1:2))
  expect_visible(flag(xNA, -1:-2))
  expect_visible(flag(xNA, 0))

  expect_visible(flag(xNA, -2:2, f))
  expect_visible(flag(xNA, 1:2, f))
  expect_visible(flag(xNA, -1:-2, f))
  expect_visible(flag(xNA, 0, f))

  expect_visible(flag(xNA, -2:2, f, t))
  expect_visible(flag(xNA, 1:2, f, t))
  expect_visible(flag(xNA, -1:-2, f, t))
  expect_visible(flag(xNA, 0, f, t))
})

test_that("flag performs (panel-) matrix lags and leads without errors", {

  expect_visible(flag(m, -2:2))
  expect_visible(flag(m, 1:2))
  expect_visible(flag(m, -1:-2))
  expect_visible(flag(m, 0))

  expect_visible(flag(m, -2:2, g))
  expect_visible(flag(m, 1:2, g))
  expect_visible(flag(m, -1:-2, g))
  expect_visible(flag(m, 0, g))

  expect_visible(flag(m, -2:2, g, td))
  expect_visible(flag(m, 1:2, g, td))
  expect_visible(flag(m, -1:-2, g, td))
  expect_visible(flag(m, 0, g, td))
})

test_that("flag performs (panel-) data.frame lags and leads without errors", {

  expect_visible(flag(data, -2:2))
  expect_visible(flag(data, 1:2))
  expect_visible(flag(data, -1:-2))
  expect_visible(flag(data, 0))

  expect_visible(flag(data, -2:2, g))
  expect_visible(flag(data, 1:2, g))
  expect_visible(flag(data, -1:-2, g))
  expect_visible(flag(data, 0, g))

  expect_visible(flag(data, -2:2, g, td))
  expect_visible(flag(data, 1:2, g, td))
  expect_visible(flag(data, -1:-2, g, td))
  expect_visible(flag(data, 0, g, td))
})

test_that("flag correctly handles unordered time-series and panel-series computations", {
  expect_equal(flag(x, -2:2, t = 1:100), flag(x, -2:2))
  expect_equal(flag(xNA, -2:2, t = 1:100), flag(xNA, -2:2))
  expect_equal(flag(m, -2:2, t = seq_along(td)), flag(m, -2:2))
  expect_equal(flag(data, -2:2, t = seq_along(td)), flag(data, -2:2))

  expect_equal(flag(xuo, -2:2, t = t2uo)[o,], unclass(flag(x, -2:2)))
  expect_equal(flag(xNAuo, -2:2, t = t2uo)[o,], unclass(flag(xNA, -2:2)))
  expect_equal(flag(muo, -2:2, t = t2duo)[od,], unclass(flag(m, -2:2)))
  expect_equal(flag(datauo, -2:2, t = t2duo)[od,], flag(data, -2:2))

  expect_equal(flag(xuo, -2:2, fuo, tuo)[o,], unclass(flag(x, -2:2, f, t)))
  expect_equal(flag(xNAuo, -2:2, fuo, tuo)[o,], unclass(flag(xNA, -2:2, f, t)))
  expect_equal(flag(muo, -2:2, guo, tduo)[od,], unclass(flag(m, -2:2, g, td)))
  expect_equal(flag(datauo, -2:2, guo, tduo)[od,], flag(data, -2:2, g, td))
})

test_that("flag performs numerically stable in ordered computations", {
  expect_true(all_obj_equal(replicate(50, flag(x), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(xNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(m), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(mNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(data), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(dataNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(x, 1, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(x, -2:2, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(xNA, 1, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(xNA, -2:2, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(m, 1, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(m, -2:2, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(mNA, 1, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(mNA, -2:2, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(data, 1, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(data, -2:2, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(dataNA, 1, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(dataNA, -2:2, g), simplify = FALSE)))
})

test_that("flag performs numerically stable in unordered computations", {
  expect_true(all_obj_equal(replicate(50, flag(xuo, t = t2uo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(xNAuo, t = t2uo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(muo, t = t2duo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(datauo, t = t2duo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(xuo, 1, fuo, tuo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(xuo, -2:2, fuo, tuo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(muo, 1, guo, tduo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(muo, -2:2, guo, tduo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(datauo, 1, guo, tduo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, flag(datauo, -2:2, guo, tduo), simplify = FALSE)))
})

test_that("flag handles special values in the right way", {
  # zero
  expect_equal(flag(c("a","b"),0), c("a","b"))
  expect_equal(flag(c(NaN,NaN),0), c(NaN,NaN))
  expect_equal(flag(c(Inf,Inf),0), c(Inf,Inf))
  expect_equal(flag(c(FALSE,TRUE),0), c(FALSE,TRUE))
  expect_equal(flag(c(TRUE,FALSE),0), c(TRUE,FALSE))
  # lags
  expect_equal(flag(c("a","b")), c(NA,"a"))
  expect_equal(flag(c(1,NA)), c(NA_real_,1))
  expect_equal(flag(c(NA,1)), c(NA_real_,NA_real_))
  expect_equal(flag(c(NaN,1)), c(NA_real_,NaN))
  expect_equal(flag(c(1,NaN)), c(NA_real_,1))
  expect_equal(flag(c(Inf,1)), c(NA,Inf))
  expect_equal(flag(c(1,Inf)), c(NA,1))
  expect_equal(flag(c(Inf,NA)), c(NA_real_,Inf))
  expect_equal(flag(c(NA,Inf)), c(NA_real_,NA_real_))
  expect_equal(flag(c(Inf,-Inf)), c(NA,Inf))
  expect_equal(flag(c(-Inf,Inf)), c(NA,-Inf))
  expect_equal(flag(c(Inf,Inf)), c(NA,Inf))
  expect_equal(flag(c(TRUE,TRUE)), c(NA,TRUE))
  expect_equal(flag(c(TRUE,FALSE)), c(NA,TRUE))
  expect_equal(flag(c(FALSE,TRUE)), c(NA,FALSE))
  # leads
  expect_equal(flag(c("a","b"),-1), c("b",NA))
  expect_equal(flag(c(1,NA),-1), c(NA_real_,NA_real_))
  expect_equal(flag(c(NA,1),-1), c(1,NA_real_))
  expect_equal(flag(c(NaN,1),-1), c(1,NA_real_))
  expect_equal(flag(c(1,NaN),-1), c(NaN,NA_real_))
  expect_equal(flag(c(Inf,1),-1), c(1,NA_real_))
  expect_equal(flag(c(1,Inf),-1), c(Inf,NA_real_))
  expect_equal(flag(c(Inf,NA),-1), c(NA_real_,NA_real_))
  expect_equal(flag(c(NA,Inf),-1), c(Inf,NA_real_))
  expect_equal(flag(c(Inf,-Inf),-1), c(-Inf,NA_real_))
  expect_equal(flag(c(-Inf,Inf),-1), c(Inf,NA_real_))
  expect_equal(flag(c(Inf,Inf),-1), c(Inf,NA_real_))
  expect_equal(flag(c(TRUE,TRUE),-1), c(TRUE,NA))
  expect_equal(flag(c(TRUE,FALSE),-1), c(FALSE,NA))
  expect_equal(flag(c(FALSE,TRUE),-1), c(TRUE,NA))

})

test_that("flag produces errors for wrong input", {
  # type: normally guaranteed by C++
  expect_visible(flag(mNAc))
  expect_visible(flag(wlddev))
  expect_error(flag(mNAc, f))
  expect_error(flag(x, "1"))
  # if n exceeds length(x), should give error
  expect_error(flag(x,101))
  expect_error(flag(x,-101))
  # if n exceeds average group size, should give error
  # expect_warning(flag(x,11,f)) # Some fail on i386 ??
  # expect_warning(flag(x,11,f,t))
  # expect_warning(flag(x,-11,f))
  # expect_warning(flag(x,-11,f,t))
  # passing repeated n-values or non-positive or non-consecutive diff values should give error
  expect_error(flag(x,c(1,1)))
  expect_error(flag(x,c(-1,-1)))
  expect_visible(flag(x,2:1))
  expect_visible(flag(x,0))
  expect_error(flag(x,f))   # common source of error probably is passing the factor in a wrong slot
  expect_error(flag(x,c(1,1),f))
  expect_error(flag(x,c(1,1),f,t))
  expect_visible(flag(x,2:1,f))
  expect_visible(flag(x,2:1,f,t))
  expect_visible(flag(x,0,f))
  expect_visible(flag(x,0,f,t))
  expect_error(flag(x,1,1)) # wrong inputs: passing a non-existent difference argument..
  expect_error(flag(x,1,1))
  expect_error(flag(x,1,1,f))
  expect_error(flag(x,1,1,f,t))
  expect_error(flag(x,1,-1,f))
  expect_error(flag(x,1,-1,f,t))
  # repeated values or gaps in time-variable should give error
  expect_error(flag(1:3, t = c(1,1,2)))
  expect_error(flag(1:3, t = c(1,2,2)))
  expect_error(flag(1:3, t = c(1,2,1)))
  expect_error(flag(1:10, g = c(1,1,1,2,2,2,3,3,3,3), t = c(1,1,2,1:3,1:4))) # This is the only possible statement which does not throw a reteated timevar error because the first C++ index is 0, and omap is also initialized with 0's.
  expect_error(flag(1:10, g = c(1,1,1,2,2,2,3,3,3,3), t = c(1,1,1,1:3,1:4)))
  expect_error(flag(1:10, g = c(1,1,1,2,2,2,3,3,3,3), t = c(1:3,1:3,1,1,3,4)))
  expect_error(flag(1:10, g = c(1,1,1,2,2,2,3,3,3,3), t = c(1,2,2,1:3,1:4)))
  expect_visible(flag(1:10, g = c(1,1,1,2,2,2,3,3,3,3), t = c(1,2,4,1:3,1:4)))
  expect_error(flag(1:10, g = c(1,2,1,2,2,2,3,3,3,3), t = c(1:3,1:3,1:4)))
  expect_visible(flag(1:10, g = c(1,1,1,2,2,2,3,3,4,3), t = c(1:3,1:3,1:4)))
  # The usual stuff: Wrongly sized grouping vectors or time-variables
  expect_error(flag(1:3, t = 1:2))
  expect_error(flag(1:3, t = 1:4))
  expect_error(flag(1:3, g = 1:2))
  expect_error(flag(1:3, g = 1:4))
  expect_error(flag(1:4, g = c(1,1,2,2), t = c(1,2,1)))
  expect_error(flag(1:4, g = c(1,2,2), t = c(1,2,1,2)))
})

# L and F

F <- getNamespace("collapse")$F

test_that("F performs like baselead", {
  expect_equal(F(1:10, -1), baselead(1:10))
  expect_equal(F(1:10, -2), baselead(1:10, 2))
  expect_equal(F(-1:1, -1), baselead(-1:1))
  expect_equal(F(x, -1), baselead(x))
  expect_equal(F(x, -2), baselead(x, 2))
  expect_equal(F(xNA, -1), baselead(xNA))
  expect_equal(F(xNA, -2), baselead(xNA, 2))
  expect_equal(F(m, -1, stubs = FALSE), dapply(m, baselead))
  expect_equal(F(m, -2, stubs = FALSE), dapply(m, baselead, 2))
  expect_equal(F(mNA, -1, stubs = FALSE), dapply(mNA, baselead))
  expect_equal(F(mNA, -2, stubs = FALSE), dapply(mNA, baselead, 2))
  expect_equal(F(num_vars(data), -1, stubs = FALSE), dapply(num_vars(data), baselead))
  expect_equal(F(num_vars(data), -2, stubs = FALSE), dapply(num_vars(data), baselead, 2))
  expect_equal(F(num_vars(dataNA), -1, stubs = FALSE), dapply(num_vars(dataNA), baselead))
  expect_equal(F(num_vars(dataNA), -2, stubs = FALSE), dapply(num_vars(dataNA), baselead, 2))
  expect_equal(F(x, -1, f), BY(x, f, baselead, use.g.names = FALSE))
  expect_equal(F(x, -2, f), BY(x, f, baselead, 2, use.g.names = FALSE))
  expect_equal(F(xNA, -1, f), BY(xNA, f, baselead, use.g.names = FALSE))
  expect_equal(F(xNA, -2, f), BY(xNA, f, baselead, 2, use.g.names = FALSE))
  expect_equal(F(m, -1, g, stubs = FALSE), BY(m, g, baselead, use.g.names = FALSE))
  expect_equal(F(m, -2, g, stubs = FALSE), BY(m, g, baselead, 2, use.g.names = FALSE))
  expect_equal(F(mNA, -1, g, stubs = FALSE), BY(mNA, g, baselead, use.g.names = FALSE))
  expect_equal(F(mNA, -2, g, stubs = FALSE), BY(mNA, g, baselead, 2, use.g.names = FALSE))
  expect_equal(F(num_vars(data), -1, g, stubs = FALSE), BY(num_vars(data), g, baselead, use.g.names = FALSE))
  expect_equal(F(num_vars(data), -2, g, stubs = FALSE), BY(num_vars(data), g, baselead, 2, use.g.names = FALSE))
  expect_equal(F(num_vars(dataNA), -1, g, stubs = FALSE), BY(num_vars(dataNA), g, baselead, use.g.names = FALSE))
  expect_equal(F(num_vars(dataNA), -2, g, stubs = FALSE), BY(num_vars(dataNA), g, baselead, 2, use.g.names = FALSE))
  # Adding time-variable: Computing fully identified panel-lags !!
  expect_equal(F(x, -1, f, t), BY(x, f, baselead, use.g.names = FALSE))
  expect_equal(F(x, -2, f, t), BY(x, f, baselead, 2, use.g.names = FALSE))
  expect_equal(F(xNA, -1, f, t), BY(xNA, f, baselead, use.g.names = FALSE))
  expect_equal(F(xNA, -2, f, t), BY(xNA, f, baselead, 2, use.g.names = FALSE))
  expect_equal(F(m, -1, g, td, stubs = FALSE), BY(m, g, baselead, use.g.names = FALSE))
  expect_equal(F(m, -2, g, td, stubs = FALSE), BY(m, g, baselead, 2, use.g.names = FALSE))
  expect_equal(F(mNA, -1, g, td, stubs = FALSE), BY(mNA, g, baselead, use.g.names = FALSE))
  expect_equal(F(mNA, -2, g, td, stubs = FALSE), BY(mNA, g, baselead, 2, use.g.names = FALSE))
  expect_equal(F(num_vars(data), -1, g, td, stubs = FALSE), BY(num_vars(data), g, baselead, use.g.names = FALSE))
  expect_equal(F(num_vars(data), -2, g, td, stubs = FALSE), BY(num_vars(data), g, baselead, 2, use.g.names = FALSE))
  expect_equal(F(num_vars(dataNA), -1, g, td, stubs = FALSE), BY(num_vars(dataNA), g, baselead, use.g.names = FALSE))
  expect_equal(F(num_vars(dataNA), -2, g, td, stubs = FALSE), BY(num_vars(dataNA), g, baselead, 2, use.g.names = FALSE))
})

test_that("L and F do the opposite of one another", {
  expect_equal(L(1:10, -2:2), F(1:10, 2:-2))
  expect_equal(L(m, -2:2), F(m, 2:-2))
  expect_equal(L(data, -2:2), F(data, 2:-2))
})

test_that("L produces errors for wrong input", {
  # type: normally guaranteed by C++
  expect_visible(L(mNAc))
  expect_visible(L(wlddev))
  expect_error(L(mNAc, f))
  expect_error(L(x, "1"))
  # if n exceeds length(x), should give error
  expect_error(L(x,101))
  expect_error(L(x,-101))
  # if n exceeds average group size, should give error
  # expect_warning(L(x,11,f)) -> some fail on i336
  # expect_warning(L(x,11,f,t))
  # expect_warning(L(x,-11,f))
  # expect_warning(L(x,-11,f,t))
  # passing repeated n-values or non-positive or non-consecutive diff values should give error
  expect_error(L(x,c(1,1)))
  expect_error(L(x,c(-1,-1)))
  expect_visible(L(x,2:1))
  expect_visible(L(x,0))
  expect_error(L(x,f))   # common source of error probably is passing the factor in a wrong slot
  expect_error(L(x,c(1,1),f))
  expect_error(L(x,c(1,1),f,t))
  expect_visible(L(x,2:1,f))
  expect_visible(L(x,2:1,f,t))
  expect_visible(L(x,0,f))
  expect_visible(L(x,0,f,t))
  expect_error(L(x,1,1)) # wrong inputs: passing a non-existent difference argument..
  expect_error(L(x,1,1))
  expect_error(L(x,1,1,f))
  expect_error(L(x,1,1,f,t))
  expect_error(L(x,1,-1,f))
  expect_error(L(x,1,-1,f,t))
  # repeated values or gaps in time-variable should give error
  expect_error(L(1:3, t = c(1,1,2)))
  expect_error(L(1:3, t = c(1,2,2)))
  expect_error(L(1:3, t = c(1,2,1)))
  expect_error(L(1:10, g = c(1,1,1,2,2,2,3,3,3,3), t = c(1,1,2,1:3,1:4))) # This is the only possible statement which does not throw a reteated timevar error because the first C++ index is 0, and omap is also initialized with 0's.
  expect_error(L(1:10, g = c(1,1,1,2,2,2,3,3,3,3), t = c(1,1,1,1:3,1:4)))
  expect_error(L(1:10, g = c(1,1,1,2,2,2,3,3,3,3), t = c(1:3,1:3,1,1,3,4)))
  expect_error(L(1:10, g = c(1,1,1,2,2,2,3,3,3,3), t = c(1,2,2,1:3,1:4)))
  expect_visible(L(1:10, g = c(1,1,1,2,2,2,3,3,3,3), t = c(1,2,4,1:3,1:4)))
  expect_error(L(1:10, g = c(1,2,1,2,2,2,3,3,3,3), t = c(1:3,1:3,1:4)))
  expect_visible(L(1:10, g = c(1,1,1,2,2,2,3,3,4,3), t = c(1:3,1:3,1:4)))
  # The usual stuff: Wrongly sized grouping vectors or time-variables
  expect_error(L(1:3, t = 1:2))
  expect_error(L(1:3, t = 1:4))
  expect_error(L(1:3, g = 1:2))
  expect_error(L(1:3, g = 1:4))
  expect_error(L(1:4, g = c(1,1,2,2), t = c(1,2,1)))
  expect_error(L(1:4, g = c(1,2,2), t = c(1,2,1,2)))
})

test_that("L.data.frame method is foolproof", {
  expect_visible(L(wlddev))
  expect_visible(L(wlddev, by = wlddev$iso3c))
  expect_error(L(wlddev, t = ~year))
  expect_visible(L(wlddev, 1, wlddev$iso3c))
  expect_visible(L(wlddev, -2:2, ~iso3c))
  expect_visible(L(wlddev, 1, ~iso3c + region))
  expect_visible(L(wlddev, -2:2, wlddev$iso3c, wlddev$year))
  expect_visible(L(wlddev, -2:2, ~iso3c, ~year))
  expect_visible(L(wlddev, cols = 9:12))
  expect_visible(L(wlddev, -1:1,~iso3c, cols = 9:12))
  expect_visible(L(wlddev, -1:1,wlddev$iso3c, wlddev$year, cols = 9:12))
  expect_visible(L(wlddev, -1:1,~iso3c, ~year, cols = 9:12))
  expect_visible(L(wlddev, cols = c("PCGDP","LIFEEX")))
  expect_visible(L(wlddev, -1:1,wlddev$iso3c, cols = c("PCGDP","LIFEEX")))
  expect_visible(L(wlddev, -1:1,~iso3c, cols = c("PCGDP","LIFEEX")))
  expect_visible(L(wlddev, -1:1,wlddev$iso3c, wlddev$year, cols = c("PCGDP","LIFEEX")))
  expect_visible(L(wlddev, -1:1,~iso3c, ~year, cols = c("PCGDP","LIFEEX")))

  expect_visible(L(wlddev, cols = NULL))
  expect_visible(L(wlddev, -1:1,wlddev$iso3c, cols = NULL))
  expect_visible(L(wlddev, -1:1,~iso3c, cols = NULL))
  expect_visible(L(wlddev, -1:1,wlddev$iso3c, wlddev$year, cols = NULL))
  expect_visible(L(wlddev, -1:1,~iso3c, ~year, cols = NULL))
  expect_error(L(wlddev, cols = 9:14))
  expect_error(L(wlddev, -1:1,~iso3c, ~year, cols = 9:14))
  expect_error(L(wlddev, cols = c("PCGDP","LIFEEX","bla")))
  expect_error(L(wlddev, -1:1,~iso3c, ~year, cols = c("PCGDP","LIFEEX","bla")))

  expect_warning(L(wlddev, w = 4))
  expect_warning(L(wlddev, g = 4))
  expect_error(L(wlddev, t = "year"))
  expect_error(L(wlddev, by = ~year2))
  expect_error(L(wlddev, t = ~year + region))
  expect_error(L(wlddev, data))
  expect_error(L(wlddev, -1:1,"iso3c"))
  expect_error(L(wlddev, -1:1,~iso3c2))
  expect_error(L(wlddev, -1:1,~iso3c + bla))
  expect_error(L(wlddev, -1:1,t = rnorm(30)))
  expect_error(L(wlddev, -1:1,by = rnorm(30)))
  expect_error(L(wlddev, -1:1,mtcars$mpg, 1:29))
  expect_error(L(wlddev, -1:1,mtcars$mpg, mtcars$cyl)) # this gives a repeated values error first because length(g) == length(t)
  expect_error(L(wlddev,-1:1, ~iso3c2, ~year2))
  expect_error(L(wlddev, cols = ~bla))
  expect_visible(L(wlddev, -1:1,wlddev$iso3c, ~year, cols = 9:12))
  expect_visible(L(wlddev, -1:1,~iso3c, wlddev$year, cols = 9:12))
  expect_error(L(wlddev, -1:1,wlddev$iso3c, ~year + bla, cols = 9:12))
  expect_error(L(wlddev, -1:1,~iso3c3, ~year, cols = 9:12))
  expect_error(L(wlddev, cols = c("PC3GDP","LIFEEX")))
})




