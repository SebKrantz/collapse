context("fcumsum")

# rm(liso = ls())
set.seed(101)
x <- abs(10*rnorm(100))
xNA <- x
xNA[sample.int(100, 20)] <- NA
xNA[1L] <- NA
f <- as.factor(rep(1:10, each = 10))
t <- as.factor(rep(1:100))

data <- wlddev[wlddev$iso3c %in% c("BLZ","IND","USA","SRB","GRL"), ]
g <- GRP(droplevels(data$iso3c))
td <- as.factor(data$year)
dataNA <- na_insert(data)
m <- as.matrix(data)
suppressWarnings(storage.mode(m) <- "numeric")
mNAc <- as.matrix(dataNA)
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
mNAuo <- mNA[od, ]
datauo = data[od, ]
dataNAuo = dataNA[od, ]
guo = as_factor_GRP(g)[od]
tduo = td[od]
t2duo = seq_along(od)[od]
od = order(od)

basecumsum <- function(x, na.rm = TRUE, fill = FALSE) {
  ax <- attributes(x)
  if(!na.rm || !anyNA(x)) return(`attributes<-`(cumsum(x), ax))
  cc <- which(!is.na(x))
  x[cc] <- cumsum(x[cc])
  if(!fill) return(x)
  if(is.na(x[1L])) x[1L] <- 0L
  data.table::nafill(x, type = "locf")
}


test_that("fcumsum performs like basecumsum", {
  # No groups, no ordering
  expect_equal(fcumsum(-10:10), basecumsum(-10:10))
  expect_equal(fcumsum(-10:10, na.rm = FALSE), basecumsum(-10:10, na.rm = FALSE))
  expect_equal(fcumsum(-10:10, fill = TRUE), basecumsum(-10:10, fill = TRUE))
  expect_equal(fcumsum(x), basecumsum(x))
  expect_equal(fcumsum(x, na.rm = FALSE), basecumsum(x, na.rm = FALSE))
  expect_equal(fcumsum(x, fill = TRUE), basecumsum(x, fill = TRUE))
  expect_equal(fcumsum(xNA), basecumsum(xNA))
  expect_equal(fcumsum(xNA, na.rm = FALSE), basecumsum(xNA, na.rm = FALSE))
  expect_equal(fcumsum(xNA, fill = TRUE), basecumsum(xNA, fill = TRUE))
  expect_equal(fcumsum(m), dapply(m, basecumsum))
  expect_equal(fcumsum(m, na.rm = FALSE), dapply(m, basecumsum, na.rm = FALSE))
  expect_equal(fcumsum(m, fill = TRUE), dapply(m, basecumsum, fill = TRUE))
  expect_equal(fcumsum(mNA), dapply(mNA, basecumsum))
  expect_equal(fcumsum(mNA, na.rm = FALSE), dapply(mNA, basecumsum, na.rm = FALSE))
  expect_equal(fcumsum(mNA, fill = TRUE), dapply(mNA, basecumsum, fill = TRUE))
  expect_equal(fcumsum(num_vars(data)), dapply(num_vars(data), basecumsum))
  expect_equal(fcumsum(num_vars(data), na.rm = FALSE), dapply(num_vars(data), basecumsum, na.rm = FALSE))
  expect_equal(fcumsum(num_vars(data), fill = TRUE), dapply(num_vars(data), basecumsum, fill = TRUE))
  expect_equal(fcumsum(num_vars(dataNA)), dapply(num_vars(dataNA), basecumsum))
  #expect_equal(fcumsum(num_vars(dataNA), na.rm = FALSE), dapply(num_vars(dataNA), basecumsum, na.rm = FALSE))
  expect_equal(fcumsum(num_vars(dataNA), fill = TRUE), dapply(num_vars(dataNA), basecumsum, fill = TRUE))
  # With groups, no ordering
  expect_equal(fcumsum(x, f), BY(x, f, basecumsum, use.g.names = FALSE))
  expect_equal(fcumsum(x, na.rm = FALSE, f), BY(x, f, basecumsum, na.rm = FALSE, use.g.names = FALSE))
  expect_equal(fcumsum(x, f, fill = TRUE), BY(x, f, basecumsum, fill = TRUE, use.g.names = FALSE))
  expect_equal(fcumsum(xNA, f), BY(xNA, f, basecumsum, use.g.names = FALSE))
  expect_equal(fcumsum(xNA, na.rm = FALSE, f), BY(xNA, f, basecumsum, na.rm = FALSE, use.g.names = FALSE))
  expect_equal(fcumsum(xNA, f, fill = TRUE), BY(xNA, f, basecumsum, fill = TRUE, use.g.names = FALSE))
  expect_equal(fcumsum(m, g), BY(m, g, basecumsum, use.g.names = FALSE))
  expect_equal(fcumsum(m, na.rm = FALSE, g), BY(m, g, basecumsum, na.rm = FALSE, use.g.names = FALSE))
  expect_equal(fcumsum(m, g, fill = TRUE), BY(m, g, basecumsum, fill = TRUE, use.g.names = FALSE))
  expect_equal(fcumsum(mNA, g), BY(mNA, g, basecumsum, use.g.names = FALSE))
  expect_equal(fcumsum(mNA, na.rm = FALSE, g), BY(mNA, g, basecumsum, na.rm = FALSE, use.g.names = FALSE))
  expect_equal(fcumsum(mNA, g, fill = TRUE), BY(mNA, g, basecumsum, fill = TRUE, use.g.names = FALSE))
  expect_equal(fcumsum(num_vars(data), g), BY(num_vars(data), g, basecumsum, use.g.names = FALSE))
  expect_equal(fcumsum(num_vars(data), na.rm = FALSE, g), BY(num_vars(data), g, basecumsum, na.rm = FALSE, use.g.names = FALSE))
  expect_equal(fcumsum(num_vars(data), g, fill = TRUE), BY(num_vars(data), g, basecumsum, fill = TRUE, use.g.names = FALSE))
  expect_equal(fcumsum(num_vars(dataNA), g), BY(num_vars(dataNA), g, basecumsum, use.g.names = FALSE))
  # expect_equal(fcumsum(num_vars(dataNA), na.rm = FALSE, g), BY(num_vars(dataNA), g, basecumsum, na.rm = FALSE, use.g.names = FALSE))
  expect_equal(fcumsum(num_vars(dataNA), g, fill = TRUE), BY(num_vars(dataNA), g, basecumsum, fill = TRUE, use.g.names = FALSE))
})


test_that("fcumsum correctly handles unordered time-series and panel-series computations", {
  # With ordering, no groups: 1
  expect_equal(fcumsum(x, o = 1:100), fcumsum(x))
  expect_equal(fcumsum(x, o = 1:100, na.rm = FALSE), fcumsum(x, na.rm = FALSE))
  expect_equal(fcumsum(x, o = 1:100, fill = TRUE), fcumsum(x, fill = TRUE))
  expect_equal(fcumsum(xNA, o = 1:100), fcumsum(xNA))
  expect_equal(fcumsum(xNA, o = 1:100, na.rm = FALSE), fcumsum(xNA, na.rm = FALSE))
  expect_equal(fcumsum(xNA, o = 1:100, fill = TRUE), fcumsum(xNA, fill = TRUE))
  expect_equal(fcumsum(m, o = seq_row(m)), fcumsum(m))
  expect_equal(fcumsum(m, o = seq_row(m), na.rm = FALSE), fcumsum(m, na.rm = FALSE))
  expect_equal(fcumsum(m, o = seq_row(m), fill = TRUE), fcumsum(m, fill = TRUE))
  expect_equal(fcumsum(mNA, o = seq_row(m)), fcumsum(mNA))
  expect_equal(fcumsum(mNA, o = seq_row(m), na.rm = FALSE), fcumsum(mNA, na.rm = FALSE))
  expect_equal(fcumsum(mNA, o = seq_row(m), fill = TRUE), fcumsum(mNA, fill = TRUE))
  expect_equal(fcumsum(num_vars(data), o = seq_row(data)), fcumsum(num_vars(data)))
  expect_equal(fcumsum(num_vars(data), o = seq_row(data), na.rm = FALSE), fcumsum(num_vars(data), na.rm = FALSE))
  expect_equal(fcumsum(num_vars(data), o = seq_row(data), fill = TRUE), fcumsum(num_vars(data), fill = TRUE))
  expect_equal(fcumsum(num_vars(dataNA), o = seq_row(data)), fcumsum(num_vars(dataNA)))
  #expect_equal(fcumsum(num_vars(dataNA), o = seq_row(data), na.rm = FALSE), fcumsum(num_vars(dataNA), na.rm = FALSE))
  expect_equal(fcumsum(num_vars(dataNA), o = seq_row(data), fill = TRUE), fcumsum(num_vars(dataNA), fill = TRUE))
  # With ordering, no groups: 2
  expect_equal(fcumsum(xuo, o = t2uo)[o], fcumsum(x))
  expect_equal(fcumsum(xuo, o = t2uo, na.rm = FALSE)[o], fcumsum(x, na.rm = FALSE))
  expect_equal(fcumsum(xuo, o = t2uo, fill = TRUE)[o], fcumsum(x, fill = TRUE))
  expect_equal(fcumsum(xNAuo, o = t2uo)[o], fcumsum(xNA))
  expect_equal(fcumsum(xNAuo, o = t2uo, na.rm = FALSE)[o], fcumsum(xNA, na.rm = FALSE))
  expect_equal(fcumsum(xNAuo, o = t2uo, fill = TRUE)[o], fcumsum(xNA, fill = TRUE))
  expect_equal(fcumsum(muo, o = t2duo)[od, ], fcumsum(m))
  expect_equal(fcumsum(muo, o = t2duo, na.rm = FALSE)[od, ], fcumsum(m, na.rm = FALSE))
  expect_equal(fcumsum(muo, o = t2duo, fill = TRUE)[od, ], fcumsum(m, fill = TRUE))
  expect_equal(fcumsum(mNAuo, o = t2duo)[od, ], fcumsum(mNA))
  expect_equal(fcumsum(mNAuo, o = t2duo, na.rm = FALSE)[od, ], fcumsum(mNA, na.rm = FALSE))
  expect_equal(fcumsum(mNAuo, o = t2duo, fill = TRUE)[od, ], fcumsum(mNA, fill = TRUE))
  expect_equal(fcumsum(num_vars(datauo), o = t2duo)[od, ], fcumsum(num_vars(data)))
  expect_equal(fcumsum(num_vars(datauo), o = t2duo, na.rm = FALSE)[od, ], fcumsum(num_vars(data), na.rm = FALSE))
  expect_equal(fcumsum(num_vars(datauo), o = t2duo, fill = TRUE)[od, ], fcumsum(num_vars(data), fill = TRUE))
  expect_equal(fcumsum(num_vars(dataNAuo), o = t2duo)[od, ], fcumsum(num_vars(dataNA)))
  #expect_equal(fcumsum(num_vars(dataNAuo), o = t2duo, na.rm = FALSE)[od, ], fcumsum(num_vars(dataNA), na.rm = FALSE))
  expect_equal(fcumsum(num_vars(dataNAuo), o = t2duo, fill = TRUE)[od, ], fcumsum(num_vars(dataNA), fill = TRUE))
  # With ordering and groups
  expect_equal(fcumsum(xuo, fuo, tuo)[o], fcumsum(x, f, t))
  expect_equal(fcumsum(xuo, fuo, tuo, na.rm = FALSE)[o], fcumsum(x, f, t, na.rm = FALSE))
  expect_equal(fcumsum(xuo, fuo, tuo, fill = TRUE)[o], fcumsum(x, f, t, fill = TRUE))
  expect_equal(fcumsum(xNAuo, fuo, tuo)[o], fcumsum(xNA, f, t))
  # expect_equal(fcumsum(xNAuo, fuo, tuo, na.rm = FALSE)[o], fcumsum(xNA, fuo, tuo, na.rm = FALSE))
  expect_equal(fcumsum(xNAuo, fuo, tuo, fill = TRUE)[o], fcumsum(xNA, f, t, fill = TRUE))
  expect_equal(fcumsum(muo, guo, tduo)[od, ], fcumsum(m, g, td))
  expect_equal(fcumsum(muo, guo, tduo, na.rm = FALSE)[od, ], fcumsum(m, g, td, na.rm = FALSE))
  expect_equal(fcumsum(muo, guo, tduo, fill = TRUE)[od, ], fcumsum(m, g, td, fill = TRUE))
  expect_equal(fcumsum(mNAuo, guo, tduo)[od, ], fcumsum(mNA, g, td))
  expect_equal(fcumsum(mNAuo, guo, tduo, na.rm = FALSE)[od, ], fcumsum(mNA, g, td, na.rm = FALSE))
  expect_equal(fcumsum(mNAuo, guo, tduo, fill = TRUE)[od, ], fcumsum(mNA, g, td, fill = TRUE))
  expect_equal(fcumsum(num_vars(datauo), guo, tduo)[od, ], fcumsum(num_vars(data), g, td))
  expect_equal(fcumsum(num_vars(datauo), guo, tduo, na.rm = FALSE)[od, ], fcumsum(num_vars(data), g, td, na.rm = FALSE))
  expect_equal(fcumsum(num_vars(datauo), guo, tduo, fill = TRUE)[od, ], fcumsum(num_vars(data), g, td, fill = TRUE))
  expect_equal(fcumsum(num_vars(dataNAuo), guo, tduo)[od, ], fcumsum(num_vars(dataNA), g, td))
  expect_equal(fcumsum(num_vars(dataNAuo), guo, tduo, na.rm = FALSE)[od, ], fcumsum(num_vars(dataNA), g, td, na.rm = FALSE))
  expect_equal(fcumsum(num_vars(dataNAuo), guo, tduo, fill = TRUE)[od, ], fcumsum(num_vars(dataNA), g, td, fill = TRUE))
})

test_that("fcumsum performs numerically stable in ordered computations", {
  expect_true(all_obj_equal(replicate(50, fcumsum(x), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(xNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(m), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(mNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(nv(data)), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(nv(dataNA)), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(x, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(x, f, fill = TRUE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(xNA, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(xNA, f, fill = TRUE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(m, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(m, g, fill = TRUE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(mNA, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(mNA, g, fill = TRUE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(nv(data), g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(nv(data), g, fill = TRUE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(nv(dataNA), g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(nv(dataNA), g, fill = TRUE), simplify = FALSE)))
})

test_that("fcumsum performs numerically stable in unordered computations", {
  expect_true(all_obj_equal(replicate(50, fcumsum(xuo, o = t2uo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(xNAuo, o = t2uo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(muo, o = t2duo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(nv(datauo), o = t2duo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(xuo, fuo, tuo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(xuo, fuo, tuo, fill = TRUE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(muo, guo, tduo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(muo, guo, tduo, fill = TRUE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(nv(datauo), guo, tduo), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fcumsum(nv(datauo), guo, tduo, fill = TRUE), simplify = FALSE)))
})

# test_that("fcumsum handles special values in the right way", {
#   # zero
#   expect_equal(fcumsum(c("a","b"),0), c("a","b"))
#   expect_equal(fcumsum(c(NaN,NaN),0), c(NaN,NaN))
#   expect_equal(fcumsum(c(Inf,Inf),0), c(Inf,Inf))
#   expect_equal(fcumsum(c(FALSE,TRUE),0), c(FALSE,TRUE))
#   expect_equal(fcumsum(c(TRUE,FALSE),0), c(TRUE,FALSE))
#   # cumsums
#   expect_equal(fcumsum(c("a","b")), c(NA,"a"))
#   expect_equal(fcumsum(c(1,NA)), c(NA_real_,1))
#   expect_equal(fcumsum(c(NA,1)), c(NA_real_,NA_real_))
#   expect_equal(fcumsum(c(NaN,1)), c(NA_real_,NaN))
#   expect_equal(fcumsum(c(1,NaN)), c(NA_real_,1))
#   expect_equal(fcumsum(c(Inf,1)), c(NA,Inf))
#   expect_equal(fcumsum(c(1,Inf)), c(NA,1))
#   expect_equal(fcumsum(c(Inf,NA)), c(NA_real_,Inf))
#   expect_equal(fcumsum(c(NA,Inf)), c(NA_real_,NA_real_))
#   expect_equal(fcumsum(c(Inf,-Inf)), c(NA,Inf))
#   expect_equal(fcumsum(c(-Inf,Inf)), c(NA,-Inf))
#   expect_equal(fcumsum(c(Inf,Inf)), c(NA,Inf))
#   expect_equal(fcumsum(c(TRUE,TRUE)), c(NA,TRUE))
#   expect_equal(fcumsum(c(TRUE,FALSE)), c(NA,TRUE))
#   expect_equal(fcumsum(c(FALSE,TRUE)), c(NA,FALSE))
#   # leads
#   expect_equal(fcumsum(c("a","b"),-1), c("b",NA))
#   expect_equal(fcumsum(c(1,NA),-1), c(NA_real_,NA_real_))
#   expect_equal(fcumsum(c(NA,1),-1), c(1,NA_real_))
#   expect_equal(fcumsum(c(NaN,1),-1), c(1,NA_real_))
#   expect_equal(fcumsum(c(1,NaN),-1), c(NaN,NA_real_))
#   expect_equal(fcumsum(c(Inf,1),-1), c(1,NA_real_))
#   expect_equal(fcumsum(c(1,Inf),-1), c(Inf,NA_real_))
#   expect_equal(fcumsum(c(Inf,NA),-1), c(NA_real_,NA_real_))
#   expect_equal(fcumsum(c(NA,Inf),-1), c(Inf,NA_real_))
#   expect_equal(fcumsum(c(Inf,-Inf),-1), c(-Inf,NA_real_))
#   expect_equal(fcumsum(c(-Inf,Inf),-1), c(Inf,NA_real_))
#   expect_equal(fcumsum(c(Inf,Inf),-1), c(Inf,NA_real_))
#   expect_equal(fcumsum(c(TRUE,TRUE),-1), c(TRUE,NA))
#   expect_equal(fcumsum(c(TRUE,FALSE),-1), c(FALSE,NA))
#   expect_equal(fcumsum(c(FALSE,TRUE),-1), c(TRUE,NA))
#
# })
#
# test_that("fcumsum produces errors for wrong input", {
#   # type: normally guaranteed by C++
#   expect_visible(fcumsum(mNAc))
#   expect_visible(fcumsum(wlddev))
#   expect_error(fcumsum(mNAc, f))
#   expect_error(fcumsum(x, "1"))
#   # if n exceeds length(x), should give error
#   expect_error(fcumsum(x,101))
#   expect_error(fcumsum(x,-101))
#   # if n exceeds average group size, should give error
#   # expect_warning(fcumsum(x,11,f)) # Some fail on i386 ??
#   # expect_warning(fcumsum(x,11,f,t))
#   # expect_warning(fcumsum(x,-11,f))
#   # expect_warning(fcumsum(x,-11,f,t))
#   # passing repeated n-values or non-positive or non-consecutive diff values should give error
#   expect_error(fcumsum(x,c(1,1)))
#   expect_error(fcumsum(x,c(-1,-1)))
#   expect_visible(fcumsum(x,2:1))
#   expect_visible(fcumsum(x,0))
#   expect_error(fcumsum(x,f))   # common source of error probably is passing the factor in a wrong slot
#   expect_error(fcumsum(x,c(1,1),f))
#   expect_error(fcumsum(x,c(1,1),f,t))
#   expect_visible(fcumsum(x,2:1,f))
#   expect_visible(fcumsum(x,2:1,f,t))
#   expect_visible(fcumsum(x,0,f))
#   expect_visible(fcumsum(x,0,f,t))
#   expect_error(fcumsum(x,1,1)) # wrong inputs: passing a non-existent difference argument..
#   expect_error(fcumsum(x,1,1))
#   expect_error(fcumsum(x,1,1,f))
#   expect_error(fcumsum(x,1,1,f,t))
#   expect_error(fcumsum(x,1,-1,f))
#   expect_error(fcumsum(x,1,-1,f,t))
#   # repeated values or gaps in time-variable should give error
#   expect_error(fcumsum(1:3, o = c(1,1,2)))
#   expect_error(fcumsum(1:3, o = c(1,2,2)))
#   expect_error(fcumsum(1:3, o = c(1,2,1)))
#   expect_error(fcumsum(1:10, g = c(1,1,1,2,2,2,3,3,3,3), o = c(1,1,2,1:3,1:4))) # This is the only possible statement which does not throw a reteated timevar error because the first C++ index is 0, and omap is also initialized with 0's.
#   expect_error(fcumsum(1:10, g = c(1,1,1,2,2,2,3,3,3,3), o = c(1,1,1,1:3,1:4)))
#   expect_error(fcumsum(1:10, g = c(1,1,1,2,2,2,3,3,3,3), o = c(1:3,1:3,1,1,3,4)))
#   expect_error(fcumsum(1:10, g = c(1,1,1,2,2,2,3,3,3,3), o = c(1,2,2,1:3,1:4)))
#   expect_visible(fcumsum(1:10, g = c(1,1,1,2,2,2,3,3,3,3), o = c(1,2,4,1:3,1:4)))
#   expect_error(fcumsum(1:10, g = c(1,2,1,2,2,2,3,3,3,3), o = c(1:3,1:3,1:4)))
#   expect_visible(fcumsum(1:10, g = c(1,1,1,2,2,2,3,3,4,3), o = c(1:3,1:3,1:4)))
#   # The usual stuff: Wrongly sized grouping vectors or time-variables
#   expect_error(fcumsum(1:3, o = 1:2))
#   expect_error(fcumsum(1:3, o = 1:4))
#   expect_error(fcumsum(1:3, g = 1:2))
#   expect_error(fcumsum(1:3, g = 1:4))
#   expect_error(fcumsum(1:4, g = c(1,1,2,2), o = c(1,2,1)))
#   expect_error(fcumsum(1:4, g = c(1,2,2), o = c(1,2,1,2)))
# })






