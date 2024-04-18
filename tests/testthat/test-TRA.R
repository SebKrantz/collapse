context("TRA")



bmean <- base::mean

# rm(list = ls())
set.seed(101)
d <- na_insert(iris[1:4])
v <- d$Sepal.Length
m <- as.matrix(d)
f <- iris$Species

# For sweep
replace_NA <- function(x, y) `[<-`(x, is.na(x), value = y[is.na(x)])
replace <- function(x, y) `[<-`(y, is.na(x), value = NA)          # `[<-`(x, !is.na(x), value = y)
replace_fill <- function(x, y) y                                  # rep(y, length(x))
"%" <- function(x, y) x * (100 / y)
"-%%" <- function(x, y) x - (x %% y)
# "-+" <- function(x, y) x - y + bmean(x, na.rm = TRUE)


test_that("TRA performs like sweep", {
  ops <- c("replace_NA","replace_fill", "replace", "-", "+", "*", "/", "%", "%%", "-%%")
  for(i in ops) {
    expect_equal(drop(sweep(qM(v), 2L, bmean(v, na.rm = TRUE), i)), TRA(v, bmean(v, na.rm = TRUE), i))
    expect_equal(`attributes<-`(sweep(qM(m), 2L, colMeans(m, na.rm = TRUE), i), attributes(m)), TRA(m, colMeans(m, na.rm = TRUE), i))
    expect_equal(setNames(qDF(sweep(d, 2L, colMeans(qM(d), na.rm = TRUE), i)), names(d)), TRA(d, colMeans(qM(d), na.rm = TRUE), i))
  }
  for(i in ops) {
    expect_equal(unlist(Map(function(x, y) drop(sweep(qM(x), 2L, y, i)), rsplit(v, f), as.list(fmean(v, f))), use.names = FALSE),
                 TRA(v, fmean(v, f), i, f))
    expect_equal(`attributes<-`(do.call(rbind, Map(function(x, y) sweep(qM(x), 2L, y, i), lapply(rsplit(qDF(m), f), qM), mrtl(fmean(m, f)))), attributes(m)),
                 TRA(m, fmean(m, f), i, f))
    expect_equal(`attributes<-`(unlist2d(Map(function(x, y) sweep(x, 2L, y, i), rsplit(d, f), mrtl(qM(fmean(d, f)))), idcols = FALSE), attributes(d)),
                 TRA(d, fmean(d, f), i, f))
  }
})

test_that("TRA performs like built-in version", {
  for(i in c(0L, seq_len(10)[-4])) {
   expect_equal(TRA(v, fmean(v), i), fmean(v, TRA = i))
   expect_equal(TRA(m, fmean(m), i), fmean(m, TRA = i))
   expect_equal(TRA(d, fmean(d), i), fmean(d, TRA = i))
  }
 for(i in c(0L, seq_len(10))) {
   expect_equal(TRA(v, fmean(v, f), i, f), fmean(v, f, TRA = i))
   expect_equal(TRA(m, fmean(m, f), i, f), fmean(m, f, TRA = i))
   expect_equal(TRA(d, fmean(d, f), i, f), fmean(d, f, TRA = i))
 }
})

test_that("TRA performs like fbetween and fwithin", {
    expect_equal(TRA(v, fmean(v), 1L), fbetween(v, fill = TRUE))
    expect_equal(TRA(v, fmean(v), 2L), fbetween(v))
    expect_equal(TRA(v, fmean(v), 3L), fwithin(v))
    expect_equal(TRA(m, fmean(m), 1L), fbetween(m, fill = TRUE))
    expect_equal(TRA(m, fmean(m), 2L), fbetween(m))
    expect_equal(TRA(m, fmean(m), 3L), fwithin(m))
    expect_equal(TRA(d, fmean(d), 1L), fbetween(d, fill = TRUE))
    expect_equal(TRA(d, fmean(d), 2L), fbetween(d))
    expect_equal(TRA(d, fmean(d), 3L), fwithin(d))

    expect_equal(TRA(v, fmean(v, f), 1L, f), fbetween(v, f, fill = TRUE))
    expect_equal(TRA(v, fmean(v, f), 2L, f), fbetween(v, f))
    expect_equal(TRA(v, fmean(v, f), 3L, f), fwithin(v, f))
    expect_equal(TRA(v, fmean(v, f), 4L, f), fwithin(v, f, mean = "overall.mean"))
    expect_equal(TRA(m, fmean(m, f), 1L, f), fbetween(m, f, fill = TRUE))
    expect_equal(TRA(m, fmean(m, f), 2L, f), fbetween(m, f))
    expect_equal(TRA(m, fmean(m, f), 3L, f), fwithin(m, f))
    expect_equal(TRA(m, fmean(m, f), 4L, f), fwithin(m, f, mean = "overall.mean"))
    expect_equal(TRA(d, fmean(d, f), 1L, f), fbetween(d, f, fill = TRUE))
    expect_equal(TRA(d, fmean(d, f), 2L, f), fbetween(d, f))
    expect_equal(TRA(d, fmean(d, f), 3L, f), fwithin(d, f))
    expect_equal(TRA(d, fmean(d, f), 4L, f), fwithin(d, f, mean = "overall.mean"))
})

test_that("TRA gives errors for wrong input", {
  expect_warning(TRA(v, fmean(v), bla = 1))
  expect_warning(TRA(m, fmean(m), bla = 1))
  expect_warning(TRA(d, fmean(d), bla = 1))
  expect_error(TRA(v, 1:2))
  expect_error(TRA(m, 1:2))
  expect_error(TRA(d, 1:2))
  expect_error(TRA(v, as.character(fmean(v))))
  expect_error(TRA(m, as.character(fmean(m))))
  expect_error(TRA(d, as.character(fmean(d))))
  expect_error(TRA(v, fmean(v, f), "-", f[-1]))
  expect_error(TRA(m, fmean(m, f), "-", f[-1]))
  expect_error(TRA(d, fmean(d, f), "-", f[-1]))
  expect_error(TRA(v, fmean(v), 19L))
  expect_error(TRA(m, fmean(m), 19L))
  expect_error(TRA(d, fmean(d), 19L))
  expect_error(TRA(v, fmean(v), "bla"))
  expect_error(TRA(m, fmean(m), "bla"))
  expect_error(TRA(d, fmean(d), "bla"))
})


test_that("TRA handles different data types as intended", {
  # Vector & Matrix: Simple
  expect_true(is.character(fnobs(na_insert(letters), TRA = "replace_NA")))
  expect_true(is.integer(fnobs(letters, TRA = "replace_fill")))
  expect_true(is.integer(fnobs(na_insert(letters), TRA = "replace")))
  for(i in c("-", "+", "*", "/", "%", "%%", "-%%"))  expect_error(fnobs(letters, TRA = i))
  expect_true(is.double(fnobs(na_insert(AirPassengers), TRA = "replace_NA")))
  expect_true(is.integer(fnobs(AirPassengers, TRA = "replace_fill")))
  expect_true(is.integer(fnobs(AirPassengers, TRA = "replace")))
  for(i in c("-", "+", "*", "/", "%", "%%", "-%%"))  expect_true(is.numeric(fnobs(AirPassengers, TRA = i)))
  expect_true(is.double(fnobs(na_insert(EuStockMarkets), TRA = "replace_NA")))
  expect_true(is.integer(fnobs(EuStockMarkets, TRA = "replace_fill")))
  expect_true(is.integer(fnobs(EuStockMarkets, TRA = "replace")))
  for(i in c("-", "+", "*", "/", "%", "%%", "-%%"))  expect_true(is.numeric(fnobs(EuStockMarkets, TRA = i)))
  # Vector & Matrix: Grouped
  set.seed(101)
  expect_error(fnobs(letters, sample.int(3, length(letters), TRUE), TRA = "replace_NA"))
  expect_true(is.integer(fnobs(letters, sample.int(3, length(letters), TRUE), TRA = "replace_fill")))
  expect_true(is.integer(fnobs(letters, sample.int(3, length(letters), TRUE), TRA = "replace")))
  for(i in c("-", "-+", "+", "*", "/", "%", "%%", "-%%"))  expect_error(fnobs(letters, sample.int(3, length(letters), TRUE), TRA = i))
  expect_true(is.double(fnobs(AirPassengers, sample.int(3, length(AirPassengers), TRUE), TRA = "replace_NA")))
  expect_true(is.integer(fnobs(AirPassengers, sample.int(3, length(AirPassengers), TRUE), TRA = "replace_fill")))
  expect_true(is.integer(fnobs(AirPassengers, sample.int(3, length(AirPassengers), TRUE), TRA = "replace")))
  for(i in c("-", "-+", "+", "*", "/", "%", "%%", "-%%"))  expect_true(is.numeric(fnobs(AirPassengers, sample.int(3, length(AirPassengers), TRUE), TRA = i)))
  expect_true(is.double(fnobs(EuStockMarkets, sample.int(3, nrow(EuStockMarkets), TRUE), TRA = "replace_NA")))
  expect_true(is.integer(fnobs(EuStockMarkets, sample.int(3, nrow(EuStockMarkets), TRUE), TRA = "replace_fill")))
  expect_true(is.integer(fnobs(EuStockMarkets, sample.int(3, nrow(EuStockMarkets), TRUE), TRA = "replace")))
  for(i in c("-", "-+", "+", "*", "/", "%", "%%", "-%%"))  expect_true(is.numeric(fnobs(EuStockMarkets, sample.int(3, nrow(EuStockMarkets), TRUE),  TRA = i)))

  # Date Frame: Simple
  expect_equal(vtypes(fndistinct(wlddev, TRA = "replace_NA")), vtypes(wlddev))
  expect_equal(unname(vtypes(fndistinct(wlddev, TRA = "replace"))), rep("integer", 13))
  expect_equal(unname(vtypes(fndistinct(wlddev, TRA = "replace_fill"))), rep("integer", 13))
  expect_equal(vtypes(fmode(wlddev, TRA = "replace_NA")), vtypes(wlddev))
  expect_equal(vtypes(ffirst(wlddev, TRA = "replace_NA")), vtypes(wlddev))
  expect_equal(vtypes(flast(wlddev, TRA = "replace_NA")), vtypes(wlddev))
  expect_equal(vtypes(fmode(wlddev, TRA = "replace_fill")), vtypes(wlddev))
  expect_equal(vtypes(ffirst(wlddev, TRA = "replace_fill")), vtypes(wlddev))
  expect_equal(vtypes(flast(wlddev, TRA = "replace_fill")), vtypes(wlddev))
  expect_equal(vtypes(fmode(wlddev, TRA = "replace")), vtypes(wlddev))
  expect_equal(vtypes(ffirst(wlddev, TRA = "replace")), vtypes(wlddev))
  expect_equal(vtypes(flast(wlddev, TRA = "replace")), vtypes(wlddev))
  for(i in c("-", "+", "*", "/", "%", "%%", "-%%")) expect_equal(unname(vtypes(fnobs(nv(wlddev), TRA = i))), rep("double", 7))
  for(i in c("-", "+", "*", "/", "%", "%%", "-%%")) expect_error(fnobs(wlddev, TRA = i))
  # Date Frame: Grouped
  expect_equal(unname(vtypes(fndistinct(wlddev, wlddev$iso3c, TRA = "replace"))), rep("integer", 13))
  expect_equal(unname(vtypes(fndistinct(wlddev, wlddev$iso3c, TRA = "replace_fill"))), rep("integer", 13))
  expect_error(fndistinct(wlddev, wlddev$iso3c, TRA = "replace_NA"))
  expect_equal(vtypes(fmode(wlddev, wlddev$iso3c, TRA = "replace_NA")), vtypes(wlddev))
  expect_equal(vtypes(ffirst(wlddev, wlddev$iso3c, TRA = "replace_NA")), vtypes(wlddev))
  expect_equal(vtypes(flast(wlddev, wlddev$iso3c, TRA = "replace_NA")), vtypes(wlddev))
  expect_equal(vtypes(fmode(wlddev, wlddev$iso3c, TRA = "replace_fill")), vtypes(wlddev))
  expect_equal(vtypes(ffirst(wlddev, wlddev$iso3c, TRA = "replace_fill")), vtypes(wlddev))
  expect_equal(vtypes(flast(wlddev, wlddev$iso3c, TRA = "replace_fill")), vtypes(wlddev))
  expect_equal(vtypes(fmode(wlddev, wlddev$iso3c, TRA = "replace")), vtypes(wlddev))
  expect_equal(vtypes(ffirst(wlddev, wlddev$iso3c, TRA = "replace")), vtypes(wlddev))
  expect_equal(vtypes(flast(wlddev, wlddev$iso3c, TRA = "replace")), vtypes(wlddev))
  for(i in c("-", "-+", "+", "*", "/", "%", "%%", "-%%")) expect_equal(unname(vtypes(fnobs(nv(wlddev), wlddev$iso3c, TRA = i))), rep("double", 7))
  for(i in c("-", "-+", "+", "*", "/", "%", "%%", "-%%")) expect_error(fnobs(wlddev, wlddev$iso3c, TRA = i))

})
