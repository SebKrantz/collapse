context("TRA")

d <- na_insert(iris[1:4])
v <- d$Sepal.Length
m <- as.matrix(d)
f <- iris$Species

test_that("TRA performs like built-in version", {
  for(i in seq_len(8)[-4]) {
   expect_equal(TRA(v, fmean(v), i), fmean(v, TRA = i))
   expect_equal(TRA(m, fmean(m), i), fmean(m, TRA = i))
   expect_equal(TRA(d, fmean(d), i), fmean(d, TRA = i))
  }
 for(i in seq_len(8)) {
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
  expect_error(TRA(v, fmean(v), bla = 1))
  expect_error(TRA(m, fmean(m), bla = 1))
  expect_error(TRA(d, fmean(d), bla = 1))
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
