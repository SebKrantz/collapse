context("setop")

if(!is.null(attributes(identical(FALSE, TRUE)))) stop("OECD label issue")

d <- mtcars$mpg
dc <- copyv(d, 0, 0)
i <- as.integer(mtcars$cyl)
ic <- copyv(i, 0, 0)
dm <- as.matrix(mtcars) + 1
dmc <- copyv(dm, 0, 0)
im <- dm
storage.mode(im) <- "integer"
imc <- copyv(im, 0, 0)
dr <- dm[nrow(dm), ]
ir <- im[nrow(im), ]
ddf <- mtcars %c+% 1
idf <- dapply(ddf, as.integer)
ddfc <- copyv(ddf, 0, 0)
idfc <- copyv(idf, 0, 0)
ops <- c("+", "-", "*", "/")


test_that("setop works in scalar-vector operations", {
  expect_equal(i %+=% 2 %-=% 2, ic)
  expect_equal(i %+=% 2L %-=% 2L, ic)
  expect_equal(i %*=% 2 %/=% 2, ic)
  expect_equal(i %*=% 2L %/=% 2L, ic)
  expect_equal(d %+=% 2 %-=% 2, dc)
  expect_equal(d %+=% 2L %-=% 2L, dc)
  expect_equal(d %*=% 2 %/=% 2, dc)
  expect_equal(d %*=% 2L %/=% 2L, dc)
  expect_equal(i %+=% dc %-=% trunc(dc), ic) # Problem: The computation creates a decimal which is then rounded down...
  expect_equal(i %+=% ic %-=% ic, ic)
  expect_equal(i %*=% dc %/=% trunc(dc), ic)
  expect_equal(i %*=% ic %/=% ic, ic)
  expect_equal(d %+=% dc %-=% dc, dc)
  expect_equal(d %+=% ic %-=% ic, dc)
  expect_equal(d %*=% dc %/=% dc, dc)
  expect_equal(d %*=% ic %/=% ic, dc)
  expect_identical(i, ic)
  expect_equal(d, dc)

  # Same with setop function
  for(o in ops) setop(i, o, 2); expect_identical(i, ic)
  for(o in ops) setop(d, o, 2); expect_equal(d, dc)
  for(o in ops) setop(i, o, 2L); expect_identical(i, ic)
  for(o in ops) setop(d, o, 2L); expect_equal(d, dc)
  for(o in ops) setop(i, o, trunc(dc)); expect_identical(i, ic)
  for(o in ops) setop(d, o, dc); expect_equal(d, dc)
  for(o in ops) setop(i, o, ic); expect_identical(i, ic)
  for(o in ops) setop(d, o, ic); expect_equal(d, dc)

})

test_that("setop works in scalar-vector-matrix operations", {
  # Matrix & Scalar
  expect_equal(im %+=% 2 %-=% 2, imc)
  expect_equal(im %+=% 2L %-=% 2L, imc)
  expect_equal(im %*=% 2 %/=% 2, imc)
  expect_equal(im %*=% 2L %/=% 2L, imc)
  expect_equal(dm %+=% 2 %-=% 2, dmc)
  expect_equal(dm %+=% 2L %-=% 2L, dmc)
  expect_equal(dm %*=% 2 %/=% 2, dmc)
  expect_equal(dm %*=% 2L %/=% 2L, dmc)
  # Matrix & Vector
  expect_equal(im %+=% trunc(dc) %-=% trunc(dc), imc)
  expect_equal(im %+=% ic %-=% ic, imc)
  expect_equal(im %*=% trunc(dc) %/=% trunc(dc), imc)
  expect_equal(im %*=% ic %/=% ic, imc)
  expect_equal(dm %+=% dc %-=% dc, dmc)
  expect_equal(dm %+=% ic %-=% ic, dmc)
  expect_equal(dm %*=% dc %/=% dc, dmc)
  expect_equal(dm %*=% ic %/=% ic, dmc)
  # Matrix & Matrix
  expect_equal(im %+=% trunc(dmc) %-=% trunc(dmc), imc)
  expect_equal(im %+=% imc %-=% imc, imc)
  expect_equal(im %*=% trunc(dmc) %/=% trunc(dmc), imc)
  expect_equal(im %*=% imc %/=% imc, imc)
  expect_equal(dm %+=% dmc %-=% dmc, dmc)
  expect_equal(dm %+=% imc %-=% imc, dmc)
  expect_equal(dm %*=% dmc %/=% dmc, dmc)
  expect_equal(dm %*=% imc %/=% imc, dmc)

  expect_identical(im, imc)
  expect_equal(dm, dmc)

  # Same with setop function
  # Matrix & Scalar
  for(o in ops) setop(im, o, 2); expect_identical(im, imc)
  for(o in ops) setop(dm, o, 2); expect_equal(dm, dmc)
  for(o in ops) setop(im, o, 2L); expect_identical(im, imc)
  for(o in ops) setop(dm, o, 2L); expect_equal(dm, dmc)
  # Matrix & Vector
  for(o in ops) setop(im, o, trunc(dc)); expect_identical(im, imc)
  for(o in ops) setop(dm, o, dc); expect_equal(dm, dmc)
  for(o in ops) setop(im, o, ic); expect_identical(im, imc)
  for(o in ops) setop(dm, o, ic); expect_equal(dm, dmc)
  # Matrix & Matrix
  for(o in ops) setop(im, o, trunc(dmc)); expect_identical(im, imc)
  for(o in ops) setop(dm, o, dmc); expect_equal(dm, dmc)
  for(o in ops) setop(im, o, imc); expect_identical(im, imc)
  for(o in ops) setop(dm, o, imc); expect_equal(dm, dmc)
  # Row-wise Matrix & Vector
  for(o in ops) setop(im, o, trunc(dr), rowwise = TRUE); expect_identical(im, imc)
  for(o in ops) setop(dm, o, dr, rowwise = TRUE); expect_equal(dm, dmc)
  for(o in ops) setop(im, o, ir, rowwise = TRUE); expect_identical(im, imc)
  for(o in ops) setop(dm, o, ir, rowwise = TRUE); expect_equal(dm, dmc)
  # Comparison with TRA (only for doubles)
  if(requireNamespace("data.table", quietly = TRUE)) {
  for(o in ops) {
    expect_equal(setop(dm, o, dr, rowwise = TRUE), TRA(dmc, dr, o))
    dm <- data.table::copy(dmc)
    expect_equal(setop(dm, o, ir, rowwise = TRUE), TRA(dmc, ir, o))
    dm <- data.table::copy(dmc)
  }
  }
})

test_that("setop works in operations involving data frames", {
  # DF & Scalar
  expect_equal(idf %+=% 2 %-=% 2, idfc)
  expect_equal(idf %+=% 2L %-=% 2L, idfc)
  expect_equal(idf %*=% 2 %/=% 2, idfc)
  expect_equal(idf %*=% 2L %/=% 2L, idfc)
  expect_equal(ddf %+=% 2 %-=% 2, ddfc)
  expect_equal(ddf %+=% 2L %-=% 2L, ddfc)
  expect_equal(ddf %*=% 2 %/=% 2, ddfc)
  expect_equal(ddf %*=% 2L %/=% 2L, ddfc)
  # DF & Vector
  expect_equal(idf %+=% trunc(dc) %-=% trunc(dc), idfc)
  expect_equal(idf %+=% ic %-=% ic, idfc)
  expect_equal(idf %*=% trunc(dc) %/=% trunc(dc), idfc)
  expect_equal(idf %*=% ic %/=% ic, idfc)
  expect_equal(ddf %+=% dc %-=% dc, ddfc)
  expect_equal(ddf %+=% ic %-=% ic, ddfc)
  expect_equal(ddf %*=% dc %/=% dc, ddfc)
  expect_equal(ddf %*=% ic %/=% ic, ddfc)
  # DF & DF
  expect_equal(idf %+=% trunc(ddfc) %-=% trunc(ddfc), idfc)
  expect_equal(idf %+=% idfc %-=% idfc, idfc)
  expect_equal(idf %*=% trunc(ddfc) %/=% trunc(ddfc), idfc)
  expect_equal(idf %*=% idfc %/=% idfc, idfc)
  expect_equal(ddf %+=% ddfc %-=% ddfc, ddfc)
  expect_equal(ddf %+=% idfc %-=% idfc, ddfc)
  expect_equal(ddf %*=% ddfc %/=% ddfc, ddfc)
  expect_equal(ddf %*=% idfc %/=% idfc, ddfc)

  expect_identical(idf, idfc)
  expect_equal(ddf, ddfc)

  # Same with setop function
  # DF & Scalar
  for(o in ops) setop(idf, o, 2); expect_identical(idf, idfc)
  for(o in ops) setop(ddf, o, 2); expect_equal(ddf, ddfc)
  for(o in ops) setop(idf, o, 2L); expect_identical(idf, idfc)
  for(o in ops) setop(ddf, o, 2L); expect_equal(ddf, ddfc)
  # DF & Vector
  for(o in ops) setop(idf, o, trunc(dc)); expect_identical(idf, idfc)
  for(o in ops) setop(ddf, o, dc); expect_equal(ddf, ddfc)
  for(o in ops) setop(idf, o, ic); expect_identical(idf, idfc)
  for(o in ops) setop(ddf, o, ic); expect_equal(ddf, ddfc)
  # DF & DF
  for(o in ops) setop(idf, o, trunc(ddfc)); expect_identical(idf, idfc)
  for(o in ops) setop(ddf, o, ddfc); expect_equal(ddf, ddfc)
  for(o in ops) setop(idf, o, idfc); expect_identical(idf, idfc)
  for(o in ops) setop(ddf, o, idfc); expect_equal(ddf, ddfc)
  # Row-wise DF & Vector
  for(o in ops) setop(idf, o, trunc(dr), rowwise = TRUE); expect_identical(idf, idfc)
  for(o in ops) setop(ddf, o, dr, rowwise = TRUE); expect_equal(ddf, ddfc)
  for(o in ops) setop(idf, o, ir, rowwise = TRUE); expect_identical(idf, idfc)
  for(o in ops) setop(ddf, o, ir, rowwise = TRUE); expect_equal(ddf, ddfc)
  # Comparison with TRA (only for doubles)
  if(requireNamespace("data.table", quietly = TRUE)) {
  for(o in ops) {
    expect_equal(setop(ddf, o, dr, rowwise = TRUE), TRA(ddfc, dr, o))
    ddf <- data.table::copy(ddfc)
    expect_equal(setop(ddf, o, ir, rowwise = TRUE), TRA(ddfc, ir, o))
    ddf <- data.table::copy(ddfc)
  }
  }
})
