context("anyv, allv, whichv, setv, copyv etc.")



# d <- replace_NA(wlddev, cols = 9:13)

test_that("whichv works well", {
  expect_identical(whichv(wlddev$country, "Chad"), which(wlddev$country == "Chad"))
  expect_identical(whichv(wlddev$country, "Chad", invert = TRUE), which(wlddev$country != "Chad"))
  expect_identical(whichNA(wlddev$PCGDP), which(is.na(wlddev$PCGDP)))
  expect_identical(whichNA(wlddev$PCGDP, invert = TRUE), which(!is.na(wlddev$PCGDP)))
  expect_identical(whichv(is.na(wlddev$PCGDP), FALSE), which(!is.na(wlddev$PCGDP)))
})


test_that("anyv, allv and whichv work properly", {
  for(i in seq_along(wlddev)) {
    vec <- .subset2(wlddev, i)
    v <- vec[trunc(runif(1L, 1L, length(vec)))]
    if(is.na(v)) v <- flast(vec)
    expect_identical(which(vec == v), whichv(vec, v))
    if(!anyNA(vec)) expect_identical(which(vec != v), whichv(vec, v, TRUE))
    expect_identical(all(vec == v), allv(vec, v))
    expect_identical(any(vec == v), anyv(vec, v))
    vecNA <- is.na(vec)
    expect_identical(which(vecNA), whichNA(vec))
    expect_identical(which(!vecNA), whichNA(vec, TRUE))
    expect_identical(all(vecNA), allNA(vec))
    expect_identical(any(vecNA), anyNA(vec))
  }
  if(identical(Sys.getenv("NCRAN"), "TRUE")) {
  expect_true(allv(rep(0.0004, 1000), 0.0004))
  expect_false(allv(rep(0.0004, 1000), 0.0005))
  }
})

if(requireNamespace("data.table", quietly = TRUE)) {

wldcopy <- data.table::copy(wlddev)
mtccopy <- data.table::copy(mtcars)

test_that("setv and copyv work properly", {
  for(FUN in list(copyv, setv)) {
    for(i in seq_along(wlddev)) {
      # print(i)
      vec <- .subset2(wlddev, i)
      v <- vec[trunc(runif(1L, 1L, length(vec)))]
      r <- vec[trunc(runif(1L, 1L, length(vec)))]
      if(is.na(v)) v <- flast(vec)
      vl <- vec == v
      nvl <- vec != v
      vna <- is.na(vec)
      expect_identical(FUN(vec, v, r), replace(vec, vl, r))
      expect_identical(FUN(vec, which(vl), r, vind1 = TRUE), replace(vec, which(vl), r))
      expect_identical(FUN(vec, 10:1000, r), replace(vec, 10:1000, r))
      expect_identical(FUN(vec, NA, r), replace(vec, vna, r))
      expect_identical(FUN(vec, vl, r), replace(vec, vl, r))
      expect_identical(FUN(vec, 258L, r, vind1 = TRUE), replace(vec, 258L, r))
      expect_identical(FUN(vec, vl, r, invert = TRUE), replace(vec, !vl, r))
      expect_identical(FUN(vec, which(nvl), r), replace(vec, which(nvl), r))
      expect_error(FUN(vec, which(vl), r, invert = TRUE, vind1 = TRUE))
      # expect_error(FUN(vec, which(nvl), r, invert = TRUE))
      if(anyNA(vl)) {
        setv(vl, NA, FALSE)
        setv(nvl, NA, FALSE)
      }
      expect_identical(FUN(vec, v, vec), replace(vec, vl, vec[vl]))
      expect_identical(FUN(vec, NA, vec), replace(vec, vna, vec[vna]))
      expect_identical(FUN(vec, vl, vec), replace(vec, vl, vec[vl]))
      expect_identical(FUN(vec, vl, vec, invert = TRUE), replace(vec, nvl, vec[nvl]))
      expect_identical(FUN(vec, which(vl), vec), replace(vec, vl, vec[vl]))
      expect_identical(FUN(vec, which(nvl), vec), replace(vec, nvl, vec[nvl]))
      # expect_error(FUN(vec, which(nvl), vec, invert = TRUE))
    }
    replr <- function(x, i, v) {
      x[i, ] <- v
      x
    }
    expect_identical(FUN(mtcars, 1, 2), replace(mtcars, mtcars == 1, 2))
    expect_identical(FUN(mtcars, 1, 2, invert = TRUE), replace(mtcars, mtcars != 1, 2))
    if(identical(FUN, copyv)) expect_visible(FUN(mtcars, 1, mtcars$mpg, invert = TRUE)) else
      expect_invisible(FUN(mtcars, 1, mtcars$mpg, invert = TRUE))
    expect_identical(FUN(mtcars, 23L, mtcars$mpg, vind1 = TRUE), replr(mtcars, 23L, mtcars$mpg[23L]))
    expect_identical(FUN(mtcars, 3:6, mtcars$mpg), replr(mtcars, 3:6, mtcars$mpg[3:6]))
    expect_identical(FUN(mtcars, 23L, mtcars, vind1 = TRUE), replr(mtcars, 23L, mtcars[23L, ]))
    expect_identical(FUN(mtcars, 3:6, mtcars), replr(mtcars, 3:6, mtcars[3:6, ]))
    expect_error(FUN(mtcars, 23, mtcars$mpg[4:10]))
    expect_warning(FUN(mtcars, 23, mtcars[4:10]))
    expect_error(FUN(mtcars, 23L, mtcars$mpg[4:10], vind1 = TRUE))
    expect_warning(FUN(mtcars, 23L, mtcars[4:10], vind1 = TRUE))
    expect_error(FUN(mtcars, 3:6, mtcars$mpg[4:10]))
    expect_warning(FUN(mtcars, 3:6, mtcars[4:10]))
    if(identical(FUN, copyv)) {
    expect_identical(wlddev, wldcopy)
    expect_identical(mtcars, mtccopy)
    }
  }
})

wlddev <- wldcopy
mtcars <- mtccopy

}
