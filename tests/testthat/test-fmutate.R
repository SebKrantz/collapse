context("fsummarise and fmutate")


expect_equal(1, 1)

if(requireNamespace("magrittr", quietly = TRUE) && requireNamespace("dplyr", quietly = TRUE)) {
library(magrittr)

bmean <- base::mean
bsum <- base::sum
bsd <- stats::sd
bmin <- base::min
bmax <- base::max

NCRAN <- identical(Sys.getenv("NCRAN"), "TRUE")

mtc <- dplyr::as_tibble(mtcars)
gmtc <- dplyr::group_by(mtc, cyl, vs, am)

expect_equal(gsplit(mtcars$mpg, GRP(gmtc), TRUE), split(mtcars$mpg, as_factor_GRP(GRP(gmtc))))

if(NCRAN) {

test_that("fsummarise works like dplyr::summarise for tagged vector expressions", {

 # Simple computations
 expect_equal(smr(mtc, mu = bmean(mpg), sigma = bsd(mpg)), dplyr::summarise(mtc, mu = bmean(mpg), sigma = bsd(mpg)))
 # TODO: Could expand like this as well... but who needs this?
 # expect_false(all_obj_equal(smr(mtc, mu = bmean(mpg), sigma = bsd(mpg), q = quantile(mpg)),
 #                            dplyr::summarise(mtc, mu = bmean(mpg), sigma = bsd(mpg), q = quantile(mpg))))

 expect_equal(smr(mtc, mu = bmean(mpg) + bsd(mpg)), dplyr::summarise(mtc, mu = bmean(mpg) + bsd(mpg)))
 expect_equal(smr(mtc, mu = bmean(mpg) + 3), dplyr::summarise(mtc, mu = bmean(mpg) + 3))
 q <- 5
 expect_equal(smr(mtc, mu = bmean(mpg) + q), dplyr::summarise(mtc, mu = bmean(mpg) + q))
 v <- mtcars$disp
 expect_equal(smr(mtc, mu = bmean(mpg) + bmean(v)), dplyr::summarise(mtc, mu = bmean(mpg) + bmean(v)))

 # Grouped computations
 expect_equal(smr(gmtc, mpg = fmean(mpg)), dplyr::summarise(gmtc, mpg = bmean(mpg), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = bmean(mpg)), dplyr::summarise(gmtc, mpg = bmean(mpg), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = fmean(mpg), carb = fmax(carb)),
              dplyr::summarise(gmtc, mpg = bmean(mpg), carb = bmax(carb), .groups = "drop"))
 expect_equal(smr(gmtc, mpg = fmean(mpg), carb = bmax(carb)),
              dplyr::summarise(gmtc, mpg = bmean(mpg), carb = bmax(carb), .groups = "drop"))
 expect_equal(smr(gmtc, mpg = bmean(mpg), carb = bmax(carb)),
              dplyr::summarise(gmtc, mpg = bmean(mpg), carb = bmax(carb), .groups = "drop"))
 expect_equal(smr(gmtc, mpg = bmean(mpg), carb = fmax(carb)),
              dplyr::summarise(gmtc, mpg = bmean(mpg), carb = bmax(carb), .groups = "drop"))

 expect_equal(fsummarise(gmtc, mpg = bmean(mpg), carb = fmax(carb), keep.group_vars = FALSE),
              fsummarise(gmtc, mpg = bmean(mpg), carb = fmax(carb)) %>% slt(-cyl,-vs,-am))

 # Multi-return values
 expect_equal(smr(gmtc, mpg = quantile(mpg)),
              dplyr::summarise(gmtc, mpg = quantile(mpg), .groups = "drop") %>% tfm(mpg = unname(mpg)))

 # More complex expressions
 expect_equal(smr(gmtc, mpg = bmean(mpg) + 1),
              dplyr::summarise(gmtc, mpg = bmean(mpg) + 1, .groups = "drop"))

  expect_equal(smr(gmtc, mpg = bmean(mpg) + q),
              dplyr::summarise(gmtc, mpg = bmean(mpg) + q, .groups = "drop"))
 expect_equal(smr(gmtc, mpg = quantile(mpg) + q),
              dplyr::summarise(gmtc, mpg = quantile(mpg) + q, .groups = "drop") %>% tfm(mpg = unname(mpg)))

 expect_equal(smr(gmtc, mpg = bmean(mpg) + bmax(v)),
              dplyr::summarise(gmtc, mpg = bmean(mpg) + bmax(v), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = quantile(mpg) + bmax(v)),
              dplyr::summarise(gmtc, mpg = quantile(mpg) + bmax(v), .groups = "drop") %>% tfm(mpg = unname(mpg)))

 expect_equal(smr(gmtc, mpg = bmean(log(mpg))),
              dplyr::summarise(gmtc, mpg = bmean(log(mpg)), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = bmean(log(mpg)) + bmax(qsec)),
              dplyr::summarise(gmtc, mpg = bmean(log(mpg)) + bmax(qsec), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = quantile(mpg) + bmax(qsec)),
              dplyr::summarise(gmtc, mpg = quantile(mpg) + bmax(qsec), .groups = "drop") %>% tfm(mpg = unname(mpg)))

 expect_equal(smr(gmtc, mpg = fmean(log(mpg)) + fmax(qsec)),
              dplyr::summarise(gmtc, mpg = bmean(log(mpg)) + bmax(qsec), .groups = "drop"))

 expect_false(all_obj_equal(smr(gmtc, mpg = fmean(log(mpg)) + bmax(qsec)),
              dplyr::summarise(gmtc, mpg = bmean(log(mpg)) + bmax(qsec), .groups = "drop")))

 # Testing expressions turned into functions:
 mid_fun <- function(x) (bmin(x) + bmax(x)) / 2
 expect_true(all_obj_equal(smr(gmtc, mid_mpg = (bmin(mpg) + bmax(mpg)) / 2),
                           smr(gmtc, mid_mpg = (fmin(mpg) + fmax(mpg)) / 2),
                           smr(gmtc, mid_mpg = mid_fun(mpg)),
              dplyr::summarise(gmtc, mid_mpg = (bmin(mpg) + bmax(mpg)) / 2, .groups = "drop")))

 # Adding global variable:
 expect_true(all_obj_equal(smr(gmtc, mid_mpg = (bmin(mpg) + bmax(mpg)) / 2 + q),
                           smr(gmtc, mid_mpg = (fmin(mpg) + fmax(mpg)) / 2 + q),
                           smr(gmtc, mid_mpg = mid_fun(mpg) + q),
                           dplyr::summarise(gmtc, mid_mpg = (bmin(mpg) + bmax(mpg)) / 2 + q, .groups = "drop")))

 # Weighted computations
 expect_equal(smr(gmtc, mpg = weighted.mean(mpg, wt)), dplyr::summarise(gmtc, mpg = weighted.mean(mpg, wt), .groups = "drop"))
 expect_equal(smr(gmtc, mpg = fmean(mpg, wt)), dplyr::summarise(gmtc, mpg = fmean(mpg, w = wt), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = weighted.mean(mpg, wt) + 5.5), dplyr::summarise(gmtc, mpg = weighted.mean(mpg, wt) + 5.5, .groups = "drop"))
 expect_equal(smr(gmtc, mpg = fmean(mpg, wt) + 5.5), dplyr::summarise(gmtc, mpg = fmean(mpg, w = wt) + 5.5, .groups = "drop"))

 expect_equal(smr(gmtc, mpg = weighted.mean(mpg, wt) + q), dplyr::summarise(gmtc, mpg = weighted.mean(mpg, wt) + q, .groups = "drop"))
 expect_equal(smr(gmtc, mpg = fmean(mpg, wt) + q), dplyr::summarise(gmtc, mpg = fmean(mpg, w = wt) + q, .groups = "drop"))

 expect_equal(smr(gmtc, mpg = weighted.mean(mpg, wt) + bmax(v)), dplyr::summarise(gmtc, mpg = weighted.mean(mpg, wt) + bmax(v), .groups = "drop"))
 expect_equal(smr(gmtc, mpg = fmean(mpg, wt) + bmax(v)), dplyr::summarise(gmtc, mpg = fmean(mpg, w = wt) + bmax(v), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = weighted.mean(mpg, wt) + bmax(qsec)),
              dplyr::summarise(gmtc, mpg = weighted.mean(mpg, wt) + bmax(qsec), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = fmean(mpg, wt) + fmax(qsec)),
              dplyr::summarise(gmtc, mpg = fmean(mpg, w = wt) + bmax(qsec), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = quantile(mpg) + weighted.mean(mpg, wt)),
              dplyr::summarise(gmtc, mpg = quantile(mpg) + weighted.mean(mpg, wt), .groups = "drop") %>% tfm(mpg = unname(mpg)))

 expect_warning(smr(gmtc, mpg = quantile(mpg) + fmean(mpg, wt)))

})

}

wld <- dplyr::as_tibble(wlddev)
gwld <- dplyr::group_by(wlddev, iso3c)

if(NCRAN) {

test_that("fsummarise works like dplyr::summarise with across and simple usage", {

  # Simple usage
  expect_true(all_obj_equal(fsummarise(mtc, across(cyl:drat, bsum)),
                            fsummarise(mtc, across(cyl:drat, fsum)),
                            dplyr::summarise(mtc, dplyr::across(cyl:drat, bsum))))

  expect_true(all_obj_equal(fsummarise(mtc, across(5:8, bsum)),
                            fsummarise(mtc, across(5:8, fsum)),
                            dplyr::summarise(mtc, dplyr::across(5:8, bsum))))

  expect_true(all_obj_equal(fsummarise(mtc, across(-(5:8), bsum)),
                            fsummarise(mtc, across(-(5:8), fsum, .apply = FALSE)),
                            dplyr::summarise(mtc, dplyr::across(-(5:8), bsum))))

  expect_true(all_obj_equal(fsummarise(wld, across(is.numeric, bsum, na.rm = TRUE)),
                            fsummarise(wld, across(is.numeric, fsum)) %>% dapply(unattrib, drop = FALSE),
                            dplyr::summarise(wld, dplyr::across(where(is.numeric), bsum, na.rm = TRUE))))

  expect_true(all_obj_equal(fsummarise(mtc, across(NULL, bsum, na.rm = TRUE)),
                            fsummarise(mtc, across(NULL, fsum)),
                            dplyr::summarise(mtc, dplyr::across(everything(), bsum, na.rm = TRUE))))

  expect_equal(fsummarise(mtc, across(cyl:vs, bsum)),
               fsummarise(mtc, cyl = bsum(cyl), across(disp:qsec, fsum), vs = fsum(vs)))

  # Simple programming use
  vlist <- list(mtc %>% fselect(cyl:drat, return = "names"), 5:8, NULL) # -(5:8), is.numeric
  for(i in seq_along(vlist)) {
    expect_true(all_obj_equal(fsummarise(mtc, across(vlist[[i]], bsum)),
                              fsummarise(mtc, across(vlist[[i]], fsum)),
                              dplyr::summarise(mtc, dplyr::across(if(is.null(vlist[[i]])) everything() else vlist[[i]], bsum))))
    v <- vlist[[i]]
    expect_true(all_obj_equal(fsummarise(mtc, across(v, bsum)),
                              fsummarise(mtc, across(v, fsum)),
                              dplyr::summarise(mtc, dplyr::across(if(is.null(v)) everything() else v, bsum))))
  }

  # Simple usage and multiple functions
  expect_true(all_obj_equal(fsummarise(mtc, across(cyl:drat, list(bmean, bsum))),
                            fsummarise(mtc, across(cyl:drat, list(bmean = fmean, bsum = fsum))),
                            dplyr::summarise(mtc, dplyr::across(cyl:drat, list(bmean = bmean, bsum = bsum)))))

  expect_true(all_obj_equal(fsummarise(mtc, across(NULL, list(bmean, bsum))),
                            fsummarise(mtc, across(NULL, list(bmean = fmean, bsum = fsum))),
                            dplyr::summarise(mtc, dplyr::across(everything(), list(bmean = bmean, bsum = bsum)))))

  # Passing additional arguments
  expect_true(all_obj_equal(fsummarise(mtc, across(cyl:drat, bsum, na.rm = FALSE)),
                            fsummarise(mtc, across(cyl:drat, fsum, na.rm = FALSE)),
                            dplyr::summarise(mtc, dplyr::across(cyl:drat, bsum, na.rm = FALSE))))

  expect_true(all_obj_equal(fsummarise(mtc, across(cyl:drat, weighted.mean, w = wt)),
                            fsummarise(mtc, across(cyl:drat, fmean, w = wt)),
                            dplyr::summarise(mtc, dplyr::across(cyl:drat, weighted.mean, w = wt))))

  expect_true(all_obj_equal(fsummarise(mtc, across(cyl:drat, list(mean = weighted.mean, sum = fsum), w = wt)),
                            fsummarise(mtc, across(cyl:drat, list(mean = fmean, sum = fsum), w = wt)),
                            dplyr::summarise(mtc, dplyr::across(cyl:drat, list(mean = weighted.mean, sum = fsum), w = wt))))

  # Simple programming use
  flist <- list(bsum, list(bmean = bmean, bsum = bsum), list(bmean, bsum)) # c("bmean", "bsum"), c(mean = "fmean", sum = "fsum")
  for(i in seq_along(flist)) {
    expect_equal(fsummarise(mtc, across(cyl:drat, flist[[i]])),
                 dplyr::summarise(mtc, dplyr::across(cyl:drat, flist[[i]])))
    f <- flist[[i]]
    expect_equal(fsummarise(mtc, across(cyl:drat, f)),
                 dplyr::summarise(mtc, dplyr::across(cyl:drat, f)))
  }

})


test_that("fsummarise works like dplyr::summarise with across and grouped usage", {

  # Simple usage
  expect_true(all_obj_equal(fsummarise(gmtc, across(hp:drat, bsum)),
                            fsummarise(gmtc, across(hp:drat, fsum)),
                            dplyr::summarise(gmtc, dplyr::across(hp:drat, bsum), .groups = "drop")))

  expect_true(all_obj_equal(fsummarise(gmtc, across(5:7, bsum)),
                            fsummarise(gmtc, across(5:7, fsum)),
                            dplyr::summarise(gmtc, dplyr::across(4:6, bsum), .groups = "drop")))

  expect_true(all_obj_equal(fsummarise(gwld, across(is.numeric, bsum, na.rm = TRUE)) %>% setLabels(NULL),
                            fsummarise(gwld, across(is.numeric, fsum)) %>% replace_NA() %>% setLabels(NULL),
                            dplyr::summarise(gwld, dplyr::across(where(is.numeric), bsum, na.rm = TRUE))))

  expect_true(all_obj_equal(fsummarise(gmtc, across(NULL, bsum, na.rm = TRUE)) %>% setLabels(NULL),
                            fsummarise(gmtc, across(NULL, fsum)) %>% setLabels(NULL),
                            dplyr::summarise(gmtc, dplyr::across(everything(), bsum, na.rm = TRUE), .groups = "drop")))

  expect_equal(fsummarise(gmtc, across(NULL, bsum, na.rm = TRUE), keep.group_vars = FALSE),
               fsummarise(gmtc, across(NULL, bsum, na.rm = TRUE)) %>% slt(-cyl,-vs,-am))

  expect_equal(fsummarise(gmtc, across(cyl:vs, bsum)),
               fsummarise(gmtc, cyl = bsum(cyl), across(disp:qsec, fsum), vs = fsum(vs)))

  # Simple programming use
  vlist <- list(mtc %>% fselect(hp:drat, return = "names"), NULL) # -(5:8), is.numeric
  for(i in seq_along(vlist)) {
    expect_true(all_obj_equal(fsummarise(gmtc, across(vlist[[i]], bsum)),
                              fsummarise(gmtc, across(vlist[[i]], fsum)),
                              dplyr::summarise(gmtc, dplyr::across(if(is.null(vlist[[i]])) everything() else vlist[[i]], bsum), .groups = "drop")))
    v <- vlist[[i]]
    expect_true(all_obj_equal(fsummarise(gmtc, across(v, bsum)),
                              fsummarise(gmtc, across(v, fsum)),
                              dplyr::summarise(gmtc, dplyr::across(if(is.null(v)) everything() else v, bsum), .groups = "drop")))
  }

  # Simple usage and multiple functions
  expect_true(all_obj_equal(fsummarise(gmtc, across(hp:drat, list(bmean, bsum))),
                            fsummarise(gmtc, across(hp:drat, list(bmean = fmean, bsum = fsum))),
                            dplyr::summarise(gmtc, dplyr::across(hp:drat, list(bmean = bmean, bsum = bsum)), .groups = "drop")))

  expect_true(all_obj_equal(fsummarise(gmtc, across(NULL, list(bmean, bsum))),
                            fsummarise(gmtc, across(NULL, list(bmean = fmean, bsum = fsum))),
                            dplyr::summarise(gmtc, dplyr::across(everything(), list(bmean = bmean, bsum = bsum)), .groups = "drop")))

  # Passing additional arguments
  expect_true(all_obj_equal(fsummarise(gwld, across(c("PCGDP", "LIFEEX"), bsum, na.rm = TRUE))  %>% setLabels(NULL),
                            fsummarise(gwld, across(c("PCGDP", "LIFEEX"), fsum, na.rm = TRUE))  %>% setLabels(NULL) %>% replace_NA(),
                            dplyr::summarise(gwld, dplyr::across(c("PCGDP", "LIFEEX"), bsum, na.rm = TRUE), .groups = "drop")))

  expect_true(all_obj_equal(fsummarise(gmtc, across(hp:drat, weighted.mean, w = wt)),
                            fsummarise(gmtc, across(hp:drat, fmean, w = wt)),
                            dplyr::summarise(gmtc, dplyr::across(hp:drat, weighted.mean, w = wt), .groups = "drop")))

  expect_equal(fsummarise(gmtc, across(cyl:vs, weighted.mean, w = wt)),
               fsummarise(gmtc, cyl = weighted.mean(cyl, wt), across(disp:qsec, fmean, w = wt), vs = fmean(vs, wt)))

  expect_true(all_obj_equal(fsummarise(gmtc, across(hp:drat, list(mean = weighted.mean, sum = fsum), w = wt)),
                            fsummarise(gmtc, across(hp:drat, list(mean = fmean, sum = fsum), w = wt)),
                            dplyr::summarise(gmtc, dplyr::across(hp:drat, list(mean = weighted.mean, sum = fsum), w = wt), .groups = "drop")))

  # Simple programming use
  flist <- list(bsum, list(bmean = bmean, bsum = bsum), list(bmean, bsum)) # c("bmean", "bsum"), c(mean = "fmean", sum = "fsum")
  for(i in seq_along(flist)) {
    expect_equal(fsummarise(gmtc, across(hp:drat, flist[[i]])),
                 dplyr::summarise(gmtc, dplyr::across(hp:drat, flist[[i]]), .groups = "drop"))
    f <- flist[[i]]
    expect_equal(fsummarise(gmtc, across(hp:drat, f)),
                 dplyr::summarise(gmtc, dplyr::across(hp:drat, f), .groups = "drop"))
  }

})

}

test_that("fsummarise miscellaneous things", {

  expect_equal(smr(gmtc, acr(disp:hp, c("bmean", "bsd"))),
               fsummarise(gmtc, across(disp:hp, c("bmean", "bsd"), .transpose = FALSE)) %>%
                 colorderv(c(4,6,5,7), pos = "exchange"))

  expect_identical(names(smr(gmtc, acr(disp:hp, fmean, .names = TRUE)))[4:5], c("disp_fmean", "hp_fmean"))
  expect_identical(names(smr(gmtc, acr(disp:hp, bmean, .names = TRUE)))[4:5], c("disp_bmean", "hp_bmean"))

  pwcorDF <- function(x, w = NULL) qDF(pwcor(x, w = w), "var")
  expect_equal(
    mtcars %>% gby(cyl) %>% smr(acr(disp:hp, pwcorDF, .apply = FALSE)),
    rsplit(mtcars, disp + hp ~ cyl) %>% lapply(pwcorDF) %>% unlist2d("cyl", "var") %>% tfm(cyl = as.numeric(cyl))
  )

  if(identical(Sys.getenv("LOCAL"), "TRUE")) # No tests depending on suggested package (except for major ones).
  expect_equal(
    mtcars %>% gby(cyl) %>% smr(acr(disp:hp, pwcorDF, w = wt, .apply = FALSE)),
    rsplit(mtcars, disp + hp + wt ~ cyl) %>% lapply(function(x) pwcorDF(gv(x, 1:2), w = x$wt)) %>%
      unlist2d("cyl", "var") %>% tfm(cyl = as.numeric(cyl))
  )

  if(requireNamespace("data.table", quietly = TRUE)) {
  lmest <- function(x) list(Mods = list(lm(disp~., x)))
  expect_equal(
    qDT(mtcars) %>% gby(cyl) %>% smr(acr(disp:hp, lmest, .apply = FALSE)),
    qDT(mtcars) %>% rsplit(disp + hp ~ cyl) %>% lapply(lmest) %>% data.table::rbindlist(idcol = "cyl") %>%
      tfm(cyl = as.numeric(cyl))
  )
  }
})

if(NCRAN) {

test_that("fmutate works as intended for simple usage", {

  expect_equal(fmutate(mtc, bla = 1), dplyr::mutate(mtc, bla = 1))
  expect_equal(fmutate(mtc, bla = list(1)), dplyr::mutate(mtc, bla = list(1)))
  expect_equal(fmutate(mtc, bla = as.list(mpg)), dplyr::mutate(mtc, bla = as.list(mpg)))
  expect_equal(fmutate(mtc, mu = bmean(mpg)), dplyr::mutate(mtc, mu = bmean(mpg)))
  expect_equal(fmutate(mtc, mu = bmean(mpg), mpg = NULL), dplyr::mutate(mtc, mu = bmean(mpg), mpg = NULL))
  expect_equal(fmutate(mtc, mu = bmean(mpg), dmu = mpg - mu), dplyr::mutate(mtc, mu = bmean(mpg), dmu = mpg - mu))
  expect_equal(fmutate(mtc, mu = log(mpg)), dplyr::mutate(mtc, mu = log(mpg)))
  expect_equal(fmutate(mtc, mu = log(mpg), dmu = mpg - mu), dplyr::mutate(mtc, mu = log(mpg), dmu = mpg - mu))

  expect_true(all_obj_equal(
    dplyr::mutate(mtc, dmu = mpg - bmean(mpg)),
    fmutate(mtc, dmu = mpg - bmean(mpg)),
    fmutate(mtc, dmu = mpg - fmean(mpg)),
    fmutate(mtc, dmu = fmean(mpg, TRA = "-")),
    fmutate(mtc, dmu = fwithin(mpg))
  ))

  # With groups:
  expect_equal(fmutate(gmtc, bla = 1), dplyr::mutate(gmtc, bla = 1))
  expect_equal(fmutate(gmtc, mu = bmean(mpg)), dplyr::mutate(gmtc, mu = bmean(mpg)))
  expect_equal(fmutate(gmtc, mu = bmean(mpg), mpg = NULL), dplyr::mutate(gmtc, mu = bmean(mpg), mpg = NULL))
  expect_equal(fmutate(gmtc, mu = bmean(mpg), dmu = mpg - mu), dplyr::mutate(gmtc, mu = bmean(mpg), dmu = mpg - mu))
  expect_equal(fmutate(gmtc, mu = log(mpg)), dplyr::mutate(gmtc, mu = log(mpg)))
  expect_equal(fmutate(gmtc, mu = log(mpg), dmu = mpg - mu), dplyr::mutate(gmtc, mu = log(mpg), dmu = mpg - mu))

  expect_true(all_obj_equal(
    dplyr::mutate(gmtc, dmu = mpg - bmean(mpg)),
    fmutate(gmtc, dmu = mpg - bmean(mpg)),
    fmutate(gmtc, dmu = mpg - fmean(mpg)),
    fmutate(gmtc, dmu = fmean(mpg, TRA = "-")),
    fmutate(gmtc, dmu = fwithin(mpg))
  ))

})

}

test_that("fmutate with across works like ftransformv", {

   expect_true(all_obj_equal(

     ftransformv(mtcars, cyl:vs, fwithin, w = wt, apply = TRUE),
     ftransformv(mtcars, cyl:vs, fwithin, w = wt, apply = FALSE),
     fmutate(mtcars, across(cyl:vs, fwithin, w = wt, .apply = TRUE)),
     fmutate(mtcars, across(cyl:vs, fwithin, w = wt, .apply = FALSE))
     # fmutate(mtcars, fwithin(.data, w = .data[["wt"]]), .cols = slt(., cyl:vs, return = "names"))

   ))


  expect_true(all_obj_equal(

    ftransformv(mtcars, cyl:vs, fwithin, g = list(cyl, vs, am), apply = TRUE) %>% setRownames(),
    ftransformv(mtcars, cyl:vs, fwithin, g = list(cyl, vs, am), apply = FALSE) %>% setRownames(),
    fmutate(gmtc, across(cyl:vs, fwithin, .apply = TRUE)) %>% qDF(),
    fmutate(gmtc, across(cyl:vs, fwithin, .apply = FALSE)) %>% qDF(),
    fmutate(gmtc, across(cyl:vs, fmean, TRA = "-", .apply = TRUE)) %>% qDF(),
    fmutate(gmtc, across(cyl:vs, fmean, TRA = "-", .apply = FALSE)) %>% qDF(),
    fmutate(gmtc, across(cyl:vs, function(x) x - bmean(x), .apply = TRUE)) %>% qDF(),
    fmutate(gmtc, across(cyl:vs, function(x) lapply(x, function(y) y - bmean(y)), .apply = FALSE)) %>% qDF(),
    gmtc %>% fmutate(fwithin(.data), .cols = slt(., cyl:vs, return = "names")) %>% qDF()

  ))

  expect_true(all_obj_equal(

    ftransformv(mtcars, cyl:vs, fwithin, g = list(cyl, vs, am), w = wt, apply = TRUE) %>% setRownames(),
    ftransformv(mtcars, cyl:vs, fwithin, g = list(cyl, vs, am), w = wt, apply = FALSE) %>% setRownames(),
    fmutate(gmtc, across(cyl:vs, fwithin, w = wt, .apply = TRUE)) %>% qDF(),
    fmutate(gmtc, across(cyl:vs, fwithin, w = wt, .apply = FALSE)) %>% qDF(),
    fmutate(gmtc, across(cyl:vs, fmean, TRA = "-", w = wt, .apply = TRUE)) %>% qDF(),
    fmutate(gmtc, across(cyl:vs, fmean, TRA = "-", w = wt, .apply = FALSE)) %>% qDF(),
    fmutate(gmtc, across(cyl:vs, function(x, w) x - weighted.mean(x, w), w = wt, .apply = TRUE)) %>% qDF()

  ))

})

test_that("fmutate with across reorders correctly", {

  for(i in seq_col(wlddev)) {
    gdf <- fgroup_by(wlddev, i)
    expect_true(all_identical(
      wlddev,
      fungroup(fmutate(gdf, across(c(PCGDP, LIFEEX), identity))),
      fungroup(fmutate(gdf, across(.fns = identity))),
      fungroup(fmutate(gdf, list(PCGDP = PCGDP, LIFEEX = LIFEEX))),
      fungroup(fmutate(gdf, (.data), .cols = .c(PCGDP, LIFEEX))),
      fungroup(fmutate(gdf, (.data)))
    ))
  }

})

test_that("fsummarise and fmutate with arbitrary expressions", {

  expect_true(
    all_obj_equal(
    fsummarise(gmtc, qDF(cor(cbind(mpg, wt, hp))), mpg_qs = quantile(mpg, c(0.25, 0.5, 0.75))),
    fsummarise(gmtc, acr(c(mpg, wt, hp), function(x) qDF(cor(x)), .apply = FALSE), mpg_qs = quantile(mpg, c(0.25, 0.5, 0.75))),
    fsummarise(gmtc, qDF(cor(.data)), .cols = .c(mpg, wt, hp), mpg_qs = quantile(mpg, c(0.25, 0.5, 0.75))),
    fsummarise(gmtc, qDF(cor(slt(.data, mpg, wt, hp))), mpg_qs = quantile(mpg, c(0.25, 0.5, 0.75))))
  )

  expect_true(
    all_obj_equal(
      fmutate(gmtc, fscale(list(mpg = mpg, wt = wt, hp = hp)), bla = 1, mu = fmean(mpg), su = sum(hp)),
      fmutate(gmtc, acr(c(mpg, wt, hp), fscale), bla = 1, mu = fmean(mpg), su = sum(hp)),
      fmutate(gmtc, acr(c(mpg, wt, hp), function(x) fscale(x), .apply = FALSE), bla = 1, mu = fmean(mpg), su = sum(hp)),
      fmutate(gmtc, fscale(.data), .cols = .c(mpg, wt, hp), bla = 1, mu = fmean(mpg), su = sum(hp)),
      fmutate(gmtc, fscale(slt(.data, mpg, wt, hp)), bla = 1, mu = fmean(mpg), su = sum(hp)))
  )

  expect_equal(fmutate(gmtc, acr(NULL, fscale)), fmutate(gmtc, fscale(.data)))
  expect_equal(fmutate(gmtc, acr(mpg:carb, fscale)), fmutate(gmtc, fscale(.data), .cols = seq_col(gmtc)))

})

if(NCRAN) {

test_that("fmutate miscellaneous", {

  expect_true(length(fmutate(mtcars, across(cyl:vs, W, w = wt, .names = NULL))) > 15)
  expect_true(length(fmutate(mtcars, across(cyl:vs, list(D, W), .names = FALSE, .transpose = FALSE))) > 15)

  expect_equal(  fmutate(mtcars, across(cyl:vs, L, stubs = FALSE)),
                 fmutate(mtcars, across(cyl:vs, flag))
  )

  expect_true(length(fmutate(mtcars, across(cyl:vs, L))) > length(fmutate(mtcars, across(cyl:vs, flag))))

  expect_equal(
    fmutate(mtcars, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "used"),
    dplyr::mutate(mtcars, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "used")
  )

  expect_equal(
    fmutate(mtcars, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "unused"),
    dplyr::mutate(mtcars, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "unused")
  )

  expect_equal(
    fmutate(mtcars, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "none"),
    dplyr::transmute(mtcars, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb)
  )

  expect_identical(names(fmutate(mtcars, a = mpg, b = a, c = cyl, hp = wt, .keep = "unused")), c(setdiff(names(mtcars), .c(mpg, cyl, wt)), letters[1:3]))
  expect_identical(names(fmutate(mtcars, a = mpg, b = a, c = cyl, hp = wt, .keep = "none")), c("a", "b", "c", "hp"))

  expect_equal(
    fmutate(gmtc, a = fmax(mpg), b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "used"),
    dplyr::mutate(gmtc, a = bmax(mpg), b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "used")
  )

  expect_equal(
    fmutate(gmtc, a = fmax(mpg), b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "unused"),
    dplyr::mutate(gmtc, a = bmax(mpg), b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "unused")
  )

  expect_equal(
    fmutate(gmtc, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "none"),
    dplyr::transmute(gmtc, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb)
  )

  # Inconsistent with the above and also inefficient...
  # expect_equal(
  #   fmutate(gmtc, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, cyl = cyl, .keep = "none"),
  #   dplyr::mutate(gmtc, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, cyl = cyl, .keep = "none")
  # )

  expect_equal(flast(names(fmutate(mtcars, across(cyl:vs, function(x) list(ps = kit::psum(x)), .apply = FALSE)))), "ps")

  expect_equal(
    fmutate(mtcars, across(cyl:vs, data.table::shift, .apply = FALSE, .names = FALSE)),
    fmutate(mtcars, across(cyl:vs, data.table::shift))
  )

  # Testing expressions turned into functions:
  bcumsum = base::cumsum
  lorentz_fun <- function(x) bcumsum(x) / bsum(x)
  gmtc = mtc %>% roworder(mpg) %>% dplyr::group_by(cyl, vs, am)
  expect_true(all_obj_equal(mtt(gmtc, lorentz_mpg = bcumsum(mpg) / bsum(mpg)),
                            mtt(gmtc, lorentz_mpg = lorentz_fun(mpg)),
                            mtt(gmtc, lorentz_mpg = fcumsum(mpg) / fsum(mpg)), # doesn't work because of global sorting...
                            dplyr::mutate(gmtc, lorentz_mpg = bcumsum(mpg) / bsum(mpg))))

  # Adding global variable:
  q = 5
  expect_true(all_obj_equal(mtt(gmtc, lorentz_mpg = bcumsum(mpg) / bsum(mpg) + q),
                            mtt(gmtc, lorentz_mpg = lorentz_fun(mpg) + q),
                            mtt(gmtc, lorentz_mpg = fcumsum(mpg) / fsum(mpg) + q),
                            dplyr::mutate(gmtc, lorentz_mpg = bcumsum(mpg) / bsum(mpg) + q)))

})

}

test_that(".names works properly", {
  expect_equal(
    smr(gmtc, acr(c(hp, wt), list(sum, max, min))),
    smr(gmtc, acr(c(hp, wt), list(sum, max, min), .names = TRUE))
  )
  expect_equal(
   smr(gmtc, acr(c(hp, wt), list(sum, max, min))),
   smr(gmtc, acr(c(hp, wt), list(sum, max, min), .names = function(c, f) paste0(c, "_", f)))
  )
  expect_equal(
    smr(gmtc, acr(c(hp, wt), list(sum, max, min), .names = "flip")),
    smr(gmtc, acr(c(hp, wt), list(sum, max, min), .names = function(c, f) paste0(f, "_", c)))
  )
  expect_equal(
    smr(gmtc, acr(c(hp, wt), list(sum, max, min), .transpose = FALSE)),
    smr(gmtc, acr(c(hp, wt), list(sum, max, min), .names = function(c, f) paste0(c, "_", f), .transpose = FALSE))
  )
  expect_equal(
    smr(gmtc, acr(c(hp, wt), list(sum, max, min), .names = "flip", .transpose = FALSE)),
    smr(gmtc, acr(c(hp, wt), list(sum, max, min), .names = function(c, f) paste0(f, "_", c), .transpose = FALSE))
  )
  expect_equal(
    names(smr(gmtc, acr(c(hp, wt), list(sum, max, min), .names = FALSE))),
    .c(cyl, vs, am, hp, hp, hp, wt, wt, wt)
  )
  expect_equal(
    names(smr(gmtc, acr(c(hp, wt), list(sum, max, min), .names = FALSE, .transpose = FALSE))),
    .c(cyl, vs, am, hp, wt, hp, wt, hp, wt)
  )
  expect_equal(
    names(smr(gmtc, acr(c(hp, wt), sum, .names = FALSE))),
    .c(cyl, vs, am, hp, wt)
  )
  expect_equal(
    names(smr(gmtc, acr(c(hp, wt), sum, .names = FALSE, .transpose = FALSE))),
    .c(cyl, vs, am, hp, wt)
  )
  expect_equal(
    names(smr(gmtc, acr(c(hp, wt), sum, .names = TRUE))),
    .c(cyl, vs, am, hp_sum, wt_sum)
  )
  expect_equal(
    names(smr(gmtc, acr(c(hp, wt), sum, .names = "flip"))),
    .c(cyl, vs, am, sum_hp, sum_wt)
  )
  expect_equal(
    names(smr(gmtc, acr(c(hp, wt), sum, .names = "flip", .transpose = FALSE))),
    .c(cyl, vs, am, sum_hp, sum_wt)
  )

})


test_that("Warnings for unnamed scalar and vector-valued arguments passed", {
  tf <- function(x, ...) x
  expect_warning(mtt(gmtc, acr(hp:carb, tf, TRUE, wt)))
  expect_warning(mtt(gmtc, acr(hp:carb, tf, wt, TRUE)))
  expect_warning(mtt(gmtc, acr(hp:carb, tf, TRUE, wt, .apply = FALSE)))
  expect_warning(mtt(gmtc, acr(hp:carb, tf, wt, TRUE, .apply = FALSE)))
})

if(FALSE) {


  fmutate(mtcars, across(cyl:vs, list(D, W), .names = TRUE, .transpose = FALSE)) %>% head(3)
  fmutate(mtcars, across(cyl:vs, list(D, W), .names = TRUE, .transpose = TRUE)) %>% head(3)

  fmutate(mtcars, across(cyl:vs, list(D, W), .names = TRUE, .apply = FALSE, .transpose = FALSE)) %>% head(3)
  fmutate(mtcars, across(cyl:vs, list(D, W), .names = FALSE, .apply = FALSE, .transpose = FALSE)) %>% head(3)

  fmutate(mtcars, across(cyl:vs, list(W, kit::psum), w = wt)) %>% head(3)
  fmutate(mtcars, across(cyl:vs, kit::psum)) %>% head(3)

  fmutate(mtcars, across(cyl:vs, identity, .apply = FALSE)) # 51 microesecond median on windows
  fmutate(mtcars, across(cyl:vs, identity)) # 62 microesecond median on windows

  fmutate(mtcars, across(cyl:vs, L))

  # TODO: Test all potential issues with environemtns etc. See if there are smarter ways to
  # incorporate internal functions, data and objects in the global environment.


}

}
