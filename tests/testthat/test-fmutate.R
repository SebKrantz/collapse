context("fsummarise and fmutate")

library(dplyr)

mtc <- as_tibble(mtcars)
gmtc <- group_by(mtc, cyl, vs, am)

expect_equal(gsplit(mtcars$mpg, GRP(gmtc), TRUE), split(mtcars$mpg, as_factor_GRP(GRP(gmtc))))


test_that("fsummarise works like summarise for tagged vector expressions", {

 # Simple computations
 expect_equal(smr(mtc, mu = mean(mpg), sigma = sd(mpg)), summarise(mtc, mu = mean(mpg), sigma = sd(mpg)))
 # TODO: Could expand like this as well... but who needs this?
 # expect_false(all_obj_equal(smr(mtc, mu = mean(mpg), sigma = sd(mpg), q = quantile(mpg)),
 #                            summarise(mtc, mu = mean(mpg), sigma = sd(mpg), q = quantile(mpg))))

 expect_equal(smr(mtc, mu = mean(mpg) + sd(mpg)), summarise(mtc, mu = mean(mpg) + sd(mpg)))
 expect_equal(smr(mtc, mu = mean(mpg) + 3), summarise(mtc, mu = mean(mpg) + 3))
 q <- 5
 expect_equal(smr(mtc, mu = mean(mpg) + q), summarise(mtc, mu = mean(mpg) + q))
 v <- mtcars$disp
 expect_equal(smr(mtc, mu = mean(mpg) + mean(v)), summarise(mtc, mu = mean(mpg) + mean(v)))

 # Grouped computations
 expect_equal(smr(gmtc, mpg = fmean(mpg)), summarise(gmtc, mpg = mean(mpg), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = mean(mpg)), summarise(gmtc, mpg = mean(mpg), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = fmean(mpg), carb = fmax(carb)),
              summarise(gmtc, mpg = mean(mpg), carb = max(carb), .groups = "drop"))
 expect_equal(smr(gmtc, mpg = fmean(mpg), carb = max(carb)),
              summarise(gmtc, mpg = mean(mpg), carb = max(carb), .groups = "drop"))
 expect_equal(smr(gmtc, mpg = mean(mpg), carb = max(carb)),
              summarise(gmtc, mpg = mean(mpg), carb = max(carb), .groups = "drop"))
 expect_equal(smr(gmtc, mpg = mean(mpg), carb = fmax(carb)),
              summarise(gmtc, mpg = mean(mpg), carb = max(carb), .groups = "drop"))

 expect_equal(fsummarise(gmtc, mpg = mean(mpg), carb = fmax(carb), keep.group_vars = FALSE),
              fsummarise(gmtc, mpg = mean(mpg), carb = fmax(carb)) |> slt(-cyl,-vs,-am))

 # Multi-return values
 expect_equal(smr(gmtc, mpg = quantile(mpg)),
              summarise(gmtc, mpg = quantile(mpg), .groups = "drop") |> tfm(mpg = unname(mpg)))

 # More complex expressions
 expect_equal(smr(gmtc, mpg = mean(mpg) + 1),
              summarise(gmtc, mpg = mean(mpg) + 1, .groups = "drop"))

  expect_equal(smr(gmtc, mpg = mean(mpg) + q),
              summarise(gmtc, mpg = mean(mpg) + q, .groups = "drop"))
 expect_equal(smr(gmtc, mpg = quantile(mpg) + q),
              summarise(gmtc, mpg = quantile(mpg) + q, .groups = "drop") |> tfm(mpg = unname(mpg)))

 expect_equal(smr(gmtc, mpg = mean(mpg) + max(v)),
              summarise(gmtc, mpg = mean(mpg) + max(v), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = quantile(mpg) + max(v)),
              summarise(gmtc, mpg = quantile(mpg) + max(v), .groups = "drop") |> tfm(mpg = unname(mpg)))

 expect_equal(smr(gmtc, mpg = mean(log(mpg))),
              summarise(gmtc, mpg = mean(log(mpg)), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = mean(log(mpg)) + max(qsec)),
              summarise(gmtc, mpg = mean(log(mpg)) + max(qsec), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = quantile(mpg) + max(qsec)),
              summarise(gmtc, mpg = quantile(mpg) + max(qsec), .groups = "drop") |> tfm(mpg = unname(mpg)))

 expect_equal(smr(gmtc, mpg = fmean(log(mpg)) + fmax(qsec)),
              summarise(gmtc, mpg = mean(log(mpg)) + max(qsec), .groups = "drop"))

 expect_false(all_obj_equal(smr(gmtc, mpg = fmean(log(mpg)) + max(qsec)),
              summarise(gmtc, mpg = mean(log(mpg)) + max(qsec), .groups = "drop")))

 # Weighted computations
 expect_equal(smr(gmtc, mpg = weighted.mean(mpg, wt)), summarise(gmtc, mpg = weighted.mean(mpg, wt), .groups = "drop"))
 expect_equal(smr(gmtc, mpg = fmean(mpg, wt)), summarise(gmtc, mpg = fmean(mpg, w = wt), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = weighted.mean(mpg, wt) + 5.5), summarise(gmtc, mpg = weighted.mean(mpg, wt) + 5.5, .groups = "drop"))
 expect_equal(smr(gmtc, mpg = fmean(mpg, wt) + 5.5), summarise(gmtc, mpg = fmean(mpg, w = wt) + 5.5, .groups = "drop"))

 expect_equal(smr(gmtc, mpg = weighted.mean(mpg, wt) + q), summarise(gmtc, mpg = weighted.mean(mpg, wt) + q, .groups = "drop"))
 expect_equal(smr(gmtc, mpg = fmean(mpg, wt) + q), summarise(gmtc, mpg = fmean(mpg, w = wt) + q, .groups = "drop"))

 expect_equal(smr(gmtc, mpg = weighted.mean(mpg, wt) + max(v)), summarise(gmtc, mpg = weighted.mean(mpg, wt) + max(v), .groups = "drop"))
 expect_equal(smr(gmtc, mpg = fmean(mpg, wt) + max(v)), summarise(gmtc, mpg = fmean(mpg, w = wt) + max(v), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = weighted.mean(mpg, wt) + max(qsec)),
              summarise(gmtc, mpg = weighted.mean(mpg, wt) + max(qsec), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = fmean(mpg, wt) + fmax(qsec)),
              summarise(gmtc, mpg = fmean(mpg, w = wt) + max(qsec), .groups = "drop"))

 expect_equal(smr(gmtc, mpg = quantile(mpg) + weighted.mean(mpg, wt)),
              summarise(gmtc, mpg = quantile(mpg) + weighted.mean(mpg, wt), .groups = "drop") |> tfm(mpg = unname(mpg)))

 expect_warning(smr(gmtc, mpg = quantile(mpg) + fmean(mpg, wt)))

})

wld <- as_tibble(wlddev)
gwld <- group_by(wlddev, country)

test_that("fsummarise works like summarise with across and simple usage", {

  # Simple usage
  expect_true(all_obj_equal(fsummarise(mtc, across(cyl:drat, sum)),
                            fsummarise(mtc, across(cyl:drat, fsum)),
                            summarise(mtc, across(cyl:drat, sum))))

  expect_true(all_obj_equal(fsummarise(mtc, across(5:8, sum)),
                            fsummarise(mtc, across(5:8, fsum)),
                            summarise(mtc, across(5:8, sum))))

  expect_true(all_obj_equal(fsummarise(mtc, across(-(5:8), sum)),
                            fsummarise(mtc, across(-(5:8), fsum, .apply = FALSE)),
                            summarise(mtc, across(-(5:8), sum))))

  expect_true(all_obj_equal(fsummarise(wld, across(is.numeric, sum, na.rm = TRUE)),
                            fsummarise(wld, across(is.numeric, fsum)) |> dapply(unattrib, drop = FALSE),
                            summarise(wld, across(where(is.numeric), sum, na.rm = TRUE))))

  expect_true(all_obj_equal(fsummarise(mtc, across(NULL, sum, na.rm = TRUE)),
                            fsummarise(mtc, across(NULL, fsum)),
                            summarise(mtc, across(everything(), sum, na.rm = TRUE))))

  expect_equal(fsummarise(mtc, across(cyl:vs, sum)),
               fsummarise(mtc, cyl = sum(cyl), across(disp:qsec, fsum), vs = fsum(vs)))

  # Simple programming use
  vlist <- list(mtc |> fselect(cyl:drat, return = "names"), 5:8, NULL) # -(5:8), is.numeric
  for(i in seq_along(vlist)) {
    expect_true(all_obj_equal(fsummarise(mtc, across(vlist[[i]], sum)),
                              fsummarise(mtc, across(vlist[[i]], fsum)),
                              summarise(mtc, across(if(is.null(vlist[[i]])) everything() else vlist[[i]], sum))))
    v <- vlist[[i]]
    expect_true(all_obj_equal(fsummarise(mtc, across(v, sum)),
                              fsummarise(mtc, across(v, fsum)),
                              summarise(mtc, across(if(is.null(v)) everything() else v, sum))))
  }

  # Simple usage and multiple functions
  expect_true(all_obj_equal(fsummarise(mtc, across(cyl:drat, list(mean, sum))),
                            fsummarise(mtc, across(cyl:drat, list(mean = fmean, sum = fsum))),
                            summarise(mtc, across(cyl:drat, list(mean = mean, sum = sum)))))

  expect_true(all_obj_equal(fsummarise(mtc, across(NULL, list(mean, sum))),
                            fsummarise(mtc, across(NULL, list(mean = fmean, sum = fsum))),
                            summarise(mtc, across(everything(), list(mean = mean, sum = sum)))))

  # Passing additional arguments
  expect_true(all_obj_equal(fsummarise(mtc, across(cyl:drat, sum, na.rm = FALSE)),
                            fsummarise(mtc, across(cyl:drat, fsum, na.rm = FALSE)),
                            summarise(mtc, across(cyl:drat, sum, na.rm = FALSE))))

  expect_true(all_obj_equal(fsummarise(mtc, across(cyl:drat, weighted.mean, w = wt)),
                            fsummarise(mtc, across(cyl:drat, fmean, w = wt)),
                            summarise(mtc, across(cyl:drat, weighted.mean, w = wt))))

  expect_true(all_obj_equal(fsummarise(mtc, across(cyl:drat, list(mean = weighted.mean, sum = fsum), w = wt)),
                            fsummarise(mtc, across(cyl:drat, list(mean = fmean, sum = fsum), w = wt)),
                            summarise(mtc, across(cyl:drat, list(mean = weighted.mean, sum = fsum), w = wt))))

  # Simple programming use
  flist <- list(sum, list(mean = mean, sum = sum), list(mean, sum)) # c("mean", "sum"), c(mean = "fmean", sum = "fsum")
  for(i in seq_along(flist)) {
    expect_equal(fsummarise(mtc, across(cyl:drat, flist[[i]])),
                 summarise(mtc, across(cyl:drat, flist[[i]])))
    f <- flist[[i]]
    expect_equal(fsummarise(mtc, across(cyl:drat, f)),
                 summarise(mtc, across(cyl:drat, f)))
  }

})


test_that("fsummarise works like summarise with across and grouped usage", {

  # Simple usage
  expect_true(all_obj_equal(fsummarise(gmtc, across(hp:drat, sum)),
                            fsummarise(gmtc, across(hp:drat, fsum)),
                            summarise(gmtc, across(hp:drat, sum), .groups = "drop")))

  expect_true(all_obj_equal(fsummarise(gmtc, across(5:7, sum)),
                            fsummarise(gmtc, across(5:7, fsum)),
                            summarise(gmtc, across(4:6, sum), .groups = "drop")))

  expect_true(all_obj_equal(fsummarise(gwld, across(is.numeric, sum, na.rm = TRUE)) |> setLabels(NULL),
                            fsummarise(gwld, across(is.numeric, fsum)) |> replace_NA() |> setLabels(NULL),
                            summarise(gwld, across(where(is.numeric), sum, na.rm = TRUE))))

  expect_true(all_obj_equal(fsummarise(gmtc, across(NULL, sum, na.rm = TRUE)) |> setLabels(NULL),
                            fsummarise(gmtc, across(NULL, fsum)) |> setLabels(NULL),
                            summarise(gmtc, across(everything(), sum, na.rm = TRUE), .groups = "drop")))

  expect_equal(fsummarise(gmtc, across(NULL, sum, na.rm = TRUE), keep.group_vars = FALSE),
               fsummarise(gmtc, across(NULL, sum, na.rm = TRUE)) |> slt(-cyl,-vs,-am))

  expect_equal(fsummarise(gmtc, across(cyl:vs, sum)),
               fsummarise(gmtc, cyl = sum(cyl), across(disp:qsec, fsum), vs = fsum(vs)))

  # Simple programming use
  vlist <- list(mtc |> fselect(hp:drat, return = "names"), NULL) # -(5:8), is.numeric
  for(i in seq_along(vlist)) {
    expect_true(all_obj_equal(fsummarise(gmtc, across(vlist[[i]], sum)),
                              fsummarise(gmtc, across(vlist[[i]], fsum)),
                              summarise(gmtc, across(if(is.null(vlist[[i]])) everything() else vlist[[i]], sum), .groups = "drop")))
    v <- vlist[[i]]
    expect_true(all_obj_equal(fsummarise(gmtc, across(v, sum)),
                              fsummarise(gmtc, across(v, fsum)),
                              summarise(gmtc, across(if(is.null(v)) everything() else v, sum), .groups = "drop")))
  }

  # Simple usage and multiple functions
  expect_true(all_obj_equal(fsummarise(gmtc, across(hp:drat, list(mean, sum))),
                            fsummarise(gmtc, across(hp:drat, list(mean = fmean, sum = fsum))),
                            summarise(gmtc, across(hp:drat, list(mean = mean, sum = sum)), .groups = "drop")))

  expect_true(all_obj_equal(fsummarise(gmtc, across(NULL, list(mean, sum))),
                            fsummarise(gmtc, across(NULL, list(mean = fmean, sum = fsum))),
                            summarise(gmtc, across(everything(), list(mean = mean, sum = sum)), .groups = "drop")))

  # Passing additional arguments
  expect_true(all_obj_equal(fsummarise(gwld, across(c("PCGDP", "LIFEEX"), sum, na.rm = TRUE))  |> setLabels(NULL),
                            fsummarise(gwld, across(c("PCGDP", "LIFEEX"), fsum, na.rm = TRUE))  |> setLabels(NULL) |> replace_NA(),
                            summarise(gwld, across(c("PCGDP", "LIFEEX"), sum, na.rm = TRUE), .groups = "drop")))

  expect_true(all_obj_equal(fsummarise(gmtc, across(hp:drat, weighted.mean, w = wt)),
                            fsummarise(gmtc, across(hp:drat, fmean, w = wt)),
                            summarise(gmtc, across(hp:drat, weighted.mean, w = wt), .groups = "drop")))

  expect_equal(fsummarise(gmtc, across(cyl:vs, weighted.mean, w = wt)),
               fsummarise(gmtc, cyl = weighted.mean(cyl, wt), across(disp:qsec, fmean, w = wt), vs = fmean(vs, wt)))

  expect_true(all_obj_equal(fsummarise(gmtc, across(hp:drat, list(mean = weighted.mean, sum = fsum), w = wt)),
                            fsummarise(gmtc, across(hp:drat, list(mean = fmean, sum = fsum), w = wt)),
                            summarise(gmtc, across(hp:drat, list(mean = weighted.mean, sum = fsum), w = wt), .groups = "drop")))

  # Simple programming use
  flist <- list(sum, list(mean = mean, sum = sum), list(mean, sum)) # c("mean", "sum"), c(mean = "fmean", sum = "fsum")
  for(i in seq_along(flist)) {
    expect_equal(fsummarise(gmtc, across(hp:drat, flist[[i]])),
                 summarise(gmtc, across(hp:drat, flist[[i]]), .groups = "drop"))
    f <- flist[[i]]
    expect_equal(fsummarise(gmtc, across(hp:drat, f)),
                 summarise(gmtc, across(hp:drat, f), .groups = "drop"))
  }

})


test_that("fsummarise miscellaneous things", {

  expect_equal(smr(gmtc, acr(disp:hp, c("mean", "sd"))),
               fsummarise(gmtc, across(disp:hp, c("mean", "sd"), .transpose = FALSE)) |>
                 colorderv(c(4,6,5,7), pos = "exchange"))

  expect_identical(names(smr(gmtc, acr(disp:hp, fmean, .names = TRUE)))[4:5], c("disp_fmean", "hp_fmean"))
  expect_identical(names(smr(gmtc, acr(disp:hp, mean, .names = TRUE)))[4:5], c("disp_mean", "hp_mean"))

  pwcorDF <- function(x, w = NULL) qDF(pwcor(x, w = w), "var")
  expect_equal(
    mtcars |> gby(cyl) |> smr(acr(disp:hp, pwcorDF, .apply = FALSE)),
    rsplit(mtcars, disp + hp ~ cyl) %>% lapply(pwcorDF) %>% unlist2d("cyl", "var") %>% tfm(cyl = as.numeric(cyl))
  )

  expect_equal(
    mtcars |> gby(cyl) |> smr(acr(disp:hp, pwcorDF, w = wt, .apply = FALSE)),
    rsplit(mtcars, disp + hp + wt ~ cyl) %>% lapply(function(x) pwcorDF(gv(x, 1:2), w = x$wt)) %>%
      unlist2d("cyl", "var") %>% tfm(cyl = as.numeric(cyl))
  )

  lmest <- function(x) list(Mods = list(lm(disp~., x)))
  expect_equal(
    qDT(mtcars) |> gby(cyl) |> smr(acr(disp:hp, lmest, .apply = FALSE)),
    qDT(mtcars) |> rsplit(disp + hp ~ cyl) %>% lapply(lmest) %>% data.table::rbindlist(idcol = "cyl") %>%
      tfm(cyl = as.numeric(cyl))
  )

})



test_that("fmutate works as intended for simple usage", {

  expect_equal(fmutate(mtc, mu = mean(mpg)), mutate(mtc, mu = mean(mpg)))
  expect_equal(fmutate(mtc, mu = mean(mpg), mpg = NULL), mutate(mtc, mu = mean(mpg), mpg = NULL))
  expect_equal(fmutate(mtc, mu = mean(mpg), dmu = mpg - mu), mutate(mtc, mu = mean(mpg), dmu = mpg - mu))
  expect_equal(fmutate(mtc, mu = log(mpg)), mutate(mtc, mu = log(mpg)))
  expect_equal(fmutate(mtc, mu = log(mpg), dmu = mpg - mu), mutate(mtc, mu = log(mpg), dmu = mpg - mu))

  expect_true(all_obj_equal(
    mutate(mtc, dmu = mpg - mean(mpg)),
    fmutate(mtc, dmu = mpg - mean(mpg)),
    fmutate(mtc, dmu = mpg - fmean(mpg)),
    fmutate(mtc, dmu = fmean(mpg, TRA = "-")),
    fmutate(mtc, dmu = fwithin(mpg))
  ))

  # With groups:
  expect_equal(fmutate(gmtc, mu = mean(mpg)), mutate(gmtc, mu = mean(mpg)))
  expect_equal(fmutate(gmtc, mu = mean(mpg), mpg = NULL), mutate(gmtc, mu = mean(mpg), mpg = NULL))
  expect_equal(fmutate(gmtc, mu = mean(mpg), dmu = mpg - mu), mutate(gmtc, mu = mean(mpg), dmu = mpg - mu))
  expect_equal(fmutate(gmtc, mu = log(mpg)), mutate(gmtc, mu = log(mpg)))
  expect_equal(fmutate(gmtc, mu = log(mpg), dmu = mpg - mu), mutate(gmtc, mu = log(mpg), dmu = mpg - mu))

  expect_true(all_obj_equal(
    mutate(gmtc, dmu = mpg - mean(mpg)),
    fmutate(gmtc, dmu = mpg - mean(mpg)),
    fmutate(gmtc, dmu = mpg - fmean(mpg)),
    fmutate(gmtc, dmu = fmean(mpg, TRA = "-")),
    fmutate(gmtc, dmu = fwithin(mpg))
  ))

})

test_that("fmutate with across works like ftransformv", {

   expect_true(all_obj_equal(

     ftransformv(mtcars, cyl:vs, fwithin, w = wt, apply = TRUE),
     ftransformv(mtcars, cyl:vs, fwithin, w = wt, apply = FALSE),
     fmutate(mtcars, across(cyl:vs, fwithin, w = wt, .apply = TRUE)),
     fmutate(mtcars, across(cyl:vs, fwithin, w = wt, .apply = FALSE))

   ))


  expect_true(all_obj_equal(

    ftransformv(mtcars, cyl:vs, fwithin, g = list(cyl, vs, am), apply = TRUE) |> setRownames(),
    ftransformv(mtcars, cyl:vs, fwithin, g = list(cyl, vs, am), apply = FALSE) |> setRownames(),
    fmutate(gmtc, across(cyl:vs, fwithin, .apply = TRUE)) |> qDF(),
    fmutate(gmtc, across(cyl:vs, fwithin, .apply = FALSE)) |> qDF(),
    fmutate(gmtc, across(cyl:vs, fmean, TRA = "-", .apply = TRUE)) |> qDF(),
    fmutate(gmtc, across(cyl:vs, fmean, TRA = "-", .apply = FALSE)) |> qDF(),
    fmutate(gmtc, across(cyl:vs, function(x) x - mean(x), .apply = TRUE)) |> qDF(),
    fmutate(gmtc, across(cyl:vs, function(x) lapply(x, function(y) y - mean(y)), .apply = FALSE)) |> qDF()

  ))

  expect_true(all_obj_equal(

    ftransformv(mtcars, cyl:vs, fwithin, g = list(cyl, vs, am), w = wt, apply = TRUE) |> setRownames(),
    ftransformv(mtcars, cyl:vs, fwithin, g = list(cyl, vs, am), w = wt, apply = FALSE) |> setRownames(),
    fmutate(gmtc, across(cyl:vs, fwithin, w = wt, .apply = TRUE)) |> qDF(),
    fmutate(gmtc, across(cyl:vs, fwithin, w = wt, .apply = FALSE)) |> qDF(),
    fmutate(gmtc, across(cyl:vs, fmean, TRA = "-", w = wt, .apply = TRUE)) |> qDF(),
    fmutate(gmtc, across(cyl:vs, fmean, TRA = "-", w = wt, .apply = FALSE)) |> qDF(),
    fmutate(gmtc, across(cyl:vs, function(x, w) x - weighted.mean(x, w), w = wt, .apply = TRUE)) |> qDF()

  ))

})


test_that("fmutate miscellaneous", {

  expect_true(length(fmutate(mtcars, across(cyl:vs, W, w = wt, .names = NULL))) > 15)
  expect_true(length(fmutate(mtcars, across(cyl:vs, list(D, W), .names = FALSE, .transpose = FALSE))) > 15)

  expect_equal(  fmutate(mtcars, across(cyl:vs, L, stubs = FALSE)),
                 fmutate(mtcars, across(cyl:vs, flag))
  )

  expect_true(length(fmutate(mtcars, across(cyl:vs, L))) > length(fmutate(mtcars, across(cyl:vs, flag))))

  expect_equal(
    fmutate(mtcars, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "used"),
    mutate(mtcars, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "used")
  )

  expect_equal(
    fmutate(mtcars, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "unused"),
    mutate(mtcars, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "unused")
  )

  expect_equal(
    fmutate(mtcars, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "none"),
    mutate(mtcars, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "none")
  )

  expect_identical(names(fmutate(mtcars, a = mpg, b = a, c = cyl, hp = wt, .keep = "unused")), c(setdiff(names(mtcars), .c(mpg, cyl, wt)), letters[1:3]))
  expect_identical(names(fmutate(mtcars, a = mpg, b = a, c = cyl, hp = wt, .keep = "none")), c("a", "b", "c", "hp"))

  expect_equal(
    fmutate(gmtc, a = fmax(mpg), b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "used"),
    mutate(gmtc, a = max(mpg), b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "used")
  )

  expect_equal(
    fmutate(gmtc, a = fmax(mpg), b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "unused"),
    mutate(gmtc, a = max(mpg), b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "unused")
  )

  expect_equal(
    fmutate(gmtc, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "none"),
    mutate(gmtc, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, .keep = "none")
  )

  # Inconsistent with the above and also inefficient...
  # expect_equal(
  #   fmutate(gmtc, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, cyl = cyl, .keep = "none"),
  #   mutate(gmtc, a = mpg, b = a + hp + disp, c = cyl, hp = wt + carb, cyl = cyl, .keep = "none")
  # )

  expect_equal(flast(names(fmutate(mtcars, across(cyl:vs, function(x) list(ps = kit::psum(x)), .apply = FALSE)))), "ps")

  expect_equal(
    fmutate(mtcars, across(cyl:vs, data.table::shift, .apply = FALSE, .names = FALSE)),
    fmutate(mtcars, across(cyl:vs, data.table::shift))
  )

})

if(FALSE) {


  fmutate(mtcars, across(cyl:vs, list(D, W), .names = TRUE, .transpose = FALSE)) |> head(3)
  fmutate(mtcars, across(cyl:vs, list(D, W), .names = TRUE, .transpose = TRUE)) |> head(3)

  fmutate(mtcars, across(cyl:vs, list(D, W), .names = TRUE, .apply = FALSE, .transpose = FALSE)) |> head(3)
  fmutate(mtcars, across(cyl:vs, list(D, W), .names = FALSE, .apply = FALSE, .transpose = FALSE)) |> head(3)

  fmutate(mtcars, across(cyl:vs, list(W, kit::psum), w = wt)) |> head(3)
  fmutate(mtcars, across(cyl:vs, kit::psum)) |> head(3)

  fmutate(mtcars, across(cyl:vs, identity, .apply = FALSE)) # 51 microesecond median on windows
  fmutate(mtcars, across(cyl:vs, identity)) # 62 microesecond median on windows

  fmutate(mtcars, across(cyl:vs, L))

  # TODO: Test all potential issues with environemtns etc. See if there are smarter ways to
  # incorporate internal functions, data and objects in the global environment.


}
