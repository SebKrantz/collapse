context("collapse and data.table integration")

if(!is.null(attributes(identical(FALSE, TRUE)))) stop("OECD label issue")

bmean <- base::mean

# TODO: Check memory allocation, particularly where names<- and attr<- are used.
# Also check attribute handling helpers with atomic and S4 objects !!
expect_equal(1, 1)

if(requireNamespace("data.table", quietly = TRUE) && requireNamespace("magrittr", quietly = TRUE)) {

options(warn = -1L)
library(data.table)
library(magrittr)
mtcDT <- qDT(roworderv(mtcars))
irisDT <- qDT(ss(iris, 1:100))
n <- 5L
# copy <- identity

# assignInNamespace("cedta.override", c(data.table:::cedta.override, "collapse"), "data.table")
assignInNamespace("cedta.override", "collapse", "data.table")

options(warn = 1L)

test_that("creating columns and printing works after passing a data.table through collapse functions", {

  expect_true(is.data.table(mtcDT))
  expect_true(is.data.table(irisDT))

  expect_output(print(mtcDT))
  expect_identical(names(mtcDT), names(mtcars))
  expect_silent(mtcDT[, col := 1])
  expect_output(print(mtcDT))
  expect_silent(mtcDT[, col := NULL])
  expect_identical(names(mtcDT), names(mtcars))
  expect_output(print(mtcDT))
  expect_silent(irisDT[, col := 1])
  expect_silent(irisDT[, col := NULL])

  # Statistical functions give warning
  dt <- fscale(copy(mtcDT))
  expect_warning(dt[, new := 1])
  expect_output(print(dt))

  dt <- fsum(copy(mtcDT), TRA = 1)
  expect_warning(dt[, new := 1])
  expect_output(print(dt))

  dt <- fsum(copy(mtcDT), drop = FALSE)
  expect_warning(dt[, new := 1])
  expect_output(print(dt))

  for(i in 1:n) {
    if(!identical(copy, identity)) mtcDT <- qDT(mtcDT)
    expect_silent(mtcDT[, col := 1])
    expect_silent(mtcDT[, col := NULL])
    expect_identical(names(mtcDT), names(mtcars))
    expect_identical(length(mtcDT), length(mtcars))
  }

  # Other functions should work:
  for(i in 1:n) {
  dt <- fgroup_by(mtcDT, cyl)
  expect_identical(names(dt), names(mtcars))
  # print(ltl(dt))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  # print(ltl(dt))
  }

  for(i in 1:n) {
  dt2 <- fgroup_vars(dt)
  expect_silent(dt2[, new := 1])
  expect_output(print(dt2))
  }

  for(i in 1:n) {
    dt <- fungroup(fgroup_by(mtcDT, c(2,8:9)))
    expect_identical(names(dt), names(mtcars))
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- funique(copy(mtcDT))
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- funique(copy(mtcDT), cols = "cyl")
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- fselect(copy(mtcDT), -mpg, -hp)
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- fselect(copy(mtcDT), col2 = disp, wt:carb)
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  fselect(dt, col2, new) <- NULL
  expect_silent(dt[, ncol := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- fsubset(copy(mtcDT), cyl == 4)
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- fsubset(copy(mtcDT), cyl == 4, bla = mpg, vs:am)
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- copy(mtcDT) %>% smr(mean_mpg = fmean(mpg))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- copy(mtcDT) %>% smr(mean_mpg = bmean(mpg))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- copy(mtcDT) %>% gby(cyl) %>% smr(mean_mpg = fmean(mpg))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- copy(mtcDT) %>% gby(cyl) %>% smr(mean_mpg = bmean(mpg))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- ftransform(copy(mtcDT), bla = 1)
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  settransform(dt, bla2 = 1)
  expect_silent(dt[, new2 := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  ftransform(dt) <- list(sds = mtcDT$qsec)
  expect_silent(dt[, new3 := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- fcompute(copy(mtcDT), bla = mpg + cyl, df = 1, keep = 7:10)
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- roworderv(copy(mtcDT))
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- roworder(copy(mtcDT), cyl, -vs)
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- roworderv(copy(mtcDT), cols = 1:2)
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- colorderv(copy(mtcDT))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- colorder(copy(mtcDT), vs, cyl, am)
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- frename(copy(mtcDT), carb = bla, mpg = x)
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- frename(copy(mtcDT), toupper)
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  setrename(dt, MPG = ABC, new = NEW)
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- get_vars(copy(irisDT), 1:3)
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  get_vars(dt, 1) <- irisDT$Species
  expect_silent(dt[, new2 := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  get_vars(dt, 1) <- NULL
  expect_silent(dt[, new3 := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- get_vars(irisDT, 1:3) %>% add_vars(gv(irisDT, 4))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  add_vars(dt) <- list(Sp = irisDT$Species)
  expect_silent(dt[, new2 := 1])
  expect_output(print(dt))
  }

  wldDT <- qDT(wlddev)
  for(i in .c(num_vars, nv, cat_vars, char_vars, fact_vars, logi_vars, date_vars)) {
    # print(i)
    # Iris data
    FUN <- match.fun(i)
    dt <- FUN(irisDT)
    expect_identical(names(dt), FUN(iris, "names"))
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
    rm(dt)

    dt <- irisDT
    eval(substitute(FUN(dt) <- NULL, list(FUN = as.name(i))))
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
    rm(dt)

    # wlddev data
    dt <- FUN(wldDT)
    expect_identical(names(dt), FUN(wlddev, "names"))
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
    rm(dt)

    dt <- wldDT
    eval(substitute(FUN(dt) <- NULL, list(FUN = as.name(i))))
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
    rm(dt)
  }

  for(i in 1:n) {
    dt <- relabel(copy(wldDT), toupper)
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
    setrelabel(dt, PCGDP = "GRP per cap", LIFEEX = "LE")
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
  }


  for(i in 1:n) {
  dt <- qDT(qTBL(qDF(qDT(GGDC10S))))
  expect_identical(names(dt), names(GGDC10S))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
    dt <- fdroplevels(copy(wldDT))
    expect_identical(names(dt), names(wlddev))
    expect_true(!anyNA(vlabels(dt)))
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
  }

  for(i in 1:n) {
  m <- qM(mtcars)
  dt <- qDT(m)
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  expect_output(print(mtcDT[, qDT(pwcor(.SD)), by = cyl, .SDcols = c("mpg", "hp", "carb")]))
  expect_output(print(melt(qDT(GGDC10S)[, qDT(pwcor(.SD)), by = .(Variable, Country), .SDcols = 6:15], 1:2)))

  for(i in 1:n) {
  dt <- as_character_factor(wldDT)
  expect_identical(names(dt), names(wlddev))
  expect_true(!anyNA(vlabels(dt)))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- as_character_factor(wldDT, keep.attr = FALSE)
  expect_identical(names(dt), names(wlddev))
  expect_true(anyNA(vlabels(dt)))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  options(warn = -1L)
  for(i in 1:n) {
  dt <- as_numeric_factor(wldDT)
  expect_identical(names(dt), names(wlddev))
  expect_true(!anyNA(vlabels(dt)))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- as_numeric_factor(wldDT, keep.attr = FALSE)
  expect_identical(names(dt), names(wlddev))
  expect_true(anyNA(vlabels(dt)))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }
  options(warn = 1L)

  for(i in 1:n) {
  dt <- collap(wldDT, ~ iso3c)
  expect_identical(names(dt), names(wlddev))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- collapv(wldDT, 1)
  expect_identical(names(dt), names(wlddev))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- collapg(gby(wldDT, 1))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- dapply(copy(mtcDT), log)
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- dapply(copy(mtcDT), log, return = "data.frame")
  expect_identical(names(dt), names(mtcars))
  expect_error(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  l <- rsplit(copy(mtcDT), ~cyl)
  expect_silent(for(i in seq_along(l)) l[[i]][, new := 1])
  expect_output(print(l))
  expect_output(print(l[[1]]))
  }

  for(i in 1:n) {
  dt <- unlist2d(l, DT = TRUE)
  expect_silent(dt[, new45 := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- na_omit(copy(mtcDT), cols = 1:2)
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- na_omit(copy(mtcDT))
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- na_insert(copy(mtcDT))
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- copy(wldDT)
  vlabels(wldDT) <- NULL
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- copy(mtcDT) %>% add_stub("B")
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- copy(mtcDT) %>% add_stub("B") %>% rm_stub("B")
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- copy(mtcDT) %>% setRownames
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
  dt <- copy(mtcDT) %>% frename(toupper) %>% setColnames(names(mtcars))
  expect_identical(names(dt), names(mtcars))
  expect_silent(dt[, new := 1])
  expect_output(print(dt))
  }

  for(i in 1:n) {
    dt <- replace_NA(copy(wldDT), cols = is.numeric)
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
  }

  for(i in 1:n) {
    dt <- replace_NA(copy(mtcDT), set = TRUE, cols = is.numeric)
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
  }

  for(i in 1:n) {
    dt <- replace_Inf(copy(wldDT))
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
  }

  for(i in 1:n) {
    dt <- replace_outliers(copy(wldDT), 3)
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
  }

  for(i in 1:n) {
    dt <- recode_num(copy(wldDT), `1` = 2)
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
  }

  for(i in 1:n) {
    dt <- recode_char(copy(wldDT), Uganda = "UGA")
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
  }

  for(i in 1:n) {
    dt <- pad(copy(mtcDT), 1:3)
    expect_silent(dt[, new := 1])
    expect_output(print(dt))
  }

})

}

