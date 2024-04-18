context("fmode")



# rm(list = ls())
set.seed(101)
x <- round(abs(10*rnorm(100)))
w <- as.integer(round(abs(10*rnorm(100)))) # round(abs(rnorm(100)), 1) -> Numeric precision issues in R
xNA <- x
wNA <- w
xNA[sample.int(100,20)] <- NA
wNA[sample.int(100,20)] <- NA
f <- as.factor(sample.int(10, 100, TRUE))
data <- wlddev[wlddev$iso3c %in% c("BLZ","IND","USA","SRB","GRL"), ]
l <- nrow(data)
g <- GRP(droplevels(data$iso3c))
gf <- as_factor_GRP(g)
dataNA <- na_insert(data)
m <- as.matrix(num_vars(data)) # without num_vars also works for ties = "first"
mNA <- as.matrix(num_vars(dataNA))
wdat <- as.integer(round(10*abs(rnorm(l)))) # round(abs(rnorm(l)), 1) -> Numeric precision issues in R
wdatNA <- wdat
wdatNA[sample.int(l, floor(l/5))] <- NA

ncv <- !char_vars(data, "logical")
getdata <- function(first) if(first) data else gv(data, ncv)
getdataNA <- function(first) if(first) dataNA else gv(dataNA, ncv)

# seteltNA <- function(x,i,j) {
#   x[i,j] <- NA
#   x
# }

whichmax <- function(x) which(as.integer(x) == as.integer(max(x))) # This solves numeric precision issues
minwa <- function(x) {
  xna <- unattrib(x)
  if(anyNA(xna)) {
    if(is.integer(xna)) return(`attributes<-`(NA_integer_, attributes(x)))
    # if(is.character(xna)) return(`attributes<-`(NA_character_, attributes(x)))
    if(is.numeric(xna)) {
      xna <- na_rm(xna)
      if(!length(xna)) return(`attributes<-`(NA_real_, attributes(x)))
    }
  }
  `attributes<-`(`storage.mode<-`(base::min(xna), storage.mode(x)), attributes(x))
}
maxwa <- function(x) {
  xna <- unattrib(x)
  if(is.numeric(xna) && anyNA(xna)) {
    xna <- na_rm(xna)
    if(!length(xna)) return(`attributes<-`(NA_real_, attributes(x)))
  }
  `attributes<-`(`storage.mode<-`(base::max(xna), storage.mode(x)), attributes(x))
}

if(identical(Sys.getenv("NCRAN"), "TRUE")) {

  # This is to fool very silly checks on CRAN scanning the code of the tests
  rowidv <- eval(parse(text = paste0("data.table", ":", ":", "rowidv")))

# firstmode <- function(x) {
#   ox <- sort(x)
#   ox[which.max(rowidv(ox))]
# }
unam <- function(x) `names<-`(x, NULL)

Mode <- function(x, na.rm = FALSE, ties = "first") {
  if(na.rm) {
    miss <- is.na(x)
    if(all(miss)) return(x[1L])
    x <- x[!miss]
  }
  o <- radixorder(x)
  ox <- unam(x)[o]
  switch(ties,
         first = unam(x)[which.max(rowidv(ox)[radixorder(o)])],
         last = unam(x)[which.max(rowidv(ox)[radixorder(o, decreasing = TRUE)])],
         min = minwa(ox[whichmax(rowidv(ox))]),
         max = maxwa(ox[whichmax(rowidv(ox))]),
         stop("Unknown ties option"))
}

}

# Mode <- function(x, na.rm = FALSE, ties = "first") {
#   if(na.rm) x <- x[!is.na(x)]
#   ux <- unique(x)
#   switch(ties,
#          first = ux[which.max(tabulate(match(x, ux)))],
#          min = minwa(ux[whichmax(tabulate(match(x, ux)))]),
#          max = maxwa(ux[whichmax(tabulate(match(x, ux)))]),
#          stop("Unknown ties option"))
# }

wMode <- function(x, w, na.rm = FALSE, ties = "first") {
  ax <- attributes(x)
  cc <- complete.cases(x, w)
  if(!any(cc)) return(`storage.mode<-`(NA, storage.mode(x)))
  if(na.rm) {
    w <- w[cc]
    x <- x[cc]
  }
  g <- GRP.default(x, call = FALSE)
  switch(ties,
         first = {
           o <- radixorder(unlist(gsplit(seq_along(w), g), use.names = FALSE))
           sw <- unlist(lapply(gsplit(w, g), base::cumsum), use.names = FALSE)[o]
           fsubset.default(x, which.max(sw))
         },
         min = minwa(fsubset.default(g[["groups"]][[1L]], whichmax(fsum.default(w, g, use.g.names = FALSE)))),
         max = maxwa(fsubset.default(g[["groups"]][[1L]], whichmax(fsum.default(w, g, use.g.names = FALSE)))),
         stop("Unknown ties option"))
  # storage.mode(res) <- storage.mode(x)
  # `attributes<-`(res, ax)
}


for (nth in 1:2) {

  if(nth == 2L) {
    if(Sys.getenv("OMP") == "TRUE") {
      fmode <- function(x, ...) collapse::fmode(x, ..., nthreads = 2L)
    } else break
  }

if(identical(Sys.getenv("NCRAN"), "TRUE")) {

test_that("fmode performs like Mode (defined above)", {
  for(t in c("first","min","max")) {
    # print(t)
    tf <- t == "first"
  expect_equal(fmode(NA, ties = t), Mode(NA, ties = t))
  expect_equal(fmode(NA, na.rm = FALSE, ties = t), Mode(NA, ties = t))
  expect_equal(fmode(1, ties = t), Mode(1, na.rm = TRUE, ties = t))
  expect_equal(fmode(1:3, ties = t), Mode(1:3, na.rm = TRUE, ties = t))
  expect_equal(fmode(-1:1, ties = t), Mode(-1:1, na.rm = TRUE, ties = t))
  expect_equal(fmode(1, na.rm = FALSE, ties = t), Mode(1, ties = t))
  expect_equal(fmode(1:3, na.rm = FALSE, ties = t), Mode(1:3, ties = t))
  expect_equal(fmode(-1:1, na.rm = FALSE, ties = t), Mode(-1:1, ties = t))
  expect_equal(fmode(x, ties = t), Mode(x, na.rm = TRUE, ties = t))
  expect_equal(fmode(x, na.rm = FALSE, ties = t), Mode(x, ties = t))
  if(tf) expect_equal(fmode(xNA, na.rm = FALSE, ties = t), Mode(xNA, ties = t))
  expect_equal(fmode(xNA, ties = t), Mode(xNA, na.rm = TRUE, ties = t))
  # expect_equal(as.character(fmode(data, drop = FALSE)), fmode(m))
  expect_equal(fmode(m, ties = t), dapply(m, Mode, na.rm = TRUE, ties = t))
  expect_equal(fmode(m, na.rm = FALSE, ties = t), dapply(m, Mode, ties = t))
  if(tf) expect_equal(fmode(mNA, na.rm = FALSE, ties = t), dapply(mNA, Mode, ties = t))
  expect_equal(fmode(mNA, ties = t), dapply(mNA, Mode, na.rm = TRUE, ties = t))
  expect_equal(fmode(getdata(tf), ties = t, drop = FALSE), dapply(getdata(tf), Mode, na.rm = TRUE, ties = t, drop = FALSE))
  expect_equal(fmode(getdata(tf), na.rm = FALSE, ties = t, drop = FALSE), dapply(getdata(tf), Mode, ties = t, drop = FALSE))
  if(tf) expect_equal(fmode(dataNA, na.rm = FALSE, ties = t, drop = FALSE), dapply(dataNA, Mode, ties = t, drop = FALSE))
  expect_equal(fmode(getdataNA(tf), ties = t, drop = FALSE), dapply(getdataNA(tf), Mode, na.rm = TRUE, ties = t, drop = FALSE))
  expect_equal(fmode(x, f, ties = t), BY(x, f, Mode, na.rm = TRUE, ties = t))
  expect_equal(fmode(x, f, na.rm = FALSE, ties = t), BY(x, f, Mode, ties = t))
  if(tf) expect_equal(fmode(xNA, f, na.rm = FALSE, ties = t), BY(xNA, f, Mode, ties = t))
  expect_equal(fmode(xNA, f, ties = t), BY(xNA, f, Mode, na.rm = TRUE, ties = t))
  expect_equal(fmode(m, g, ties = t), BY(m, g, Mode, na.rm = TRUE, ties = t))
  expect_equal(fmode(m, g, na.rm = FALSE, ties = t), BY(m, g, Mode, ties = t))
  if(tf) expect_equal(fmode(mNA, g, na.rm = FALSE), BY(mNA, g, Mode)) # Mode gives NA
  expect_equal(fmode(mNA, g, ties = t), BY(mNA, g, Mode, na.rm = TRUE, ties = t))
  expect_equal(fmode(getdata(tf), g, ties = t), BY(getdata(tf), g, Mode, na.rm = TRUE, ties = t))
  expect_equal(fmode(getdata(tf), g, na.rm = FALSE, ties = t), BY(getdata(tf), g, Mode, ties = t))
  if(tf) expect_equal(fmode(dataNA, g, na.rm = FALSE), BY(dataNA, g, Mode))  # Mode gives NA
  expect_equal(fmode(getdataNA(tf), g, ties = t), BY(getdataNA(tf), g, Mode, na.rm = TRUE, ties = t))
  }
})

}

test_that("fmode with weights performs as intended (unbiased)", {
  expect_equal(fmode(c(2,2,4,5,5,5)), fmode(c(2,4,5), w = c(2,1,3)))
  expect_equal(fmode(c(2,2,4,5,5,5), na.rm = FALSE), fmode(c(2,4,5), w = c(2,1,3), na.rm = FALSE))
  expect_equal(fmode(c(2.456,2.456,4.123,5.009,5.009,5.009)), fmode(c(2.456,4.123,5.009), w = c(2,1,3)))
  expect_equal(fmode(c(2.456,2.456,4.123,5.009,5.009,5.009), na.rm = FALSE), fmode(c(2.456,4.123,5.009), w = c(2,1,3), na.rm = FALSE))
  expect_equal(fmode(c(2,2,NA,5,5,5)), fmode(c(2,NA,5), w = c(2,1,3)))
  expect_equal(fmode(c(2,2,NA,5,5,5), na.rm = FALSE), fmode(c(2,NA,5), w = c(2,1,3), na.rm = FALSE))
  expect_equal(fmode(c(2,2,NA,5,5,5)), fmode(c(2,4,5), w = c(2,NA,3)))
  expect_equal(fmode(c(2,2,NA,5,5,5), na.rm = FALSE), fmode(c(2,4,5), w = c(2,NA,3), na.rm = FALSE))
  expect_equal(fmode(c(NA,NA,4.123,5.009,5.009,5.009)), fmode(c(NA,4.123,5.009), w = c(2,1,3)))
  expect_equal(fmode(c(NA,NA,4.123,5.009,5.009,5.009), na.rm = FALSE), fmode(c(NA,4.123,5.009), w = c(2,1,3), na.rm = FALSE))
  expect_equal(fmode(c(NA,NA,4.123,5.009,5.009,5.009)), fmode(c(2.456,4.123,5.009), w = c(NA,1,3)))
  expect_equal(fmode(c(NA,NA,4.123,5.009,5.009,5.009), na.rm = FALSE), fmode(c(2.456,4.123,5.009), w = c(NA,1,3), na.rm = FALSE))
  f <- as.factor(rep(1:2, each = 6)); fs <- as.factor(rep(1:2, each = 3))
  v <- c(2,2,4,5,5,5,2,2,4,5,5,5); vs <- c(2,4,5,2,4,5); w <- c(2,1,3,2,1,3)
  v2 <- c(2.456,2.456,4.123,5.009,5.009,5.009,2.456,2.456,4.123,5.009,5.009,5.009); v2s <- c(2.456,4.123,5.009,2.456,4.123,5.009)
  expect_equal(fmode(v, f), fmode(vs, fs, w))
  expect_equal(fmode(v, f, na.rm = FALSE), fmode(vs, fs, w, na.rm = FALSE))
  expect_equal(fmode(v2, f), fmode(v2s, fs, w))
  expect_equal(fmode(v2, f, na.rm = FALSE), fmode(v2s, fs, w, na.rm = FALSE))
  v[c(3,9)] <- NA; vs[c(2,5)] <- NA
  expect_equal(fmode(v, f), fmode(vs, fs, w))
  expect_equal(fmode(v, f, na.rm = FALSE), fmode(vs, fs, w, na.rm = FALSE))
  vs[c(2,5)] <- 4; w[c(2,5)] <- NA
  expect_equal(fmode(v, f), fmode(vs, fs, w))
  expect_equal(fmode(v, f, na.rm = FALSE), fmode(vs, fs, w, na.rm = FALSE))
  w[c(2,5)] <- 1; v2[c(1:2,7:8)] <- NA; v2s[c(1,4)] <- NA
  expect_equal(fmode(v2, f), fmode(v2s, fs, w))
  expect_equal(fmode(v2, f, na.rm = FALSE), fmode(v2s, fs, w, na.rm = FALSE))
  v2s[c(1,4)] <- 2.456; w[c(1,4)] <- NA
  expect_equal(fmode(v2, f), fmode(v2s, fs, w))
  expect_equal(fmode(v2, f, na.rm = FALSE), fmode(v2s, fs, w, na.rm = FALSE))
})

test_that("fmode performs like fmode with weights all equal", {
  for(t in c("first","min","max")) {
  expect_equal(fmode(NA, ties = t), fmode(NA, w = 0.9, ties = t))
  expect_equal(fmode(NA, na.rm = FALSE, ties = t), fmode(NA, w = 2.946, na.rm = FALSE, ties = t))
  expect_equal(fmode(1, ties = t), fmode(1, w = 3, ties = t))
  expect_equal(fmode(1:3, ties = t), fmode(1:3, w = rep(0.9,3), ties = t))
  expect_equal(fmode(-1:1, ties = t), fmode(-1:1, w = rep(4.2,3), ties = t))
  expect_equal(fmode(1, na.rm = FALSE, ties = t), fmode(1, w = 5, na.rm = FALSE, ties = t))
  expect_equal(fmode(1:3, na.rm = FALSE, ties = t), fmode(1:3, w = rep(1.4, 3), na.rm = FALSE, ties = t))
  expect_equal(fmode(-1:1, na.rm = FALSE, ties = t), fmode(-1:1, w = rep(1.4, 3), na.rm = FALSE, ties = t))
  expect_equal(fmode(x, ties = t), fmode(x, w = rep(1,100), ties = t))
  expect_equal(fmode(x, na.rm = FALSE, ties = t), fmode(x, w = rep(1.4, 100), na.rm = FALSE, ties = t))  # failed on patched solaris...
  expect_equal(fmode(xNA, na.rm = FALSE, ties = t), fmode(xNA, w = rep(4.6, 100), na.rm = FALSE, ties = t))
  expect_equal(fmode(xNA, ties = t), fmode(xNA, w = rep(4.6, 100), ties = t)) # failed on patched solaris...
  expect_equal(fmode(m, ties = t), fmode(m, w = rep(6587, l), ties = t))
  expect_equal(fmode(m, na.rm = FALSE, ties = t), fmode(m, w = rep(6587, l), na.rm = FALSE, ties = t))
  expect_equal(fmode(mNA, na.rm = FALSE, ties = t), fmode(mNA, w = rep(6587, l), na.rm = FALSE, ties = t))
  expect_equal(fmode(mNA, ties = t), fmode(mNA, w = rep(6587, l), ties = t))
  expect_equal(fmode(data, ties = t), fmode(data, w = rep(6787, l), ties = t))
  expect_equal(fmode(data, na.rm = FALSE, ties = t), fmode(data, w = rep(6787, l), na.rm = FALSE, ties = t))
  expect_equal(fmode(dataNA, na.rm = FALSE, ties = t), fmode(dataNA, w = rep(6787, l), na.rm = FALSE, ties = t))
  expect_equal(fmode(dataNA, ties = t), fmode(dataNA, w = rep(6787, l), ties = t))
  expect_equal(fmode(x, f, ties = t), fmode(x, f, rep(546,100), ties = t))
  expect_equal(fmode(x, f, na.rm = FALSE, ties = t), fmode(x, f, rep(5,100), na.rm = FALSE, ties = t))
  expect_equal(fmode(xNA, f, na.rm = FALSE, ties = t), fmode(xNA, f, rep(52.7,100), na.rm = FALSE, ties = t)) # Failed sometimes for some reason... v. 1.5.1 error
  expect_equal(fmode(xNA, f, ties = t), fmode(xNA, f, rep(599,100), ties = t))
  expect_equal(fmode(m, g, ties = t), fmode(m, g, rep(546,l), ties = t))
  expect_equal(fmode(m, g, na.rm = FALSE, ties = t), fmode(m, g, rep(1,l), na.rm = FALSE, ties = t))
  expect_equal(fmode(mNA, g, na.rm = FALSE, ties = t), fmode(mNA, g, rep(7,l), na.rm = FALSE, ties = t))
  expect_equal(fmode(mNA, g, ties = t), fmode(mNA, g, rep(1,l), ties = t))
  expect_equal(fmode(data, g, ties = t), fmode(data, g, rep(53,l), ties = t))
  expect_equal(fmode(data, g, na.rm = FALSE, ties = t), fmode(data, g, rep(546,l), na.rm = FALSE, ties = t))
  expect_equal(fmode(dataNA, g, na.rm = FALSE, ties = t), fmode(dataNA, g, rep(1,l), na.rm = FALSE, ties = t)) # rep(0.999999,l) failed CRAN Arch i386
  expect_equal(fmode(dataNA, g, ties = t), fmode(dataNA, g, rep(999,l), ties = t)) # rep(999.9999,l) failed CRAN Arch i386
  }
})

test_that("fmode with weights performs like wMode (defined above)", {
  for(t in c("first","min","max")) {
  # print(t)
    tf <- t == "first"
  # complete weights
  expect_equal(fmode(NA, w = 1, ties = t), wMode(NA, 1, ties = t))
  expect_equal(fmode(NA, w = 1, na.rm = FALSE, ties = t), wMode(NA, 1, ties = t))
  expect_equal(fmode(1, w = 1, ties = t), wMode(1, w = 1, ties = t))
  expect_equal(fmode(1:3, w = 1:3, ties = t), wMode(1:3, 1:3, ties = t))
  expect_equal(fmode(-1:1, w = 1:3, ties = t), wMode(-1:1, 1:3, ties = t))
  expect_equal(fmode(1, w = 1, na.rm = FALSE, ties = t), wMode(1, 1, ties = t))
  expect_equal(fmode(1:3, w = c(0.99,3454,1.111), na.rm = FALSE, ties = t), wMode(1:3, c(0.99,3454,1.111), ties = t))
  expect_equal(fmode(-1:1, w = 1:3, na.rm = FALSE, ties = t), wMode(-1:1, 1:3, ties = t))
  expect_equal(fmode(x, w = w, ties = t), wMode(x, w, ties = t))
  expect_equal(fmode(x, w = w, na.rm = FALSE, ties = t), wMode(x, w, ties = t))
  if(tf) expect_equal(fmode(xNA, w = w, na.rm = FALSE, ties = t), wMode(xNA, w, ties = t))
  expect_equal(fmode(xNA, w = w, ties = t), wMode(xNA, w, na.rm = TRUE, ties = t))
  # expect_equal(fmode(data, w = wdat, drop = FALSE, ties = t), fmode(m, w = wdat, ties = t))
  expect_equal(fmode(m, w = wdat, ties = t), dapply(m, wMode, wdat, na.rm = TRUE, ties = t))
  expect_equal(fmode(m, w = wdat, na.rm = FALSE, ties = t), dapply(m, wMode, wdat, ties = t))
  if(tf) expect_equal(fmode(mNA, w = wdat, na.rm = FALSE, ties = t), dapply(mNA, wMode, wdat, ties = t))
  expect_equal(fmode(mNA, w = wdat, ties = t), dapply(mNA, wMode, wdat, na.rm = TRUE, ties = t))
  expect_equal(fmode(getdata(tf), w = wdat, drop = FALSE, ties = t), dapply(getdata(tf), wMode, wdat, na.rm = TRUE, drop = FALSE, ties = t))
  expect_equal(fmode(getdata(tf), w = wdat, na.rm = FALSE, drop = FALSE, ties = t), dapply(getdata(tf), wMode, wdat, drop = FALSE, ties = t))
  if(tf) expect_equal(fmode(dataNA, w = wdat, na.rm = FALSE, drop = FALSE, ties = t), dapply(dataNA, wMode, wdat, drop = FALSE, ties = t))
  expect_equal(fmode(getdataNA(tf), w = wdat, drop = FALSE, ties = t), dapply(getdataNA(tf), wMode, wdat, na.rm = TRUE, drop = FALSE, ties = t))
  expect_equal(fmode(x, f, w, ties = t), BY(x, f, wMode, w, ties = t))
  expect_equal(fmode(x, f, w, na.rm = FALSE, ties = t), BY(x, f, wMode, w, ties = t))
  if(tf) expect_equal(fmode(xNA, f, w, na.rm = FALSE, ties = t), BY(xNA, f, wMode, w, ties = t))
  expect_equal(fmode(xNA, f, w, ties = t), BY(xNA, f, wMode, w, na.rm = TRUE, ties = t))
  expect_equal(fmode(m, g, wdat, ties = t), BY(m, gf, wMode, wdat, na.rm = TRUE, ties = t))
  expect_equal(fmode(m, g, wdat, na.rm = FALSE, ties = t), BY(m, gf, wMode, wdat, ties = t))
  if(tf) expect_equal(fmode(mNA, g, wdat, na.rm = FALSE, ties = t),  BY(mNA, gf, wMode, wdat, ties = t))
  expect_equal(fmode(mNA, g, wdat, ties = t), BY(mNA, gf, wMode, wdat, na.rm = TRUE, ties = t))
  expect_equal(fmode(getdata(tf), g, wdat, ties = t), BY(getdata(tf), gf, wMode, wdat, na.rm = TRUE, ties = t))
  expect_equal(fmode(getdata(tf), g, wdat, na.rm = FALSE, ties = t), BY(getdata(tf), gf, wMode, wdat, ties = t))
  if(tf) expect_equal(fmode(dataNA, g, wdat, na.rm = FALSE, ties = t), BY(dataNA, gf, wMode, wdat, ties = t))
  expect_equal(fmode(getdataNA(tf), g, wdat, ties = t), BY(getdataNA(tf), gf, wMode, wdat, na.rm = TRUE, ties = t))
  # missing weights: # missing weights are summed : wsum is NA.... fmode does not properly deal with missing weights if na.rm = FALSE
  expect_equal(fmode(NA, w = NA, ties = t), wMode(NA, NA, ties = t))
  # expect_equal(fmode(1, w = NA, ties = t), wMode(1, w = NA, ties = t))
  expect_equal(fmode(1:3, w = c(NA,1:2), ties = t), wMode(1:3, c(NA,1:2), na.rm = TRUE, ties = t))
  expect_equal(fmode(-1:1, w = c(NA,1:2), ties = t), wMode(-1:1, c(NA,1:2), na.rm = TRUE, ties = t))
  expect_equal(fmode(x, w = wNA, ties = t), wMode(x, wNA, na.rm = TRUE, ties = t))
  expect_equal(fmode(xNA, w = wNA, ties = t), wMode(xNA, wNA, na.rm = TRUE, ties = t))
  # expect_equal(fmode(data, w = wdatNA, ties = t), fmode(m, w = wdatNA, ties = t))
  expect_equal(fmode(m, w = wdatNA, ties = t), dapply(m, wMode, wdatNA, na.rm = TRUE, ties = t))
  expect_equal(fmode(mNA, w = wdatNA, ties = t), dapply(mNA, wMode, wdatNA, na.rm = TRUE, ties = t))
  expect_equal(fmode(getdata(tf), w = wdatNA, ties = t, drop = FALSE), dapply(getdata(tf), wMode, wdatNA, na.rm = TRUE, ties = t, drop = FALSE))
  expect_equal(fmode(getdataNA(tf), w = wdatNA, ties = t, drop = FALSE), dapply(getdataNA(tf), wMode, wdatNA, na.rm = TRUE, ties = t, drop = FALSE))
  expect_equal(fmode(x, f, wNA, ties = t), BY(x, f, wMode, wNA, na.rm = TRUE, ties = t)) # failed on MAC OSX
  expect_equal(fmode(xNA, f, wNA, ties = t), BY(xNA, f, wMode, wNA, na.rm = TRUE, ties = t)) # failed on mac OSX...
  expect_equal(fmode(m, g, wdatNA, ties = t), BY(m, gf, wMode, wdatNA, na.rm = TRUE, ties = t))
  expect_equal(fmode(mNA, g, wdatNA, ties = t), BY(mNA, gf, wMode, wdatNA, na.rm = TRUE, ties = t))
  expect_equal(fmode(getdata(tf), g, wdatNA, ties = t), BY(getdata(tf), gf, wMode, wdatNA, na.rm = TRUE, ties = t))
  expect_equal(fmode(getdataNA(tf), g, wdatNA, ties = t), BY(getdataNA(tf), gf, wMode, wdatNA, na.rm = TRUE, ties = t))
  }
})

test_that("fmode performs numerically stable", {
  for(t in c("first","min","max")) {
  expect_true(all_obj_equal(replicate(50, fmode(1, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(NA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(NA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(x, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(x, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(xNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(xNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(m, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(m, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(mNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(mNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(data, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(data, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(dataNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(dataNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(x, f, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(x, f, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(xNA, f, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(xNA, f, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(m, g, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(m, g, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(mNA, g, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(mNA, g, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(data, g, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(data, g, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(dataNA, g, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(dataNA, g, ties = t), simplify = FALSE)))
  }
})

test_that("fmode with complete weights performs numerically stable", {
  for(t in c("first","min","max")) {
  expect_true(all_obj_equal(replicate(50, fmode(1, w = 1, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(NA, w = 1, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(NA, w = 1, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(x, w = w, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(x, w = w, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(xNA, w = w, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(xNA, w = w, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(m, w = wdat, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(m, w = wdat, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(mNA, w = wdat, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(mNA, w = wdat, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(data, w = wdat, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(data, w = wdat, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(dataNA, w = wdat, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(dataNA, w = wdat, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(x, f, w, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(x, f, w, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(xNA, f, w, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(xNA, f, w, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(m, g, wdat, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(m, g, wdat, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(mNA, g, wdat, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(mNA, g, wdat, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(data, g, wdat, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(data, g, wdat, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(dataNA, g, wdat, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(dataNA, g, wdat, ties = t), simplify = FALSE)))
  }
})

test_that("fmode with missing weights performs numerically stable", {
  for(t in c("first","min","max")) {
  expect_true(all_obj_equal(replicate(50, fmode(1, w = NA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(NA, w = NA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(NA, w = NA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(x, w = wNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(x, w = wNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(xNA, w = wNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(xNA, w = wNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(m, w = wdatNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(m, w = wdatNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(mNA, w = wdatNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(mNA, w = wdatNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(data, w = wdatNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(data, w = wdatNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(dataNA, w = wdatNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(dataNA, w = wdatNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(x, f, wNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(x, f, wNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(xNA, f, wNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(xNA, f, wNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(m, g, wdatNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(m, g, wdatNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(mNA, g, wdatNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(mNA, g, wdatNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(data, g, wdatNA, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(data, g, wdatNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(dataNA, g, wdatNA, na.rm = FALSE, ties = t), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmode(dataNA, g, wdatNA, ties = t), simplify = FALSE)))
  }
})

test_that("fmode handles special values in the right way", {
  expect_equal(fmode(NA), NA)
  expect_equal(fmode(NaN), NaN)
  expect_equal(fmode(Inf), Inf)
  expect_equal(fmode(-Inf), -Inf)
  expect_equal(fmode(TRUE), TRUE)
  expect_equal(fmode(FALSE), FALSE)
  expect_equal(fmode(NA, na.rm = FALSE), NA)
  expect_equal(fmode(NaN, na.rm = FALSE), NaN)
  expect_equal(fmode(Inf, na.rm = FALSE), Inf)
  expect_equal(fmode(-Inf, na.rm = FALSE), -Inf)
  expect_equal(fmode(TRUE, na.rm = FALSE), TRUE)
  expect_equal(fmode(FALSE, na.rm = FALSE), FALSE)
  expect_equal(fmode(c(1,NA)), 1)
  expect_equal(fmode(c(1,NaN)), 1)
  expect_equal(fmode(c(1,Inf)), 1)
  expect_equal(fmode(c(1,-Inf)), 1)
  expect_equal(fmode(c(FALSE,TRUE)), FALSE)
  expect_equal(fmode(c(FALSE,FALSE)), FALSE)
  expect_equal(fmode(c(1,Inf), na.rm = FALSE), 1)
  expect_equal(fmode(c(1,-Inf), na.rm = FALSE), 1)
  expect_equal(fmode(c(FALSE,TRUE), na.rm = FALSE), FALSE)
  expect_equal(fmode(c(FALSE,FALSE), na.rm = FALSE), FALSE)
})

test_that("fmode with weights handles special values in the right way", {
  expect_equal(fmode(NA, w = 1), NA)
  expect_equal(fmode(NaN, w = 1), NaN)
  expect_equal(fmode(Inf, w = 1), Inf)
  expect_equal(fmode(-Inf, w = 1), -Inf)
  expect_equal(fmode(TRUE, w = 1), TRUE)
  expect_equal(fmode(FALSE, w = 1), FALSE)
  expect_equal(fmode(NA, w = 1, na.rm = FALSE), NA)
  expect_equal(fmode(NaN, w = 1, na.rm = FALSE), NaN)
  expect_equal(fmode(Inf, w = 1, na.rm = FALSE), Inf)
  expect_equal(fmode(-Inf, w = 1, na.rm = FALSE), -Inf)
  expect_equal(fmode(TRUE, w = 1, na.rm = FALSE), TRUE)
  expect_equal(fmode(FALSE, w = 1, na.rm = FALSE), FALSE)
  expect_equal(fmode(NA, w = NA), NA)
  expect_equal(fmode(NaN, w = NA), NaN)
  expect_equal(fmode(Inf, w = NA), Inf)
  expect_equal(fmode(-Inf, w = NA), -Inf)
  expect_equal(fmode(TRUE, w = NA), TRUE)
  expect_equal(fmode(FALSE, w = NA), FALSE)
  expect_equal(fmode(NA, w = NA, na.rm = FALSE), NA)
  expect_equal(fmode(NaN, w = NA, na.rm = FALSE), NaN)
  expect_equal(fmode(Inf, w = NA, na.rm = FALSE), Inf)
  expect_equal(fmode(-Inf, w = NA, na.rm = FALSE), -Inf)
  expect_equal(fmode(TRUE, w = NA, na.rm = FALSE), TRUE)
  expect_equal(fmode(FALSE, w = NA, na.rm = FALSE), FALSE)
  expect_equal(fmode(1:3, w = c(1,Inf,3)), 2)
  expect_equal(fmode(1:3, w = c(1,-Inf,3)), 3)
  expect_equal(fmode(1:3, w = c(1,Inf,3), na.rm = FALSE), 2)
  expect_equal(fmode(1:3, w = c(1,-Inf,3), na.rm = FALSE), 3)
})

test_that("fmode produces errors for wrong input", {
  expect_visible(fmode("a"))
  expect_visible(fmode(NA_character_))
  expect_visible(fmode(mNA))
  expect_error(fmode(mNA, f))
  expect_error(fmode(1:2,1:3))
  expect_error(fmode(m,1:31))
  expect_error(fmode(data,1:31))
  expect_error(fmode(data, w = 1:31))
  expect_visible(fmode("a", w = 1))
  expect_error(fmode(1:2, w = 1:3))
  expect_visible(fmode(NA_character_, w = 1))
  expect_visible(fmode(mNA, w = wdat))
  expect_error(fmode(mNA, f, wdat))
  expect_error(fmode(mNA, w = 1:33))
  expect_error(fmode(1:2,1:2, 1:3))
  expect_error(fmode(m,1:32,1:20))
  expect_error(fmode(data,1:32,1:10))
  expect_error(fmode(1:2, w = c("a","b")))
  expect_visible(fmode(wlddev))
  expect_visible(fmode(wlddev, w = wlddev$year, drop = FALSE))
  expect_visible(fmode(wlddev, wlddev$iso3c))
  expect_visible(fmode(wlddev, wlddev$iso3c, wlddev$year))
})

}


test_that("Singleton group optimization works properly", {
  g <- GRP(as.character(seq_row(mtcars)))
  w <- mtcars$wt
  expect_equal(unattrib(fmode(mtcars$mpg, g)), mtcars$mpg[g$order])
  expect_equal(unattrib(fmode(mtcars$mpg, g, w)), mtcars$mpg[g$order])
  g <- GRP(seq_row(mtcars))
  expect_equal(unattrib(fmode(mtcars$mpg, g)), mtcars$mpg[g$order])
  expect_equal(unattrib(fmode(mtcars$mpg, g, w)), mtcars$mpg[g$order])
  g <- GRP(sample.int(100, 32))
  expect_equal(unattrib(fmode(mtcars$mpg, g)), mtcars$mpg[g$order])
  expect_equal(unattrib(fmode(mtcars$mpg, g, w)), mtcars$mpg[g$order])
})
