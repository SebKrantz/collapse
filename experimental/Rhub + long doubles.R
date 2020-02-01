
library(rhub)
rhubcheck <- check_for_cran(path = "C:/Users/Sebastian Krantz/Documents/R/collapse_1.0.0.tar.gz",
  platforms = unique(c(rhub:::default_cran_check_platforms("."),"macos-elcapitan-release")))


check_for_cran(path = "C:/Users/Sebastian Krantz/Documents/R/collapse_1.0.0.tar.gz",
               platforms = c("macos-elcapitan-release","ubuntu-gcc-release",
                             "linux-x86_64-rocker-gcc-san","fedora-clang-devel"))


## Let's run some benchmarks and compare fsum against data.table and base::rowsum
# Starting with small data
mtcDT <- qDT(mtcars)
f <- qF(mtcars$cyl)

library(microbenchmark)
microbenchmark(mtcDT[, lapply(.SD, sum), by = f],
               rowsum(mtcDT, f, reorder = FALSE),
               fsum(mtcDT, f, na.rm = FALSE), unit = "relative")
# Now larger data
tdata <- qDT(replicate(100, rnorm(1e5), simplify = FALSE)) # 100 columns with 100.000 obs
f <- qF(sample.int(1e4, 1e5, TRUE))                        # A factor with 10.000 groups

microbenchmark(tdata[, lapply(.SD, sum), by = f],
               rowsum(tdata, f, reorder = FALSE),
               fsum(tdata, f), unit = "relative")
# My results:
expr      min       lq     mean   median       uq       max neval cld
tdata[, lapply(.SD, sum), by = f] 2.646992 2.975489 2.834771 3.081313 3.120070 1.2766475   100   c
rowsum(tdata, f, reorder = FALSE) 1.747567 1.753313 1.629036 1.758043 1.839348 0.2720937   100  b
fsum(tdata, f, na.rm = FALSE) 1.000000 1.000000 1.000000 1.000000 1.000000 1.0000000   100 a

# long vs not long

PM = qM(num_vars(PRIO))

> system.time(fsum(PM))
User      System verstrichen
1.26        0.00        1.28
> system.time(fsum(PM))
User      System verstrichen
1.26        0.00        1.27
> system.time(fmean(PM))
User      System verstrichen
1.44        0.00        1.46
> system.time(fsd(PM))
User      System verstrichen
0.97        0.00        1.02
> system.time(fsd(PM))
User      System verstrichen
0.82        0.00        0.84
> system.time(fsd(PM))
User      System verstrichen
0.88        0.00        0.92
> rsd = fsd(PM)
> fm = fmean(PM)
> fs = fsum(PM)
> system.time(fscale(PM))
User      System verstrichen
1.78        2.30       11.08
> system.time(fscale(PM))
User      System verstrichen
1.68        0.91        8.58
> gc()
used   (Mb) gc trigger   (Mb)  max used   (Mb)
Ncells   3708319  198.1    5640972  301.3   5640972  301.3
Vcells 288953017 2204.6  703047783 5363.9 584230291 4457.4
> gc()
used   (Mb) gc trigger   (Mb)  max used   (Mb)
Ncells   3707812  198.1    5640972  301.3   5640972  301.3
Vcells 288952179 2204.6  703047783 5363.9 584230291 4457.4
> system.time(qsu(PM))
User      System verstrichen
1.31        1.16        8.45
> system.time(qsu(PM))
User      System verstrichen
1.49        0.60        2.09
> su =  qsu(PM)
> gc()

# They are pretty much the same with or without long doubles !! all.equal tolerance 1e-14 about.
# Note that qsu seems not much faster after removing long doubles !!! (see about higher moment precision !!)
