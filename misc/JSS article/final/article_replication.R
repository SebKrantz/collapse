#####################################################
### Replication Script for the JSS Article:
### collapse: Advanced and Fast Statistical Computing
###           and Data Transformation in R
### -------------------------------------------------
### By: Sebastian Krantz, IfW Kiel
### E-Mail: sebastian.krantz@ifw-kiel.de
#####################################################

###################################################
### code chunk number 1: preliminaries
###################################################

rm(list = ls()); gc()
r_opts <- options(prompt = "R> ", continue = "+  ", width = 77, digits = 4, useFancyQuotes = FALSE, warn = 1)

# Loading libraries and installing if unavailable
if(!requireNamespace("fastverse", quietly = TRUE)) install.packages("fastverse")
options(fastverse.styling = FALSE)
library(fastverse) # loads data.table, collapse, magrittr and kit (not used)
# Package versions used in the article:
# fastverse 0.3.4, collapse 2.1.2, data.table 1.17.0, nycflights23 0.2.0,
# bench 1.1.4, dplyr 1.1.4, tidyr 1.3.1, matrixStats 1.5.0, janitor 2.2.0
pkg <- c("bench", "dplyr", "tidyr", "matrixStats", "janitor", "nycflights23")
if(!all(avail <- is_installed(pkg))) install.packages(pkg[!avail])

# Reset collapse options (if set)
oldopts <- set_collapse(nthreads = 1L, remove = NULL, stable.algo = TRUE, sort = TRUE,
             digits = 2L, stub = TRUE, verbose = 1L, mask = NULL, na.rm = TRUE)

# For presentation: removing whitespace around labels
vlabels(GGDC10S) <- trimws(vlabels(GGDC10S))

# This is the benchmark function
bmark <- function(...) {
  bench::mark(..., min_time = 2, check = FALSE) |> janitor::clean_names() |>
  fselect(expression, min, median, mem_alloc, n_itr, n_gc, total_time) |>
  fmutate(expression = names(expression)) |> dapply(as.character) |> qDF()
}


###################################################
### code chunk number 2: collapse_topics
###################################################
.COLLAPSE_TOPICS
help("collapse-documentation")


###################################################
### code chunk number 3: faststatfun
###################################################
fmean(iris$Sepal.Length)
fmean(num_vars(iris))
identical(fmean(num_vars(iris)), fmean(as.matrix(num_vars(iris))))
fmean(iris$Sepal.Length, g = iris$Species)
fmean(num_vars(iris), g = iris$Species, nthreads = 4)


###################################################
### code chunk number 4: TRA
###################################################
fmean(iris$Sepal.Length, g = iris$Species, TRA = "fill")[c(1:5, 51:55)]


###################################################
### code chunk number 5: BY
###################################################
BY(num_vars(iris), g = iris$Species, FUN = mad) |> head(1)


###################################################
### code chunk number 6: fscale
###################################################
fscale(iris$Sepal.Length, g = iris$Species) |> fsd(iris$Species)


###################################################
### code chunk number 7: wlddev_distinct
###################################################
fndistinct(wlddev)


###################################################
### code chunk number 8: fsummarise
###################################################
wlddev |> fgroup_by(country, income) |>
  fsummarise(n = fnobs(LIFEEX),
             diff = fmax(LIFEEX) %-=% fmin(LIFEEX),
             mean_diff = fmean(fdiff(LIFEEX)),
             cor_PCGDP = pwcor(Dlog(LIFEEX), Dlog(PCGDP)),
             across(c(LIFEEX, PCGDP, POP), flast)) |>
  collap( ~ income, w = ~ POP, keep.w = FALSE) |> print(digits = 2)


###################################################
### code chunk number 9: fsummarise_2
###################################################
wlddev |> fsubset(year >= 2015, income, PCGDP:POP) |>
  fgroup_by(income) |> fmean(POP)


###################################################
### code chunk number 10: glinreg
###################################################
wlddev |> collap(POP ~ region + year, FUN = fsum) |>
  fmutate(POP = POP / 1e6) |> fgroup_by(region) |>
  fmutate(dmy = fmean(year, TRA = "-")) |>
  fsummarise(beta = fsum(POP, dmy) %/=% fsum(dmy, dmy),
             POP20 = flast(POP)) |>
  fmutate(POP21 = POP20 + beta, POP22 = POP21 + beta,
          POP23 = POP22 + beta, POP24 = POP23 + beta)


###################################################
### code chunk number 11: wstats
###################################################
wlddev |> fsubset(is.finite(POP)) |> fgroup_by(region) |>
  fmutate(o = radixorder(GRPid(), LIFEEX)) |>
  fsummarise(min = fmin(LIFEEX),
             Q1 = fnth(LIFEEX, 0.25, POP, o = o, ties = "q8"),
             mean = fmean(LIFEEX, POP),
             median = fmedian(LIFEEX, POP, o = o),
             Q3 = fnth(LIFEEX, 0.75, POP, o = o, ties = "q8"),
             max = fmax(LIFEEX))


###################################################
### code chunk number 12: GRP
###################################################
str(g <- GRP(wlddev, ~ income + OECD))


###################################################
### code chunk number 13: gcomp
###################################################
add_vars(g$groups,
         get_vars(wlddev, "country") |> fmode(g, wlddev$POP, use = FALSE),
         get_vars(wlddev, c("PCGDP", "LIFEEX")) |> fmean(g, wlddev$POP, use = F),
         get_vars(wlddev, "POP") |> fsum(g, use = FALSE))


###################################################
### code chunk number 14: gtrans
###################################################
add_vars(wlddev) <- get_vars(wlddev, c("PCGDP", "LIFEEX")) |>
  fwithin(g, mean = "overall.mean") |> add_stub("center_")


###################################################
### code chunk number 15: fsumm_bench
###################################################
bmark(collapse = collapse::fsummarise(mtcars, sum = sum(mpg)),
      dplyr = dplyr::summarise(mtcars, sum = sum(mpg)))


###################################################
### code chunk number 16: ts_example
###################################################
fgrowth(airmiles) |> round(2)


###################################################
### code chunk number 17: exports
###################################################
set.seed(101)
exports <- expand.grid(y = 2001:2010, c = paste0("c", 1:10),
                       s = paste0("s", 1:10)) |> fmutate(v = abs(rnorm(1e3))) |>
  colorder(c, s) |> fsubset(-sample.int(1e3, 500))


###################################################
### code chunk number 18: ts_example_2
###################################################
.c(y, v) %=% fsubset(exports, c == "c1" & s == "s7", -c, -s)
print(y)
fgrowth(v, t = y) |> round(2)
fgrowth(v, -1:3, t = y) |> head(4)


###################################################
### code chunk number 19: example_growth
###################################################
G(exports, -1:2, by = v ~ c + s, t = ~ y) |> head(3)
tfm(exports, fgrowth(list(v = v), -1:2, g = list(c, s), t = y)) |> head(3)


###################################################
### code chunk number 20: example_growth_continued
###################################################
A <- exports |> fgroup_by(c, s) |> fmutate(gv = G(v, t = y)) |> fungroup()
head(B <- exports |> fmutate(gv = G(v, g = list(c, s), t = y)), 4)
identical(A, B)


###################################################
### code chunk number 21: indexing
###################################################
exportsi <- exports |> findex_by(c, s, y)
exportsi |> G() |> head(5)
exportsi |> findex() |> print(2)


###################################################
### code chunk number 22: indexing_2
###################################################
vi <- exportsi$v; str(vi, width = 70, strict = "cut")
is_irregular(vi)
vi |> psmat() |> head(3)
fdiff(vi) |> psmat() |> head(3)


###################################################
### code chunk number 23: indexing_3
###################################################
settransform(exportsi, v_ld = Dlog(v))
lm(v_ld ~ L(v_ld, 1:2), exportsi) |> summary() |> coef() |> round(3)


###################################################
### code chunk number 24: rollmean
###################################################
BY(vi, ix(vi)$c.s, data.table::frollmean, 5) |> head(10)


###################################################
### code chunk number 25: joins
###################################################
teacher <- data.frame(id = 1:4, names = c("John", "Jane", "Bob", "Carl"),
                      age = c(35, 32, 42, 67), subject = c("Math", "Econ", "Stats", "Trade"))
course <- data.frame(id = c(1, 2, 2, 3, 5), semester = c(1, 1, 2, 1, 2),
                     course = c("Math I", "Microecon", "Macroecon", "Stats I", "History"))
join(teacher, course, on = "id")


###################################################
### code chunk number 26: joins_2
###################################################
join(teacher, course, how = "full", multiple = TRUE, column = TRUE)


###################################################
### code chunk number 27: joins_3
###################################################
join(teacher, course, multiple = TRUE, attr = "jn") |> attr("jn") |> str(strict.width = "cut", width = 70)


###################################################
### code chunk number 28: joins_4
###################################################
join(teacher, course, on = "id", validate = "1:1") |>
  tryCatch(error = function(e) strwrap(e) |> cat(sep = "\n"))


###################################################
### code chunk number 29: joins_4.5
###################################################
join(teacher, course, on = "id", require = list(x = 0.8)) |>
  tryCatch(error = function(e) substr(e, 102, 200) |> cat())


###################################################
### code chunk number 30: joins_5
###################################################
for (h in c("semi", "anti")) join(teacher, course, how = h) |> print()


###################################################
### code chunk number 31: joins_6
###################################################
course$names <- teacher$names[course$id]
join(teacher, course, on = "id", how = "inner", multiple = TRUE)


###################################################
### code chunk number 32: joins_7
###################################################
join(teacher, course, on = "id", multiple = TRUE, drop.dup.cols = "y")


###################################################
### code chunk number 33: pivots
###################################################
data <- GGDC10S |>
  fmutate(Label = ifelse(Variable == "VA", "Value Added", "Employment")) |>
  fsubset(is.finite(AGR), Country, Variable, Label, Year, AGR:MAN)
namlab(data, N = TRUE, Ndistinct = TRUE, class = TRUE)


###################################################
### code chunk number 34: pivot_longer
###################################################
head(dl <- pivot(data, ids = 1:4, names = list("Sectorcode", "Value"),
                 labels = "Sector", how = "longer"))


###################################################
### code chunk number 35: pivot_wider
###################################################
head(dw <- pivot(data, c("Country", "Year"), names = "Variable",
                 labels = "Label", how = "w"))
namlab(dw)


###################################################
### code chunk number 36: pivot_recast
###################################################
head(dr <- pivot(data, c("Country", "Year"),
                 names = list(from = "Variable", to = "Sectorcode"),
                 labels = list(from = "Label", to = "Sector"), how = "r"))
vlabels(dr)[3:6]


###################################################
### code chunk number 37: pivot_agg
###################################################
head(dr_agg <- pivot(data, "Country", c("AGR", "MIN", "MAN"), how = "r",
                     names = list(from = "Variable", to = "Sectorcode"),
                     labels = list(from = "Label", to = "Sector"), FUN = "mean"))


###################################################
### code chunk number 38: list_proc
###################################################
dl <- GGDC10S |> rsplit( ~ Country + Variable)
dl$ARG$EMP$AGR[1:12]


###################################################
### code chunk number 39: list_proc_2
###################################################
result <- list()
for (country in c("ARG", "BRA", "CHL")) {
  for (variable in c("EMP", "VA")) {
    m <- lm(log(AGR+1) ~ log(MIN+1) + log(MAN+1) + Year,
            data = dl[[country]][[variable]])
    result[[country]][[variable]] <- list(model = m, BIC = BIC(m),
                                          summary = summary(m))
  }
}


###################################################
### code chunk number 40: list_proc_3
###################################################
str(r_sq_l <- result |> get_elem("r.squared"))
rowbind(r_sq_l, idcol = "Country", return = "data.frame")


###################################################
### code chunk number 41: list_proc_3.5
###################################################
r_sq_l |> t_list() |> rowbind(idcol = "Variable", return = "data.frame")


###################################################
### code chunk number 42: list_proc_4
###################################################
result$ARG$EMP$summary$coefficients


###################################################
### code chunk number 43: list_proc_5
###################################################
result |> get_elem("coefficients") |> get_elem(is.matrix) |>
  unlist2d(idcols = c("Country", "Variable"),
           row.names = "Covariate") |> head(3)


###################################################
### code chunk number 44: varying
###################################################
varying(wlddev, ~ iso3c)


###################################################
### code chunk number 45: pdec
###################################################
LIFEEXi <- reindex(wlddev$LIFEEX, wlddev$iso3c)
all.equal(fvar(LIFEEXi), fvar(fbetween(LIFEEXi)) + fvar(fwithin(LIFEEXi)))


###################################################
### code chunk number 46: qsu_1
###################################################
qsu(LIFEEXi)


###################################################
### code chunk number 47: qsu_2
###################################################
qsu(wlddev, by = LIFEEX ~ income, w = ~ POP)


###################################################
### code chunk number 48: descr
###################################################
wlda15 <- wlddev |> fsubset(year >= 2015) |> fgroup_by(iso3c) |> flast()
wlda15 |> descr(by = income + LIFEEX ~ OECD)


###################################################
### code chunk number 49: qtab
###################################################
wlda15 |> with(qtab(OECD, income))


###################################################
### code chunk number 50: qtab_2
###################################################
wlda15 |> with(qtab(OECD, income, w = POP) / 1e6)


###################################################
### code chunk number 51: qtab_3
###################################################
wlda15 |> with(qtab(OECD, income, w = LIFEEX, wFUN = fmean))


###################################################
### code chunk number 52: qtab_4
###################################################
wlda15 |> with(qtab(OECD, income, w = LIFEEX, wFUN = fmean,
                    wFUN.args = list(w = POP)))


###################################################
### code chunk number 53: bench_1
###################################################
set.seed(101)
int <- 1:1000; g_int <- sample.int(1000, 1e7, replace = TRUE)
char <- c(letters, LETTERS, month.abb, month.name)
g_char <- sample(char <- outer(char, char, paste0), 1e7, TRUE)
bmark(base_int = unique(g_int), collapse_int = funique(g_int))
bmark(base_char = unique(g_char), collapse_char = funique(g_char))
bmark(base_int = match(g_int, int), collapse_int = fmatch(g_int, int))
bmark(base_char = match(g_char, char), data.table_char =
        chmatch(g_char, char), collapse_char = fmatch(g_char, char))


###################################################
### code chunk number 54: bench_2
###################################################
set_collapse(na.rm = FALSE, sort = FALSE, nthreads = 4)
m <- matrix(rnorm(1e7), ncol = 1000)
bmark(R = colSums(m), collapse = fsum(m))
bmark(R = colMeans(m), collapse = fmean(m))
bmark(MS = matrixStats::colMedians(m), collapse = fmedian(m))


###################################################
### code chunk number 55: bench_2.1
###################################################
g <- sample.int(1e3, 1e4, TRUE)
bmark(R = rowsum(m, g), collapse = fsum(m, g))


###################################################
### code chunk number 56: bench_flights_setup
###################################################
fastverse_extend(nycflights23, dplyr, data.table); setDTthreads(4)
list(flights, airports, airlines, planes, weather) |> sapply(nrow)
flights |> fselect(month, day, origin, dest) |> fnunique()


###################################################
### code chunk number 57: bench_flights_summarise
###################################################
vars <- .c(dep_delay, arr_delay, air_time, distance, hour, minute)
bmark(dplyr = flights |> group_by(month, day, origin, dest) |>
        summarise(across(all_of(vars), sum), .groups = "drop"),
      data.table = qDT(flights)[, lapply(.SD, sum), .SDcols = vars,
                                by = .(month, day, origin, dest)],
      collapse = flights |> fgroup_by(month, day, origin, dest) |>
        get_vars(vars) |> fsum())


###################################################
### code chunk number 58: bench_flights_summarise_2
###################################################
bmark(dplyr_mean = flights |> group_by(month, day, origin, dest) |>
        summarise(across(all_of(vars), mean), .groups = "drop"),
      data.table_mean = qDT(flights)[, lapply(.SD, mean), .SDcols = vars,
                                     by = .(month, day, origin, dest)],
      collapse_mean = flights |> fgroup_by(month, day, origin, dest) |>
        get_vars(vars) |> fmean())
bmark(dplyr_median = flights |> group_by(month, day, origin, dest) |>
        summarise(across(all_of(vars), median), .groups = "drop"),
      data.table_median = qDT(flights)[, lapply(.SD, median), .SDcols = vars,
                                       by = .(month, day, origin, dest)],
      collapse_median = flights |> fgroup_by(month, day, origin, dest) |>
        get_vars(vars) |> fmedian())


###################################################
### code chunk number 59: bench_flights_summarise_3
###################################################
bmark(dplyr = flights |> group_by(month, day, origin, dest) |>
        summarise(rng = max(arr_delay) - min(arr_delay), .groups = "drop"),
      data.table = qDT(flights)[, .(rng = max(arr_delay) - min(arr_delay)),
                                by = .(month, day, origin, dest)],
      collapse = flights |> fgroup_by(month, day, origin, dest) |>
        fsummarise(rng = fmax(arr_delay) - fmin(arr_delay)))


###################################################
### code chunk number 60: bench_flights_join_setup
###################################################
flights |> join(weather, on = c("origin", "time_hour")) |>
  join(planes, on = "tailnum") |> join(airports, on = c(dest = "faa")) |>
  join(airlines, on = "carrier") |> dim() |>
  capture.output() |> substr(1, 76) |> cat(sep = "\n")


###################################################
### code chunk number 61: bench_flights_join
###################################################
bmark(
  dplyr_joins = flights |>
    left_join(weather, by = c("origin", "time_hour"), multiple = "first") |>
    left_join(planes, by = "tailnum", multiple = "first") |>
    left_join(airports, by = c(dest = "faa"), multiple = "first") |>
    left_join(airlines, by = "carrier", multiple = "first"),

  data.table_joins = qDT(airlines)[qDT(airports)[qDT(planes)[qDT(weather)[qDT(flights), on = c("origin", "time_hour"), mult = "first"], on = "tailnum", mult = "first"], on = c("faa" = "dest"), mult = "first"], on = "carrier", mult = "first"],

  collapse_joins = flights |>
    join(weather, on = c("origin", "time_hour"), verbose = 0) |>
    join(planes, on = "tailnum", verbose = 0) |>
    join(airports, on = c("dest" = "faa"), verbose = 0) |>
    join(airlines, on = "carrier", verbose = 0)
)


###################################################
### code chunk number 62: bench_pivot_long
###################################################
bmark(tidyr = tidyr::pivot_longer(flights, cols = vars),
      data.table = qDT(flights) |> melt(measure = vars),
      collapse = pivot(flights, values = vars))


###################################################
### code chunk number 63: bench_pivot_wide
###################################################
bmark(tidyr = tidyr::pivot_wider(flights, id_cols = .c(month, day, dest),
                                 names_from = "origin", values_from = vars, values_fn = sum),
      data.table = dcast(qDT(flights), month + day + dest ~ origin,
                         value.var = vars, fun = sum),
      collapse_fsum = pivot(flights, .c(month, day, dest), vars,
                            "origin", how = "wider", FUN = fsum),
      collapse_itnl = pivot(flights, .c(month, day, dest), vars,
                            "origin", how = "wider", FUN = "sum"))

###################################################
### code chunk number 64: reset_options
###################################################
options(r_opts)
set_collapse(oldopts)

###################################################
### Print Session Information
###################################################

sessionInfo()

