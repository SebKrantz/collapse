#####################################################
### Replication Script for the JSS Article:
### collapse: Advanced and Fast Statistical Computing
###           and Data Transformation in R
### -------------------------------------------------
### By: Sebastian Krantz, IfW Kiel
### E-Mail: sebastian.krantz@ifw-kiel.de
#####################################################

###################################################
### code chunk number 0: Preliminaries
###################################################

options(prompt = "R> ", continue = "+  ", width = 77, digits = 4, useFancyQuotes = FALSE, warn = 1)

# Loading libraries and installing if unavailable
if(!requireNamespace("fastverse", quietly = TRUE)) install.packages("fastverse")
options(fastverse.styling = FALSE)
library(fastverse) # loads data.table, collapse, magrittr and kit (not used)
fastverse_extend(microbenchmark, Rfast, fixest, install = TRUE) # loads and installs if unavailable
# Package versions used in the article:
# fastverse 0.3.2, collapse 2.0.10, data.table 1.15.0, magrittr 2.0.3,
# microbenchmark 1.4.10, Rfast 2.1.0, and fixest 0.11.3

###################################################
### code chunk number 1: collapse Topics and Documentation
###################################################
.COLLAPSE_TOPICS
help("collapse-documentation")


###################################################
### code chunk number 2: Fast Statistical Functions: Basic Examples
###################################################
fmean(mtcars$mpg)
fmean(EuStockMarkets)
fmean(mtcars[5:10])
fmean(mtcars$mpg, w = mtcars$wt)
fmean(mtcars$mpg, g = mtcars$cyl)
fmean(mtcars$mpg, g = mtcars$cyl, w = mtcars$wt)
fmean(mtcars[5:10], g = mtcars$cyl, w = mtcars$wt)
fmean(mtcars$mpg, g = mtcars$cyl, TRA = "fill") |> head(20)


###################################################
### code chunk number 3: Airquality Dataset
###################################################
fnobs(airquality)


###################################################
### code chunk number 4: Imputation by Reference
###################################################
fmedian(airquality[1:2], airquality$Month, TRA = "replace_na", set = TRUE)


###################################################
### code chunk number 5: Transformation Example
###################################################
airquality |> fmutate(
  rad_day = fsum(as.double(Solar.R), Day, TRA = "/"),
  ozone_deg = Ozone / Temp,
  ozone_amed = Ozone > fmedian(Ozone, Month, TRA = "fill"),
  ozone_resid = fmean(Ozone, list(Month, ozone_amed), ozone_deg, "-")
) |> head(3)


###################################################
### code chunk number 6: GRP Objects
###################################################
str(g <- GRP(mtcars, ~ cyl + vs + am))


###################################################
### code chunk number 7: Aggregation with GRP Objects
###################################################
dat <- get_vars(mtcars, c("mpg", "disp")); w <- mtcars$wt
add_vars(g$groups,
  fmean(dat, g, w, use.g.names = FALSE) |> add_stub("w_mean_"),
  fsd(dat, g, w, use.g.names = FALSE) |> add_stub("w_sd_")) |> head(2)


###################################################
### code chunk number 8: Transformation with GRP Objects
###################################################
mtcars |> add_vars(fmean(dat, g, w, "-") |> add_stub("w_demean_"),
                   fscale(dat, g, w) |> add_stub("w_scale_")) |> head(2)


###################################################
### code chunk number 9: fsummarise Integration
###################################################
mtcars |>
  fsubset(mpg > 11) |>
  fgroup_by(cyl, vs, am) |>
  fsummarise(across(c(mpg, carb, hp), fmean),
             qsec_w_med = fmean(qsec, wt)) |> head(2)


###################################################
### code chunk number 10: grouped_df Methods for Fast Statistical Functions
###################################################
mtcars |>
  fsubset(mpg > 11, cyl, vs, am, mpg, carb, hp, wt) |>
  fgroup_by(cyl, vs, am) |>
  fmean(wt) |> head(2)


###################################################
### code chunk number 11: Vectorized Grouped Linear Regression
###################################################
mtcars |>
 fgroup_by(vs) |>
 fmutate(dm_carb = fmean(carb, TRA = "-")) |>
 fsummarise(slope = fsum(mpg, dm_carb) %/=% fsum(dm_carb^2))


###################################################
### code chunk number 12: Advanced Weighted Group Statistics
###################################################
mtcars |>
    fgroup_by(cyl, vs, am) |>
    fmutate(o = radixorder(GRPid(), mpg)) |>
    fsummarise(mpg_min = fmin(mpg),
               mpg_Q1 = fnth(mpg, 0.25, wt, o = o, ties = "q8"),
               mpg_mean = fmean(mpg, wt),
               mpg_median = fmedian(mpg, wt, o = o, ties = "q8"),
               mpg_mode = fmode(mpg, wt, ties = "max"),
               mpg_Q3 = fnth(mpg, 0.75, wt, o = o, ties = "q8"),
               mpg_max = fmax(mpg)) |> head(3)


###################################################
### code chunk number 13: Data Aggregation with collap()
###################################################
collap(wlddev, country + PCGDP + LIFEEX ~ year + income, w = ~ POP) |>
  head(4)


###################################################
### code chunk number 14: Growth Rate of Airmiles Time Series
###################################################
fgrowth(airmiles) |> round(2)


###################################################
### code chunk number 15: Creating an Irregular Series and Demonstrating Indexation
###################################################
am_ir <- airmiles[-c(3, 15)]
t <- time(airmiles)[-c(3, 15)]
fgrowth(am_ir, t = t) |> round(2)
fgrowth(am_ir, -1:3, t = t) |> head(4)


###################################################
### code chunk number 16: Ad-Hoc Transformations on World Bank Panel Data Supplied with collapse
###################################################
G(wlddev, c(1, 10), by = POP + LIFEEX ~ iso3c, t = ~ year) |> head(3)
settransform(wlddev, POP_growth = G(POP, g = iso3c, t = year))


###################################################
### code chunk number 17: Integration with Data Manipualtion Functions
###################################################
wlddev |> fgroup_by(iso3c) |> fselect(iso3c, year, POP, LIFEEX) |>
  fmutate(across(c(POP, LIFEEX), G, t = year)) |> head(2)


###################################################
### code chunk number 18: Two Solutions for Grouped Scaling
###################################################
iris |> fgroup_by(Species) |> fscale() |> head(2)
STD(iris, ~ Species) |> head(2)


###################################################
### code chunk number 19: Fixed Effects Regression a la Mundlak (1978)
###################################################
lm(mpg ~ carb + B(carb, cyl), data = mtcars) |> coef()


###################################################
### code chunk number 20: Detrending with Country-Level Cubic Polynomials: Requires {fixest}
###################################################
HDW(wlddev, PCGDP + LIFEEX ~ iso3c * poly(year, 3), stub = F) |> head(2)


###################################################
### code chunk number 21: Indexed Frame
###################################################
wldi <- wlddev |> findex_by(iso3c, year)
wldi |> fsubset(-3, iso3c, year, PCGDP:POP) |> G() |> head(4)


###################################################
### code chunk number 22: Indexed Series
###################################################
LIFEEXi = wldi$LIFEEX
str(LIFEEXi, width = 70, strict = "cut")
c(is_irregular(LIFEEXi), is_irregular(LIFEEXi[-5]))
G(LIFEEXi[c(1:5, 7:10)])


###################################################
### code chunk number 23: Demonstrating Deep Indexation
###################################################
settransform(wldi, PCGDP_ld = Dlog(PCGDP))
lm(D(LIFEEX) ~ L(PCGDP_ld, 0:5) + B(PCGDP_ld), wldi) |>
  summary() |> coef() |> round(3)


###################################################
### code chunk number 24: Using 3rd Party Functions: Rolling Average
###################################################
BY(LIFEEXi, findex(LIFEEXi)$iso3c, data.table::frollmean, 5) |> head(10)


###################################################
### code chunk number 25: Joins: Adding Join Column
###################################################
df1 <- data.frame(id1 = c(1, 1, 2, 3),
                  id2 = c("a", "b", "b", "c"),
                  name = c("John", "Jane", "Bob", "Carl"),
                  age = c(35, 28, 42, 50))
df2 <- data.frame(id1 = c(1, 2, 3, 3),
                  id2 = c("a", "b", "c", "e"),
                  salary = c(60000, 55000, 70000, 80000),
                  dept = c("IT", "Marketing", "Sales", "IT"))

join(df1, df2, on = c("id1", "id2"), how = "full", column = TRUE)


###################################################
### code chunk number 26: Validation + Join Attribute
###################################################
join(df1, df2, on = c("id1", "id2"), validate = "1:1", attr = "join") |>
  attr("join") |> str(width = 70, strict = "cut")


###################################################
### code chunk number 27: Overidentification Warning
###################################################
df2$name = df1$name
join(df1, df2) |> capture.output(type="m") |> strwrap(77) |> cat(sep="\n")


###################################################
### code chunk number 28: Automatic Renaming
###################################################
join(df1, df2, on = c("id1", "id2"))


###################################################
### code chunk number 29: Data for Pivots
###################################################
data <- data.frame(type = rep(c("A", "B"), each = 2),
            type_name = rep(c("Apples", "Bananas"), each = 2),
            id = rep(1:2, 2), r = abs(rnorm(4)), h = abs(rnorm(4)*2))
setrelabel(data, id = "Fruit Id", r = "Fruit Radius", h = "Fruit Height")
print(data)
vlabels(data)


###################################################
### code chunk number 30: Pivot Longer
###################################################
(dl <- pivot(data, ids = c("type", "type_name", "id"), labels = "label"))
vlabels(dl)


###################################################
### code chunk number 31: Pivot Wider
###################################################
(dw <- pivot(data, "id", names = "type", labels = "type_name", how = "w"))
namlab(dw)


###################################################
### code chunk number 32: Pivot Recast
###################################################
(dr <- pivot(data, ids = "id", names = list(from = "type"),
             labels = list(from = "type_name", to = "label"), how = "r"))
vlabels(dr)


###################################################
### code chunk number 33: Recursive Splitting: Creates Nested List of Data Frames
###################################################
(dl <- mtcars |> rsplit(mpg + hp + carb ~ vs + am)) |> str(max.level = 2)


###################################################
### code chunk number 34: Fitting Linear Models and Obtaining Coefficient Matrices
###################################################
nest_lm_coef <- dl |>
  rapply2d(lm, formula = mpg ~ .) |>
  rapply2d(summary, classes = "lm") |>
  get_elem("coefficients")

nest_lm_coef |> str(give.attr = FALSE, strict = "cut")


###################################################
### code chunk number 35: Unlisting to Data Frame
###################################################
nest_lm_coef |> unlist2d(c("vs", "am"), row.names = "variable") |> head(2)


###################################################
### code chunk number 36: Removing Generated Series (Hidden)
###################################################
wldi <- wldi[1:13]

###################################################
### code chunk number 37: Which Columns/Countries have Time Varying Information?
###################################################
varying(wldi)
varying(wldi, any_group = FALSE) |> head(3)


###################################################
### code chunk number 38: Demonstrating Panel Decomposition
###################################################
all.equal(fvar(W(LIFEEXi)) + fvar(B(LIFEEXi)), fvar(LIFEEXi))


###################################################
### code chunk number 39: Panel Summary Statistics
###################################################
qsu(LIFEEXi)


###################################################
### code chunk number 40: Weighted Panel Summary Statistics by Groups
###################################################
qsu(LIFEEXi, g = wlddev$OECD, w = wlddev$POP, higher = TRUE) |> aperm()


###################################################
### code chunk number 41: Detailed (Grouped, Weighted) Statistical Description
###################################################
descr(wlddev, LIFEEX ~ OECD, w = ~ replace_na(POP))


###################################################
### code chunk number 42: qtab: Basic Usage
###################################################
library(magrittr) # World after 2015 (latest country data)
wlda15 <- wlddev |> fsubset(year >= 2015) |> fgroup_by(iso3c) |> flast()
wlda15 %$% qtab(OECD, income)


###################################################
### code chunk number 43: qtab: Population Counts
###################################################
wlda15 %$% qtab(OECD, income, w = POP) %>% divide_by(1e6)


###################################################
### code chunk number 44: qtab: Average Life Expectancy
###################################################
wlda15 %$% qtab(OECD, income, w = LIFEEX, wFUN = fmean) %>% replace_na(0)


###################################################
### code chunk number 45: qtab: Population Weighted Average Life Expectancy
###################################################
wlda15 %$% qtab(OECD, income, w = LIFEEX, wFUN = fmean,
                wFUN.args = list(w = POP)) %>% replace_na(0)


###################################################
### code chunk number 46: Benchmark: Statistics and Data Manipulation
###################################################
setDTthreads(4)
set_collapse(na.rm = FALSE, sort = FALSE, nthreads = 4)
set.seed(101)
m <- matrix(rnorm(1e7), ncol = 1000)
data <- qDT(replicate(100, rnorm(1e5), simplify = FALSE))
g <- sample.int(1e4, 1e5, TRUE)

microbenchmark(R = colMeans(m),
               Rfast = Rfast::colmeans(m, parallel = TRUE, cores = 4),
               collapse = fmean(m))
microbenchmark(R = rowsum(data, g, reorder = FALSE),
               data.table = data[, lapply(.SD, sum), by = g],
               collapse = fsum(data, g))
add_vars(data) <- g
microbenchmark(data.table = data[, lapply(.SD, median), by = g],
               collapse = data |> fgroup_by(g) |> fmedian())
d <- data.table(g = unique(g), x = 1, y = 2, z = 3)
microbenchmark(data.table = d[data, on = "g"],
               collapse = join(data, d, on = "g", verbose = 0))
microbenchmark(data.table = melt(data, "g"),
               collapse = pivot(data, "g"))
settransform(data, id = rowid(g))
cols <- grep("^V", names(data), value = TRUE)
microbenchmark(data.table = dcast(data, g ~ id, value.var = cols),
          collapse = pivot(data, ids = "g", names = "id", how = "w"))


###################################################
### code chunk number 47: Benchmark: Unique Values and Matching
###################################################
set.seed(101)
g_int <- sample.int(1e3, 1e7, replace = TRUE)
char <- c(letters, LETTERS, month.abb, month.name)
char <- outer(char, char, paste0)
g_char <- sample(char, 1e7, replace = TRUE)
microbenchmark(base_int = unique(g_int), collapse_int = funique(g_int),
            base_char = unique(g_char), collapse_char = funique(g_char))
microbenchmark(base_int = match(g_int, 1:1000),
               collapse_int = fmatch(g_int, 1:1000),
               base_char = match(g_char, char),
               data.table_char = chmatch(g_char, char),
               collapse_char = fmatch(g_char, char), times = 10)


###################################################
### Print Session Information
###################################################

sessionInfo()
