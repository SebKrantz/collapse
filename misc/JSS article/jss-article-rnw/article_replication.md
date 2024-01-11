

```r
#####################################################
### Replication Script for the JSS Article:
### collapse: Advanced and Fast Statistical Computing
###           and Data Transformation in R
```

```r
### By: Sebastian Krantz, IfW Kiel
### E-Mail: sebastian.krantz@ifw-kiel.de
#####################################################

###################################################
### code chunk number 1: Preliminaries
###################################################

options(prompt = "R> ", continue = "+  ", width = 80, digits = 4, useFancyQuotes = FALSE, warn = 1)
library(data.table)     # v1.14.10
```

```
## Warning: package 'data.table' was built under R version 4.3.1
```

```
## data.table 1.14.10 using 1 threads (see ?getDTthreads).  Latest news: r-datatable.com
```

```
## **********
## This installation of data.table has not detected OpenMP support. It should still work but in single-threaded mode.
## This is a Mac. Please read https://mac.r-project.org/openmp/. Please engage with Apple and ask them for support. Check r-datatable.com for updates, and our Mac instructions here: https://github.com/Rdatatable/data.table/wiki/Installation. After several years of many reports of installation problems on Mac, it's time to gingerly point out that there have been no similar problems on Windows or Linux.
## **********
```

```r
library(magrittr)       # v2.0.3
library(microbenchmark) # v1.4.10
library(collapse)       # v2.0.8
```

```
## collapse 2.0.8, see ?`collapse-package` or ?`collapse-documentation`
```

```
## 
## Attaching package: 'collapse'
```

```
## The following object is masked from 'package:stats':
## 
##     D
```

```r
# Also used: {Rfast} v2.1.0 and {fixest} v0.11.2

###################################################
### code chunk number 2: collapse Topics and Documentation
###################################################
.COLLAPSE_TOPICS
```

```
##  [1] "collapse-documentation"     "fast-statistical-functions"
##  [3] "fast-grouping-ordering"     "fast-data-manipulation"    
##  [5] "quick-conversion"           "advanced-aggregation"      
##  [7] "data-transformations"       "time-series-panel-series"  
##  [9] "list-processing"            "summary-statistics"        
## [11] "recode-replace"             "efficient-programming"     
## [13] "small-helpers"              "collapse-options"
```

```r
help("collapse-documentation")


###################################################
### code chunk number 3: Airquality Dataset
###################################################
fnobs(airquality)
```

```
##   Ozone Solar.R    Wind    Temp   Month     Day 
##     116     146     153     153     153     153
```

```r
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
```

```
##   Ozone Solar.R Wind Temp Month Day rad_day ozone_deg ozone_amed ozone_resid
## 1    41     190  7.4   67     5   1   0.191    0.6119       TRUE     -10.279
## 2    36     118  8.0   72     5   2   0.135    0.5000       TRUE     -15.279
## 3    12     149 12.6   74     5   3   0.168    0.1622      FALSE      -3.035
```

```r
###################################################
### code chunk number 6: GRP Objects
###################################################
str(g <- GRP(mtcars, ~ cyl + vs + am))
```

```
## Class 'GRP'  hidden list of 9
##  $ N.groups    : int 7
##  $ group.id    : int [1:32] 4 4 3 5 6 5 6 2 2 5 ...
##  $ group.sizes : int [1:7] 1 3 7 3 4 12 2
##  $ groups      :'data.frame':	7 obs. of  3 variables:
##   ..$ cyl: num [1:7] 4 4 4 6 6 8 8
##   ..$ vs : num [1:7] 0 1 1 0 1 0 0
##   ..$ am : num [1:7] 1 0 1 1 0 0 1
##  $ group.vars  : chr [1:3] "cyl" "vs" "am"
##  $ ordered     : Named logi [1:2] TRUE FALSE
##   ..- attr(*, "names")= chr [1:2] "ordered" "sorted"
##  $ order       : int [1:32] 27 8 9 21 3 18 19 20 26 28 ...
##   ..- attr(*, "starts")= int [1:7] 1 2 5 12 15 19 31
##   ..- attr(*, "maxgrpn")= int 12
##   ..- attr(*, "sorted")= logi FALSE
##  $ group.starts: int [1:7] 27 8 3 1 4 5 29
##  $ call        : language GRP.default(X = mtcars, by = ~cyl + vs + am)
```

```r
###################################################
### code chunk number 7: Aggregation with GRP Objects
###################################################
dat <- get_vars(mtcars, c("mpg", "disp")); w <- mtcars$wt
add_vars(g$groups,
  fmean(dat, g, w, use.g.names = FALSE) |> add_stub("w_mean_"),
  fsd(dat, g, w, use.g.names = FALSE) |> add_stub("w_sd_")) |> head(2)
```

```
##   cyl vs am w_mean_mpg w_mean_disp w_sd_mpg w_sd_disp
## 1   4  0  1      26.00       120.3    0.000       0.0
## 2   4  1  0      23.02       137.1    1.236      11.6
```

```r
###################################################
### code chunk number 8: Transformation with GRP Objects
###################################################
mtcars |> add_vars(fmean(dat, g, w, "-") |> add_stub("w_demean_"),
                   fscale(dat, g, w) |> add_stub("w_scale_")) |> head(2)
```

```
##               mpg cyl disp  hp drat    wt  qsec vs am gear carb w_demean_mpg
## Mazda RX4      21   6  160 110  3.9 2.620 16.46  0  1    4    4       0.4357
## Mazda RX4 Wag  21   6  160 110  3.9 2.875 17.02  0  1    4    4       0.4357
##               w_demean_disp w_scale_mpg w_scale_disp
## Mazda RX4             5.027      0.6657       0.6657
## Mazda RX4 Wag         5.027      0.6657       0.6657
```

```r
###################################################
### code chunk number 9: fsummarise Integration
###################################################
mtcars |>
  fsubset(mpg > 11) |>
  fgroup_by(cyl, vs, am) |>
  fsummarise(across(c(mpg, carb, hp), fmean),
             qsec_w_med = fmean(qsec, wt)) |> head(2)
```

```
##   cyl vs am  mpg  carb    hp qsec_w_med
## 1   4  0  1 26.0 2.000 91.00      16.70
## 2   4  1  0 22.9 1.667 84.67      21.04
```

```r
###################################################
### code chunk number 10: grouped_df Methods for Fast Statistical Functions
###################################################
mtcars |>
  fsubset(mpg > 11, cyl, vs, am, mpg, carb, hp, wt) |>
  fgroup_by(cyl, vs, am) |>
  fmean(wt) |> head(2)
```

```
##   cyl vs am sum.wt   mpg carb   hp
## 1   4  0  1  2.140 26.00 2.00 91.0
## 2   4  1  0  8.805 23.02 1.72 83.6
```

```r
###################################################
### code chunk number 11: Vectorized Grouped Linear Regression
###################################################
mtcars |>
 fgroup_by(vs) |>
 fmutate(dm_carb = fmean(carb, TRA = "-")) |>
 fsummarise(slope = fsum(mpg, dm_carb) %/=% fsum(dm_carb^2))
```

```
##   vs   slope
## 1  0 -0.5557
## 2  1 -2.0706
```

```r
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
```

```
##   cyl vs am mpg_min mpg_Q1 mpg_mean mpg_median mpg_mode mpg_Q3 mpg_max
## 1   4  0  1    26.0  26.00    26.00      26.00     26.0  26.00    26.0
## 2   4  1  0    21.5  22.10    23.02      23.17     24.4  24.38    24.4
## 3   4  1  1    21.4  22.29    27.74      27.85     30.4  31.79    33.9
```

```r
###################################################
### code chunk number 13: Data Aggregation with collap()
###################################################
collap(wlddev, country + PCGDP + LIFEEX ~ year + income, w = ~ POP) |>
  head(4)
```

```
##         country year              income   PCGDP LIFEEX       POP
## 1 United States 1960         High income 12768.7  68.59 7.495e+08
## 2      Ethiopia 1960          Low income   658.5  38.33 1.474e+08
## 3         India 1960 Lower middle income   500.8  45.27 9.280e+08
## 4         China 1960 Upper middle income  1166.1  49.86 1.184e+09
```

```r
###################################################
### code chunk number 14: Growth Rate of Airmiles Time Series
###################################################
fgrowth(airmiles) |> round(2)
```

```
## Time Series:
## Start = 1937 
## End = 1960 
## Frequency = 1 
##  [1]    NA 16.50 42.29 54.03 31.65  2.38 15.23 33.29 54.36 76.92  2.71 -2.10
## [13] 12.91 18.51 32.03 18.57 17.82 13.61 18.19 12.83 13.32  0.01 15.49  4.25
```

```r
###################################################
### code chunk number 15: Creating an Irregular Series and Demonstrating Indexation
###################################################
am_ir <- airmiles[-c(3, 15)]
t <- time(airmiles)[-c(3, 15)]
fgrowth(am_ir, t = t) |> round(2)
```

```
##  [1]    NA 16.50    NA 31.65  2.38 15.23 33.29 54.36 76.92  2.71 -2.10 12.91
## [13] 18.51    NA 17.82 13.61 18.19 12.83 13.32  0.01 15.49  4.25
```

```r
fgrowth(am_ir, -1:3, t = t) |> head(4)
```

```
##          FG1   --    G1  L2G1  L3G1
## [1,] -14.167  412    NA    NA    NA
## [2,]      NA  480 16.50    NA    NA
## [3,] -24.043 1052    NA 119.2 155.3
## [4,]  -2.327 1385 31.65    NA 188.5
```

```r
###################################################
### code chunk number 16: Ad-Hoc Transformations on World Bank Panel Data Supplied with collapse
###################################################
G(wlddev, c(1, 10), by = POP + LIFEEX ~ iso3c, t = ~ year) |> head(3)
```

```
##   iso3c year G1.POP L10G1.POP G1.LIFEEX L10G1.LIFEEX
## 1   AFG 1960     NA        NA        NA           NA
## 2   AFG 1961  1.917        NA     1.590           NA
## 3   AFG 1962  1.985        NA     1.544           NA
```

```r
settransform(wlddev, POP_growth = G(POP, g = iso3c, t = year))


###################################################
### code chunk number 17: Integration with Data Manipualtion Functions
###################################################
wlddev |> fgroup_by(iso3c) |> fselect(iso3c, year, POP, LIFEEX) |>
  fmutate(across(c(POP, LIFEEX), G, t = year)) |> head(2)
```

```
##   iso3c year     POP LIFEEX G1.POP G1.LIFEEX
## 1   AFG 1960 8996973  32.45     NA        NA
## 2   AFG 1961 9169410  32.96  1.917      1.59
```

```r
###################################################
### code chunk number 18: Two Solutions for Grouped Scaling
###################################################
iris |> fgroup_by(Species) |> fscale() |> head(2)
```

```
##   Species Sepal.Length Sepal.Width Petal.Length Petal.Width
## 1  setosa       0.2667      0.1899       -0.357     -0.4365
## 2  setosa      -0.3007     -1.1291       -0.357     -0.4365
```

```r
STD(iris, ~ Species) |> head(2)
```

```
##   Species STD.Sepal.Length STD.Sepal.Width STD.Petal.Length STD.Petal.Width
## 1  setosa           0.2667          0.1899           -0.357         -0.4365
## 2  setosa          -0.3007         -1.1291           -0.357         -0.4365
```

```r
###################################################
### code chunk number 19: Fixed Effects Regression a la Mundlak (1978)
###################################################
lm(mpg ~ carb + B(carb, cyl), data = mtcars) |> coef()
```

```
##  (Intercept)         carb B(carb, cyl) 
##      34.8297      -0.4655      -4.7750
```

```r
###################################################
### code chunk number 20: Detrending with Country-Level Cubic Polynomials: Requires {fixest}
###################################################
HDW(wlddev, PCGDP + LIFEEX ~ iso3c * poly(year, 3), stub = F) |> head(2)
```

```
##    PCGDP   LIFEEX
## 1  9.964 0.023670
## 2 14.045 0.006743
```

```r
###################################################
### code chunk number 21: Indexed Frame
###################################################
wldi <- wlddev |> findex_by(iso3c, year)
wldi |> fsubset(-3, iso3c, year, PCGDP:POP) |> G() |> head(4)
```

```
##   iso3c year G1.PCGDP G1.LIFEEX G1.GINI G1.ODA G1.POP
## 1   AFG 1960       NA        NA      NA     NA     NA
## 2   AFG 1961       NA     1.590      NA  98.75  1.917
## 3   AFG 1963       NA        NA      NA     NA     NA
## 4   AFG 1964       NA     1.448      NA  24.48  2.112
## 
## Indexed by:  iso3c [1] | year [4 (61)]
```

```r
###################################################
### code chunk number 22: Indexed Series
###################################################
LIFEEXi = wldi$LIFEEX
str(LIFEEXi, width = 70, strict = "cut")
```

```
##  'indexed_series' num [1:13176] 32.4 33 33.5 34 34.5 ...
##  - attr(*, "label")= chr "Life expectancy at birth, total (years)"
##  - attr(*, "index_df")=Classes 'index_df', 'pindex' and 'data.frame'..
##   ..$ iso3c: Factor w/ 216 levels "ABW","AFG","AGO",..: 2 2 2 2 2 2 ..
##   .. ..- attr(*, "label")= chr "Country Code"
##   ..$ year : Ord.factor w/ 61 levels "1960"<"1961"<..: 1 2 3 4 5 6 7..
##   .. ..- attr(*, "label")= chr "Year"
```

```r
c(is_irregular(LIFEEXi), is_irregular(LIFEEXi[-5]))
```

```
## [1] FALSE  TRUE
```

```r
G(LIFEEXi[c(1:5, 7:10)])
```

```
## [1]    NA 1.590 1.544 1.494 1.448    NA 1.366 1.362 1.365
## 
## Indexed by:  iso3c [1] | year [9 (61)]
```

```r
###################################################
### code chunk number 23: Demonstrating Deep Indexation
###################################################
settransform(wldi, PCGDP_ld = Dlog(PCGDP))
lm(D(LIFEEX) ~ L(PCGDP_ld, 0:5) + B(PCGDP_ld), wldi) |>
  summary() |> coef() |> round(3)
```

```
##                    Estimate Std. Error t value Pr(>|t|)
## (Intercept)           0.299      0.007  44.412    0.000
## L(PCGDP_ld, 0:5)--    0.300      0.080   3.735    0.000
## L(PCGDP_ld, 0:5)L1    0.269      0.081   3.332    0.001
## L(PCGDP_ld, 0:5)L2    0.227      0.079   2.854    0.004
## L(PCGDP_ld, 0:5)L3    0.200      0.078   2.563    0.010
## L(PCGDP_ld, 0:5)L4    0.143      0.076   1.871    0.061
## L(PCGDP_ld, 0:5)L5    0.095      0.073   1.301    0.193
## B(PCGDP_ld)          -1.021      0.316  -3.234    0.001
```

```r
###################################################
### code chunk number 24: Using 3rd Party Functions: Rolling Average
###################################################
BY(LIFEEXi, findex(LIFEEXi)$iso3c, data.table::frollmean, 5) |> head(10)
```

```
##  [1]    NA    NA    NA    NA 33.46 33.96 34.46 34.95 35.43 35.92
## 
## Indexed by:  iso3c [1] | year [10 (61)]
```

```r
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
```

```
## full join: df1[id1, id2] 3/4 (75%) <m:m> df2[id1, id2] 3/4 (75%)
```

```
##   id1 id2 name age salary      dept   .join
## 1   1   a John  35  60000        IT matched
## 2   1   b Jane  28     NA      <NA>     df1
## 3   2   b  Bob  42  55000 Marketing matched
## 4   3   c Carl  50  70000     Sales matched
## 5   3   e <NA>  NA  80000        IT     df2
```

```r
###################################################
### code chunk number 26: Validation + Join Attribute
###################################################
join(df1, df2, on = c("id1", "id2"), validate = "1:1", attr = "join") |>
  attr("join") |> str(width = 70, strict = "cut")
```

```
## left join: df1[id1, id2] 3/4 (75%) <1:1> df2[id1, id2] 3/4 (75%)
## List of 3
##  $ call   : language join(x = df1, y = df2, on = c("id1", "id2"), v"..
##  $ on.cols:List of 2
##   ..$ x: chr [1:2] "id1" "id2"
##   ..$ y: chr [1:2] "id1" "id2"
##  $ match  : 'qG' int [1:4] 1 NA 2 3
##   ..- attr(*, "N.nomatch")= int 1
##   ..- attr(*, "N.groups")= int 4
##   ..- attr(*, "N.distinct")= int 3
```

```r
###################################################
### code chunk number 27: Overidentification Warning
###################################################
df2$name = df1$name
join(df1, df2) |> capture.output(type="m") |> strwrap(77) |> cat(sep="\n")
```

```
## Warning in fmatch(x[ixon], y[iyon], nomatch = NA_integer_, count = count, :
## Overidentified match/join: the first 2 of 3 columns uniquely match the records.
## With overid > 0, fmatch() continues to match columns. Consider removing columns
## or setting overid = 0 to terminate the algorithm after 2 columns (the results
## may differ, see ?fmatch). Alternatively set overid = 2 to silence this warning.
```

```
## left join: df1[id1, id2, name] 1/4 (25%) <m:m> df2[id1, id2, name] 1/4 (25%)
##   id1 id2 name age salary dept
## 1   1   a John  35  60000   IT
## 2   1   b Jane  28     NA <NA>
## 3   2   b  Bob  42     NA <NA>
## 4   3   c Carl  50     NA <NA>
```

```r
###################################################
### code chunk number 28: Automatic Renaming
###################################################
join(df1, df2, on = c("id1", "id2"))
```

```
## left join: df1[id1, id2] 3/4 (75%) <m:m> df2[id1, id2] 3/4 (75%)
## duplicate columns: name => renamed using suffix '_df2' for y
```

```
##   id1 id2 name age salary      dept name_df2
## 1   1   a John  35  60000        IT     John
## 2   1   b Jane  28     NA      <NA>     <NA>
## 3   2   b  Bob  42  55000 Marketing     Jane
## 4   3   c Carl  50  70000     Sales      Bob
```

```r
###################################################
### code chunk number 29: Data for Pivots
###################################################
data <- data.frame(type = rep(c("A", "B"), each = 2),
            type_name = rep(c("Apples", "Bananas"), each = 2),
            id = rep(1:2, 2), r = abs(rnorm(4)), h = abs(rnorm(4)*2))
setrelabel(data, id = "Fruit Id", r = "Fruit Radius", h = "Fruit Height")
print(data)
```

```
##   type type_name id       r      h
## 1    A    Apples  1 0.26103 0.4104
## 2    A    Apples  2 0.36128 1.0518
## 3    B   Bananas  1 1.04000 0.5072
## 4    B   Bananas  2 0.05624 2.5290
```

```r
vlabels(data)
```

```
##           type      type_name             id              r              h 
##             NA             NA     "Fruit Id" "Fruit Radius" "Fruit Height"
```

```r
###################################################
### code chunk number 30: Pivot Longer
###################################################
(dl <- pivot(data, ids = c("type", "type_name", "id"), labels = "label"))
```

```
##   type type_name id variable        label   value
## 1    A    Apples  1        r Fruit Radius 0.26103
## 2    A    Apples  2        r Fruit Radius 0.36128
## 3    B   Bananas  1        r Fruit Radius 1.04000
## 4    B   Bananas  2        r Fruit Radius 0.05624
## 5    A    Apples  1        h Fruit Height 0.41043
## 6    A    Apples  2        h Fruit Height 1.05181
## 7    B   Bananas  1        h Fruit Height 0.50718
## 8    B   Bananas  2        h Fruit Height 2.52899
```

```r
vlabels(dl)
```

```
##       type  type_name         id   variable      label      value 
##         NA         NA "Fruit Id"         NA         NA         NA
```

```r
###################################################
### code chunk number 31: Pivot Wider
###################################################
(dw <- pivot(data, "id", names = "type", labels = "type_name", how = "w"))
```

```
##   id    r_A     r_B    h_A    h_B
## 1  1 0.2610 1.04000 0.4104 0.5072
## 2  2 0.3613 0.05624 1.0518 2.5290
```

```r
namlab(dw)
```

```
##   Variable                  Label
## 1       id               Fruit Id
## 2      r_A  Fruit Radius - Apples
## 3      r_B Fruit Radius - Bananas
## 4      h_A  Fruit Height - Apples
## 5      h_B Fruit Height - Bananas
```

```r
###################################################
### code chunk number 32: Pivot Recast
###################################################
(dr <- pivot(data, ids = "id", names = list(from = "type"),
             labels = list(from = "type_name", to = "label"), how = "r"))
```

```
##   id variable        label      A       B
## 1  1        r Fruit Radius 0.2610 1.04000
## 2  2        r Fruit Radius 0.3613 0.05624
## 3  1        h Fruit Height 0.4104 0.50718
## 4  2        h Fruit Height 1.0518 2.52899
```

```r
vlabels(dr)
```

```
##         id   variable      label          A          B 
## "Fruit Id"         NA         NA   "Apples"  "Bananas"
```

```r
###################################################
### code chunk number 33: Recursive Splitting: Creates Nested List of Data Frames
###################################################
(dl <- mtcars |> rsplit(mpg + hp + carb ~ vs + am)) |> str(max.level = 2)
```

```
## List of 2
##  $ 0:List of 2
##   ..$ 0:'data.frame':	12 obs. of  3 variables:
##   ..$ 1:'data.frame':	6 obs. of  3 variables:
##  $ 1:List of 2
##   ..$ 0:'data.frame':	7 obs. of  3 variables:
##   ..$ 1:'data.frame':	7 obs. of  3 variables:
```

```r
###################################################
### code chunk number 34: Fitting Linear Models and Obtaining Coefficient Matrices
###################################################
nest_lm_coef <- dl |>
  rapply2d(lm, formula = mpg ~ .) |>
  rapply2d(summary, classes = "lm") |>
  get_elem("coefficients")

nest_lm_coef |> str(give.attr = FALSE, strict = "cut")
```

```
## List of 2
##  $ 0:List of 2
##   ..$ 0: num [1:3, 1:4] 15.8791 0.0683 -4.5715 3.655 0.0345 ...
##   ..$ 1: num [1:3, 1:4] 26.9556 -0.0319 -0.308 2.293 0.0149 ...
##  $ 1:List of 2
##   ..$ 0: num [1:3, 1:4] 30.896903 -0.099403 -0.000332 3.346033 0.03587 ...
##   ..$ 1: num [1:3, 1:4] 37.0012 -0.1155 0.4762 7.3316 0.0894 ...
```

```r
###################################################
### code chunk number 35: Unlisting to Data Frame
###################################################
nest_lm_coef |> unlist2d(c("vs", "am"), row.names = "variable") |> head(2)
```

```
##   vs am    variable Estimate Std. Error t value Pr(>|t|)
## 1  0  0 (Intercept) 15.87915    3.65495   4.345 0.001865
## 2  0  0          hp  0.06832    0.03449   1.981 0.078938
```

```r
###################################################
### code chunk number 36: Removing Generated Series (Hidden)
###################################################
wldi <- wldi[1:13]

###################################################
### code chunk number 37: Which Columns/Countries have Time Varying Information?
###################################################
varying(wldi)
```

```
## country    date    year  decade  region  income    OECD   PCGDP  LIFEEX    GINI 
##   FALSE    TRUE    TRUE    TRUE   FALSE   FALSE   FALSE    TRUE    TRUE    TRUE 
##     ODA     POP 
##    TRUE    TRUE
```

```r
varying(wldi, any_group = FALSE) |> head(3)
```

```
##     country date year decade region income  OECD PCGDP LIFEEX GINI  ODA  POP
## ABW   FALSE TRUE TRUE   TRUE  FALSE  FALSE FALSE  TRUE   TRUE   NA TRUE TRUE
## AFG   FALSE TRUE TRUE   TRUE  FALSE  FALSE FALSE  TRUE   TRUE   NA TRUE TRUE
## AGO   FALSE TRUE TRUE   TRUE  FALSE  FALSE FALSE  TRUE   TRUE TRUE TRUE TRUE
```

```r
###################################################
### code chunk number 38: Demonstrating Panel Decomposition
###################################################
all.equal(fvar(W(LIFEEXi)) + fvar(B(LIFEEXi)), fvar(LIFEEXi))
```

```
## [1] TRUE
```

```r
###################################################
### code chunk number 39: Panel Summary Statistics
###################################################
qsu(LIFEEXi)
```

```
##              N/T     Mean       SD      Min      Max
## Overall    11670  64.2963  11.4764   18.907  85.4171
## Between      207  64.9537   9.8936  40.9663  85.4171
## Within   56.3768  64.2963   6.0842  32.9068  84.4198
```

```r
###################################################
### code chunk number 40: Weighted Panel Summary Statistics by Groups
###################################################
qsu(LIFEEXi, g = wlddev$OECD, w = wlddev$POP, higher = TRUE) |> aperm()
```

```
## , , FALSE
## 
##              N/T     Mean      SD      Min      Max     Skew    Kurt
## Overall     9503  63.5476  9.2368   18.907  85.4171  -0.7394  2.7961
## Between      171  63.5476  6.0788  43.0905  85.4171  -0.8041   3.082
## Within   55.5731  65.8807  6.9545  30.3388  82.8832  -1.0323  4.0998
## 
## , , TRUE
## 
##              N/T     Mean      SD      Min      Max     Skew    Kurt
## Overall     2156  74.9749  5.3627   45.369  84.3563  -1.2966  6.5505
## Between       36  74.9749  2.9256  66.2983  78.6733  -1.3534  4.5999
## Within   59.8889  65.8807  4.4944  44.9513  77.2733   -0.627  3.9839
```

```r
###################################################
### code chunk number 41: Detailed (Grouped, Weighted) Statistical Description
###################################################
descr(wlddev, LIFEEX ~ OECD, w = ~ replace_na(POP))
```

```
## Dataset: wlddev, 1 Variables, N = 13176, WeightSum = 313233706778
## Grouped by: OECD [2]
##            N   Perc       WeightSum  Perc
## FALSE  10980  83.33  2.49344474e+11  79.6
## TRUE    2196  16.67  6.38892329e+10  20.4
## --------------------------------------------------------------------------------
## LIFEEX (numeric): Life expectancy at birth, total (years)
## Statistics (N = 11659, 11.51% NAs)
##           N   Perc  Ndist   Mean    SD    Min    Max   Skew  Kurt
## FALSE  9503  81.51   8665  63.55  9.24  18.91  85.42  -0.74   2.8
## TRUE   2156  18.49   2016  74.97  5.36  45.37  84.36   -1.3  6.55
## 
## Quantiles
##           1%     5%    10%    25%    50%    75%    90%    95%    99%
## FALSE  41.32  37.63  48.98   57.5  65.87  69.68   74.1  32.39  76.18
## TRUE   56.67  65.65  69.69  71.84  75.32  78.61  81.26  81.23  83.59
## --------------------------------------------------------------------------------
```

```r
###################################################
### code chunk number 42: qtab: Basic Usage
###################################################
library(magrittr) # World after 2015 (latest country data)
wlda15 <- wlddev |> fsubset(year >= 2015) |> fgroup_by(iso3c) |> flast()
wlda15 %$% qtab(OECD, income)
```

```
##        income
## OECD    High income Low income Lower middle income Upper middle income
##   FALSE          45         30                  47                  58
##   TRUE           34          0                   0                   2
```

```r
###################################################
### code chunk number 43: qtab: Population Counts
###################################################
wlda15 %$% qtab(OECD, income, w = POP) %>% divide_by(1e6)
```

```
##        income
## OECD    High income Low income Lower middle income Upper middle income
##   FALSE       93.01     694.89             3063.54             2459.71
##   TRUE      1098.75       0.00                0.00              211.01
```

```r
###################################################
### code chunk number 44: qtab: Average Life Expectancy
###################################################
wlda15 %$% qtab(OECD, income, w = LIFEEX, wFUN = fmean) %>% replace_na(0)
```

```
##        income
## OECD    High income Low income Lower middle income Upper middle income
##   FALSE       78.75      62.81               68.30               73.81
##   TRUE        81.09       0.00                0.00               76.37
```

```r
###################################################
### code chunk number 45: qtab: Population Weighted Average Life Expectancy
###################################################
wlda15 %$% qtab(OECD, income, w = LIFEEX, wFUN = fmean,
                wFUN.args = list(w = POP)) %>% replace_na(0)
```

```
##        income
## OECD    High income Low income Lower middle income Upper middle income
##   FALSE       77.91      63.81               68.76               75.93
##   TRUE        81.13       0.00                0.00               76.10
```

```r
###################################################
### code chunk number 46: Benchmark: Statistics and Data Manipulation
###################################################
set_collapse(na.rm = FALSE, sort = FALSE, nthreads = 1)
set.seed(101)
m <- matrix(rnorm(1e7), ncol = 1000)
data <- qDT(replicate(100, rnorm(1e5), simplify = FALSE))
g <- sample.int(1e4, 1e5, TRUE)

microbenchmark(R = colMeans(m),
               Rfast = Rfast::colmeans(m),
               collapse = fmean(m))
```

```
## Warning in microbenchmark(R = colMeans(m), Rfast = Rfast::colmeans(m), collapse
## = fmean(m)): less accurate nanosecond times to avoid potential integer
## overflows
```

```
## Unit: milliseconds
##      expr   min    lq  mean median    uq    max neval
##         R 9.493 9.662 9.913  9.842 9.906 13.703   100
##     Rfast 4.738 4.801 5.364  4.919 4.979 40.632   100
##  collapse 1.402 1.569 1.869  1.677 1.884  7.335   100
```

```r
microbenchmark(R = rowsum(data, g, reorder = FALSE),
               data.table = data[, lapply(.SD, sum), by = g],
               collapse = fsum(data, g))
```

```
## Unit: milliseconds
##        expr    min     lq   mean median     uq   max neval
##           R 10.555 11.001 12.002 11.126 11.367 55.73   100
##  data.table 15.797 16.532 17.566 16.705 17.036 51.62   100
##    collapse  5.099  5.464  6.131  5.577  5.864 12.85   100
```

```r
add_vars(data) <- g
microbenchmark(data.table = data[, lapply(.SD, median), by = g],
               collapse = data |> fgroup_by(g) |> fmedian())
```

```
## Unit: milliseconds
##        expr   min    lq  mean median    uq   max neval
##  data.table 129.2 130.7 132.2  131.4 132.4 142.4   100
##    collapse 116.9 118.0 119.2  118.7 119.7 127.4   100
```

```r
d <- data.table(g = unique(g), x = 1, y = 2, z = 3)
microbenchmark(data.table = d[data, on = "g"],
               collapse = join(data, d, on = "g", verbose = 0))
```

```
## Unit: milliseconds
##        expr    min     lq   mean median     uq    max neval
##  data.table 20.956 24.659 36.358 26.198 36.002 78.567   100
##    collapse  1.225  1.348  1.384  1.378  1.412  1.798   100
```

```r
microbenchmark(data.table = melt(data, "g"),
               collapse = pivot(data, "g"))
```

```
## Unit: milliseconds
##        expr   min    lq  mean median    uq   max neval
##  data.table 11.45 14.80 21.52  16.17 17.72 76.84   100
##    collapse 11.49 13.36 21.95  16.14 17.46 65.15   100
```

```r
settransform(data, id = rowid(g))
cols = grep("^V", names(data), value = TRUE)
microbenchmark(data.table = dcast(data, g ~ id, value.var = cols),
          collapse = pivot(data, ids = "g", names = "id", how = "w"))
```

```
## Unit: milliseconds
##        expr   min     lq  mean median     uq   max neval
##  data.table 68.30 107.68 111.4 111.78 116.48 133.2   100
##    collapse 35.18  76.12  79.4  78.34  83.67 105.3   100
```

```r
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
```

```
## Unit: milliseconds
##           expr    min     lq  mean median     uq    max neval
##       base_int 60.026 61.009 63.98 63.811 65.181 113.13   100
##   collapse_int  8.494  9.153 10.04  9.311  9.787  16.91   100
##      base_char 91.912 94.418 96.91 96.118 97.607 143.70   100
##  collapse_char 21.173 22.595 23.30 22.929 23.318  30.64   100
```

```r
microbenchmark(base_int = match(g_int, 1:1000),
               collapse_int = fmatch(g_int, 1:1000),
               base_char = match(g_char, char),
               data.table_char = chmatch(g_char, char),
               collapse_char = fmatch(g_char, char), times = 10)
```

```
## Unit: milliseconds
##             expr   min     lq   mean median      uq     max neval
##         base_int 26.69 26.809 28.908 27.055  32.516  33.839    10
##     collapse_int  8.59  8.733  8.812  8.787   8.859   9.248    10
##        base_char 94.55 94.624 98.492 96.228 101.662 108.488    10
##  data.table_char 41.42 41.516 43.104 41.709  42.128  48.874    10
##    collapse_char 36.57 36.654 37.747 36.792  37.865  44.645    10
```

```r
###################################################
### Print Session Information
###################################################

sessionInfo()
```

```
## R version 4.3.0 (2023-04-21)
## Platform: aarch64-apple-darwin20 (64-bit)
## Running under: macOS Ventura 13.4.1
## 
## Matrix products: default
## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: Europe/Berlin
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] collapse_2.0.8        microbenchmark_1.4.10 magrittr_2.0.3       
## [4] data.table_1.14.10   
## 
## loaded via a namespace (and not attached):
##  [1] Formula_1.2-5       numDeriv_2016.8-1.1 xfun_0.39          
##  [4] fixest_0.11.2       lattice_0.21-8      RcppZiggurat_0.1.6 
##  [7] zoo_1.8-12          knitr_1.43          parallel_4.3.0     
## [10] RcppParallel_5.1.7  dreamerr_1.3.0      sandwich_3.1-0     
## [13] grid_4.3.0          compiler_4.3.0      rstudioapi_0.14    
## [16] tools_4.3.0         Rfast_2.1.0         nlme_3.1-162       
## [19] evaluate_0.21       Rcpp_1.0.11
```

