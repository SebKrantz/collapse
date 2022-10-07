
  ``` r
# System: M1 MAC (2020), 16GB, CRAN Binaries
library(fastverse)
#> -- Attaching packages --------------------------------------- fastverse 0.2.4 --
#> collapse       1.8.7      testthat       3.1.4
#> magrittr       2.0.3      microbenchmark 1.4.9
#> data.table     1.14.2     kit            0.0.12
#> -- Conflicts ------------------------------------------ fastverse_conflicts() --
#> testthat::equals()       masks magrittr::equals()
#> testthat::is_less_than() masks magrittr::is_less_than()
#> testthat::not()          masks magrittr::not()
fastverse_extend(dplyr, arrow, microbenchmark)
#> -- Attaching extension packages ----------------------------- fastverse 0.2.4 --
#> dplyr 1.0.9     arrow 8.0.0
#> -- Conflicts ------------------------------------------ fastverse_conflicts() --
#> dplyr::between()   masks data.table::between()
#> dplyr::filter()    masks stats::filter()
#> dplyr::first()     masks data.table::first()
#> dplyr::intersect() masks base::intersect()
#> arrow::is_in()     masks magrittr::is_in()
#> dplyr::lag()       masks stats::lag()
#> dplyr::last()      masks data.table::last()
#> arrow::matches()   masks dplyr::matches(), testthat::matches()
#> dplyr::setdiff()   masks base::setdiff()
#> dplyr::setequal()  masks base::setequal()
#> arrow::timestamp() masks utils::timestamp()
#> dplyr::union()     masks base::union()

# Data From: https://www.kaggle.com/datasets/neilclack/nyc-taxi-trip-data-google-public-data?resource=download
system.time(taxi <- fread("misc/NYC_TAXI_2019/taxi_trip_data.csv", tz = "UTC"))
#>    user  system elapsed
#>   3.578   0.498   7.767
system.time(taxi_arrow <- read_csv_arrow("misc/NYC_TAXI_2019/taxi_trip_data.csv"))
#>    user  system elapsed
#>   6.810   1.772   2.001
taxi_arrow %<>% arrow_table()

fdim(taxi)
#> [1] 10000000       17
fndistinct(taxi)
#>           vendor_id     pickup_datetime    dropoff_datetime     passenger_count
#>                   3             7824766             7823041                  10
#>       trip_distance           rate_code  store_and_fwd_flag        payment_type
#>                6502                   7                   2                   5
#>         fare_amount               extra             mta_tax          tip_amount
#>                8286                  61                  24                5919
#>        tolls_amount       imp_surcharge        total_amount  pickup_location_id
#>                2463                  12               20933                 262
#> dropoff_location_id
#>                 263

# A benchmark shared on Twitter, some manipulation and aggregation across 262 groups...
# In this case even dplyr is efficient.
microbenchmark(

  dplyr = taxi |> qTBL() |>
    select(trip_distance, passenger_count, total_amount, pickup_location_id) |>
    filter(trip_distance >= 5, passenger_count > 1) |>
    mutate(total_amount_per_head = as.numeric(total_amount) / passenger_count) |>
    group_by(pickup_location_id) |>
    summarise(mean_amount = mean(total_amount_per_head, na.rm = TRUE)),

  arrow = taxi_arrow |>
    select(trip_distance, passenger_count, total_amount, pickup_location_id) |>
    filter(trip_distance >= 5, passenger_count > 1) |>
    mutate(total_amount_per_head = as.numeric(total_amount) / passenger_count) |>
    group_by(pickup_location_id) |>
    summarise(mean_amount = mean(total_amount_per_head, na.rm = TRUE)) |>
    collect(),

  data.table = taxi[trip_distance >= 5 & passenger_count > 1, .(passenger_count, total_amount, pickup_location_id)
  ][, total_amount_per_head := as.numeric(total_amount) / passenger_count
  ][, .(mean_amount = mean(total_amount_per_head, na.rm = TRUE)), keyby = pickup_location_id],

  collapse = taxi |>
    fsubset(trip_distance >= 5 & passenger_count > 1, passenger_count, total_amount, pickup_location_id) |>
    fmutate(total_amount_per_head = as.numeric(total_amount) %/=% passenger_count) |>
    fgroup_by(pickup_location_id) |>
    fsummarise(mean_amount = fmean(total_amount_per_head))

  , times = 30)
#> Warning in microbenchmark(dplyr =
#> summarise(group_by(mutate(filter(select(qTBL(taxi), : less accurate nanosecond
#> times to avoid potential integer overflows
#> Unit: milliseconds
#>        expr       min        lq     mean    median       uq      max neval
#>       dplyr 138.78168 147.13268 166.3310 158.06453 177.0337 262.6299    30
#>       arrow  89.02100  92.12376 104.2015  97.74847 111.5167 183.6385    30
#>  data.table 145.05644 153.12032 161.8363 162.31650 167.7433 192.5772    30
#>    collapse  92.60695  99.46748 108.7406 103.27551 119.7872 129.7093    30

# ~2500 groups, just aggregation: here arrow is the fastest
taxi |> fselect(vendor_id, pickup_location_id, passenger_count) |> fnunique()
#> [1] 3498

microbenchmark(

  dplyr = taxi |> qTBL() |>
    group_by(vendor_id, pickup_location_id, passenger_count) |>
    summarise(across(c(total_amount, fare_amount, tip_amount, tolls_amount, mta_tax, imp_surcharge), sum, na.rm = TRUE),
              .groups = "drop"),

  arrow_unordered = taxi_arrow |>
    group_by(vendor_id, pickup_location_id, passenger_count) |> # Note that unlike dplyr::group_by(), arrow::group_by() is unsorted!!!
    summarise(total_amount = sum(total_amount), # No support for across()
              fare_amount = sum(fare_amount),
              tip_amount = sum(tip_amount),
              tolls_amount = sum(tolls_amount),
              mta_tax = sum(mta_tax),
              imp_surcharge = sum(imp_surcharge)) |> collect(),

  data.table = taxi[, lapply(.SD, sum, na.rm = TRUE), keyby = .(vendor_id, pickup_location_id, passenger_count),
                    .SDcols = .c(total_amount, fare_amount, tip_amount, tolls_amount, mta_tax, imp_surcharge)],

  data.table_unordered = taxi[, lapply(.SD, sum, na.rm = TRUE), by = .(vendor_id, pickup_location_id, passenger_count),
                              .SDcols = .c(total_amount, fare_amount, tip_amount, tolls_amount, mta_tax, imp_surcharge)],

  collapse = taxi |> fgroup_by(vendor_id, pickup_location_id, passenger_count) |>
    fselect(total_amount, fare_amount, tip_amount, tolls_amount, mta_tax, imp_surcharge) |> fsum(),

  collapse_unordered = taxi |> fgroup_by(vendor_id, pickup_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount, fare_amount, tip_amount, tolls_amount, mta_tax, imp_surcharge) |> fsum()

  , times = 10)
#> Unit: milliseconds
#>                  expr        min         lq       mean     median         uq
#>                 dplyr 3977.57625 4097.62905 4237.63773 4260.53583 4313.47507
#>       arrow_unordered   46.04927   46.72266   51.10864   49.55049   54.07297
#>            data.table  409.12969  417.89959  427.07873  425.71030  436.93499
#>  data.table_unordered  488.22997  494.55733  522.56507  515.74851  529.74206
#>              collapse  262.96941  263.92040  268.45868  266.35814  272.59494
#>    collapse_unordered   99.97936  100.68944  103.97961  102.56460  106.70406
#>        max neval
#>  4511.4170    10
#>    59.4933    10
#>   442.9845    10
#>   607.2070    10
#>   279.4542    10
#>   110.6758    10

# Now aggregating across 210638 groups
taxi |> fselect(vendor_id, pickup_location_id, dropoff_location_id, passenger_count) |> fnunique()
#> [1] 210638

# This is where you do not want to use dplyr
# dplyr:
system.time(taxi |>
              group_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count) |>
              summarise(across(c(total_amount, fare_amount, tip_amount, tolls_amount, mta_tax, imp_surcharge, extra), sum, na.rm = TRUE),
                        .groups = "drop"))
#>    user  system elapsed
#> 239.852   1.179 242.261

# arrow, data.table and collapse, all unordered
microbenchmark(
  taxi_arrow |>
    group_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count) |>
    summarise(total_amount = sum(total_amount),
              fare_amount = sum(fare_amount),
              tip_amount = sum(tip_amount),
              tolls_amount = sum(tolls_amount),
              mta_tax = sum(mta_tax),
              imp_surcharge = sum(imp_surcharge)) |> collect(),

  data.table = taxi[, lapply(.SD, sum), by = .(vendor_id, pickup_location_id, dropoff_location_id, passenger_count),
                    .SDcols = .c(total_amount, fare_amount, tip_amount, tolls_amount, mta_tax, imp_surcharge, extra)],

  collapse = taxi |>
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount, fare_amount, tip_amount, tolls_amount, mta_tax, imp_surcharge, extra) |> fsum(na.rm = FALSE)

  , times = 10)
#> Unit: milliseconds
#>       expr       min       lq     mean   median       uq      max neval
#>      arrow  300.7408 308.8974 376.9485 354.0625 367.5049 638.2129    10
#> data.table  618.6634 630.2733 671.7836 649.6240 659.3185 876.6953    10
#>   collapse  241.4413 248.5718 268.7940 253.6805 290.6394 327.4196    10

# Advanced operations, with many groups ----------------------------------------

# Median
microbenchmark(
  arrow = taxi_arrow |>
    group_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count) |>
    summarise(total_amount = median(total_amount)) |> collect(),

  data.table = taxi[, .(total_amount = median(total_amount)), by = .(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)],
  collapse = taxi |>
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount) |> fmedian(na.rm = FALSE)
  , times = 2)
#> Warning: median() currently returns an approximate median in Arrow
#> This warning is displayed once per session.
#> Unit: milliseconds
#>        expr       min        lq      mean    median        uq       max neval
#>       arrow 8141.4202 8141.4202 8507.1152 8507.1152 8872.8102 8872.8102     2
#>  data.table  710.5617  710.5617  713.9699  713.9699  717.3781  717.3781     2
#>    collapse  360.0937  360.0937  373.7713  373.7713  387.4490  387.4490     2

# Distinct values
microbenchmark(
  arrow = taxi_arrow |>
    group_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count) |>
    summarise(total_amount = n_distinct(total_amount)) |> collect(),

  collapse = taxi |>
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount) |> fndistinct(na.rm = FALSE)
  , times = 10)
#> Unit: milliseconds
#>      expr      min       lq     mean   median       uq      max neval
#>     arrow 488.6324 495.6366 514.8036 505.6181 532.6775 585.8691    10
#>  collapse 332.2395 335.4083 342.7865 340.0869 346.9045 368.6868    10

# Standard Deviation
microbenchmark(
  arrow = taxi_arrow |>
    group_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count) |>
    summarise(total_amount = sd(total_amount)) |> collect(),

  data.table = taxi[, .(total_amount = sd(total_amount)), by = .(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)],

  collapse = taxi |>
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount) |> fsd(na.rm = FALSE),

  collapse_unstable = taxi |> # Sum of squares algorithm using long doubles, unstable in extreme cases
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount) |> fsd(stable = FALSE, na.rm = FALSE)
  , times = 10)
#> Unit: milliseconds
#>               expr      min       lq     mean   median       uq      max neval
#>              arrow 469.0867 488.8644 505.5929 502.9060 522.7925 545.1278    10
#>         data.table 768.7327 772.1942 806.6084 797.7831 831.6306 875.3769    10
#>           collapse 217.0025 219.3173 229.1510 225.0078 235.3567 257.0124    10
#>  collapse_unstable 206.4531 212.1699 223.8805 219.4320 237.1253 254.4761    10

# Advanced operations not supported by arrow -----------------------------------

# Whenever an expression is not supported you get a message like this and dplyr does the computation
taxi_arrow |>
  group_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count) |>
  mutate(total_amount = cumsum(total_amount)) |> collect()
#> Warning: Expression cumsum(total_amount) not supported in Arrow; pulling data
#> into R
#> # A tibble: 10,000,000 × 17
#> # Groups:   vendor_id, pickup_location_id, dropoff_location_id, passenger_count
#> #   [210,638]
#>    vendor_id pickup_datetime     dropoff_datetime    passenger_count
#>        <int> <dttm>              <dttm>                        <int>
#>  1         2 2018-02-02 11:18:23 2018-02-02 11:58:09               1
#>  2         1 2018-02-02 11:43:58 2018-02-02 12:08:49               1
#>  3         1 2018-02-02 13:50:55 2018-02-02 14:18:37               2
#>  4         1 2018-02-02 13:04:12 2018-02-02 13:38:32               1
#>  5         2 2018-02-02 14:42:41 2018-02-02 14:57:12               1
#>  6         1 2018-02-02 14:36:20 2018-02-02 15:14:20               2
#>  7         2 2018-02-02 15:35:16 2018-02-02 15:37:08               1
#>  8         2 2018-09-17 10:53:42 2018-09-17 11:33:22               2
#>  9         2 2018-09-17 10:24:50 2018-09-17 10:53:54               2
#> 10         1 2018-09-17 11:28:05 2018-09-17 12:02:10               1
#> # … with 9,999,990 more rows, and 13 more variables: trip_distance <dbl>,
#> #   rate_code <int>, store_and_fwd_flag <chr>, payment_type <int>,
#> #   fare_amount <dbl>, extra <dbl>, mta_tax <dbl>, tip_amount <dbl>,
#> #   tolls_amount <dbl>, imp_surcharge <dbl>, total_amount <dbl>,
#> #   pickup_location_id <int>, dropoff_location_id <int>

# Cumsum
microbenchmark(
  data.table = taxi[, cs := cumsum(total_amount), by = .(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)],
  collapse = taxi |> fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fmutate(total_amount = fcumsum(total_amount))
  , times = 10)
#> Unit: milliseconds
#>        expr      min       lq     mean   median       uq      max neval
#>  data.table 690.2950 694.1321 707.3646 700.2476 707.9255 772.1126    10
#>    collapse 208.1338 210.3989 218.1475 211.3972 223.1234 255.8180    10

# Centering
microbenchmark(
  data.table = taxi[, center := total_amount - mean(total_amount), by = .(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)],
  collapse = settransform(taxi, center = fwithin(total_amount, list(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)))
  , times = 10)
#> Unit: milliseconds
#>        expr      min       lq     mean   median       uq       max neval
#>  data.table 967.1464 972.2606 981.1940 975.9247 989.6332 1005.4241    10
#>    collapse 200.6609 204.1434 205.9979 205.3107 208.7906  209.4111    10

# Weighted mean
microbenchmark(
  data.table = taxi[, .(total_amount = weighted.mean(total_amount, trip_distance)), by = .(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)],
  collapse = taxi |>
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount, trip_distance) |> fmean(trip_distance, keep.w = FALSE, na.rm = FALSE)
  , times = 10)
#> Unit: milliseconds
#>        expr       min       lq      mean    median       uq       max neval
#>  data.table 1131.0776 1139.699 1153.5972 1149.1631 1156.582 1210.7012    10
#>    collapse  206.2943  207.036  211.6297  209.2618  211.827  233.6085    10

# Scaling
system.time(
  settransform(taxi, scaled = fscale(total_amount, list(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)))
)
#>    user  system elapsed
#>   0.220   0.033   0.254
# Mode
system.time(
  taxi |>
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount) |> fmode(na.rm = FALSE)
)
#>    user  system elapsed
#>   0.370   0.034   0.404
# Weighted median
system.time(
  taxi |>
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount, trip_distance) |> fmedian(trip_distance, keep.w = FALSE, na.rm = FALSE)
)
#>    user  system elapsed
#>   0.685   0.041   0.725
# Weighted mode
system.time(
  taxi |>
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount, fare_amount, trip_distance) |> fmode(trip_distance, keep.w = FALSE, na.rm = FALSE)
)
#>    user  system elapsed
#>   0.612   0.086   0.699
# Growth Rate
system.time(
  settransform(taxi, growth = fgrowth(total_amount, g = list(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)))
)
#>    user  system elapsed
#>   0.185   0.032   0.216
# Rolling Avareages
system.time(
  taxi[, rollavg := frollmean(total_amount, 10), by = .(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)]
)
#>    user  system elapsed
#>   2.952   0.153   3.106
# Linear Regression
system.time(
  taxi |>
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fmutate(dmx = fwithin(total_amount)) |> fsummarise(slope = fsum(dmx, trip_distance) %/=% fsum(dmx^2))
)
#>    user  system elapsed
#>   0.219   0.037   0.256
```

<sup>Created on 2022-08-04 by the [reprex package](https://reprex.tidyverse.org) (v2.0.1)</sup>
