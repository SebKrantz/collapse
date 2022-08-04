reprex::reprex({

# System: M1 MAC (2020), 16GB, CRAN Binaries
library(fastverse)
fastverse_extend(dplyr, arrow, microbenchmark)

# Data From: https://www.kaggle.com/datasets/neilclack/nyc-taxi-trip-data-google-public-data?resource=download
system.time(taxi <- fread("misc/NYC_TAXI_2019/taxi_trip_data.csv", tz = "UTC"))
system.time(taxi_arrow <- read_csv_arrow("misc/NYC_TAXI_2019/taxi_trip_data.csv"))
taxi_arrow %<>% arrow_table()

fdim(taxi)
fndistinct(taxi)

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

# ~2500 groups, just aggregation: here arrow is the fastest
taxi |> fselect(vendor_id, pickup_location_id, passenger_count) |> fnunique()

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

# Now aggregating across 210638 groups
taxi |> fselect(vendor_id, pickup_location_id, dropoff_location_id, passenger_count) |> fnunique()

# This is where you do not want to use dplyr
# dplyr:
system.time(taxi |>
    group_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count) |>
    summarise(across(c(total_amount, fare_amount, tip_amount, tolls_amount, mta_tax, imp_surcharge, extra), sum, na.rm = TRUE),
              .groups = "drop"))

# arrow, data.table and collapse, all unordered
microbenchmark(
   arrow = taxi_arrow |>
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

# Distinct values
microbenchmark(
  arrow = taxi_arrow |>
    group_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count) |>
    summarise(total_amount = n_distinct(total_amount)) |> collect(),

  collapse = taxi |>
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount) |> fndistinct(na.rm = FALSE)
, times = 10)

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

# Advanced operations not supported by arrow -----------------------------------

# Whenever an expression is not supported you get a message like this and dplyr does the computation
taxi_arrow |>
  group_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count) |>
  mutate(total_amount = cumsum(total_amount)) |> collect()

# Cumsum
microbenchmark(
  data.table = taxi[, cs := cumsum(total_amount), by = .(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)],
  collapse = taxi |> fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
              fmutate(total_amount = fcumsum(total_amount))
, times = 10)

# Centering
microbenchmark(
  data.table = taxi[, center := total_amount - mean(total_amount), by = .(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)],
  collapse = settransform(taxi, center = fwithin(total_amount, list(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)))
, times = 10)

# Weighted mean
microbenchmark(
  data.table = taxi[, .(total_amount = weighted.mean(total_amount, trip_distance)), by = .(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)],
  collapse = taxi |>
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount, trip_distance) |> fmean(trip_distance, keep.w = FALSE, na.rm = FALSE)
, times = 10)

# Scaling
system.time(
  settransform(taxi, scaled = fscale(total_amount, list(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)))
)
# Mode
system.time(
  taxi |>
  fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
  fselect(total_amount) |> fmode(na.rm = FALSE)
)
# Weighted median
system.time(
  taxi |>
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount, trip_distance) |> fmedian(trip_distance, keep.w = FALSE, na.rm = FALSE)
)
# Weighted mode
system.time(
  taxi |>
    fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
    fselect(total_amount, fare_amount, trip_distance) |> fmode(trip_distance, keep.w = FALSE, na.rm = FALSE)
)
# Growth Rate
system.time(
  settransform(taxi, growth = fgrowth(total_amount, g = list(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)))
)
# Rolling Avareages
system.time(
  taxi[, rollavg := frollmean(total_amount, 10), by = .(vendor_id, pickup_location_id, dropoff_location_id, passenger_count)]
)
# Linear Regression
system.time(
taxi |>
  fgroup_by(vendor_id, pickup_location_id, dropoff_location_id, passenger_count, sort = FALSE) |>
  fmutate(dmx = fwithin(total_amount)) |> fsummarise(slope = fsum(dmx, trip_distance) %/=% fsum(dmx^2))
)

}, wd = getwd())

# Conclusion:
# - Arrow is an impressive library providing cutting edge performance
# - Not faster than other cutting edge C/C++ libraries, especially with many groups
# - No support for many statistical operations data.table and collapse do effectively
# - Provides some cool tools for larger-than memory data (parquet files)!!
# - Provides zero-conversion data sharing between different languages!!






