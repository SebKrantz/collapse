# System: M1 MacBook Pro (2020), 16GB Unified Memory, 1TB SSD
library(collapse)
library(data.table)
char <- c(letters, LETTERS, month.abb, month.name)
g1 <- sample(outer(char, char, paste0), 1e8, replace = TRUE)
g2 <- sample.int(1e4/2, 1e8, replace = TRUE)
fndistinct(g1)
fndistinct(g2)

data <- qDT(list(g1 = g1, g2 = g2,
                 x1 = na_insert(rnorm(1e8, mean = 10), prop = 0.1),
                 x2 = rnorm(1e8, mean = 100, sd = 10),
                 x3 = na_insert(rnorm(1e8), prop = 0.05)))

system.time(with(data, x1 %+=% x3)) # math by reference (also supports matrices and data frames)
system.time(data[, x1 := x1 + x3])  # data.table does it internally

system.time(fsubset(data, g1 %==% "hh")) # %==% returns indices directly
system.time(data[g1 == "hh"]) # data.table has optimized this internally

# TODO: mclapply
gc()
# Unsorted (first-appearance order)
system.time(data |> gby(g1, sort = FALSE) |> slt(x1:x3) |> fsum())
system.time(data[, lapply(.SD, sum, na.rm = TRUE), by = g1, .SDcols = x1:x3])
system.time(data |> gby(g2, sort = FALSE) |> smr(acr(x1:x3, fsum))) # new fsummarise(across()) interface
system.time(data[, lapply(.SD, sum, na.rm = TRUE), by = g2, .SDcols = x1:x3])
system.time(data |> gby(g2, g1, sort = FALSE) |> fsum()) # 27 million groups
system.time(data[, lapply(.SD, sum, na.rm = TRUE), by = .(g2, g1), .SDcols = x1:x3])
gc()
# Sorted
system.time(data |> gby(g1) |> slt(x1:x3) |> fsum())
system.time(data[, lapply(.SD, sum, na.rm = TRUE), keyby = g1, .SDcols = x1:x3])
system.time(data |> gby(g2) |> smr(acr(x1:x3, fsum)))
system.time(data[, lapply(.SD, sum, na.rm = TRUE), keyby = g2, .SDcols = x1:x3])
system.time(data |> gby(g2, g1) |> fsum())
system.time(data[, lapply(.SD, sum, na.rm = TRUE), keyby = .(g2, g1), .SDcols = x1:x3])
gc()
# Other stats
system.time(data |> gby(g1, sort = FALSE) |> slt(x1:x3) |> fmean())
system.time(data |> gby(g1, sort = FALSE) |> slt(x1:x3) |> fmean(x2, keep.w = FALSE))
system.time(data |> gby(g1, sort = FALSE) |> slt(x1:x3) |> fmedian())

# Notes:
# - OpenMP parallelism doesn't work on MAC, so this benchmark compares serial code speed.
# - You'll get faster data.table speeds on 4+ core windows or linux machines.
# - collapse supports equally fast grouped computations on vectors and matrices.

d <- replicate(1000, na_insert(rnorm(1e5)), simplify = FALSE) |>
     setNames(paste0("V", 1:1000)) |> qDT() |>
     tfm(g = sample.int(1e4, 1e5, replace = TRUE))

system.time(d |> gby(g) |> fsum())
system.time(d[, lapply(.SD, sum, na.rm = TRUE), by = g])


system.time(d |> gby(g) |> fmean(V1))
system.time(d[, lapply(.SD, weighted.mean, w = V1, na.rm = TRUE), by = g])
