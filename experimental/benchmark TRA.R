library(microbenchmark)

# Solutions: https://stackoverflow.com/questions/8910268/general-lag-in-time-series-panel-data/60764643#60764643

# Default modulus
stats <- fmin(NWDI)
microbenchmark(TRA(NWDI, stats, "%%")) # 92 milliseconds
microbenchmark(TRA(NWDI, stats, "-%%")) # 92 milliseconds
gstats <- fmin(NWDI, g)
microbenchmark(TRA(NWDI, gstats, "%%", g)) # 99 milliseconds
microbenchmark(TRA(NWDI, gstats, "-%%", g)) # 100 milliseconds

# cmath fmod: / std::fmod: slow !! 5 seconds !!

# Adding cmath template function yourself:
