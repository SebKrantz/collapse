# rm(list = ls())
# Sys.setenv(R_TESTS = "")
library(testthat)
options(collapse_export_F = TRUE)
# library(collapse)

test_check("collapse")

