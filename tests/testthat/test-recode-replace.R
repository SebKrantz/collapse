context("recode, replace")


test_that("replace_NA and replace_Inf work well", {

expect_equal(replace_NA(airquality, 0), `[<-`(airquality, is.na(airquality), value = 0))
expect_equal(replace_NA(flag(EuStockMarkets), 0), `[<-`(flag(EuStockMarkets), is.na(flag(EuStockMarkets)), value = 0))

expect_equal(replace_Inf(dapply(mtcars, log)), `[<-`(dapply(mtcars, log), sapply(dapply(mtcars, log), is.infinite), value = NA))
expect_equal(replace_Inf(log(EuStockMarkets)), `[<-`(log(EuStockMarkets), is.infinite(log(EuStockMarkets)), value = NA))

expect_equal(replace_Inf(dapply(mtcars, log), replace.nan = TRUE), `[<-`(dapply(mtcars, log), sapply(dapply(mtcars, log), is.infinite), value = NA))
expect_equal(replace_Inf(log(EuStockMarkets), replace.nan = TRUE), `[<-`(log(EuStockMarkets), is.infinite(log(EuStockMarkets)), value = NA))

})
