context("recode, replace")

gmtc <- fgroup_by(mtcars, cyl)

test_that("replace_NA and replace_Inf work well", {

expect_equal(replace_NA(airquality, 0), `[<-`(airquality, is.na(airquality), value = 0))
expect_equal(replace_NA(airquality, 0, cols = 1:2), `[<-`(airquality, is.na(airquality), value = 0))
expect_equal(replace_NA(airquality, 0, cols = is.numeric), `[<-`(airquality, is.na(airquality), value = 0))
expect_equal(replace_NA(flag(EuStockMarkets), 0), `[<-`(flag(EuStockMarkets), is.na(flag(EuStockMarkets)), value = 0))

expect_equal(replace_Inf(dapply(mtcars, log)), `[<-`(dapply(mtcars, log), sapply(dapply(mtcars, log), is.infinite), value = NA))
expect_equal(replace_Inf(log(EuStockMarkets)), `[<-`(log(EuStockMarkets), is.infinite(log(EuStockMarkets)), value = NA))

expect_equal(replace_Inf(dapply(mtcars, log), replace.nan = TRUE), `[<-`(dapply(mtcars, log), sapply(dapply(mtcars, log), is.infinite), value = NA))
expect_equal(replace_Inf(log(EuStockMarkets), replace.nan = TRUE), `[<-`(log(EuStockMarkets), is.infinite(log(EuStockMarkets)), value = NA))

})


test_that("replace_outliers works well.", {

  expect_equal(replace_outliers(mtcars, 2), replace(mtcars, fscale(mtcars) > 2, NA))
  expect_equal(replace_outliers(gmtc, 2, single.limit = "overall_SDs"), replace(gmtc, dapply(mtcars, fscale) > 2, NA))
  expect_equal(replace_outliers(mtcars, 2, single.limit = "min"), replace(mtcars, mtcars < 2, NA))
  expect_equal(replace_outliers(mtcars, 2, single.limit = "max"), replace(mtcars, mtcars > 2, NA))

  expect_equal(replace_outliers(EuStockMarkets, 2), replace(EuStockMarkets, fscale(EuStockMarkets) > 2, NA))
  expect_equal(replace_outliers(EuStockMarkets, 2, single.limit = "overall_SDs"), replace(EuStockMarkets, dapply(EuStockMarkets, fscale) > 2, NA))
  expect_equal(replace_outliers(EuStockMarkets, 2, single.limit = "min"), replace(EuStockMarkets, EuStockMarkets < 2, NA))
  expect_equal(replace_outliers(EuStockMarkets, 2, single.limit = "max"), replace(EuStockMarkets, EuStockMarkets > 2, NA))


})

set.seed(101)
lmiss <- na_insert(letters)
month.miss <- na_insert(month.name)
char_dat <- na_insert(char_vars(GGDC10S))

options(warn = -1)

test_that("recode_char works well", {

  expect_equal(recode_char(lmiss, a = "b"), replace(lmiss, lmiss == "a", "b"))
  expect_visible(recode_char(lmiss, a = "b", missing = "a"))  # continue here to write proper tests!!..
  expect_visible(recode_char(lmiss, a = "b", missing = "c"))
  expect_visible(recode_char(lmiss, a = "b", default = "n"))
  expect_visible(recode_char(lmiss, a = "b", default = "n", missing = "c"))

  expect_visible(recode_char(month.miss, ber = NA, regex = TRUE))
  expect_visible(recode_char(month.miss, ber = NA, missing = "c", regex = TRUE))
  expect_visible(recode_char(lmiss, ber = NA, default = "n", regex = TRUE))
  expect_visible(recode_char(lmiss, ber = NA, default = "n", missing = "c", regex = TRUE))

  expect_visible(recode_char(lmiss, a = "b", e = "f"))
  expect_visible(recode_char(lmiss, a = "b", e = "f", missing = "a"))
  expect_visible(recode_char(lmiss, a = "b", e = "f", missing = "c"))
  expect_visible(recode_char(lmiss, a = "b", e = "f", default = "n"))
  expect_visible(recode_char(lmiss, a = "b", e = "f", default = "n", missing = "c"))

  expect_visible(recode_char(month.miss, ber = NA, May = "a", regex = TRUE))
  expect_visible(recode_char(month.miss, ber = NA, May = "a", missing = "c", regex = TRUE))
  expect_visible(recode_char(lmiss, ber = NA, May = "a", default = "n", regex = TRUE))
  expect_visible(recode_char(lmiss, ber = NA, May = "a", default = "n", missing = "c", regex = TRUE))

  expect_visible(recode_char(char_dat, SGP = "SINGAPORE", VA = "Value Added"))
  expect_visible(recode_char(char_dat, SGP = "SINGAPORE", VA = "Value Added", missing = "c"))
  expect_visible(recode_char(char_dat, SGP = "SINGAPORE", VA = "Value Added", default = "n"))
  expect_visible(recode_char(char_dat, SGP = "SINGAPORE", VA = "Value Added", default = "n", missing = "c"))

  expect_visible(recode_char(char_dat, saharan = "SSA", regex = TRUE))
  expect_visible(recode_char(char_dat, saharan = "SSA", regex = TRUE, missing = "c"))
  expect_visible(recode_char(char_dat, saharan = "SSA", regex = TRUE, default = "n"))
  expect_visible(recode_char(char_dat, saharan = "SSA", regex = TRUE, default = "n", missing = "c"))

})

set.seed(101)
vmiss <- na_insert(mtcars$cyl)
mtcNA <- na_insert(mtcars)
test_that("recode_num works well", {

  expect_equal(recode_num(vmiss, `4` = 5), replace(vmiss, vmiss == 4, 5))
  expect_visible(recode_num(vmiss, `4` = 5, missing = 4))  # continue here to write proper tests!!..
  expect_visible(recode_num(vmiss, `4` = 5, missing = 7))
  expect_visible(recode_num(vmiss, `4` = 5, default = 8))
  expect_visible(recode_num(vmiss, `4` = 5, default = 8, missing = 7))

  expect_visible(recode_num(vmiss, `4` = 5, `6` = 10))
  expect_visible(recode_num(vmiss, `4` = 5, `6` = 10, missing = 6))
  expect_visible(recode_num(vmiss, `4` = 5, `6` = 10, missing = 7))
  expect_visible(recode_num(vmiss, `4` = 5, `6` = 10, default = 8))
  expect_visible(recode_num(vmiss, `4` = 5, `6` = 10, default = 8, missing = 7))

  expect_visible(recode_num(mtcNA, `4` = 5, `1` = 2))
  expect_visible(recode_num(mtcNA, `4` = 5, `1` = 2, missing = 6))
  expect_visible(recode_num(mtcNA, `4` = 5, `1` = 2, default = 8))
  expect_visible(recode_num(mtcNA, `4` = 5, `1` = 2, default = 8, missing = 7))

})

options(warn = 1)
