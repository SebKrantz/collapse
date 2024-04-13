context("qtab")



set.seed(101)
wldNA <- na_insert(wlddev)

qtable <- function(...) {
  r <- qtab(...)
  oldClass(r) <- "table"
  attr(r, "sorted") <- NULL
  attr(r, "weighted") <- NULL
  r
}

ones <- alloc(1L, fnrow(wlddev))

attach(wlddev)

expect_equal(table(region, income), qtable(region, income))
expect_equal(table(income, region), qtable(income, region))
expect_equal(table(region, income, OECD), qtable(region, income, OECD))
expect_equal(table(decade, region, income, OECD), qtable(decade, region, income, OECD))
expect_equal(table(decade, country), qtable(decade, country))
expect_equal(table(iso3c, country), qtable(iso3c, country))
expect_equal(table(iso3c, decade), qtable(iso3c, decade))
expect_equal(table(iso3c, OECD), qtable(iso3c, OECD))

expect_equal(table(region, income), qtable(region, income, w = ones))
expect_equal(table(income, region), qtable(income, region, w = ones))
expect_equal(table(region, income, OECD), qtable(region, income, OECD, w = ones))
expect_equal(table(decade, region, income, OECD), qtable(decade, region, income, OECD, w = ones))
expect_equal(table(decade, country), qtable(decade, country, w = ones))
expect_equal(table(iso3c, country), qtable(iso3c, country, w = ones))
expect_equal(table(iso3c, decade), qtable(iso3c, decade, w = ones))
expect_equal(table(iso3c, OECD), qtable(iso3c, OECD, w = ones))

expect_equal(qtable(region, income, w = ones), qtable(region, income, w = ones, wFUN = sum))
expect_equal(qtable(income, region, w = ones), qtable(income, region, w = ones, wFUN = sum))
expect_equal(qtable(region, income, OECD, w = ones), qtable(region, income, OECD, w = ones, wFUN = sum))
expect_equal(qtable(decade, region, income, OECD, w = ones), qtable(decade, region, income, OECD, w = ones, wFUN = sum))
expect_equal(qtable(decade, country, w = ones), qtable(decade, country, w = ones, wFUN = sum))
expect_equal(qtable(iso3c, country, w = ones), qtable(iso3c, country, w = ones, wFUN = sum))
expect_equal(qtable(iso3c, decade, w = ones), qtable(iso3c, decade, w = ones, wFUN = sum))
expect_equal(qtable(iso3c, OECD, w = ones), qtable(iso3c, OECD, w = ones, wFUN = sum))

expect_equal(qtable(region, income, w = ones),  replace_NA(qtable(region, income, w = ones, wFUN = fsum)))
expect_equal(qtable(income, region, w = ones), replace_NA(qtable(income, region, w = ones, wFUN = fsum)))
expect_equal(qtable(region, income, OECD, w = ones), replace_NA(qtable(region, income, OECD, w = ones, wFUN = fsum)))
expect_equal(qtable(decade, region, income, OECD, w = ones), replace_NA(qtable(decade, region, income, OECD, w = ones, wFUN = fsum)))
expect_equal(qtable(decade, country, w = ones), replace_NA(qtable(decade, country, w = ones, wFUN = fsum)))
expect_equal(qtable(iso3c, country, w = ones), replace_NA(qtable(iso3c, country, w = ones, wFUN = fsum)))
expect_equal(qtable(iso3c, decade, w = ones), replace_NA(qtable(iso3c, decade, w = ones, wFUN = fsum)))
expect_equal(qtable(iso3c, OECD, w = ones), replace_NA(qtable(iso3c, OECD, w = ones, wFUN = fsum)))

detach(wlddev)

attach(wldNA)

expect_equal(table(region, income), qtable(region, income))
expect_equal(table(income, region), qtable(income, region))
expect_equal(table(region, income, OECD), qtable(region, income, OECD))
expect_equal(table(decade, region, income, OECD), qtable(decade, region, income, OECD))
expect_equal(table(decade, country), qtable(decade, country))
expect_equal(table(iso3c, country), qtable(iso3c, country))
expect_equal(table(iso3c, decade), qtable(iso3c, decade))
expect_equal(table(iso3c, OECD), qtable(iso3c, OECD))

expect_equal(table(region, income), qtable(region, income, w = ones))
expect_equal(table(income, region), qtable(income, region, w = ones))
expect_equal(table(region, income, OECD), qtable(region, income, OECD, w = ones))
expect_equal(table(decade, region, income, OECD), qtable(decade, region, income, OECD, w = ones))
expect_equal(table(decade, country), qtable(decade, country, w = ones))
expect_equal(table(iso3c, country), qtable(iso3c, country, w = ones))
expect_equal(table(iso3c, decade), qtable(iso3c, decade, w = ones))
expect_equal(table(iso3c, OECD), qtable(iso3c, OECD, w = ones))

expect_equal(qtable(region, income, w = ones), qtable(region, income, w = ones, wFUN = sum))
expect_equal(qtable(income, region, w = ones), qtable(income, region, w = ones, wFUN = sum))
expect_equal(qtable(region, income, OECD, w = ones), qtable(region, income, OECD, w = ones, wFUN = sum))
expect_equal(qtable(decade, region, income, OECD, w = ones), qtable(decade, region, income, OECD, w = ones, wFUN = sum))
expect_equal(qtable(decade, country, w = ones), qtable(decade, country, w = ones, wFUN = sum))
expect_equal(qtable(iso3c, country, w = ones), qtable(iso3c, country, w = ones, wFUN = sum))
expect_equal(qtable(iso3c, decade, w = ones), qtable(iso3c, decade, w = ones, wFUN = sum))
expect_equal(qtable(iso3c, OECD, w = ones), qtable(iso3c, OECD, w = ones, wFUN = sum))

expect_equal(qtable(region, income, w = ones),  replace_NA(qtable(region, income, w = ones, wFUN = fsum)))
expect_equal(qtable(income, region, w = ones), replace_NA(qtable(income, region, w = ones, wFUN = fsum)))
expect_equal(qtable(region, income, OECD, w = ones), replace_NA(qtable(region, income, OECD, w = ones, wFUN = fsum)))
expect_equal(qtable(decade, region, income, OECD, w = ones), replace_NA(qtable(decade, region, income, OECD, w = ones, wFUN = fsum)))
expect_equal(qtable(decade, country, w = ones), replace_NA(qtable(decade, country, w = ones, wFUN = fsum)))
expect_equal(qtable(iso3c, country, w = ones), replace_NA(qtable(iso3c, country, w = ones, wFUN = fsum)))
expect_equal(qtable(iso3c, decade, w = ones), replace_NA(qtable(iso3c, decade, w = ones, wFUN = fsum)))
expect_equal(qtable(iso3c, OECD, w = ones), replace_NA(qtable(iso3c, OECD, w = ones, wFUN = fsum)))

expect_equal(table(region, income, useNA = "ifany"), qtable(region, income, na.exclude = FALSE))
expect_equal(table(income, region, useNA = "ifany"), qtable(income, region, na.exclude = FALSE))
expect_equal(table(region, income, OECD, useNA = "ifany"), qtable(region, income, OECD, na.exclude = FALSE))
expect_equal(table(decade, region, income, OECD, useNA = "ifany"), qtable(decade, region, income, OECD, na.exclude = FALSE))
expect_equal(table(decade, country, useNA = "ifany"), qtable(decade, country, na.exclude = FALSE))
expect_equal(table(iso3c, country, useNA = "ifany"), qtable(iso3c, country, na.exclude = FALSE))
expect_equal(table(iso3c, decade, useNA = "ifany"), qtable(iso3c, decade, na.exclude = FALSE))
expect_equal(table(iso3c, OECD, useNA = "ifany"), qtable(iso3c, OECD, na.exclude = FALSE))

expect_equal(table(region, income, useNA = "ifany"), qtable(region, income, w = ones, na.exclude = FALSE))
expect_equal(table(income, region, useNA = "ifany"), qtable(income, region, w = ones, na.exclude = FALSE))
expect_equal(table(region, income, OECD, useNA = "ifany"), qtable(region, income, OECD, w = ones, na.exclude = FALSE))
expect_equal(table(decade, region, income, OECD, useNA = "ifany"), qtable(decade, region, income, OECD, w = ones, na.exclude = FALSE))
expect_equal(table(decade, country, useNA = "ifany"), qtable(decade, country, w = ones, na.exclude = FALSE))
expect_equal(table(iso3c, country, useNA = "ifany"), qtable(iso3c, country, w = ones, na.exclude = FALSE))
expect_equal(table(iso3c, decade, useNA = "ifany"), qtable(iso3c, decade, w = ones, na.exclude = FALSE))
expect_equal(table(iso3c, OECD, useNA = "ifany"), qtable(iso3c, OECD, w = ones, na.exclude = FALSE))

expect_equal(table(region, income, useNA = "ifany"), qtable(region, income, w = ones, wFUN = sum, na.exclude = FALSE))
expect_equal(table(income, region, useNA = "ifany"), qtable(income, region, w = ones, wFUN = sum, na.exclude = FALSE))
expect_equal(table(region, income, OECD, useNA = "ifany"), qtable(region, income, OECD, w = ones, wFUN = sum, na.exclude = FALSE))
expect_equal(table(decade, region, income, OECD, useNA = "ifany"), qtable(decade, region, income, OECD, w = ones, wFUN = sum, na.exclude = FALSE))
expect_equal(table(decade, country, useNA = "ifany"), qtable(decade, country, w = ones, wFUN = sum, na.exclude = FALSE))
expect_equal(table(iso3c, country, useNA = "ifany"), qtable(iso3c, country, w = ones, wFUN = sum, na.exclude = FALSE))
expect_equal(table(iso3c, decade, useNA = "ifany"), qtable(iso3c, decade, w = ones, wFUN = sum, na.exclude = FALSE))
expect_equal(table(iso3c, OECD, useNA = "ifany"), qtable(iso3c, OECD, w = ones, wFUN = sum, na.exclude = FALSE))

expect_equal(table(region, income, useNA = "ifany"), replace_NA(qtable(region, income, w = ones, wFUN = fsum, na.exclude = FALSE)))
expect_equal(table(income, region, useNA = "ifany"), replace_NA(qtable(income, region, w = ones, wFUN = fsum, na.exclude = FALSE)))
expect_equal(table(region, income, OECD, useNA = "ifany"), replace_NA(qtable(region, income, OECD, w = ones, wFUN = fsum, na.exclude = FALSE)))
expect_equal(table(decade, region, income, OECD, useNA = "ifany"), replace_NA(qtable(decade, region, income, OECD, w = ones, wFUN = fsum, na.exclude = FALSE)))
expect_equal(table(decade, country, useNA = "ifany"), replace_NA(qtable(decade, country, w = ones, wFUN = fsum, na.exclude = FALSE)))
expect_equal(table(iso3c, country, useNA = "ifany"), replace_NA(qtable(iso3c, country, w = ones, wFUN = fsum, na.exclude = FALSE)))
expect_equal(table(iso3c, decade, useNA = "ifany"), replace_NA(qtable(iso3c, decade, w = ones, wFUN = fsum, na.exclude = FALSE)))
expect_equal(table(iso3c, OECD, useNA = "ifany"), replace_NA(qtable(iso3c, OECD, w = ones, wFUN = fsum, na.exclude = FALSE)))

detach(wldNA)
rm(qtable)
