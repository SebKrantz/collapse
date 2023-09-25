context("fmatch")

test_that("fmatch works well", {
  expect_identical(wlddev$iso3c %iin% "DEU", which(wlddev$iso3c %in% "DEU"))
  expect_identical(fsubset(wlddev, iso3c %in% c("DEU", "ITA")), fsubset(wlddev, iso3c %iin% c("DEU", "ITA")))
  expect_identical(qF(1:10+0.1) %iin% 1.1, 1L) # qF(1:10+0.1) %in% 1.1 works
  # what about integers?
})

###########################
# Proper Systematic Testing
###########################

library(kit) # count()
fmatch_base <- function(x, table, nomatch = NA_integer_, count = FALSE) {
  if(is.list(x)) {
    x <- do.call(paste0, x)
    table <- do.call(paste0, table)
  }
  res <- match(x, table, nomatch)
  if(count) {
    attr(res, "N.nomatch") <- count(res, nomatch)
    attr(res, "N.groups") <- length(table)
    attr(res, "N.distinct") <- if(is.na(nomatch))
        fndistinct.default(res) else fndistinct.default(res) - anyv(res, nomatch)
    oldClass(res) <- "qG"
  }
  res
}

random_vector_pair <- function(df, replace = FALSE, max.cols = 1) {
  d <- dim(df)
  cols <- sample.int(d[2L], if(is.na(max.cols)) as.integer(1 + d[2L] * runif(1)) else max.cols, replace)
  rows_x <- sample.int(d[1L], as.integer(1 + d[1L] * runif(1)), replace)
  rows_table <- sample.int(d[1L], as.integer(1 + d[1L] * runif(1)), replace)
  list(df[rows_x, cols], df[rows_table, cols])
}

match_identcal <- function(df, replace = FALSE, max.cols = 1, nomatch = NA_integer_, count = FALSE) {
  data <- random_vector_pair(df, replace, max.cols)
  x <- data[[1]]
  table <- data[[2]]
  id <- identical(fmatch(x, table, nomatch, count, overid = 2L),
                  fmatch_base(x, table, nomatch, count))
  if(id) TRUE else data
}

wldna <- na_insert(wlddev)

test_that("fmatch works well with atomic vectors", {
  for (r in c(FALSE, TRUE)) { # r = replace
    expect_true(all(replicate(100, match_identcal(wlddev, r))))
    expect_true(all(replicate(100, match_identcal(wlddev, r, nomatch = 0L))))
    expect_true(all(replicate(100, match_identcal(wlddev, r, count = TRUE))))
    expect_true(all(replicate(100, match_identcal(wlddev, r, nomatch = 0L, count = TRUE))))
    expect_true(all(replicate(100, match_identcal(wldna, r))))
    expect_true(all(replicate(100, match_identcal(wldna, r, nomatch = 0L))))
    expect_true(all(replicate(100, match_identcal(wldna, r, count = TRUE))))
    expect_true(all(replicate(100, match_identcal(wldna, r, nomatch = 0L, count = TRUE))))
  }
})

test_that("fmatch works well with data frames / lists", {
  for (r in c(FALSE, TRUE)) { # r = replace
    expect_true(all(replicate(20, match_identcal(wlddev, r, max.cols = NA))))
    expect_true(all(replicate(20, match_identcal(wlddev, r, max.cols = NA, nomatch = 0L))))
    expect_true(all(replicate(20, match_identcal(wlddev, r, max.cols = NA, count = TRUE))))
    expect_true(all(replicate(20, match_identcal(wlddev, r, max.cols = NA, nomatch = 0L, count = TRUE))))
    expect_true(all(replicate(20, match_identcal(wldna, r, max.cols = NA))))
    expect_true(all(replicate(20, match_identcal(wldna, r, max.cols = NA, nomatch = 0L))))
    expect_true(all(replicate(20, match_identcal(wldna, r, max.cols = NA, count = TRUE))))
    expect_true(all(replicate(20, match_identcal(wldna, r, max.cols = NA, nomatch = 0L, count = TRUE))))
  }
})


wld <- wlddev |> slt(iso3c, year = PCGDP) |> roworderv()
wld <- na_insert(wld)
x <- ss(wld, sample.int(10000, replace = TRUE))
table <- ss(wld, sample.int(1000, replace = TRUE))

expect_identical(fmatch(x$year, table$year), match(x$year, table$year))
expect_identical(fmatch(x, table), fmatch_base(x, table))


########################
# AI Generated Tests
########################

test_that("fmatch returns expected results", {

  # Test with vector input
  x <- c("a", "b", "c")
  table <- c("a", "b", "d")
  expect_equal(fmatch(x, table), fmatch_base(x, table))

  # Test with list input
  tab <- wlddev[sample.int(10000, 1000), ]
  expect_equal(fmatch(wlddev, tab, overid = 2L), fmatch_base(wlddev, tab))

  # Test with nomatch argument
  expect_equal(fmatch(x, table, nomatch = 0), fmatch_base(x, table, nomatch = 0))

  # Test with count argument
  expect_equal(fmatch(x, table, count = TRUE),
               fmatch_base(x, table, count = TRUE))

})

test_that("fmatch handles NA matching correctly", {

  x <- c("a", NA, "c")
  table <- c("a", "b")

  expect_equal(fmatch(x, table), fmatch_base(x, table))
  expect_equal(fmatch(x, table, nomatch = 0),
              fmatch_base(x, table, nomatch = 0))

})

test_that("fmatch returns correct index positions", {
  x <- c("a", "b", "c", "d")
  expect_equal(fmatch("a", x), 1L)
  expect_equal(fmatch("d", x), 4L)
  expect_equal(fmatch(c("a", "c"), x), c(1L, 3L))
  expect_equal(fmatch("e", x), NA_integer_)
})

test_that("fmatch works with nomatch argument", {
  x <- c("a", "b", "c", "d")
  expect_equal(fmatch("a", x, nomatch = 0L), 1L)
  expect_equal(fmatch("e", x, nomatch = 0L), 0L)

})

test_that("fmatch works with incomparables", {
  x <- c("a", NA, "c", "d")
  expect_equal(fmatch("a", x), 1L)
  expect_equal(fmatch(NA, x), 2L)
  expect_equal(fmatch("c", x), 3L)

})

test_that("fmatch works with duplicates", {
  x <- c("a", "b", "c", "c", "d")
  expect_equal(fmatch("c", x), 3L)
})

test_that("fmatch works with integer data", {
  x <- c(1L, 2L, 3L, 4L)
  expect_equal(fmatch(1L, x), 1L)
  expect_equal(fmatch(4L, x), 4L)
  expect_equal(fmatch(c(1L, 3L), x), c(1L, 3L))
  expect_equal(fmatch(5L, x), NA_integer_)

})

test_that("fmatch works with double data", {
  x <- c(1.1, 2.2, 3.3, 4.4)
  expect_equal(fmatch(1.1, x), 1L)
  expect_equal(fmatch(4.4, x), 4L)
  expect_equal(fmatch(c(1.1, 3.3), x), c(1L, 3L))
  expect_equal(fmatch(5.5, x), NA_integer_)
})

test_that("fmatch works with factor data", {
  x <- factor(c("a", "b", "c", "d"))
  expect_equal(fmatch("a", x), 1L)
  expect_equal(fmatch("d", x), 4L)
  expect_equal(fmatch(c("a", "c"), x), c(1L, 3L))
  expect_equal(fmatch("e", x), NA_integer_)

})

test_that("fmatch works with logical data", {
  x <- c(TRUE, FALSE, TRUE, FALSE)
  expect_equal(fmatch(TRUE, x), 1L)
  expect_equal(fmatch(FALSE, x), 2L)
})

