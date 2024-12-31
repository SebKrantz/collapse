context("join")

df1 <- data.frame(
  id1 = c(1, 1, 2, 3),
  id2 = c("a", "b", "b", "c"),
  name = c("John", "Jane", "Bob", "Carl"),
  age = c(35, 28, 42, 50)
)
df2 <- data.frame(
  id1 = c(1, 2, 3, 3),
  id2 = c("a", "b", "c", "e"),
  salary = c(60000, 55000, 70000, 80000),
  dept = c("IT", "Marketing", "Sales", "IT")
)

opts <- set_collapse(verbose = 0)

for (sort in c(FALSE, TRUE)) {
  expect_identical(join(df1, df2, how = "inner", sort = sort), merge(df1, df2))
  expect_identical(join(df1, df2, how = "left", sort = sort), merge(df1, df2, all.x = TRUE))
  expect_identical(join(df1, df2, how = "right", sort = sort), merge(df1, df2, all.y = TRUE))
  expect_identical(join(df1, df2, how = "full", sort = sort), merge(df1, df2, all = TRUE))
}

expect_identical(names(join(df1, df2, on = "id2", how = "full", keep.col.order = FALSE, column = TRUE))[1:2], c("id2", ".join"))
expect_identical(names(join(df1, df2, on = "id2", how = "full", keep.col.order = FALSE, column = TRUE, multiple = TRUE))[1:2], c("id2", ".join"))
expect_identical(names(join(df1, df2, on = "id2", how = "right", keep.col.order = FALSE, column = TRUE))[1:2], c("id2", ".join"))
expect_identical(names(join(df1, df2, on = "id2", how = "right", keep.col.order = FALSE, column = TRUE, multiple = TRUE))[1:2], c("id2", ".join"))

# Different types of joins
# https://github.com/SebKrantz/collapse/issues/503
x1 = data.frame(
  id = c(1L, 1L, 2L, 3L, NA_integer_),
  t  = c(1L, 2L, 1L, 2L, NA_integer_),
  x  = 11:15
)
y1 = data.frame(
  id = c(1,2, 4),
  y  = c(11L, 15L, 16)
)

for(i in c("l","i","r","f","s","a")) {
  expect_identical(capture.output(join(df1, df2, how = i, verbose = 1))[-1], capture.output(join(df1, df2, how = i, verbose = 0)))
  expect_identical(capture.output(join(x1, y1, how = i, verbose = 1))[-1], capture.output(join(x1, y1, how = i, verbose = 0)))
}

df1 = na_insert(df1, 0.3)
df2 = na_insert(df2, 0.3)

for(i in c("l","i","r","f","s","a")) {
  expect_identical(capture.output(join(df1, df2, how = i, verbose = 1))[-1], capture.output(join(df1, df2, how = i, verbose = 0)))
}


sort_merge <- function(..., sort = FALSE) {
  res = merge(...)
  if(sort) return(roworder(res, id1, id2))
  res
}

expect_identical(join(df1, df2, how = "inner", sort = TRUE), sort_merge(df1, df2, sort = TRUE))
expect_identical(join(df1, df2, how = "left", sort = TRUE), sort_merge(df1, df2, all.x = TRUE, sort = TRUE))
expect_identical(join(df1, df2, how = "right", sort = TRUE), sort_merge(df1, df2, all.y = TRUE, sort = TRUE))

######################################
# Rigorous Testing Sort-Merge-Join
######################################

sort_join <- function(x, y, on, ...) {
  res = join(x, y, on, ...)
  roworderv(res, on)
}

random_df_pair <- function(df, replace = FALSE, max.cols = 1) {
  d <- dim(df)
  cols <- sample.int(d[2L], if(is.na(max.cols)) as.integer(1 + d[2L] * 0.75 * runif(1)) else max.cols)
  rows_x <- sample.int(d[1L], as.integer(1 + d[1L] * runif(1)), replace)
  rows_table <- sample.int(d[1L], as.integer(1 + d[1L] * runif(1)), replace)
  list(ss(df, rows_x, cols), ss(df, rows_table, cols), rows_x, rows_table, cols)
}

join_identical <- function(df, replace = FALSE, max.cols = 1, sort = TRUE, ...) {
  data <- random_df_pair(df, replace, max.cols)
  x <- data[[1]]
  y <- data[[2]]
  cols <- data[[5]]
  nam <- names(df)
  rem <- nam[-cols]
  if(length(rem) > 2L) {
    rem_x <- sample(rem, as.integer(length(rem)/2))
    rem_y <- setdiff(rem, rem_x)
    av(x) <- ss(df, data[[3]], rem_x)
    av(y) <- ss(df, data[[4]], rem_y)
  }
  if(sort) {
    id <- tryCatch(identical(join(x, y, on = nam[cols], sort = TRUE, ...),
                     sort_join(x, y, on = nam[cols], overid = 2L, ...)), error = function(e) FALSE)
  } else {
    id <- identical(join(x, y, on = nam[cols], sort = FALSE, overid = 2L, ...),
                    merge(x, y, by = nam[cols], all.x = TRUE, ...))
  }
  if(id) TRUE else list(x, y, nam[cols])
}

# (d <- join_identical(wlddev))

wldna <- na_insert(wlddev)
wldcc <- replace_NA(wlddev)

test_that("sort merge join works well with single vectors", {
  for (h in c("l","i","r","f","s","a")) {
    for (r in c(FALSE, TRUE)) { # r = replace
      expect_true(all(replicate(100, join_identical(wlddev, r, how = h))))
      expect_true(all(replicate(100, join_identical(wldna, r, how = h))))
      expect_true(all(replicate(100, join_identical(wldcc, r, how = h))))
    }
  }
})


#  (d <- join_identical(wlddev[1:8], FALSE, max.cols = 4))

wldna <- na_insert(wlddev)
wldcc <- replace_NA(wlddev)
NCRAN <- Sys.getenv("NCRAN") == "TRUE"
test_that("sort merge join works well with multiple vectors", {
  for (h in c("l", if(NCRAN) c("i","r","f","s","a") else NULL)) {
    for (r in c(FALSE, TRUE)) { # r = replace
      expect_true(all(replicate(100, join_identical(wlddev, r, max.cols = NA, how = h))))
      expect_true(all(replicate(100, join_identical(wldna, r, max.cols = NA, how = h))))
      expect_true(all(replicate(100, join_identical(wldcc, r, max.cols = NA, how = h))))
    }
  }
})


# Testing misc. issues: factors with integers and doubles
d1 = mtcars |> fcompute(v1 = mpg, g = qF(seq_len(32)+100))
d2 = mtcars |> fcompute(v2 = mpg, g = seq_len(32)+100L)

expect_true(all_identical(with(join(d1, d2, verbose = 0), list(v1, v2))))
expect_true(all_identical(with(join(d1, d2, verbose = 0, sort = TRUE), list(v1, v2))))

d2 = mtcars |> fcompute(v2 = mpg, g = seq_len(32)+100)

expect_true(all_identical(with(join(d1, d2, verbose = 0), list(v1, v2))))
expect_true(all_identical(with(join(d1, d2, verbose = 0, sort = TRUE), list(v1, v2))))


set_collapse(opts)
