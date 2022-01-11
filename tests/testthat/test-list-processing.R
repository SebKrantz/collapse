context("list-processing")

l <- lm(mpg ~cyl +  vs + am, mtcars)
# str(l, give.attr = FALSE)

is.regular <- function(x) is.atomic(x) || is.list(x)

test_that("atomic_elem and list_elem work well", {
  expect_equal(atomic_elem(l), unclass(l)[sapply(l, is.atomic)])
  expect_equal(list_elem(l), unclass(l)[sapply(l, is.list)])
  expect_equal(atomic_elem(l, keep.class = TRUE), `oldClass<-`(unclass(l)[sapply(l, is.atomic)], oldClass(l)))
  expect_equal(list_elem(l, keep.class = TRUE), `oldClass<-`(unclass(l)[sapply(l, is.list)], oldClass(l)))

  for(i in 1:6) expect_equal(atomic_elem(l, keep.class = TRUE, return = i), get_vars(l, is.atomic, return = i))

  for(i in 1:6) expect_equal(list_elem(l, keep.class = TRUE, return = i), get_vars(l, is.list, return = i))

  expect_identical(`atomic_elem<-`(l, atomic_elem(l)), l)
  expect_identical(`list_elem<-`(l, list_elem(l)), l)
  expect_error(`atomic_elem<-`(l, list_elem(l)))
  expect_error(`list_elem<-`(l, atomic_elem(l)))

})

test_that("ldepth works well", {
  expect_identical(ldepth(list(mtcars), DF.as.list = FALSE), 1L)
  expect_identical(ldepth(list(mtcars), DF.as.list = TRUE), 2L)
  expect_identical(ldepth(list(mtcars, l), DF.as.list = FALSE), 3L)
  expect_identical(ldepth(list(mtcars, l), DF.as.list = TRUE), 3L)
  expect_identical(ldepth(list(list(list(mtcars)), l), DF.as.list = FALSE), 3L)
  expect_identical(ldepth(list(list(list(mtcars)), l), DF.as.list = TRUE), 4L)
})

test_that("rapply2d works well", {
  l2 <- list(qM(mtcars), list(qM(mtcars), as.matrix(mtcars)))
  expect_equal(rapply2d(l2, fmean), rapply(l2, fmean, how = "list"))
  expect_equal(rapply2d(l[-length(l)], is.regular), rapply(l[-length(l)], is.regular, how = "list"))
})

test_that("get_elem works well", { # Could still add more tests..

  expect_true(is.matrix(get_elem(list(list(list(l))), is.matrix)))
  expect_false(is.matrix(get_elem(list(list(list(l))), is.matrix, keep.tree = TRUE)))

  l2 <- list(list(2,list("a",1)),list(1,list("b",2)))
  expect_identical(get_elem(l2, is.character), list("a", "b"))
  expect_identical(get_elem(l2, is.character, keep.tree = TRUE), list(list(list("a")),list(list("b"))))

  expect_identical(get_elem(l, "residuals"), resid(l))
  expect_identical(get_elem(l, "fit", regex = TRUE), fitted(l))
  expect_equal(get_elem(l, "tol"), 1e-7)
  expect_identical(get_elem(mtcars, 1), mtcars[[1]])
  expect_identical(get_elem(mtcars, 1, DF.as.list = TRUE), as.list(ss(mtcars, 1)))

})

test_that("reg_elem and irreg_elem work well", {
  expect_true(is_unlistable(reg_elem(l)))
  expect_false(is_unlistable(irreg_elem(l)))
  expect_true(is_unlistable(reg_elem(list(l), keep.tree = FALSE)))
  expect_true(is_unlistable(reg_elem(list(l), keep.tree = TRUE)))
  expect_false(is_unlistable(irreg_elem(list(l), keep.tree = FALSE)))
  expect_false(is_unlistable(irreg_elem(list(l), keep.tree = TRUE)))

})

test_that("has_elem works well", {
  expect_true(has_elem(l, is.matrix))
  expect_true(has_elem(l, is.data.frame))
  expect_false(has_elem(l, is.data.frame, DF.as.list = TRUE))
  expect_true(has_elem(l, is_categorical))
  expect_false(has_elem(l, is_date))
  expect_false(has_elem(l, is_qG))

  expect_false(has_elem(l, "am", recursive = FALSE))
  expect_false(has_elem(l, "pivot", recursive = FALSE))
  expect_true(has_elem(l, "pivot"))
  expect_true(has_elem(l, "am", DF.as.list = TRUE))
  expect_false(has_elem(l, "am"))
  expect_true(has_elem(l, "tol"))
  expect_false(has_elem(l, "mod"))
  expect_true(has_elem(l, "mod", regex = TRUE))
  expect_true(has_elem(l, "vot", regex = TRUE))
  expect_false(has_elem(l, "piv", regex = TRUE, recursive = FALSE))
})


