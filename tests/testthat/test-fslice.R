context("fslice")
data("iris")

test_that("fslice works with integers and no grouping", {
  N <- c(1, 5, 17)
  for (n in N) {
    # first
    expect_equal(
      dplyr::slice_head(iris, n = n),
      fslice(iris, n = n)
    )
    expect_equal(
      dplyr::slice_head(iris, n = n),
      fslice(iris, n = n, how = "first")
    )
    # last
    expect_equal(
      setRownames(dplyr::slice_tail(iris, n = n)),
      fslice(iris, n = n, how = "last")
    )
    # min
    expect_equal(
      iris |> dplyr::slice_min(Petal.Length, n = n, with_ties = FALSE),
      fslice(iris, n = n, how = "min", order.by = "Petal.Length")
    )
    # max
    expect_equal(
      iris |> dplyr::slice_max(Petal.Length, n = n, with_ties = FALSE),
      fslice(iris, n = n, how = "max", order.by = "Petal.Length")
    )
  }
})


test_that("fslice works with proportions and no grouping", {
  N <- c(0.5, 0.75)
  for (n in N) {
    # first
    expect_equal(
      dplyr::slice_head(iris, prop = n),
      fslice(iris, n = n)
    )
    expect_equal(
      dplyr::slice_head(iris, prop = n),
      fslice(iris, n = n, how = "first")
    )
    # last
    expect_equal(
      setRownames(dplyr::slice_tail(iris, prop = n)),
      fslice(iris, n = n, how = "last")
    )
    # min
    expect_equal(
      iris |> dplyr::slice_min(Petal.Length, prop = n, with_ties = FALSE),
      fslice(iris, n = n, how = "min", order.by = "Petal.Length")
    )
    # max
    expect_equal(
      iris |> dplyr::slice_max(Petal.Length, prop = n, with_ties = FALSE),
      fslice(iris, n = n, how = "max", order.by = "Petal.Length")
    )
  }
})


test_that("fslice works with grouping", {
  N <- c(1, 5, 17)
  for (n in N) {
    # first
    expect_equal(
      iris |> dplyr::group_by(Species) |> dplyr::slice_head(n = n) |> qDF(),
      fslice(iris, "Species", n = n, how = "first")
    )
    # last
    expect_equal(
      iris |> dplyr::group_by(Species) |> dplyr::slice_tail(n = n) |> qDF(),
      fslice(iris, "Species", n = n, how = "last")
    )
    # min
    expect_equal(
      iris |> dplyr::group_by(Species) |> dplyr::slice_min(Petal.Length, n = n, with_ties = FALSE) |> qDF(),
      fslice(iris, "Species", n = n, how = "min", order.by = "Petal.Length")
    )
    # max
    expect_equal(
      iris |> dplyr::group_by(Species) |> dplyr::slice_max(Petal.Length, n = n, with_ties = FALSE) |> qDF(),
      fslice(iris, "Species", n = n, how = "max", order.by = "Petal.Length")
    )
  }
})

test_that("fslice works with ties", {
  N <- 1 # c(1, 5, 17)
  for (n in N) {
    # min
    expect_equal(
      iris |> dplyr::group_by(Species) |> dplyr::slice_min(Petal.Length, n = n, with_ties = TRUE) |> qDF(),
      fslice(iris, "Species", n = n, how = "min", order.by = "Petal.Length", with.ties = TRUE)
    )
    # max
    expect_equal(
      iris |> dplyr::group_by(Species) |> dplyr::slice_max(Petal.Length, n = n, with_ties = TRUE) |> qDF(),
      fslice(iris, "Species", n = n, how = "max", order.by = "Petal.Length", with.ties = TRUE)
    )
  }
})
