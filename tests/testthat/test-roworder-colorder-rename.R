context("roworder, colorder, rename")

test_that("roworder works as intended", {

  expect_identical(roworder(mtcars, cyl, -hp), mtcars[with(mtcars, order(cyl, -hp)), ])
  expect_identical(roworder(airquality, Month, -Ozone), setRownames(airquality[with(airquality, order(Month, -Ozone)), ]))
  expect_identical(fnrow(roworder(airquality, Month, -Ozone, na.last = NA)), 116L)  # Removes the missing values in Ozone

  ## Same in standard evaluation
  expect_identical(roworderv(airquality, c("Month", "Ozone"), decreasing = c(FALSE, TRUE)), roworder(airquality, Month, -Ozone))

  ## Custom reordering
  expect_identical(roworderv(mtcars, neworder = 3:4), rbind(mtcars[3:4, ], mtcars[-(3:4), ]))               # Bring rows 3 and 4 to the front
  expect_identical(roworderv(mtcars, neworder = 3:4, pos = "end"), rbind(mtcars[-(3:4), ], mtcars[3:4, ]))  # Bring them to the end
  expect_identical(roworderv(mtcars, neworder = mtcars$vs == 1), rbind(mtcars[mtcars$vs == 1, ], mtcars[mtcars$vs != 1, ]))    # Bring rows with vs == 1 to the top
  expect_identical(ss(roworderv(mtcars, neworder = c(8, 2), pos = "exchange"), c(2,8)), ss(mtcars, c(8,2)))

})
