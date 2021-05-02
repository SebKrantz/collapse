context("fsubset and ftransform")

# rm(list = ls())

v <- na_insert(mtcars$mpg)
m <- na_insert(as.matrix(mtcars))

test_that("fsubset works like base::subset for vectors and matrices", {
  expect_equal(fsubset(v, 1:3), v[1:3])
  expect_equal(fsubset(v, 4:8), v[4:8])
  expect_error(fsubset(v, -(1:3))) # This does not work !!
  expect_equal(fsubset(v, v > 16), v[v > 16 & !is.na(v)])
  expect_equal(fsubset(m, 1:3), m[1:3, ])
  expect_equal(fsubset(m, v > 16), m[v > 16, ])
  expect_equal(fsubset(m, -(4:8)), m[-(4:8), ])

  expect_equal(fsubset(m, -(4:8), 1:5), m[-(4:8), 1:5])
  expect_equal(fsubset(m, v > 16 & !is.na(v), mpg:vs), subset(m, v > 16 & !is.na(v), mpg:vs))
  expect_equal(fsubset(m, v > 16 & !is.na(v), mpg, cyl:vs), subset(m, v > 16 & !is.na(v), c(mpg, cyl:vs)))
  expect_equal(fsubset(m, v > 16 & !is.na(v), -mpg), subset(m, v > 16 & !is.na(v), -mpg))
  expect_equal(fsubset(m, v > 16 & !is.na(v), -(mpg:vs)), subset(m, v > 16 & !is.na(v), -(mpg:vs)))
})

test_that("fsubset works like base::subset for data frames", {
  expect_equal(unattrib(fsubset(airquality, Ozone > 42)),
                   unattrib(subset(airquality, Ozone > 42)))
  expect_equal(unattrib(fsubset(airquality, Temp > 80, Ozone, Temp)),
                   unattrib(subset(airquality, Temp > 80, select = c(Ozone, Temp))))
  expect_equal(unattrib(fsubset(airquality, Day == 1, -Temp)),
                   unattrib(subset(airquality, Day == 1, select = -Temp)))
  expect_equal(unattrib(fsubset(airquality, Day == 1, -(Day:Temp))),
                   unattrib(subset(airquality, Day == 1, -(Day:Temp))))
  expect_equal(unattrib(fsubset(airquality, Day == 1, Ozone:Wind)),
                   unattrib(subset(airquality, Day == 1, Ozone:Wind)))
  expect_equal(unattrib(fsubset(airquality, Day == 1 & !is.na(Ozone), Ozone:Wind, Month)),
                   unattrib(subset(airquality, Day == 1 & !is.na(Ozone), c(Ozone:Wind, Month))))
})

test_that("fsubset column renaming", {

  expect_equal(names(fsubset(airquality, Temp > 90, OZ = Ozone, Temp)), .c(OZ, Temp))

})

test_that("ss works like an improved version of [", { # replaced setRownames wit unattrib because of unexplained test failures on some systems
  expect_equal(ss(airquality, 1:100, 1:3), airquality[1:100, 1:3])
  expect_equal(unattrib(ss(airquality, -(1:100), 1:3)), unattrib(airquality[-(1:100), 1:3]))
  expect_equal(ss(airquality, 1:100, -(1:3)), airquality[1:100, -(1:3)])
  expect_equal(unattrib(ss(airquality, -(1:100), -(1:3))), unattrib(airquality[-(1:100), -(1:3)]))
  nam <- names(airquality)[2:5]
  v <- sample.int(fnrow(airquality), 100)
  expect_equal(unattrib(ss(airquality, v, nam)), unattrib(airquality[v, nam, drop = FALSE]))
  expect_equal(unattrib(ss(airquality, -v, nam)), unattrib(airquality[-v, nam, drop = FALSE]))
  vl <- sample(c(TRUE, FALSE), fnrow(airquality), replace = TRUE)
  cl <- sample(c(TRUE, FALSE), fncol(airquality), replace = TRUE)
  expect_equal(unattrib(ss(airquality, vl, nam)), unattrib(airquality[vl, nam, drop = FALSE]))
  expect_equal(unattrib(ss(airquality, vl, cl)), unattrib(airquality[vl, cl, drop = FALSE]))
  vl <- na_insert(vl)
  cl[4L] <- NA
  expect_equal(unattrib(ss(airquality, vl, nam)), unattrib(airquality[vl & !is.na(vl), nam, drop = FALSE]))
  expect_equal(unattrib(ss(airquality, vl, cl)), unattrib(airquality[vl & !is.na(vl), cl & !is.na(cl), drop = FALSE]))
})


test_that("ftransform works like base::transform", {

  expect_equal(ftransform(airquality, Ozone = -Ozone), transform(airquality, Ozone = -Ozone))
  expect_equal(ftransform(airquality, new = Ozone / Wind * 100), transform(airquality, new = Ozone / Wind * 100))
  expect_equal(ftransform(airquality, new = -Ozone, Temp = (Temp-32)/1.8), transform(airquality, new = -Ozone, Temp = (Temp-32)/1.8))
  expect_equal(ftransform(airquality, new = -Ozone, new2 = 1, Temp = NULL), transform(airquality, new = -Ozone, new2 = 1, Temp = NULL))
  expect_equal(ftransform(airquality, Ozone = NULL, Temp = NULL), transform(airquality, Ozone = NULL, Temp = NULL))

})


test_that("fcompute works well", {

  expect_equal(fcompute(airquality, new = -Ozone, new2 = 1, keep = 1:3), ftransform(airquality[1:3], new = -Ozone, new2 = 1))
  expect_equal(names(fcompute(airquality, new = -Ozone, new2 = 1, keep = 1:3)), .c(Ozone, Solar.R, Wind, new, new2))
  expect_equal(names(fcompute(airquality, new = -Ozone, new2 = 1)), .c(new, new2))
  expect_equal(names(fcompute(airquality, Ozone = -Ozone, new = 1, keep = 1:3)), .c(Ozone, Solar.R, Wind, new))

})



# Still do wrong input...

test_that("fsubset error for wrong input", {
  # expect_error(fsubset(mtcars, mpg))
  # expect_error(fsubset(mtcars, mpg:cyl))
  expect_error(fsubset(mtcars, "mpg"))
  expect_error(fsubset(mtcars, TRUE))
  expect_error(fsubset(mtcars, mpg > 15, cyl < 4))
  expect_error(fsubset(mtcars, mpg > 15, TRUE))
  expect_error(fsubset(mtcars, mpg > 15, 35))
  expect_error(fsubset(mtcars, mpg > 15, ~mpg))
})
