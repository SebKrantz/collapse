context("subset and transform data")

v <- na_insert(mtcars$mpg)
m <- na_insert(as.matrix(mtcars))

test_that("fsubset works like base::subset for vectors and matrices", {
  expect_identical(fsubset(v, 1:3), v[1:3])
  expect_identical(fsubset(v, 4:8), v[4:8])
  expect_error(fsubset(v, -(1:3))) # This does not work !!
  expect_identical(fsubset(v, v > 16), v[v > 16 & !is.na(v)])
  expect_identical(fsubset(m, 1:3), m[1:3, ])
  expect_identical(fsubset(m, v > 16), m[v > 16, ])
  expect_identical(fsubset(m, -(4:8)), m[-(4:8), ])

  expect_identical(fsubset(m, -(4:8), 1:5), m[-(4:8), 1:5])
  expect_identical(fsubset(m, v > 16 & !is.na(v), mpg:vs), subset(m, v > 16 & !is.na(v), mpg:vs))
  expect_identical(fsubset(m, v > 16 & !is.na(v), mpg, cyl:vs), subset(m, v > 16 & !is.na(v), c(mpg, cyl:vs)))
  expect_identical(fsubset(m, v > 16 & !is.na(v), -mpg), subset(m, v > 16 & !is.na(v), -mpg))
  expect_identical(fsubset(m, v > 16 & !is.na(v), -(mpg:vs)), subset(m, v > 16 & !is.na(v), -(mpg:vs)))
})

test_that("fsubset works like base::subset for data frames", {
  expect_identical(setRownames(fsubset(airquality, Ozone > 42)),
                   setRownames(subset(airquality, Ozone > 42)))
  expect_identical(setRownames(fsubset(airquality, Temp > 80, Ozone, Temp)),
                   setRownames(subset(airquality, Temp > 80, select = c(Ozone, Temp))))
  expect_identical(setRownames(fsubset(airquality, Day == 1, -Temp)),
                   setRownames(subset(airquality, Day == 1, select = -Temp)))
  expect_identical(setRownames(fsubset(airquality, Day == 1, -(Day:Temp))),
                   setRownames(subset(airquality, Day == 1, -(Day:Temp))))
  expect_identical(setRownames(fsubset(airquality, Day == 1, Ozone:Wind)),
                   setRownames(subset(airquality, Day == 1, Ozone:Wind)))
  expect_identical(setRownames(fsubset(airquality, Day == 1 & !is.na(Ozone), Ozone:Wind, Month)),
                   setRownames(subset(airquality, Day == 1 & !is.na(Ozone), c(Ozone:Wind, Month))))
})


test_that("ss works like an improved version of [", {
  expect_identical(ss(airquality, 1:100, 1:3), airquality[1:100, 1:3])
  expect_identical(setRownames(ss(airquality, -(1:100), 1:3)), setRownames(airquality[-(1:100), 1:3]))
  expect_identical(ss(airquality, 1:100, -(1:3)), airquality[1:100, -(1:3)])
  expect_identical(setRownames(ss(airquality, -(1:100), -(1:3))), setRownames(airquality[-(1:100), -(1:3)]))
  nam <- names(airquality)[2:5]
  v <- sample.int(nrow(airquality), 100)
  expect_identical(setRownames(ss(airquality, v, nam)), setRownames(airquality[v, nam, drop = FALSE]))
  expect_identical(setRownames(ss(airquality, -v, nam)), setRownames(airquality[-v, nam, drop = FALSE]))
  vl <- sample(c(TRUE, FALSE), nrow(airquality), replace = TRUE)
  cl <- sample(c(TRUE, FALSE), ncol(airquality), replace = TRUE)
  expect_identical(setRownames(ss(airquality, vl, nam)), setRownames(airquality[vl, nam, drop = FALSE]))
  expect_identical(setRownames(ss(airquality, vl, cl)), setRownames(airquality[vl, cl, drop = FALSE]))
  vl <- na_insert(vl)
  cl[4L] <- NA
  expect_identical(setRownames(ss(airquality, vl, nam)), setRownames(airquality[vl & !is.na(vl), nam, drop = FALSE]))
  expect_identical(setRownames(ss(airquality, vl, cl)), setRownames(airquality[vl & !is.na(vl), cl & !is.na(cl), drop = FALSE]))
})


test_that("ftransform works like base::transform", {

  expect_identical(ftransform(airquality, Ozone = -Ozone), transform(airquality, Ozone = -Ozone))
  expect_identical(ftransform(airquality, new = Ozone / Wind * 100), transform(airquality, new = Ozone / Wind * 100))
  expect_identical(ftransform(airquality, new = -Ozone, Temp = (Temp-32)/1.8), transform(airquality, new = -Ozone, Temp = (Temp-32)/1.8))
  expect_identical(ftransform(airquality, new = -Ozone, new2 = 1, Temp = NULL), transform(airquality, new = -Ozone, new2 = 1, Temp = NULL))
  expect_identical(ftransform(airquality, Ozone = NULL, Temp = NULL), transform(airquality, Ozone = NULL, Temp = NULL))

})


# Still do wrong input...
