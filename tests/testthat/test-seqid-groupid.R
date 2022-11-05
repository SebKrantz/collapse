context("seqid, groupid")

if(!is.null(attributes(identical(FALSE, TRUE)))) stop("OECD label issue")

# rm(list = ls())

x <- c(1:10, 1:10)

test_that("seqid performas as expected", {

expect_identical(unattrib(seqid(x)), rep(1:2, each = 10))
expect_identical(unattrib(seqid(x)), unattrib(seqid(x, na.skip = TRUE)))
expect_identical(unattrib(seqid(c(1, NA, 3), na.skip = TRUE)), as.integer(c(1, NA, 2)))
expect_identical(unattrib(seqid(c(1, NA, 2), na.skip = TRUE)), as.integer(c(1, NA, 1)))
expect_identical(unattrib(seqid(c(1, NA, 3), na.skip = TRUE, skip.seq = TRUE)), as.integer(c(1, NA, 1)))
expect_identical(unattrib(seqid(c(1, NA, 2), na.skip = TRUE, skip.seq = TRUE)), as.integer(c(1, NA, 2)))
expect_identical(unattrib(seqid(x)), unattrib(seqid(x, na.skip = TRUE)))

set.seed(101)
xNA <- na_insert(x, prop = 0.15)
expect_true(!anyNA(seqid(xNA)))
expect_identical(is.na(seqid(xNA, na.skip = TRUE)), is.na(xNA))
xNA2 <- xNA
xNA2[c(1,20)] <- NA_integer_
expect_true(!anyNA(seqid(xNA2)))
expect_identical(is.na(seqid(xNA2, na.skip = TRUE)), is.na(xNA2))

# Start at 0
expect_equal(seqid(x, start = 0)[1], 0L)
expect_equal(seqid(x, na.skip = TRUE, start = 0)[1], 0L)
expect_identical(unclass(seqid(x, start = 0)), unclass(seqid(x, na.skip = TRUE, start = 0)))

o <- order(rnorm(20))
xuo <- x[o]
xNAuo <- xNA[o]
xNA2uo <- xNA2[o]
o <- order(o)
expect_identical(x, xuo[o])
expect_identical(xNA, xNAuo[o])
expect_identical(xNA2, xNA2uo[o])

# seqid(xuo)
# seqid(xuo, na.skip = TRUE)
# seqid(xNAuo)
# seqid(xNAuo, na.skip = TRUE)
# seqid(xNA2uo)
# seqid(xNA2uo, na.skip = TRUE)

expect_identical(seqid(xuo, o)[o], unattrib(seqid(x)))
expect_identical(seqid(xuo, o, na.skip = TRUE)[o], unattrib(seqid(x, na.skip = TRUE)))
expect_identical(seqid(xNAuo, o)[o], unattrib(seqid(xNA)))
expect_identical(seqid(xNAuo, o, na.skip = TRUE)[o], unattrib(seqid(xNA, na.skip = TRUE)))
expect_identical(seqid(xNA2uo, o)[o], unattrib(seqid(xNA2)))
expect_identical(seqid(xNA2uo, o, na.skip = TRUE)[o],  unattrib(seqid(xNA2, na.skip = TRUE)))

# Check o
expect_identical(seqid(xuo, o, check.o = FALSE)[o], unattrib(seqid(x)))
expect_identical(seqid(xuo, o, na.skip = TRUE, check.o = FALSE)[o], unattrib(seqid(x, na.skip = TRUE)))
expect_identical(seqid(xNAuo, o, check.o = FALSE)[o], unattrib(seqid(xNA)))
expect_identical(seqid(xNAuo, o, na.skip = TRUE, check.o = FALSE)[o], unattrib(seqid(xNA, na.skip = TRUE)))
expect_identical(seqid(xNA2uo, o, check.o = FALSE)[o], unattrib(seqid(xNA2)))
expect_identical(seqid(xNA2uo, o, na.skip = TRUE, check.o = FALSE)[o],  unattrib(seqid(xNA2, na.skip = TRUE)))

# Start at 0
expect_identical(seqid(xuo, o, start = 0)[o], unattrib(seqid(x, start = 0)))
expect_identical(seqid(xuo, o, na.skip = TRUE, start = 0)[o], unattrib(seqid(x, na.skip = TRUE, start = 0)))
expect_identical(seqid(xNAuo, o, start = 0)[o], unattrib(seqid(xNA, start = 0)))
expect_identical(seqid(xNAuo, o, na.skip = TRUE, start = 0)[o], unattrib(seqid(xNA, na.skip = TRUE, start = 0)))
expect_identical(seqid(xNA2uo, o, start = 0)[o], unattrib(seqid(xNA2, start = 0)))
expect_identical(seqid(xNA2uo, o, na.skip = TRUE, start = 0)[o],  unattrib(seqid(xNA2, na.skip = TRUE, start = 0)))

# Check o, start at 0
expect_identical(seqid(xuo, o, check.o = FALSE, start = 0)[o], unattrib(seqid(x, start = 0)))
expect_identical(seqid(xuo, o, na.skip = TRUE, check.o = FALSE, start = 0)[o], unattrib(seqid(x, na.skip = TRUE, start = 0)))
expect_identical(seqid(xNAuo, o, check.o = FALSE, start = 0)[o], unattrib(seqid(xNA, start = 0)))
expect_identical(seqid(xNAuo, o, na.skip = TRUE, check.o = FALSE, start = 0)[o], unattrib(seqid(xNA, na.skip = TRUE, start = 0)))
expect_identical(seqid(xNA2uo, o, check.o = FALSE, start = 0)[o], unattrib(seqid(xNA2, start = 0)))
expect_identical(seqid(xNA2uo, o, na.skip = TRUE, check.o = FALSE, start = 0)[o],  unattrib(seqid(xNA2, na.skip = TRUE, start = 0)))

})

# Testing groupid -----------------------
x <- rep(5:6, each = 10)

test_that("groupid performas as expected", {

# groupid(x)
# groupid(x, na.skip = TRUE)
set.seed(101)
xNA <- na_insert(x, prop = 0.15)
# groupid(xNA)  # desirable behavior ??
# groupid(xNA, na.skip = TRUE) # -> Yes !!
xNA2 <- xNA
xNA2[c(1,20)] <- NA_integer_
# groupid(xNA2)
# groupid(xNA2, na.skip = TRUE)

# This was an issue !!
expect_identical(groupid(c(NA,NA,1.343,NA,NA)), groupid(c(NA,NA,1L,NA,NA)))

expect_true(allNA(replicate(500, groupid(NA, na.skip = TRUE)))) #335
expect_equal(unattrib(groupid(c(NA, NA), na.skip = TRUE)), c(NA_integer_, NA_integer_))
expect_equal(unattrib(groupid(c(NA, "a"), na.skip = TRUE)), c(NA, 1L))
expect_equal(unattrib(groupid(c(NA, NA, "a"), na.skip = TRUE)), c(NA, NA, 1L))

# Start at 0
# groupid(x, start = 0)
# groupid(x, na.skip = TRUE, start = 0)
# groupid(xNA, start = 0)
# groupid(xNA, na.skip = TRUE, start = 0)
# groupid(xNA2, start = 0)
# groupid(xNA2, na.skip = TRUE, start = 0)

o <- order(rnorm(20))
xuo <- x[o]
xNAuo <- xNA[o]
xNA2uo <- xNA2[o]
o <- order(o)
expect_identical(x, xuo[o])
expect_identical(xNA, xNAuo[o])
expect_identical(xNA2, xNA2uo[o])

# groupid(xuo)
# groupid(xuo, na.skip = TRUE)
# groupid(xNAuo)
# groupid(xNAuo, na.skip = TRUE)
# groupid(xNA2uo)
# groupid(xNA2uo, na.skip = TRUE)

expect_identical(groupid(xuo, o)[o], unattrib(groupid(x)))
expect_identical(groupid(xuo, o, na.skip = TRUE)[o], unattrib(groupid(x, na.skip = TRUE)))
expect_identical(groupid(xNAuo, o)[o], unattrib(groupid(xNA)))
expect_identical(groupid(xNAuo, o, na.skip = TRUE)[o], unattrib(groupid(xNA, na.skip = TRUE)))
expect_identical(groupid(xNA2uo, o)[o], unattrib(groupid(xNA2)))
expect_identical(groupid(xNA2uo, o, na.skip = TRUE)[o],  unattrib(groupid(xNA2, na.skip = TRUE)))

# Check o
expect_identical(groupid(xuo, o, check.o = FALSE)[o], unattrib(groupid(x)))
expect_identical(groupid(xuo, o, na.skip = TRUE, check.o = FALSE)[o], unattrib(groupid(x, na.skip = TRUE)))
expect_identical(groupid(xNAuo, o, check.o = FALSE)[o], unattrib(groupid(xNA)))
expect_identical(groupid(xNAuo, o, na.skip = TRUE, check.o = FALSE)[o], unattrib(groupid(xNA, na.skip = TRUE)))
expect_identical(groupid(xNA2uo, o, check.o = FALSE)[o], unattrib(groupid(xNA2)))
expect_identical(groupid(xNA2uo, o, na.skip = TRUE, check.o = FALSE)[o],  unattrib(groupid(xNA2, na.skip = TRUE)))

# Start at 0
expect_identical(groupid(xuo, o, start = 0)[o], unattrib(groupid(x, start = 0)))
expect_identical(groupid(xuo, o, na.skip = TRUE, start = 0)[o], unattrib(groupid(x, na.skip = TRUE, start = 0)))
expect_identical(groupid(xNAuo, o, start = 0)[o], unattrib(groupid(xNA, start = 0)))
expect_identical(groupid(xNAuo, o, na.skip = TRUE, start = 0)[o], unattrib(groupid(xNA, na.skip = TRUE, start = 0)))
expect_identical(groupid(xNA2uo, o, start = 0)[o], unattrib(groupid(xNA2, start = 0)))
expect_identical(groupid(xNA2uo, o, na.skip = TRUE, start = 0)[o],  unattrib(groupid(xNA2, na.skip = TRUE, start = 0)))

# Check o, start at 0
expect_identical(groupid(xuo, o, check.o = FALSE, start = 0)[o], unattrib(groupid(x, start = 0)))
expect_identical(groupid(xuo, o, na.skip = TRUE, check.o = FALSE, start = 0)[o], unattrib(groupid(x, na.skip = TRUE, start = 0)))
expect_identical(groupid(xNAuo, o, check.o = FALSE, start = 0)[o], unattrib(groupid(xNA, start = 0)))
expect_identical(groupid(xNAuo, o, na.skip = TRUE, check.o = FALSE, start = 0)[o], unattrib(groupid(xNA, na.skip = TRUE, start = 0)))
expect_identical(groupid(xNA2uo, o, check.o = FALSE, start = 0)[o], unattrib(groupid(xNA2, start = 0)))
expect_identical(groupid(xNA2uo, o, na.skip = TRUE, check.o = FALSE, start = 0)[o],  unattrib(groupid(xNA2, na.skip = TRUE, start = 0)))

})
