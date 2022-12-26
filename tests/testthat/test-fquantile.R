context("fquantile")

probs1 <- c(0, 0.25, 0.5, 0.75, 1)
probs2 <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

for(x in mtcars) {
  for(o in list(NULL, radixorder(x))) {
    for(Qprobs in list(probs1, probs2)) {
      for(t in 5:9) {
        expect_equal(fquantile(x, Qprobs, type = t, o = o),
                      quantile(x, Qprobs, type = t))
        for(j in 1:3) {
          expect_equal(fquantile(x, Qprobs, type = t),
                       fquantile(x, Qprobs, type = t, w = rep(j + rnorm(1, sd = 0.05), 32), o = o))
        }
      }
    }
  }
}

# for(x in na_insert(airquality, 0.05)) {
#   for(o in list(NULL, radixorder(x))) {
#     for(Qprobs in list(probs1, probs2)) {
#       for(t in 5:9) {
#         expect_equal(fquantile(x, Qprobs, type = t, o = o),
#                       quantile(x, Qprobs, type = t, names = FALSE, na.rm = TRUE))
#         for(j in 1:10) {
#           w = rep(j + rnorm(1, sd = 0.05), 153)
#           expect_equal(fquantile(x, Qprobs, type = t),
#                        fquantile(x, Qprobs, type = t, w = w, o = o))
#         }
#       }
#     }
#   }
# }
