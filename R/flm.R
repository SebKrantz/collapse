
# formatcoef <- function(r, X, y) {
#   if(!is.matrix(r)) dim(r) <- c(length(r), 1L)
#   `dimnames<-`(r, list(dimnames(X)[[2L]], if(is.matrix(y)) dimnames(y)[[2L]] else NULL))
# }

# formatcoef <- function(r, y, X, drop) {
#   if(is.matrix(r)) return(`dimnames<-`(r, list(dimnames(X)[[2L]], if(is.matrix(y)) dimnames(y)[[2L]] else NULL)))
#   if(drop) return(name) # .....
#
#
#     list(dim = c(dim(X)[2L], 1L), dimnames = list(dimnames(X)[[2L]], NULL))
# }

flm <- function(...) if(is.atomic(..1)) flm.default(...) else flm.formula(...)

flm.default <- function(y, X, w = NULL, add.icpt = FALSE, #  sparse = FALSE,
                return.raw = FALSE, # only.coef
                method = c("lm", "solve", "qr", "arma", "chol", "eigen"),
                eigen.method = 3L, ...) {
  if(add.icpt) X <- cbind(`(Intercept)` = 1, X)
  n <- dim(X)[1L]
  if(n != NROW(y)) stop("NROW(y) must match nrow(X)")
  # if(sparse) X <- as(X, "dgCMatrix") # what about y ??
  if(length(w)) {
    if(length(w) != n) stop("w must be numeric and length(w) == nrow(X)")
    wts <- sqrt(w)
    if(return.raw) return(switch(method[1L],
                  lm = {
                    z <- .lm.fit(X * wts, y * wts, ...)
                    z$residuals <- z$residuals / wts # This is correct !!!
                    z
                  },
                  solve = (function(xw) solve(crossprod(xw), crossprod(xw, y * wts), ...))(X * wts),
                  qr = qr.coef(qr(X * wts, ...), y * wts),
                  arma = getenvFUN("RcppArmadillo_fastLmPure")(X * wts, y * wts), # .Call("_RcppArmadillo_fastLm_impl", X * wts, y * wts, PACKAGE = "RcppArmadillo"),
                  chol = (function(xw) chol2inv(chol(crossprod(xw), ...)) %*% crossprod(xw, y * wts))(X * wts),
                  eigen = {
                   z <- getenvFUN("RcppEigen_fastLmPure")(X * wts, y * wts, eigen.method) # .Call("RcppEigen_fastLm_Impl", X * wts, y * wts, eigen.method, PACKAGE = "RcppEigen")
                   z$residuals <- z$residuals / wts # This is correct !!!
                   z$fitted.values <- y - z$residuals
                   z
                  }, stop("Unknown method!")))

    ar <- if(is.matrix(y)) list(dim = c(dim(X)[2L], dim(y)[2L]), dimnames = list(dimnames(X)[[2L]], dimnames(y)[[2L]])) else
      list(dim = c(dim(X)[2L], 1L), dimnames = list(dimnames(X)[[2L]], NULL))

    return(`attributes<-`(switch(method[1L],
                  lm = .lm.fit(X * wts, y * wts, ...)[[2L]],
                  solve = (function(xw) solve(crossprod(xw), crossprod(xw, y * wts), ...))(X * wts),
                  qr = qr.coef(qr(`dimnames<-`(X, NULL) * wts, ...), y * wts),
                  arma = getenvFUN("RcppArmadillo_fastLmPure")(X * wts, y * wts)[[1L]], # .Call("_RcppArmadillo_fastLm_impl", X * wts, y * wts, PACKAGE = "RcppArmadillo"),
                  chol = (function(xw) chol2inv(chol(crossprod(xw), ...)) %*% crossprod(xw, y * wts))(X * wts),
                  eigen = getenvFUN("RcppEigen_fastLmPure")(X * wts, y * wts, eigen.method)[[1L]], # .Call("RcppEigen_fastLm_Impl", X * wts, y * wts, eigen.method, PACKAGE = "RcppEigen")
                  stop("Unknown method!")), ar))

  }
  if(return.raw) return(switch(method[1L],
                        lm = .lm.fit(X, y, ...),
                        solve = solve(crossprod(X), crossprod(X, y), ...),
                        qr = qr.coef(qr(X, ...), y),
                        arma = getenvFUN("RcppArmadillo_fastLmPure")(X, y),
                        chol = chol2inv(chol(crossprod(X), ...)) %*% crossprod(X, y),
                        eigen = getenvFUN("RcppEigen_fastLmPure")(X, y, eigen.method),
                        stop("Unknown method!")))

  ar <- if(is.matrix(y)) list(dim = c(dim(X)[2L], dim(y)[2L]), dimnames = list(dimnames(X)[[2L]], dimnames(y)[[2L]])) else
    list(dim = c(dim(X)[2L], 1L), dimnames = list(dimnames(X)[[2L]], NULL))

  `attributes<-`(switch(method[1L],
         lm = .lm.fit(X, y, ...)[[2L]],
         solve = solve(crossprod(X), crossprod(X, y), ...),
         qr = qr.coef(qr(`dimnames<-`(X, NULL), ...), y),
         arma = getenvFUN("RcppArmadillo_fastLmPure")(X, y)[[1L]],
         chol = chol2inv(chol(crossprod(X), ...)) %*% crossprod(X, y),
         eigen = getenvFUN("RcppEigen_fastLmPure")(X, y, eigen.method)[[1L]],
         stop("Unknown method!")), ar)

  # if(!return.raw) return(switch(method[1L], solve = formatcoef(res$coefficients, X, y), res$coefficients))
  # res
}

flm.formula <- function(formula, data = NULL, weights = NULL, add.icpt = TRUE, ...) {
  w <- substitute(weights)
  tms <- attributes(terms.formula(formula, data = data))
  pe <- tms[[".Environment"]]
  mf <- eval(tms$variables, data, pe)
  y <- mf[[1L]]
  X <- mf[-1L]
  if(length(w)) w <- eval(w, data, pe)
  names(X) <- tms$term.labels
  if(add.icpt) X <- c(list(`(Intercept)` = alloc(1, NROW(y))), X) # y could be matrix
  flm.default(y, do.call(cbind, X), w, FALSE, ...)
}

# Slower than using chol2inv (discarded)
# lmchol2 <- function(X, y) {
#   ch <- chol(crossprod(X))
#   backsolve(ch, forwardsolve(ch, crossprod(X, y), upper = TRUE, trans = TRUE))
# }
