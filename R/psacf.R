
# TODO: could use source code of C_acf and adjust for panel: https://github.com/SurajGupta/r-source/blob/a28e609e72ed7c47f6ddfbb86c85279a0750f0b7/src/library/stats/src/filter.c

psacf <- function(x, ...) UseMethod("psacf") # , x

psacf.default <- function(x, g, t = NULL, lag.max = NULL, type = c("correlation", "covariance","partial"), plot = TRUE, gscale = TRUE, ...) {
  if(!is.numeric(x)) stop("'x' must be a numeric vector")
  typei <- switch(type[1L], correlation = 1L, covariance = 2L, partial = 3L, stop("Unknown type!"))
  series <- l1orlst(as.character(substitute(x)))
  g <- G_guo(g)
  if(is.null(lag.max)) lag.max <- round(2*sqrt(length(x)/g[[1L]]))
  if(gscale) x <- fscaleCpp(x,g[[1L]],g[[2L]])
  acf <- if(typei == 2L)
    cov(x, .Call(Cpp_flaglead,x,0:lag.max,NA,g[[1L]],g[[2L]],G_t(t),FALSE), use = "pairwise.complete.obs") else
    c(1, cov(x, .Call(Cpp_flaglead,x,seq_len(lag.max),NA,g[[1L]],g[[2L]],G_t(t),FALSE), use = "pairwise.complete.obs")/fvar.default(x)) #  or complete obs ?
  d <- c(lag.max+1,1,1)
  if(typei == 3L) {
    acf <- .Call(C_pacf1, array(acf, d), lag.max)
    lag <- array(seq_len(d[1]), c(lag.max,1,1))
  } else {
    dim(acf) <- d
    lag <- array(0:lag.max, d)
  }
  acf.out <- `oldClass<-`(list(acf = acf, type = type[1L], n.used = length(x), lag = lag, series = series, snames = NULL), "acf")
  if(plot) {
    plot(acf.out, ylab = if(typei == 3L) "Panel Series Partial ACF" else "Panel Series ACF", ...)
    invisible(acf.out)
  } else {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    return(acf.out)
  }
}

psacf.data.frame <- function(x, by, t = NULL, cols = is.numeric, lag.max = NULL, type = c("correlation", "covariance","partial"), plot = TRUE, gscale = TRUE, ...) {
  typei <- switch(type[1L], correlation = 1L, covariance = 2L, partial = 3L, stop("Unknown type!"))
  series <- l1orlst(as.character(substitute(x)))
  oldClass(x) <- NULL
  if(is.call(by)) { # best way ?
    nam <- names(x)
    if(length(by) == 3L) {
      v <- ckmatch(all.vars(by[[2L]]), nam)
      by <- ckmatch(all.vars(by[[3L]]), nam)
    } else {
      by <- ckmatch(all.vars(by), nam)
      v <- if(is.null(cols)) seq_along(x)[-by] else fsetdiff(cols2int(cols, x, nam), by)
    }
    by <- if(length(by) == 1L) x[[by]] else x[by]
    if(is.call(t)) { # If time-variable supplied
      tv <- ckmatch(all.vars(t), nam, "Unknown time variable:")
      v <- fsetdiff(v, tv)
      t <- eval(if(length(tv) == 1L) t[[2L]] else attr(terms.formula(t), "variables"), x, attr(t, ".Environment")) # if(length(t) == 1L) x[[t]] else x[t]
    }
    x <- x[v]
  } else if(length(cols)) x <- x[cols2int(cols, x, names(x), FALSE)]
  lx <- length(x)
  nrx <- length(x[[1L]])
  snames <- names(x)
  attributes(x) <- NULL # already class is 0... Necessary ?
  getacf <- function(ng, g) {
    if(length(t)) t <- G_t(t)
    if(gscale) x <- fscalelCpp(x,ng,g)
    acf <- array(numeric(0), c(lag.max+1, lx, lx))
    fun <- if(typei == 2L) cov else function(x, y, ...) cov(x, y, ...)/fvar.default(x) # cor
      for(i in seq_len(lx)) {
        xim <- .Call(Cpp_flaglead,x[[i]],0:lag.max,NA,ng,g,t,FALSE)
        for(j in seq_len(lx)) acf[ , j, i] <- fun(x[[j]], xim, use = "pairwise.complete.obs") # correct !
      }
    acf
  }
  by <- G_guo(by)
  if(is.null(lag.max)) lag.max <- round(2*sqrt(nrx/by[[1L]]))
  acf <- getacf(by[[1L]], by[[2L]])
  lag <- matrix(1, lx, lx)
  lag[lower.tri(lag)] <- -1
  if(typei == 3L) {
    zvec <- double((1L+lag.max)*lx*lx)
    z <- .C(C_multi_yw, aperm(acf, 3:1), as.integer(nrx), as.integer(lag.max), as.integer(lx),
            coefs = zvec, pacf = zvec, var = zvec, aic = double(1L+lag.max), order = 0L, 1L)
    acf <- aperm(array(z$pacf, dim = c(lx, lx, lag.max + 1L)), 3:1)[-1L, , , drop = FALSE]
  }
  acf.out <- `oldClass<-`(list(acf = acf, type = type[1L], n.used = nrx,
                            lag = if(typei == 3L) outer(1L:lag.max, lag) else outer(0L:lag.max, lag),
                            series = series, snames = snames), "acf")
  if(plot) {
    plot(acf.out, ylab = if(typei == 3L) "Panel Series Partial ACF" else "Panel Series ACF",
         mar = if(lx > 2) c(3, 2.4, 2, 0.8) else par("mar"), ...)
    invisible(acf.out)
  } else {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    return(acf.out)
  }
}

psacf.pseries <- function(x, lag.max = NULL, type = c("correlation", "covariance","partial"), plot = TRUE, gscale = TRUE, ...) {
  if(!is.numeric(x)) stop("'x' must be a numeric pseries ")
  index <- uncl2pix(x)
  g <- index[[1L]]
  t <- index[[2L]]
  if(length(t) && !inherits(x, "indexed_series")) t <- plm_check_time(t)
  ng <- fnlevels(g)
  typei <- switch(type[1L], correlation = 1L, covariance = 2L, partial = 3L, stop("Unknown type!"))
  series <- l1orlst(as.character(substitute(x))) # faster ?
  if(is.null(lag.max)) lag.max <- round(2*sqrt(length(x)/ng))
  if(gscale) x <- fscaleCpp(x,ng,g)
  acf <- if(typei == 2L)
    cov(x, .Call(Cpp_flaglead,x,0:lag.max,NA,ng,g,t,FALSE), use = "pairwise.complete.obs") else
    c(1, cov(x, .Call(Cpp_flaglead,x,seq_len(lag.max),NA,ng,g,t,FALSE), use = "pairwise.complete.obs")/fvar.default(x)) # or complete obs ?
  d <- c(lag.max+1,1,1)
  if(typei == 3L) {
    acf <- .Call(C_pacf1, array(acf, d), lag.max)
    lag <- array(seq_len(d[1]), c(lag.max,1,1))
  } else {
    dim(acf) <- d
    lag <- array(0:lag.max, d)
  }
  acf.out <- `oldClass<-`(list(acf = acf, type = type[1L], n.used = length(x), lag = lag, series = series, snames = NULL), "acf")
  if (plot) {
    plot(acf.out, ylab = if(typei == 3L) "Panel Series Partial ACF" else "Panel Series ACF", ...)
    invisible(acf.out)
  } else {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    return(acf.out)
  }
}

psacf.pdata.frame <- function(x, cols = is.numeric, lag.max = NULL, type = c("correlation", "covariance","partial"), plot = TRUE, gscale = TRUE, ...) {
  typei <- switch(type[1L], correlation = 1L, covariance = 2L, partial = 3L, stop("Unknown type!"))
  series <- l1orlst(as.character(substitute(x))) # faster solution ?
  index <- uncl2pix(x)
  clx <- oldClass(x)
  oldClass(x) <- NULL
  nrx <- length(x[[1L]])
  if(length(cols)) x <- x[cols2int(cols, x, names(x), FALSE)]
  lx <- length(x)
  snames <- names(x)
  g <- index[[1L]]
  t <- index[[2L]]
  if(length(t) && !any(clx == "indexed_frame")) t <- plm_check_time(t)
  ng <- fnlevels(g)
  attributes(x) <- NULL # necessary after unclass above ?
    if(is.null(lag.max)) lag.max <- round(2*sqrt(nrx/ng))
    if(gscale) x <- fscalelCpp(x,ng,g)
    acf <- array(numeric(0), c(lag.max+1, lx, lx))
    fun <- if(typei == 2L) cov else function(x, y, ...) cov(x, y, ...)/fvar.default(x) # cor
    for(i in seq_len(lx)) {
      xim <- .Call(Cpp_flaglead,x[[i]],0:lag.max,NA,ng,g,t,FALSE)
      for(j in seq_len(lx)) acf[ , j, i] <- fun(x[[j]], xim, use = "pairwise.complete.obs") # correct !
    }
  lag <- matrix(1, lx, lx)
  lag[lower.tri(lag)] <- -1
  if(typei == 3L) {
    zvec <- double((1L+lag.max)*lx*lx)
    z <- .C(C_multi_yw, aperm(acf, 3:1), as.integer(nrx), as.integer(lag.max), as.integer(lx),
            coefs = zvec, pacf = zvec, var = zvec, aic = double(1L+lag.max), order = 0L, 1L)
    acf <- aperm(array(z$pacf, dim = c(lx, lx, lag.max + 1L)), 3:1)[-1L, , , drop = FALSE]
  }
  acf.out <- `oldClass<-`(list(acf = acf, type = type[1L], n.used = nrx,
                            lag = if(typei == 3L) outer(1L:lag.max, lag) else outer(0L:lag.max, lag),
                            series = series, snames = snames), "acf")
  if(plot) {
    plot(acf.out, ylab = if(typei == 3L) "Panel Series Partial ACF" else "Panel Series ACF",
         mar = if(lx > 2) c(3, 2.4, 2, 0.8) else par("mar"), ...)
    invisible(acf.out)
  } else {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    return(acf.out)
  }
}

pspacf <- function(x, ...) UseMethod("pspacf") # , x

pspacf.default <- function(x, g, t = NULL, lag.max = NULL, plot = TRUE, gscale = TRUE, ...) {
  if(plot)
  psacf.default(x, g, t, lag.max, "partial", plot, gscale, main = paste0("Series ",l1orlst(as.character(substitute(x)))), ...) else
  psacf.default(x, g, t, lag.max, "partial", plot, gscale, ...)
}

pspacf.pseries <- function(x, lag.max = NULL, plot = TRUE, gscale = TRUE, ...) {
  if(plot)
  psacf.pseries(x, lag.max, "partial", plot, gscale, main = paste0("Series ",l1orlst(as.character(substitute(x)))), ...) else
  psacf.pseries(x, lag.max, "partial", plot, gscale, ...)
}

pspacf.data.frame <- function(x, by, t = NULL, cols = is.numeric, lag.max = NULL, plot = TRUE, gscale = TRUE, ...) {
  psacf.data.frame(x, by, t, cols, lag.max, "partial", plot, gscale, ...)
}

pspacf.pdata.frame <- function(x, cols = is.numeric, lag.max = NULL, plot = TRUE, gscale = TRUE, ...) {
  psacf.pdata.frame(x, cols, lag.max, "partial", plot, gscale, ...)
}

psccf <- function(x, y, ...) UseMethod("psccf") # , x

psccf.default <- function(x, y, g, t = NULL, lag.max = NULL, type = c("correlation", "covariance"), plot = TRUE, gscale = TRUE, ...) {
  if(!is.numeric(x)) stop("'x' must be a numeric vector")
  if(!is.numeric(y)) stop("'y' must be a numeric vector")
  lx <- length(x)
  if(lx != length(y)) stop("length(x) must be equal to length(y)")
  typei <- switch(type[1L], correlation = 1L, covariance = 2L, partial = 3L, stop("Unknown type!"))
  snames <- paste(c(l1orlst(as.character(substitute(x))), l1orlst(as.character(substitute(x)))), collapse = " & ")
  getccf <- function(ng, g) {
    if(length(t)) t <- G_t(t)
    if(gscale) {
      x <- fscaleCpp(x,ng,g)
      y <- fscaleCpp(y,ng,g)
    }
    if(typei == 2L)
      drop(cov(x, .Call(Cpp_flaglead,y,-lag.max:lag.max,NA,ng,g,t,FALSE), use = "pairwise.complete.obs")) else
      drop(cov(x, .Call(Cpp_flaglead,y,-lag.max:lag.max,NA,ng,g,t,FALSE), use = "pairwise.complete.obs")/(fsd.default(x)*fsd.default(y))) # or complete obs ?
  }
  g <- G_guo(g)
  if(is.null(lag.max)) lag.max <- round(2*sqrt(lx/g[[1L]]))
  acf <- getccf(g[[1L]], g[[2L]])
  d <- c(2*lag.max+1,1,1)
  dim(acf) <- d
  acf.out <- `oldClass<-`(list(acf = acf, type = type[1L], n.used = lx,
                            lag = array(-lag.max:lag.max, d), series = snames, snames = snames), "acf")
  if (plot) {
    plot(acf.out, ylab = "Panel Series CCF", ...)
    invisible(acf.out)
  } else {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    return(acf.out)
  }
}

psccf.pseries <- function(x, y, lag.max = NULL, type = c("correlation", "covariance"), plot = TRUE, gscale = TRUE, ...) {
  if(!is.numeric(x)) stop("'x' must be a numeric pseries")
  if(!is.numeric(y) || !inherits(y, "pseries")) stop("'y' must be a numeric pseries")
  lx <- length(x)
  if(lx != length(y)) stop("length(x) must be equal to length(y)")
  if(!identical(findex(x), findex(y))) stop("index of x and y differs")
  index <- uncl2pix(x)
  g <- index[[1L]]
  t <- index[[2L]]
  if(length(t) && !inherits(x, "indexed_series")) t <- plm_check_time(t)
  ng <- fnlevels(g)
  typei <- switch(type[1L], correlation = 1L, covariance = 2L, partial = 3L, stop("Unknown type!"))
  snames <- paste(c(l1orlst(as.character(substitute(x))), l1orlst(as.character(substitute(x)))), collapse = " & ")
  if (gscale) {
    x <- fscaleCpp(x,ng,g)
    y <- fscaleCpp(y,ng,g)
  }
  if (is.null(lag.max)) lag.max <- round(2*sqrt(length(x)/ng))
  l_seq <- -lag.max:lag.max
  acf <- if(typei == 2L)
    drop(cov(x, .Call(Cpp_flaglead,y,l_seq,NA,ng,g,t,FALSE), use = "pairwise.complete.obs")) else
    drop(cov(x, .Call(Cpp_flaglead,y,l_seq,NA,ng,g,t,FALSE), use = "pairwise.complete.obs")/(fsd.default(x)*fsd.default(y))) # or complete obs ?
  d <- c(2*lag.max+1,1,1)
  dim(acf) <- d
  acf.out <- `oldClass<-`(list(acf = acf, type = type[1L], n.used = lx,
                            lag = array(l_seq, d), series = snames, snames = snames), "acf")
  if (plot) {
    plot(acf.out, ylab = "Panel Series CCF", ...)
    invisible(acf.out)
  } else {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    return(acf.out)
  }
}


# could do AR models also :
# psar.data.frame <- function (x, aic = TRUE, order.max = lag.max, na.action = na.fail,
#           demean = TRUE, series = NULL, var.method = 1L, ...)
# {
#   if (is.null(series))
#     series <- l1orlst(as.character(substitute(x)))
#   if (ists <- is.ts(x))
#     xtsp <- tsp(x)
#   x <- na.action(as.ts(x))
#   if (anyNA(x))
#     stop("NAs in 'x'")
#   if (ists)
#     xtsp <- tsp(x)
#   xfreq <- frequency(x)
#   x <- as.matrix(x)
#   nser <- ncol(x)
#   n.used <- nrow(x)
#   if (demean) {
#     x.mean <- colMeans(x)
#     x <- sweep(x, 2L, x.mean, check.margin = FALSE)
#   }
#   else x.mean <- rep(0, nser)
#   order.max <- if (is.null(order.max))
#     floor(10 * log10(n.used))
#   else floor(order.max)
#   if (order.max < 1L)
#     stop("'order.max' must be >= 1")
#   xacf <- acf(x, type = "cov", plot = FALSE, lag.max = order.max)$acf
#   z <- .C(stats:::C_"multi_yw",
#           aperm(xacf, 3:1),
#           as.integer(n.used),
#           as.integer(order.max),
#           as.integer(nser),
#           coefs = double((1L +order.max) * nser * nser),
#           pacf = double((1L + order.max) * nser * nser),
#           var = double((1L + order.max) * nser * nser),
#           aic = double(1L + order.max),
#           order = integer(1L),
#           as.integer(aic))
#   partialacf <- aperm(array(z$pacf, dim = c(nser, nser, order.max +
#                                               1L)), 3:1)[-1L, , , drop = FALSE]
#   var.pred <- aperm(array(z$var, dim = c(nser, nser, order.max +
#                                            1L)), 3:1)
#   xaic <- setNames(z$aic - bmin(z$aic), 0:order.max)
#   order <- z$order
#   resid <- x
#   if (order > 0) {
#     ar <- -aperm(array(z$coefs, dim = c(nser, nser, order.max +
#                                           1L)), 3:1)[2L:(order + 1L), , , drop = FALSE]
#     for (i in 1L:order) resid[-(1L:order), ] <- resid[-(1L:order),
#                                                       ] - x[(order - i + 1L):(n.used - i), ] %*% t(ar[i,
#                                                                                                       , ])
#     resid[1L:order, ] <- NA
#   }
#   else ar <- array(dim = c(0, nser, nser))
#   var.pred <- var.pred[order + 1L, , , drop = TRUE] * n.used/(n.used -
#                                                                 nser * (demean + order))
#   if (ists) {
#     attr(resid, "tsp") <- xtsp
#     attr(resid, "class") <- c("mts", "ts")
#   }
#   snames <- colnames(x)
#   colnames(resid) <- snames
#   dimnames(ar) <- list(seq_len(order), snames, snames)
#   dimnames(var.pred) <- list(snames, snames)
#   dimnames(partialacf) <- list(1L:order.max, snames, snames)
#   res <- list(order = order, ar = ar, var.pred = var.pred,
#               x.mean = x.mean, aic = xaic, n.used = n.used, order.max = order.max,
#               partialacf = partialacf, resid = resid, method = "Yule-Walker",
#               series = series, frequency = xfreq, call = match.call())
#   oldClass(res) <- "ar"
#   return(res)
# }
