
psmat <- function(x, ...) UseMethod("psmat") # , x

psmat.default <- function(x, g, t = NULL, transpose = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.matrix(x)) stop("x is already a matrix")
  if(is.atomic(g) && length(g) == 1L) {
    if(transpose) matrix(x, ncol = round(g), dimnames =
    list(seq_len(length(x)/round(g)), paste0("GRP.",seq_len(g)))) else
    matrix(x, nrow = round(g), byrow = TRUE,
    dimnames = list(paste0("GRP.",seq_len(g)), seq_len(length(x)/round(g))))
  } else {
  if(!is.nmfactor(g)) if(is.atomic(g)) g <- qF(g, na.exclude = FALSE) else if(is_GRP(g))
                    g <- as_factor_GRP(g) else g <- as_factor_GRP(GRP.default(g, return.order = FALSE, call = FALSE))
  if(is.null(t)) {
    # message("No timevar provided: Assuming Balanced Panel")
    return(.Call(Cpp_psmat,x, g, NULL, transpose))
  } else {
    if(!is.nmfactor(t)) if(is.atomic(t)) t <- qF(t, sort = TRUE, na.exclude = FALSE) else if(is_GRP(t))
                      t <- as_factor_GRP(t) else t <- as_factor_GRP(GRP.default(t, sort = TRUE, return.order = FALSE, call = FALSE))
    return(.Call(Cpp_psmat,x, g, t, transpose))
    }
  }
}

psmat.data.frame <- function(x, by, t = NULL, cols = NULL, transpose = FALSE, array = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  oldClass(x) <- NULL # Setting globally !
  if(is.atomic(by) && length(by) == 1L) {
    nr <- .Call(C_fnrow, x)
    n <- round(by)
    if(length(cols)) x <- x[cols2int(cols, x, names(x), FALSE)]
    if(transpose) {
      dn <- list(seq_len(nr/n), paste0("GRP.",seq_len(by)))
      res <- lapply(x, matrix, ncol = n, dimnames = dn)
    } else {
      dn <- list(paste0("GRP.",seq_len(by)), seq_len(nr/n))
      res <- lapply(x, matrix, nrow = n, byrow = TRUE, dimnames = dn)
    }
  } else {
    if(is.call(by)) {
      nam <- names(x)
      if(length(by) == 3L) {
        v <- ckmatch(all.vars(by[[2L]]), nam)
        by <- ckmatch(all.vars(by[[3L]]), nam)
      } else {
        by <- ckmatch(all.vars(by), nam)
        v <- if(is.null(cols)) seq_along(x)[-by] else fsetdiff(cols2int(cols, x, nam), by)
      }
      by <- if(length(by) == 1L) x[[by]] else GRP.default(x, by, return.order = FALSE, call = FALSE)
      if(is.call(t)) { # If time-variable supplied !
        tv <- ckmatch(all.vars(t), nam, "Unknown time variable:")
        v <- fsetdiff(v, tv)
        t <- eval(if(length(tv) == 1L) t[[2L]] else attr(terms.formula(t), "variables"), x, attr(t, ".Environment")) # if(length(t) == 1L) x[[t]] else GRP.default(x, t, sort = TRUE, call = FALSE)
      }
      x <- x[v]
    } else if(length(cols)) x <- x[cols2int(cols, x, names(x), FALSE)]

    if(!is.nmfactor(by)) if(is.atomic(by)) by <- qF(by, na.exclude = FALSE) else if(is_GRP(by))
                         by <- as_factor_GRP(by) else by <- as_factor_GRP(GRP.default(by, return.order = FALSE, call = FALSE))
      if(is.null(t)) {
        # message("No timevar provided: Assuming Balanced Panel")
        res <- lapply(x, psmatCpp, by, NULL, transpose)
      } else {
        if(!is.nmfactor(t)) if(is.atomic(t)) t <- qF(t, sort = TRUE, na.exclude = FALSE) else if(is_GRP(t))
                  t <- as_factor_GRP(t) else t <- as_factor_GRP(GRP.default(t, sort = TRUE, return.order = FALSE, call = FALSE))
        res <- lapply(x, psmatCpp, by, t, transpose)
      }
  }
  if(array) {
    if(length(res) == 1L) return(res[[1L]]) else
    return(addAttributes(fsimplify2array(res), list(transpose = transpose, class = c("psmat","array"))))
  } else return(res)
}

psmat.pseries <- function(x, transpose = FALSE, drop.index.levels = "none", ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  index <- droplevels_index(uncl2pix(x, interact = TRUE), drop.index.levels)
  if(is.matrix(x)) stop("x is already a matrix")
  .Call(Cpp_psmat, x, index[[1L]], index[[2L]], transpose)
}

psmat.pdata.frame <- function(x, cols = NULL, transpose = FALSE, array = TRUE, drop.index.levels = "none", ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  index <- droplevels_index(uncl2pix(x, interact = TRUE), drop.index.levels)
  oldClass(x) <- NULL
  res <- lapply(if(is.null(cols)) x else x[cols2int(cols, x, names(x), FALSE)], psmatCpp, index[[1L]], index[[2L]], transpose)
  if(array) {
    if(length(res) == 1L) return(res[[1L]]) else
    return(addAttributes(fsimplify2array(res), list(transpose = transpose, class = c("psmat","array"))))
  } else return(res)
}

plot.psmat <- function(x, legend = FALSE,
                       colours = legend,
                       labs = NULL,
                       grid = FALSE, ...) {
  d <- dim(x)
  arl <- length(d) == 3L
  if(isFALSE(attr(x, "transpose"))) {
    x <- if(arl) aperm(x, c(2L, 1L, 3L)) else t.default(x)
    d <- dim(x)
  }
  dn <- dimnames(x)
  colours <- if(isTRUE(colours)) rainbow(d[2L]) else if(isFALSE(colours)) TRUE else colours
  t <- as.numeric(dn[[1L]])
  if(!is.na(t[1L])) {
    mint <- bmin(t)
    maxt <- bmax(t)
  } else {
    mint <- 1L
    maxt <- length(t)
  }
  ns <- d[2L]
  dots <- list(...)
  if(arl) {
    vars <- if(is.null(labs)) dn[[3L]] else labs
    nv <- d[3L]
    if(nv == 2L) mfr <- c(1L, 2L + legend) else if(nv + legend <= 4L) mfr <- c(2L, 2L) else {
      sqnv <- sqrt(nv)
      fsqnv <- floor(sqnv)
      mfr <- if(sqnv == fsqnv) c(fsqnv+legend,fsqnv) else c(fsqnv + 1L, fsqnv)
    }
    oldpar <- par(mfrow = mfr, mar = c(2.5, 2.5, 2.1, 1.5), mgp = c(2.5, 1, 0))
    on.exit(par(oldpar))
    for(i in seq_along(vars)) {
      ts.plot(ts(x[, , i], mint, maxt), main = vars[i], col = colours, xlab = NULL, ...)
      if(grid) grid()
    }
    if(legend) {
      plot(1:10, type = "n", axes = FALSE, xlab = NA, ylab = NA)
      legend(x = 0, y = if(nv == 2L) 10.5 else 10.75,  # 'topleft',
             dn[[2L]], col = colours, lty = if(any(names(dots) == "lty")) dots[["lty"]] else 1L,
             cex= if(ns > 80L) 1-sqrt(ns)/sqrt(1150) else 1, bty = "n", xpd = TRUE, # y.intersp = 0.5, x.intersp = 0.5,
             ncol = if(ns <= 10L) 1L else if(nv == 2L) floor(ns^.32) else floor(ns^.39)) # .37
    }
  } else {
    ts.plot(ts(x, mint, maxt), col = colours, ...)
    if(grid) grid()
    if(legend) legend('topleft', dn[[2L]], col = colours,
                      lty = if(any(names(dots) == "lty")) dots[["lty"]] else 1L,
                      cex= if(ns > 80L) 1-sqrt(ns)/sqrt(1150) else 1, bty = "n", xpd = TRUE,
                      # y.intersp = 0.5, x.intersp = 0.5,
                      ncol = if(d[2L] <= 10L) 1L else floor(d[2L]^.39)) #.37
  }
}

# print.psmat <- print.qsu # nah, too expensive

print.psmat <- function(x, digits = 3, ...) {
  print.default(`attr<-`(unclass(x), "transpose", NULL), digits = digits, ...)
}

`[.psmat` <- function(x, i, j, ..., drop = TRUE) {
  ret <- NextMethod()
  if(length(dim(ret)) > 1L) {
    attr(ret, "transpose") <- attr(x, "transpose")
    oldClass(ret) <- oldClass(x)
  }
  ret
}

aperm.psmat <- function(a, perm = NULL, resize = TRUE, keep.class = TRUE, ...) {
  r <- aperm.default(a, perm, resize = resize)
  if(keep.class) {
    attr(r, "transpose") <- attr(a, "transpose")
    oldClass(r) <- oldClass(a)
  }
  r
}

