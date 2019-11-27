# library(Rcpp)
# sourceCpp("R/C++/psmat.cpp", rebuild = TRUE) # Todo: What about factors and dates ?? Matrix return fastest ??
# sourceCpp("R/C++/qFqG.cpp", rebuild = TRUE)
# qF <- function(x, ordered = TRUE) {
#   if(is.factor(x)) return(x)
#   qFCpp(x, ordered)
# }

# General Note : if add ... can put extra arguments that are not used. if not, gives unuse argument error !!!

psmat <- function(x, ...) { # g, t = NULL, cols = NULL, transpose = FALSE, simplify = FALSE
  UseMethod("psmat", x)
}
psmat.default <- function(x, g, t = NULL, transpose = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.matrix(x)) stop("x is already a matrix")
  if(is.atomic(g) && length(g) == 1L) {
    if(transpose) matrix(x, ncol = round(g), dimnames =
    list(seq_len(length(x)/round(g)), paste0("GRP.",seq_len(g)))) else
    matrix(x, nrow = round(g), byrow = TRUE,
    dimnames = list(paste0("GRP.",seq_len(g)), seq_len(length(x)/round(g))))
  } else {
  if(!is.factor(g)) if(is.atomic(g)) g <- qF(g) else if(is.GRP(g))
                    g <- as.factor.GRP(g) else g <- as.factor.GRP(GRP(g)) # interaction(lapply(g, qF))
  if(is.null(t)) {
    message("No timevar provided: Assuming Balanced Panel")
    psmatCpp(x, g, NULL, transpose)
  } else {
    if(!is.factor(t)) if(is.atomic(t)) t <- qF(t) else if(is.GRP(t))
                      t <- as.factor.GRP(t) else t <- as.factor.GRP(GRP(t)) # interaction(lapply(t, qF))
    psmatCpp(x, g, t, transpose)
    }
  }
}
psmat.data.frame <- function(x, by, t = NULL, cols = NULL, transpose = FALSE, array = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.atomic(by) && length(by) == 1L) {
    n <- round(by)
    if(transpose) {
      dn <- list(seq_len(length(x)/n), paste0("GRP.",seq_len(by)))
      res <- lapply(x, matrix, ncol = n, dimnames = dn)
    } else {
      dn <- list(paste0("GRP.",seq_len(by)), seq_len(length(x)/round(by)))
      res <- lapply(x, matrix, nrow = n, byrow = TRUE, dimnames = dn)
    }
  } else {
    if(is.call(by)) {
      nam <- names(x)
      if(length(by) == 3L) {
        v <- anyNAerror(match(all.vars(by[[2L]]), nam), "Unknown by columns!")
        by <- anyNAerror(match(all.vars(by[[3L]]), nam), "Unknown by columns!")
      } else {
        by <- anyNAerror(match(all.vars(by), nam), "Unknown by columns!")
        v <- if(is.null(cols)) seq_along(x)[-by] else if(is.function(cols))
             setdiff(which(vapply(x, cols, TRUE)), by) else if(is.character(cols))
             anyNAerror(match(cols, nam), "Unknown column names!") else cols
      }
      class(x) <- NULL
      by <- if(length(by) == 1L) x[[by]] else GRP(x, by) #, return.groups = FALSE)
      if(is.call(t)) { # If time-variable supplied !!
        t <- anyNAerror(match(all.vars(t), nam), "Unknown time variable!")
        v <- setdiff(v, t)
        t <- if(length(t) == 1L) x[[t]] else GRP(x, t) #, return.groups = FALSE)
      }
      x <- x[v]
    } else if(!is.null(cols)) {
      class(x) <- NULL
      x <- if(is.function(cols)) x[vapply(x, cols, TRUE)] else x[cols]
    }
    if(!is.factor(by)) if(is.atomic(by)) by <- qF(by) else if(is.GRP(by))
      by <- as.factor.GRP(by) else by <- as.factor.GRP(GRP(by)) # interaction(lapply(by, qF))
      if(is.null(t)) {
        message("No timevar provided: Assuming Balanced Panel")
        res <- lapply(x, psmatCpp, by, NULL, transpose)
      } else {
        if(!is.factor(t)) if(is.atomic(t)) t <- qF(t) else if(is.GRP(t))
          t <- as.factor.GRP(t) else t <- as.factor.GRP(GRP(t)) # interaction(lapply(t, qF))
        res <- lapply(x, psmatCpp, by, t, transpose)
      }
  }
  if(array) {
    if(length(res) == 1L) return(res[[1L]]) else
    return(addAttributes(simplify2array(res), list(transpose = transpose, class = c("psmat","array","table"))))
  } else return(res)
}
psmat.pseries <- function(x, transpose = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  index <- attr(x, "index")
  if(is.matrix(x)) stop("x is already a matrix")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  psmatCpp(x, index[[1L]], index[[2L]], transpose)
}
psmat.pdata.frame <- function(x, transpose = FALSE, array = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  index <- attr(x, "index")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  res <- lapply(x, psmatCpp, index[[1L]], index[[2L]], transpose)
  if(array) {
    if(length(res) == 1L) return(res[[1L]]) else
    return(addAttributes(simplify2array(res), list(transpose = transpose, class = c("psmat","array","table"))))
  } else return(res)
}

plot.psmat <- function(x, legend = FALSE, colours = legend, labs = NULL, ...) {
  d <- dim(x)
  arl <- length(d) == 3L
  if(isFALSE(attr(x, "transpose"))) {
    x <- if(arl) aperm(x, c(2,1,3)) else t.default(x)
    d <- dim(x)
  }
  dn <- dimnames(x)
  colours <- if(colours) rainbow(d[2L]) else TRUE
  t <- as.numeric(dn[[1L]])
  if(!is.na(t[1L])) {
    mint <- min(t)
    maxt <- max(t)
  } else {
    mint <- 1L
    maxt <- length(t)
  }
  ns <- d[2L]
  if(arl) {
    vars <- if(is.null(labs)) dn[[3L]] else labs
    nv <- d[3L]
    if(nv == 2L) mfr <- c(1,2+legend) else if(nv+legend <= 4L) mfr <- c(2,2) else {
      sqnv <- sqrt(nv)
      fsqnv <- floor(sqnv)
      mfr <- if(sqnv == fsqnv) c(fsqnv+legend,fsqnv) else c(fsqnv+1L,fsqnv)
    }
    settings <- par(c("mfrow","mar","mgp"))
    par(mfrow = mfr, mar = c(2.5,2.5,2.1,1.5), mgp = c(2.5,1,0))
    for(i in seq_along(vars)) ts.plot(ts(x[, , i], mint, maxt), main = vars[i], col = colours, xlab = NULL, ...)
    if(legend) {
      plot(1, type="n", axes=FALSE) #, xlab="", ylab="")
      legend('topleft', dn[[2L]], col = colours, lty=1, cex= if(ns > 80L) .65 else 1, bty = "n",
             ncol = if(ns <= 10L) 1L else if(nv == 2L) floor(ns^.25) else floor(ns^.37))
    }
    par(settings)
  } else {
    ts.plot(ts(x, mint, maxt), col = colours, ...)
    if(legend) legend('topleft', dn[[2L]], col = colours, lty=1,
                      cex= if(ns > 80L) .65 else 1, bty = "n",
                      ncol = if(d[2L] <= 10L) 1L else floor(d[2L]^.37))
  }
}

# is.balanced.panel()
# is.unsorted.panel <- function(x, g, t) {
# }
#
# index = attr(dGGDC,"index")
# x = index[[1L]]
# y = index[[2L]]
# if (length(x) != length(y))
#   stop("The length of the two vectors differs\n")
# x <- x[drop = TRUE]
# y <- y[drop = TRUE]
# z <- table(x, y)
# if (any(as.vector(z) == 0)) {
#   balanced <- FALSE
# }
# else {
#   balanced <- TRUE
# }
# if (any(as.vector(z) > 1))
#   warning("duplicate couples (id-time)\n")
# return(balanced)
