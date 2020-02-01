# library(Rcpp)
# sourceCpp("C++/small_helper.cpp")

# Export --------------------------------------
vlabels <- function(X, attrn = "label") {
  if(is.atomic(X)) {
    res <- attr(X, attrn)
    if(is.null(res)) NA else res
  } else {
    res <- lapply(X, attr, attrn)
    res[vapply(res, is.null, TRUE)] <- NA
    unlist(res)
  }
}
"vlabels<-" <- function(X, attrn = "label", value) {
  if(is.atomic(X)) {
    attr(X, attrn) <- value
  } else {
    for (i in seq_along(value)) attr(X[[i]], attrn) <- value[i]
  }
  X
}
namlab <- function(X, class = FALSE, attrn = "label") {
  pasteclass <- function(x) paste(class(x), collapse = " ")
  res <- if(class) list(names(X), vapply(X, pasteclass, character(1)), vlabels(X, attrn)) else list(names(X), vlabels(X, attrn))
  attributes(res) <- list(names = if(class) c("Variable","Class","Label") else c("Variable","Label"),
                          row.names = .set_row_names(length(X)),
                          class = "data.frame")
  return(res)
}
add_stub <- function(X, stub, pre = TRUE) {
  if(!is.character(stub)) return(X)
  if(is.array(X)) {
    if(length(dim(X)) > 2L) stop("Can't stub higher dimensional arrays!")
    dn <- dimnames(X)
    dimnames(X) <- list(dn[[1L]], if(pre) paste0(stub, dn[[2L]]) else paste0(dn[[2L]], stub))
  } else attr(X, "names") <- if(pre) paste0(stub, attr(X, "names")) else paste0(attr(X, "names"), stub)
  X
}
seq_row <- function(X) seq_len(nrow(X))
seq_col <- function(X) seq_len(ncol(X))
setRownames <- function(object = nm, nm = seq_row(object)) {
  rownames(object) <- nm
  object
}
setColnames <- function(object = nm, nm) {
  colnames(object) <- nm
  object
}
setDimnames <- function(object = dn, dn) {
  dimnames(object) <- dn
  object
}

# sumcc <- function(x, y)  sum(complete.cases(x,y))
# pwNobs <- function(x) qM(dapply(x, function(y) dapply(x, sumcc, y)))

pwNobs <- function(X) { # Faster !!
  dg <- fNobs(X)
  if(is.matrix(X)) {
    cn <- dimnames(X)[[2L]]
    X <- mctl(X)
  } else {
    cn <- names(X)
    class(X) <- NULL
  }
  n <- length(X)
  N.mat <- matrix(NA, n, n, dimnames = list(cn, cn))
  diag(N.mat) <- dg
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      N.mat[i, j] <- N.mat[j, i] <- sum(complete.cases(X[[i]], X[[j]]))
    }
  }
  N.mat
}

corr.p <- function(r, n) {
  if (n < 3L) return(1)
  df <- n - 2L
  t <- sqrt(df) * r/sqrt(1 - r^2)
  return(2 * min(pt(t, df), pt(t, df, lower.tail = FALSE))) # taken from corr.test
}

corr.pmat <- function(cm, nm) {
  n <- ncol(cm)
  p.mat <- matrix(NA, n, n, dimnames = dimnames(cm))
  # diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      p.mat[i, j] <- p.mat[j, i] <- corr.p(cm[i, j], nm[i, j])
    }
  }
  p.mat
}


# test against Hmisc::rcorr -> yup, same both corr, n and p !!
pwcor <- function(X, ..., N = FALSE, P = FALSE, array = TRUE) {
  r <- cor(X, ..., use = "pairwise.complete.obs")
  if(!N && !P) return(`class<-`(r, c("pwcor","matrix")))
  n <- pwNobs(X) # what if using ... to supply y ???
  if(N) {
    res <- if(P) list(r = r, N = n, P = corr.pmat(r, n)) else list(r = r, N = n)
  } else res <- list(r = r, P = corr.pmat(r, n))
  if(array) {
    res <- simplify2array(res)
    class(res) <- c("pwcor","array","table")
  } else class(res) <- "pwcor"
  return(res)
}

pwcov <- function(X, ..., N = FALSE, P = FALSE, array = TRUE) {
  r <- cov(X, ..., use = "pairwise.complete.obs")
  if(!N && !P) return(`class<-`(r, c("pwcov","matrix")))
  n <- pwNobs(X)
  if(N) {
    res <- if(P) list(cov = r, N = n, P = corr.pmat(cor(X, ..., use = "pairwise.complete.obs"), n)) else list(cov = r, N = n)
  } else res <- list(cov = r, P = corr.pmat(cor(X, ..., use = "pairwise.complete.obs"), n))
  if(array) {
    res <- simplify2array(res)
    class(res) <- c("pwcov","array","table")
  } else class(res) <- "pwcov"
  return(res)
}

print.pwcor <- function(x, digits = 2L, sig.level = 0.05, show = c("all","lower.tri","upper.tri"), spacing = 1L, ...) {
  formfun <- function(x, dg1 = FALSE) {
    xx <- format(round(x, digits)) # , digits = digits-1
    xx <- sub("(-?)0\\.", "\\1.", xx)
    if(dg1) diag(xx) <- paste0(c("  1",rep(" ",digits-1)), collapse = "") else {
     xna <- is.na(x)
     xx[xna] <- ""
     xpos <- x >= 1 & !xna
     xx[xpos] <- sub(paste0(c(".", rep("0",digits)), collapse = ""), "", xx[xpos]) # Problem: Deletes .00 also..
    }
    return(xx)
  }
  show <- switch(show[1L], all = 1L, lower.tri = 2L, upper.tri = 3L, stop("Unknown 'show' option"))
  se <- "Allowed spacing options are 0, 1 and 2!"
  if(is.array(x)) {
    sc <- TRUE
    d <- dim(x)
    ld <- length(d)
    if(ld > 2L) {
      dn <- dimnames(x)
      d3 <- dn[[3L]]
      if(all(d3 %in% c("r","N","P"))) {
        if(length(d3) == 3L) {
          sig <- matrix(" ", d[1L], d[2L])
          sig[x[,, 3L] <= sig.level] <- "*"
          res <- sprintf(switch(spacing+1L, "%s%s(%i)", "%s%s (%i)", " %s%s (%i)", stop(se)), formfun(x[,, 1L], TRUE), sig, x[,, 2L]) # paste0(formfun(x[,, 1L]),sig,"(",x[,, 2L],")")
        } else if(d3[2L] == "P") {
          sig <- matrix(" ", d[1L], d[2L])
          sig[x[,, 2L] <= sig.level] <- "*"
          res <- sprintf(switch(spacing+1L, "%s%s", " %s%s", " %s %s", stop(se)), formfun(x[,, 1L], TRUE), sig)
        } else res <- sprintf(switch(spacing+1L, "%s(%i)", "%s (%i)", " %s (%i)", stop(se)), formfun(x[,, 1L], TRUE), x[,, 2L])
      } else {
        sc <- FALSE
        res <- duplAttributes(switch(spacing+1L, formfun(x), sprintf(" %s",formfun(x)), sprintf("  %s",formfun(x)), stop(se)), x) # remove this before publishing !!!
      }
      if(sc) attributes(res) <- list(dim = d[1:2], dimnames = dn[1:2])
    } else res <- if(spacing == 0L) formfun(x, TRUE) else
                   duplAttributes(sprintf(switch(spacing," %s","  %s",stop(se)), formfun(x, TRUE)), x)
      if(sc && show != 1L) if(show == 2L) res[upper.tri(res)] <- "" else res[lower.tri(res)] <- ""
  } else {
    if(spacing == 0L) res <- lapply(x, formfun) else {
      ff <- function(i) duplAttributes(sprintf(switch(spacing," %s","  %s",stop(se)),formfun(i)), i)
      res <- lapply(x, ff)
    }
    if(show != 1L) res <- if(show == 2L) lapply(res, function(i){i[upper.tri(i)] <- ""; i}) else
                                         lapply(res, function(i){i[lower.tri(i)] <- ""; i})
  }
  print.default(unclass(res), quote = FALSE, right = TRUE, ...)
  invisible(x)
} #print.table(dapply(round(x, digits), function(j) sub("^(-?)0.", "\\1.", j)), right = TRUE, ...) # print.table(, right = TRUE)


print.pwcov <- function(x, digits = 2L, sig.level = 0.05, show = c("all","lower.tri","upper.tri"), spacing = 1L, ...) {
  formfun <- function(x, adj = FALSE) {
    xx <- format(round(x, digits), digits = 9, big.mark = ",", big.interval = 6)
    xx <- sub("(-?)0\\.", "\\1.", xx)
    if(adj) {
      xna <- is.na(x)
      xx[xna] <- ""
      xpos <- x >= 1 & !xna
      xx[xpos] <- sub(paste0(c(".", rep("0",digits)), collapse = ""), "", xx[xpos]) # Problem: Deletes .00 also..
    }
    return(xx)
  }
  show <- switch(show[1L], all = 1L, lower.tri = 2L, upper.tri = 3L, stop("Unknown 'show' option"))
  se <- "Allowed spacing options are 0, 1 and 2!"
  if(is.array(x)) {
    sc <- TRUE
    d <- dim(x)
    ld <- length(d)
    if(ld > 2L) {
      dn <- dimnames(x)
      d3 <- dn[[3L]]
      if(all(d3 %in% c("cov","N","P"))) {
        if(length(d3) == 3L) {
          sig <- matrix(" ", d[1L], d[2L])
          sig[x[,, 3L] <= sig.level] <- "*"
          res <- sprintf(switch(spacing+1L, "%s%s(%i)", "%s%s (%i)", " %s%s (%i)", stop(se)), formfun(x[,, 1L]), sig, x[,, 2L]) # paste0(formfun(x[,, 1L]),sig,"(",x[,, 2L],")")
        } else if(d3[2L] == "P") {
          sig <- matrix(" ", d[1L], d[2L])
          sig[x[,, 2L] <= sig.level] <- "*"
          res <- sprintf(switch(spacing+1L, "%s%s", " %s%s", " %s %s", stop(se)), formfun(x[,, 1L]), sig)
        } else res <- sprintf(switch(spacing+1L, "%s(%i)", "%s (%i)", " %s (%i)", stop(se)), formfun(x[,, 1L]), x[,, 2L])
      } else {
        sc <- FALSE
        res <- duplAttributes(switch(spacing+1L, formfun(x, TRUE), sprintf(" %s",formfun(x, TRUE)), sprintf("  %s",formfun(x, TRUE)), stop(se)), x) # remove this before publishing !!!
      }
      if(sc) attributes(res) <- list(dim = d[1:2], dimnames = dn[1:2])
    } else res <- if(spacing == 0L) formfun(x) else
      duplAttributes(sprintf(switch(spacing," %s","  %s",stop(se)), formfun(x)), x)
    if(sc && show != 1L) if(show == 2L) res[upper.tri(res)] <- "" else res[lower.tri(res)] <- ""
  } else {
    if(spacing == 0L) res <- lapply(x, formfun, TRUE) else {
      ff <- function(i) duplAttributes(sprintf(switch(spacing," %s","  %s",stop(se)),formfun(i, TRUE)), i)
      res <- lapply(x, ff)
    }
    if(show != 1L) res <- if(show == 2L) lapply(res, function(i){i[upper.tri(i)] <- ""; i}) else
      lapply(res, function(i){i[lower.tri(i)] <- ""; i})
  }
  print.default(unclass(res), quote = FALSE, right = TRUE, ...)
  invisible(x)
} #print.table(dapply(round(x, digits), function(j) sub("^(-?)0.", "\\1.", j)), right = TRUE, ...) # print.table(, right = TRUE)


# print.pwcov <- function(x, digits = 2, ...) print.default(formatC(round(x, digits), format = "g",
#                                                                   digits = 9, big.mark = ",", big.interval = 6), quote = FALSE, right = TRUE, ...)




getdf <- function(x) {
  if(is.matrix(x)) return(ncol(x))
  if(is.atomic(x)) if(is.factor(x)) return(fnlevels(x)-1L) else return(1L)
  sum(vapply(x, function(i) if(is.factor(i)) fnlevels(i)-1L else 1L, 1L))
}


fFtest <- function(y, exc, X = NULL, full.df = TRUE, ...) {
  if(!is.numeric(y)) stop("y needs to be a numeric vector")
  if(!is.null(X)) {
    Xn <- NCOL(X)
    data <- if(is.numeric(X) && is.numeric(exc)) na.omit(cbind(y, X, exc)) else
            na.omit(qDT(c(list(y = y), qDF(X), qDF(exc))))
    if(full.df && is.list(data) && any(fc <- vapply(data, is.factor, TRUE))) {
      cld <- class(data)
      class(data) <- NULL
      data[fc] <- lapply(data[fc], droplevels.factor)
      df <- vapply(data, function(i) if(is.factor(i)) fnlevels(i)-1L else 1L, 1L) # getdf(data)
      k <- sum(df) # 1 for intercept added with y
      p <- sum(df[-seq_len(Xn+1L)])
      y <- data[[1L]]
      oldClass(data) <- cld
    } else {
      k <- ncol(data) # 1 for intercept added with y
      p <- NCOL(exc)
      y <- data[, 1L]
    }
    kr <- k-p-1
    vy <- var(y)
    n <- nrow(data)
    r2f <- 1 - var(fHDwithin.default(y, data[, -1L], na.rm = FALSE, ...))/vy
    r2r <- 1 - var(fHDwithin.default(y, data[, 2:(Xn+1L)], na.rm = FALSE, ...))/vy # this way is data.tabel proof !!
    ndff <- k-1
    ddff <- n-k
    Fstatf <- r2f/ndff * ddff/(1-r2f)
    pf <- pf(Fstatf, ndff, ddff, lower.tail = FALSE)
    ddfr <- n-kr-1
    Fstatr <- r2r/kr * ddfr/(1-r2r)
    pr <- pf(Fstatr, kr, ddfr, lower.tail = FALSE)
    Fstate <- (r2f - r2r)/p * ddff/(1-r2f) # https://www.youtube.com/watch?v=Pz3j4Zu8BOQ
    pe <- pf(Fstate, p, ddff, lower.tail = FALSE)
    res <- matrix(c(r2f, ndff, ddff, Fstatf, pf,
                    r2r, kr, ddfr, Fstatr, pr,
                    r2f-r2r, p, ddff, Fstate, pe), nrow = 3L, ncol = 5L, byrow = TRUE,
                  dimnames = list(c("Full Model","Restricted Model","Exclusion Rest."),
                                  c("R-Sq.","DF1","DF2","F-Stat.","P-Value")))
    class(res) <- c("fFtest","matrix")
  } else {
    u <- fHDwithin.default(y, exc) # Residuals
    miss <- attr(u, "na.rm")
    if(full.df && !is.null(miss) && !is.numeric(exc)) {
      p <- if(is.factor(exc)) fnlevels(exc[-miss, drop = TRUE])-1L else if(any(vapply(exc, is.factor, TRUE)))
           getdf(droplevels.data.frame(exc[-miss, , drop = FALSE])) else length(exc)
    } else if(full.df) {
      p <- if(is.factor(exc) || (is.list(exc) && any(vapply(exc, is.factor, TRUE)))) getdf(droplevels(exc)) else NCOL(exc)
    } else p <- NCOL(exc)
    n <- length(u)
    r2 <- 1 - var(u)/var(if(is.null(miss)) y else y[-miss]) # R-Squared
    ddf <- n-p-1
    Fstat <- r2/p * ddf/(1-r2) # F statistic for the model (the constant goes unrestricted)
    Pv <- pf(Fstat, p, ddf, lower.tail = FALSE) # P-value corresponding to the F statistic
    res <- c(`R-Sq.` = r2, `DF1` = p, `DF2` = ddf, `F-Stat.` = Fstat, `P-value` = Pv)
    class(res) <- "fFtest"
  }
  return(res)
}

print.fFtest <- function(x, digits = 3, ...) {
  xx <- unclass(format(round(x, digits)))
  xpos <- x >= 1
  xx[xpos] <- sub(paste0(c(".", rep("0",digits)), collapse = ""), "", xx[xpos]) # Problem: Deletes .00 also..
  print.default(xx, quote = FALSE, right = TRUE, ...)
}

all_identical <- function(...) {
  if(length(match.call())-1L == 1L && is.list(...)) { # https://stackoverflow.com/questions/44011918/count-number-of-arguments-passed-to-function
    all(unlist(lapply(...[-1L], identical, ...[[1L]]), use.names = FALSE)) # use vapply ??
  } else {
    l <- list(...)
    all(unlist(lapply(l[-1L], identical, l[[1L]]), use.names = FALSE)) # use vapply ??
  }
}
all.identical <- all_identical
is.categorical <- function(x) !is.numeric(x)
is.Date <- function(x) inherits(x, c("Date","POSIXlt","POSIXct"))
"%!in%" <- function(x, table) match(x, table, nomatch = 0L) == 0L
na_rm <- function(x) x[!is.na(x)] # more consistent with base than na_rm !!! if not Cpp version that's fine !!
# na.rm <- function(x) { # cpp version available, but not faster !!
#   if(!is.null(attr(x, "names"))) { # gives corruped time-series !!
#     ax <- attributes(x)
#     r <- x[!is.na(x)]
#     ax[["names"]] <- names(r)
#     setAttributes(r, ax)
#   } else duplAttributes(x[!is.na(x)], x)
# }
na_insert <- function(X, prop = 0.1) {
  if(!is.null(d <- dim(X))) {
    n <- d[1L]
    p <- d[2L]
    NAloc <- rep(FALSE, n * p)
    NAloc[sample.int(n * p, floor(n * p * prop))] <- TRUE
    X[matrix(NAloc, nrow = n, ncol = p)] <- NA
  } else if(is.atomic(X)) {
    l <- length(X)
    X[sample.int(l, floor(l * prop))] <- NA
  } else stop("X must be an atomic vector, matrix or data.frame")
  return(X)
}
fnlevels <- function(x) length(attr(x, "levels")) # make cpp version ?? -> nope, slower !!
as.numeric_factor <- function(X) {
  if(is.atomic(X)) return(as.numeric(attr(X, "levels"))[X])

    fcts <- vapply(X, is.factor, TRUE, USE.NAMES = FALSE)
    # if(all(fcts)) return(dapply(X, function(x) as.numeric(attr(x, "levels"))[x]))
    clx <- class(X)
    class(X) <- NULL
    X[fcts] <- lapply(X[fcts], function(x) as.numeric(attr(x, "levels"))[x])
    return(`oldClass<-`(X, clx))
}




setRow.names <- function(df, nm) {
  attr(df, "row.names") <- nm
  df
}
remove_attributes <- function(x) {
  attributes(x) <- NULL
  x
}
addAttributes <- function(x, a) {
  ax <- attributes(x)
  .Call(Cpp_setAttributes, x, c(ax, a))
}
TRAtoInt <- function(x) # A lot faster than match based verion !!!
  switch(x, replace_fill = 1L, replace = 2L, `-` = 3L, `-+` = 4L, `/` = 5L, `%` = 6L, `+` = 7L, `*` = 8L,
            stop("Unknown transformation!"))
condsetn <- function(x, value, cond) {
  if(cond) attr(x, "names") <- value
  x
}
give_nam <- function(x, gn, stub) {
  if(!gn) return(x)
  attr(x, "names") <- paste0(stub, attr(x, "names"))
  x
}
anyNAerror <- function(x, e) if(anyNA(x)) stop(e) else x
cols2int <- function(cols, x, nam) {
 if(is.numeric(cols)) {
  if(max(abs(cols)) > length(x)) stop("Index out of range abs(1:length(x))")
  return(cols)
 } else if(is.function(cols))
  return(which(vapply(x, cols, TRUE))) else if(is.character(cols))
  return(anyNAerror(match(cols, nam), "Unknown column names!")) else if(is.logical(cols)) {
    if(length(cols) != length(x)) stop("Logical subsetting vector must match columns!")
    return(which(cols))
  } else stop("cols must be a function, character vector, numeric indices or logical vector!")
}
cols2log <- function(cols, x, nam) {
  if(is.logical(cols)) if(length(cols) == length(x)) return(cols) else stop("Logical subsetting vector must match columns!")
  if(is.function(cols)) return(vapply(x, cols, TRUE))
  r <- logical(length(x))
  if(is.character(cols)) {
    r[anyNAerror(match(cols, nam), "Unknown column names!")] <- TRUE
  } else if(is.numeric(cols)) {
    if(max(abs(cols)) > length(r)) stop("Index out of range abs(1:length(x))")
    r[cols] <- TRUE
  } else stop("cols must be a function, character vector, numeric indices or logical vector!")
  return(r)
}
colsubset <- function(x, ind) { # also works for grouped tibbles !!
  ax <- attributes(x)
  attributes(x) <- NULL # faster than unclass !! and good here since vapply on unclassed is faster !!
  if(is.numeric(ind)) {
    if(max(abs(ind)) > length(x)) stop("Index out of range abs(1:length(x))")
  } else if(is.logical(ind)) {
    if(length(ind) != length(x)) stop("Logical subsetting vector must match length(x)")
  } else {
    ind <- if(is.function(ind)) vapply(x, ind, TRUE, USE.NAMES = FALSE) else
      anyNAerror(match(ind, ax[["names"]]), "Unknown column names!")
  }
  ax[["names"]] <- ax[["names"]][ind]
  return(.Call(Cpp_setAttributes, x[ind], ax)) # return(`attributes<-`(x[ind], ax)) # This is slow on large data -> a lot of checks !!!
}

at2GRP <- function(x) {
  if(is.nmfactor(x))
  return(list(length(attr(x, "levels")), x, NULL)) else {
    res <- list(NULL, NULL, NULL)
    res[[2L]] <- qG(x, ordered = FALSE, na.exclude = FALSE)
    res[[1L]] <- attr(res[[2L]], "N.groups")
    return(res)
  }
}
# ret2int <- function(x) match(ret, c("cols","all","add")) # return = c("cols","all","add")
G_t <- function(x, m = TRUE) {
  if(is.null(x)) {
    if(m) message("Panel-lag computed without timevar: Assuming ordered data")
    return(x)
  } else if(is.atomic(x)) {
    if(is.integer(x)) return(x) else return(qG(x, na.exclude = FALSE)) # make sure it is ordered !!! qG already ckecks factor !!
  } else if(is.GRP(x)) return(x[[2L]]) else return(GRP(x, return.groups = FALSE)[[2L]])
}
rgrep <- function(exp, nam, ...) if(length(exp) > 1L) funique(unlist(lapply(exp, grep, nam, ...), use.names = FALSE)) else grep(exp, nam, ...)
NROW2 <- function(x, d) if(length(d)) d[1L] else length(x)
NCOL2 <- function(d, ilv) if(ilv) d[2L] else 1L
charorNULL <- function(x) if(is.character(x)) x else NULL
# more security here??
unique_factor <- function(x) {
  res <- seq_along(attr(x, "levels"))
  .Call(Cpp_duplAttributes, res, x)
}
dotstostr <- function(...) {
  args <- deparse(substitute(c(...)))
  nc <- nchar(args)
  substr(args, 2, nc) # 3, nc-1 for no brackets !!
}
is.nmfactor <- function(x) inherits(x, "factor") && (inherits(x, "na.included") || !anyNA(unclass(x)))
addNA2 <- function(x) {
  if(!anyNA(unclass(x))) return(x)
  ax <- attributes(x)
  if(!anyNA(ax[["levels"]])) ax[["levels"]] <- c(ax[["levels"]], NA)
  attributes(x) <- NULL
  x[is.na(x)] <- length(ax[["levels"]])
  return(setAttributes(x, ax))
}

# addNA2 <- function(x) {
#   clx <- c(class(x), "na.included")
#   if(!anyNA(unclass(x))) return(`oldClass<-`(x, clx))
#   ll <- attr(x, "levels")
#   if(!anyNA(ll)) ll <- c(ll, NA)
#   return(`oldClass<-`(factor(x, levels = ll, exclude = NULL), clx))
# }
