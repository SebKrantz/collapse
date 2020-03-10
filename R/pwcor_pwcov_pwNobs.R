# sumcc <- function(x, y)  sum(complete.cases(x,y))
# pwNobs <- function(x) qM(dapply(x, function(y) dapply(x, sumcc, y)))

pwNobs <- function(X) { # Faster !!
  dg <- fNobs(X)
  if(is.matrix(X)) {
    cn <- dimnames(X)[[2L]]
    X <- mctl(X)
  } else {
    class(X) <- NULL
    cn <- names(X)
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
  } else if(is.list(x)) {
    if(spacing == 0L) res <- lapply(x, formfun) else {
      ff <- function(i) duplAttributes(sprintf(switch(spacing," %s","  %s",stop(se)),formfun(i)), i)
      res <- lapply(x, ff)
    }
    if(show != 1L) res <- if(show == 2L) lapply(res, function(i){i[upper.tri(i)] <- ""; i}) else
      lapply(res, function(i){i[lower.tri(i)] <- ""; i})
  } else res <- formfun(x)
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
  } else if(is.list(x)) {
    if(spacing == 0L) res <- lapply(x, formfun, TRUE) else {
      ff <- function(i) duplAttributes(sprintf(switch(spacing," %s","  %s",stop(se)),formfun(i, TRUE)), i)
      res <- lapply(x, ff)
    }
    if(show != 1L) res <- if(show == 2L) lapply(res, function(i){i[upper.tri(i)] <- ""; i}) else
      lapply(res, function(i){i[lower.tri(i)] <- ""; i})
  } else res <- formfun(x)
  print.default(unclass(res), quote = FALSE, right = TRUE, ...)
  invisible(x)
} #print.table(dapply(round(x, digits), function(j) sub("^(-?)0.", "\\1.", j)), right = TRUE, ...) # print.table(, right = TRUE)


# print.pwcov <- function(x, digits = 2, ...) print.default(formatC(round(x, digits), format = "g",
#                                                                   digits = 9, big.mark = ",", big.interval = 6), quote = FALSE, right = TRUE, ...)
