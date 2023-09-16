# sumcc <- function(x, y)  bsum(complete.cases(x,y))
# pwnobs <- function(x) qM(dapply(x, function(y) dapply(x, sumcc, y)))

pwnobs <- function(X) {
  if(is.atomic(X) && is.matrix(X)) return(.Call(Cpp_pwnobsm, X)) # cn <- dimnames(X)[[2L]] # X <- mctl(X)
  if(!is.list(X)) stop("X must be a matrix or data.frame!") # -> if unequal length will warn below !!
  dg <- fnobs.data.frame(X)
  oldClass(X) <- NULL
  n <- length(X)
  nr <- .Call(C_fnrow, X)
  N.mat <- diag(dg)
  for (i in 1:(n - 1L)) {
    miss <- is.na(X[[i]]) # faster than complete.cases, also for large data ! // subsetting X[[j]] faster ?? -> NOPE !
    for (j in (i + 1L):n) N.mat[i, j] <- N.mat[j, i] <- nr - bsum(miss | is.na(X[[j]])) # bsum(complete.cases(X[[i]], X[[j]]))
  }
  dimnames(N.mat) <- list(names(dg), names(dg))
  N.mat
}

pwNobs <- function(X) {
  .Deprecated(msg = "'pwNobs' was renamed to 'pwnobs'. It will be removed end of 2023, see help('collapse-renamed').")
  pwnobs(X)
}
# corr.p <- function(r, n) {
#   if (n < 3L) return(1)
#   df <- n - 2L
#   t <- sqrt(df) * r/sqrt(1 - r^2)
#   return(2 * bmin(pt(t, df), pt(t, df, lower.tail = FALSE))) # taken from corr.test
# }

corr.pmat <- function(cm, nm) {
  df <- nm - 2L
  acm <- abs(cm)
  diag(acm) <- NA_real_ # tiny bit faster here vs below..
  `attributes<-`(2 * pt(sqrt(df) * acm/sqrt(1 - acm^2), df, lower.tail = FALSE),
                 attributes(cm))
  # n <- ncol(cm)
  # p.mat <- matrix(NA, n, n, dimnames = dimnames(cm))
  # for (i in 1:(n - 1)) {
  #   for (j in (i + 1):n) {
  #     p.mat[i, j] <- p.mat[j, i] <- corr.p(cm[i, j], nm[i, j])
  #   }
  # }
  # p.mat
}

complpwnobs <- function(X) {
  # if(is.list(X)) { # Not needed anymore because now always coercing to matrix...
  #   n <- length(unclass(X))
  #   coln <- attr(X, "names")
  # } else {
    n <- ncol(X)
    coln <- dimnames(X)[[2L]]
  # }
  matrix(bsum(complete.cases(X)), n, n, dimnames = list(coln, coln))
}

# Test:
# all.equal(Hmisc::rcorr(qM(mtcars))$P, corr.pmat(r, n))
namat <- function(X) {
  nc <- dim(X)[2L]
  cn <- dimnames(X)[[2L]]
  mat <- rep(NA_real_, nc * nc)
  dim(mat) <- c(nc, nc)
  diag(mat) <- 1
  dimnames(mat) <- list(cn, cn)
  mat
}

nmat <- function(n, X) {
  nc <- dim(X)[2L]
  cn <- dimnames(X)[[2L]]
  mat <- rep(n, nc * nc)
  dim(mat) <- c(nc, nc)
  dimnames(mat) <- list(cn, cn)
  mat
}

# Check speed of it ...
# Also check weighted cor p-value against lm() with weights -> Good !!

# -> This is good
# all.equal(unattrib(cov.wt(mtcars, w, cor = TRUE)$cor), unattrib(pwcor(mtcars, w = w)))
# all.equal(unattrib(cov.wt(mtcars, w, cor = TRUE)$cor), unattrib(pwcor(mtcars, w = w, use = "complete.obs")))
# all.equal(pwcor(mtcars, w = w), pwcor(mtcars, w = w, use = "complete.obs"))
pwcor <- function(X, ..., w = NULL, N = FALSE, P = FALSE, array = TRUE, use = "pairwise.complete.obs") {
  if(is.list(X)) X <- do.call(cbind, X)
  lcc <- FALSE
  if(is.null(w)) r <- cor(X, ..., use = use) else if(use == "pairwise.complete.obs")
  r <- getenvFUN("weights_wtd.cors")(X, ..., weight = w) else {
    if(!missing(...)) stop("y is currently not supported with weighted correlations and use != 'pairwise.complete.obs'")
    cc <- which(complete.cases(X, w))
    lcc <- length(cc)
    if(use == "all.obs" && lcc != length(w)) stop("missing observations in cov/cor")
    if(lcc) {
      if(lcc != length(w)) {
        X <- X[cc, , drop = FALSE]
        w <- w[cc]
      }
      r <- cov2cor(crossprod(sqrt(w) * BWmCpp(X, w = w, narm = FALSE))) # all.equal(cov2cor(crossprod(sqrt(w) * BWmCpp(X, w = w, narm = FALSE))), weights::wtd.cors(X, weight = w))
    } else r <- switch(use, complete.obs = stop("no complete element pairs"), namat(X))
  }
  if(!(N || P)) return(`oldClass<-`(r, c("pwcor", "matrix")))
  n <- if(lcc) nmat(lcc, X) else switch(use, pairwise.complete.obs = pwnobs(X), complpwnobs(X)) # TODO: what about weights paiwrise ? # what if using ... to supply y ???
  if(N) {
    res <- if(P) list(r = r, N = n, P = corr.pmat(r, n)) else list(r = r, N = n)
  } else res <- list(r = r, P = corr.pmat(r, n))
  if(array) {
    res <- fsimplify2array(res)
    oldClass(res) <- c("pwcor","array","table")
  } else oldClass(res) <- "pwcor"
  res
}

# Not all equal...
# all.equal(unattrib(cov.wt(mtcars, w)$cov), unattrib(pwcov(mtcars, w = w)))
# all.equal(unattrib(cov.wt(mtcars, w)$cov), unattrib(pwcov(mtcars, w = w, use = "complete.obs")))
# all.equal(pwcov(mtcars, w = w), pwcov(mtcars, w = w, use = "complete.obs")) -> Yes !

pwcov <- function(X, ..., w = NULL, N = FALSE, P = FALSE, array = TRUE, use = "pairwise.complete.obs") {
  if(is.list(X)) X <- do.call(cbind, X)
  lcc <- FALSE
  if(is.null(w)) r <- cov(X, ..., use = use) else if(use == "pairwise.complete.obs") {
    r <- getenvFUN("weights_wtd.cors")(X, ..., weight = w)
    # sw <- bsum(w, na.rm = TRUE)
    Xsd <- fsd(X, w = w) # * (sw-1) / (1 - bsum((w/sw)^2)) # cov.wt, method = "unbiased" ???
    r <- if(missing(...)) r * outer(Xsd, Xsd) else r * outer(Xsd, fsd(..., w = w))
  } else {
    if(!missing(...)) stop("y is currently not supported with weighted correlations and use != 'pairwise.complete.obs'")
    cc <- which(complete.cases(X, w))
    lcc <- length(cc)
    if(use == "all.obs" && lcc != length(w)) stop("missing observations in cov/cor")
    if(lcc) {
      if(lcc != length(w)) {
        X <- X[cc, , drop = FALSE]
        w <- w[cc]
      }
      r <- crossprod(sqrt(w) * BWmCpp(X, w = w, narm = FALSE)) / (bsum(w) - 1) # Check numeric accuracy !
      # w <- w/bsum(w) # same method as cov.wt, method = "unbiased"
      # r <- crossprod(sqrt(w) * BWmCpp(X, w = w, narm = FALSE)) / (1 - bsum(w^2))
    } else r <- switch(use, complete.obs = stop("no complete element pairs"), namat(X)) # namat correct ??
  }
  if(!(N || P)) return(`oldClass<-`(r, c("pwcov", "matrix")))
  n <- if(lcc) nmat(lcc, X) else switch(use, pairwise.complete.obs = pwnobs(X), complpwnobs(X)) # TODO: what about weights paiwrise ?
  if(N) {                                           # good ??? // cov(X) / outer(fsd(X), fsd(X))
    res <- if(P) list(cov = r, N = n, P = corr.pmat(cov2cor(r), n)) else list(cov = r, N = n) # what about x and y here ??
  } else res <- list(cov = r, P = corr.pmat(cov2cor(r), n))
  if(array) {
    res <- fsimplify2array(res)
    oldClass(res) <- c("pwcov","array","table")
  } else oldClass(res) <- "pwcov"
  res
}

print.pwcor <- function(x, digits = .op[["digits"]], sig.level = 0.05, show = c("all","lower.tri","upper.tri"), spacing = 1L, return = FALSE, ...) {
  formfun <- function(x, dg1 = FALSE) {
    xx <- format(round(x, digits)) # , digits = digits-1
    xx <- sub("(-?)0\\.", "\\1.", xx)
    if(dg1) {
      dgx <- diag(xx)
      new1 <- paste0(c("  1", rep(" ",digits-1)), collapse = "")
      if(!all(st <- startsWith(dgx, " 1") | startsWith(dgx, "1"))) { # can have positive or negative values...
        dgx[st] <- new1
        diag(xx) <- dgx
      } else diag(xx) <- new1
    } else {
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
  if(return) return(unclass(res))
  print.default(unclass(res), quote = FALSE, right = TRUE, ...)
  invisible(x)
} #print.table(dapply(round(x, digits), function(j) sub("^(-?)0.", "\\1.", j)), right = TRUE, ...) # print.table(, right = TRUE)


print.pwcov <- function(x, digits = .op[["digits"]], sig.level = 0.05, show = c("all","lower.tri","upper.tri"), spacing = 1L, return = FALSE, ...) {
  formfun <- function(x, adj = FALSE) {
    xx <- format(round(x, digits), digits = 9, big.mark = "'", big.interval = 6)
    # xx <- sub("(-?)0\\.", "\\1.", xx) # Not needed here...
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
  if(return) return(unclass(res))
  print.default(unclass(res), quote = FALSE, right = TRUE, ...)
  invisible(x)
} #print.table(dapply(round(x, digits), function(j) sub("^(-?)0.", "\\1.", j)), right = TRUE, ...) # print.table(, right = TRUE)

# print.pwcov <- function(x, digits = 2, ...) print.default(formatC(round(x, digits), format = "g",
#                         digits = 9, big.mark = "'", big.interval = 6), quote = FALSE, right = TRUE, ...)


`[.pwcor` <- `[.pwcov` <- function(x, i, j, ..., drop = TRUE) `oldClass<-`(NextMethod(), oldClass(x))


