
getdf <- function(x) {
  if(is.atomic(x)) if(is.factor(x)) return(fnlevels(x)-1L) else return(1L)
  bsum(vapply(unattrib(x), function(i) if(is.factor(i)) fnlevels(i)-1L else 1L, 1L))
}

fFtest <- function(...) if(is.call(..1) || is.call(..2)) fFtest.formula(...) else fFtest.default(...)

fFtest.default <- function(y, exc, X = NULL, w = NULL, full.df = TRUE, ...) {
  if(!is.numeric(y)) stop("y needs to be a numeric vector")
  if(!is.null(X)) {
    Xn <- fNCOL(X)
    atl <- is.atomic(X) && is.numeric(X) && is.atomic(exc) && is.numeric(exc)
    if(length(w)) {
      if(atl) {
        cc <- which(complete.cases(w, y, X, exc))
        if(length(cc) < length(w)) {
          data <- cbind(y, X, exc)[cc, , drop = FALSE]
          w <- w[cc]
        }
      } else {
        data <- na_omit(qDF(c(list(w = w), list(y = y), qDF(X), qDF(exc))))
        w <- .subset2(data, 1L)
        data[[1L]] <- NULL
      }
    } else {
      data <- if(atl) na_omit(cbind(y, X, exc)) else na_omit(qDF(c(list(y = y), qDF(X), qDF(exc))))
    }
    if(full.df && !atl && any(fc <- .Call(C_vtypes, data, 2L))) { # vapply(unattrib(data), is.factor, TRUE)
      cld <- oldClass(data)
      oldClass(data) <- NULL
      data[fc] <- lapply(data[fc], fdroplevels.factor)
      df <- vapply(unattrib(data), function(i) if(is.factor(i)) fnlevels(i)-1L else 1L, 1L) # getdf(data)
      k <- bsum(df) # 1 for intercept added with y
      p <- bsum(df[(Xn+2L):length(df)])
      y <- data[[1L]]
      oldClass(data) <- cld
    } else {
      p <- fNCOL(exc)
      if(atl) {
        k <- ncol(data) # 1 for intercept added with y
        y <- data[, 1L]
      } else {
        k <- length(unclass(data)) # 1 for intercept added with y
        y <- .subset2(data, 1L)
      }
    }
    kr <- k-p-1L
    vy <- fvar.default(y, w = w)
    if(atl) {
      n <- nrow(data)
      r2f <- 1 - fvar.default(fhdwithin.default(y, data[, -1L], w, na.rm = FALSE, ...), w = w)/vy
      r2r <- 1 - fvar.default(fhdwithin.default(y, data[, 2:(Xn+1L)], w, na.rm = FALSE, ...), w = w)/vy
    } else {
      n <- fnrow2(data)
      r2f <- 1 - fvar.default(fhdwithin.default(y, fcolsubset(data, -1L), w, na.rm = FALSE, ...), w = w)/vy
      r2r <- 1 - fvar.default(fhdwithin.default(y, fcolsubset(data, 2:(Xn+1L)), w, na.rm = FALSE, ...), w = w)/vy
    }
    ndff <- k-1L
    ddff <- n-k
    Fstatf <- r2f/ndff * ddff/(1-r2f)
    pf <- pf(Fstatf, ndff, ddff, lower.tail = FALSE)
    ddfr <- n-kr-1L
    Fstatr <- r2r/kr * ddfr/(1-r2r)
    pr <- pf(Fstatr, kr, ddfr, lower.tail = FALSE)
    Fstate <- (r2f - r2r)/p * ddff/(1-r2f)
    pe <- pf(Fstate, p, ddff, lower.tail = FALSE)
    res <- matrix(c(r2f, ndff, ddff, Fstatf, pf,
                    r2r, kr, ddfr, Fstatr, pr,
                    r2f-r2r, p, ddff, Fstate, pe), nrow = 3L, ncol = 5L, byrow = TRUE,
                  dimnames = list(c("Full Model","Restricted Model","Exclusion Rest."),
                                  c("R-Sq.","DF1","DF2","F-Stat.","P-Value")))
    oldClass(res) <- c("fFtest","matrix")
  } else {
    u <- fhdwithin.default(y, exc, w, na.rm = TRUE, ...) # Residuals
    miss <- attr(u, "na.rm")
    if(!is.null(miss)) w <- w[-miss]
    if(full.df && length(miss) && !is.atomic(exc) && !is.numeric(exc)) {
      p <- if(is.factor(exc)) fnlevels(exc[-miss, drop = TRUE])-1L else if(any(.Call(C_vtypes, exc, 2L))) # vapply(unattrib(exc), is.factor, TRUE)
        getdf(fdroplevels.data.frame(ss(exc, -miss))) else length(unclass(exc))
    } else if(full.df) {
      p <- if(is.factor(exc) || (is.list(exc) && any(.Call(C_vtypes, exc, 2L)))) getdf(fdroplevels(exc)) else fNCOL(exc) # vapply(unattrib(exc), is.factor, TRUE)
    } else p <- fNCOL(exc)
    n <- length(u)
    r2 <- 1 - fvar.default(u, w = w)/fvar.default(if(is.null(miss)) y else y[-miss], w = w) # R-Squared
    ddf <- n-p-1L
    Fstat <- r2/p * ddf/(1-r2) # F statistic for the model (the constant goes unrestricted)
    Pv <- pf(Fstat, p, ddf, lower.tail = FALSE) # P-value corresponding to the F statistic
    res <- c(`R-Sq.` = r2, `DF1` = p, `DF2` = ddf, `F-Stat.` = Fstat, `P-value` = Pv)
    oldClass(res) <- "fFtest"
  }
  res
}

fFtest.formula <- function(formula, data = NULL, weights = NULL, ...) {
  w <- substitute(weights)
  pe <- parent.frame()
  if(length(w)) w <- eval(w, data, pe)
  if(!any(all.names(formula) == "|")) { # Standard formula (no X term)
    tms <- attributes(terms.formula(formula, data = data))
    mf <- eval(tms$variables, data, pe)
    exc <- mf[-1L]
    names(exc) <- tms$term.labels
    return(fFtest.default(mf[[1L]], exc, NULL, w, ...))
  }
  y <- eval(formula[[2L]], data, pe)
  fml <- formula[[3L]]
  exc <- attributes(terms.formula(call("~", fml[[2L]]), data = data))
  exc <- eval(exc$variables, data, pe)
  X <- attributes(terms.formula(call("~", fml[[3L]]), data = data))
  X <- eval(X$variables, data, pe)
  fFtest.default(y, exc, X, w, ...)
}

print.fFtest <- function(x, digits = 3, ...) {
  xx <- unclass(format(round(x, digits)))
  xpos <- x >= 1
  xx[xpos] <- sub(paste0(c(".", rep("0",digits)), collapse = ""), "", xx[xpos]) # Problem: Deletes .00 also..
  print.default(xx, quote = FALSE, right = TRUE, ...)
}
