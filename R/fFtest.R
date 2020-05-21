
getdf <- function(x) {
  if(is.atomic(x)) if(is.factor(x)) return(fnlevels(x)-1L) else return(1L)
  sum(vapply(unattrib(x), function(i) if(is.factor(i)) fnlevels(i)-1L else 1L, 1L))
}


fFtest <- function(y, exc, X = NULL, full.df = TRUE, ...) {
  if(!is.numeric(y)) stop("y needs to be a numeric vector")
  if(!is.null(X)) {
    Xn <- fNCOL(X)
    atl <- is.atomic(X) && is.numeric(X) && is.atomic(exc) && is.numeric(exc)
    data <- if(atl) na_omit(cbind(y, X, exc)) else na_omit(qDF(c(list(y = y), qDF(X), qDF(exc))))
    if(full.df && !atl && any(fc <- vapply(unattrib(data), is.factor, TRUE))) {
      cld <- class(data)
      class(data) <- NULL
      data[fc] <- lapply(data[fc], droplevels.factor)
      df <- vapply(unattrib(data), function(i) if(is.factor(i)) fnlevels(i)-1L else 1L, 1L) # getdf(data)
      k <- sum(df) # 1 for intercept added with y
      p <- sum(df[(Xn+2L):length(df)])
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
    kr <- k-p-1
    vy <- var(y)
    if(atl) {
      n <- nrow(data)
      r2f <- 1 - var(fHDwithin.default(y, data[, -1L], na.rm = FALSE, ...))/vy
      r2r <- 1 - var(fHDwithin.default(y, data[, 2:(Xn+1L)], na.rm = FALSE, ...))/vy
    } else {
      n <- fnrow2(data)
      r2f <- 1 - var(fHDwithin.default(y, fcolsubset(data, -1L), na.rm = FALSE, ...))/vy
      r2r <- 1 - var(fHDwithin.default(y, fcolsubset(data, 2:(Xn+1L)), na.rm = FALSE, ...))/vy
    }
    ndff <- k-1
    ddff <- n-k
    Fstatf <- r2f/ndff * ddff/(1-r2f)
    pf <- pf(Fstatf, ndff, ddff, lower.tail = FALSE)
    ddfr <- n-kr-1
    Fstatr <- r2r/kr * ddfr/(1-r2r)
    pr <- pf(Fstatr, kr, ddfr, lower.tail = FALSE)
    Fstate <- (r2f - r2r)/p * ddff/(1-r2f)
    pe <- pf(Fstate, p, ddff, lower.tail = FALSE)
    res <- matrix(c(r2f, ndff, ddff, Fstatf, pf,
                    r2r, kr, ddfr, Fstatr, pr,
                    r2f-r2r, p, ddff, Fstate, pe), nrow = 3L, ncol = 5L, byrow = TRUE,
                  dimnames = list(c("Full Model","Restricted Model","Exclusion Rest."),
                                  c("R-Sq.","DF1","DF2","F-Stat.","P-Value")))
    class(res) <- c("fFtest","matrix")
  } else {
    u <- fHDwithin.default(y, exc, na.rm = TRUE) # Residuals
    miss <- attr(u, "na.rm")
    if(full.df && !is.null(miss) && !is.atomic(exc) && !is.numeric(exc)) {
      p <- if(is.factor(exc)) fnlevels(exc[-miss, drop = TRUE])-1L else if(any(vapply(unattrib(exc), is.factor, TRUE)))
        getdf(droplevels.data.frame(ss(exc, -miss))) else length(unclass(exc))
    } else if(full.df) {
      p <- if(is.factor(exc) || (is.list(exc) && any(vapply(unattrib(exc), is.factor, TRUE)))) getdf(droplevels(exc)) else fNCOL(exc)
    } else p <- fNCOL(exc)
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
