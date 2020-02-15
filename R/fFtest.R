
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
