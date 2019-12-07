# library(Rcpp)
# sourceCpp("src/mrtl_mctl.cpp")
# sourceCpp("src/qF_qG.cpp", rebuild = TRUE) # https://gallery.rcpp.org/articles/fast-factor-generation/
qF <- function(x, ordered = TRUE, na.exclude = TRUE) {
  if(is.factor(x)) return(x)
  .Call(Cpp_qF, x, ordered, na.exclude)
}
qG <- function(x, ordered = TRUE, na.exclude = TRUE) {
  if(is.factor(x)) return(x)
  .Call(Cpp_qG, x, ordered, na.exclude)
}
# what about attribute preervation ??
# -> I think it is not good having all kinds of stuff attached to a matrix ??
# also dapply and BY convert known data types (matrix or data.frame. These functions can convert anything, so we need to be more careful !!)
qDF <- function(X) {
  if(is.atomic(X)) {
    d <- dim(X)
    ld <- length(d)
    if(ld == 2L)
    return(.Call(Cpp_mctl, X, TRUE, 1L)) else if (ld > 2L) {
      dn <- dimnames(X)
      dim(X) <- c(d[1L], prod(d[-1L]))
      if(!is.null(dn)) {
        for (i in 2L:ld) if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
        dimnames(X) <- list(dn[[1L]], interaction(expand.grid(dn[-1L]))) # Good??
      }
      return(.Call(Cpp_mctl, X, TRUE, 1L))
    } else {
      lx <- length(X)
      X <- `names<-`(list(X), deparse(substitute(X)))
      attr(X, "row.names") <- .set_row_names(lx)
      class(X) <- "data.frame"
      return(X)
    }
  } else {
    # if(inherits(X, "data.frame")) return(X)
    if(is.null(attr(X, "names"))) attr(X, "names") <- paste0("V", seq_along(X))
    if(is.null(attr(X, "row.names"))) attr(X, "row.names") <- .set_row_names(length(X[[1L]]))
    class(X) <- "data.frame"
    return(X)
  }
}
qDT <- function(X) { # what if already DT ?? return ??
  if(is.atomic(X)) {
    d <- dim(X)
    ld <- length(d)
    if(ld == 2L)
    return(.Call(Cpp_mctl, X, TRUE, 2L)) else if(ld > 2L) {
      dn <- dimnames(X)
      dim(X) <- c(d[1L], prod(d[-1L]))
      if(!is.null(dn)) {
        for (i in 2L:ld) if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
        dimnames(X) <- list(dn[[1L]], interaction(expand.grid(dn[-1L]))) # Good??
      }
      return(.Call(Cpp_mctl, X, TRUE, 2L))
    } else {
      lx <- length(X)
      X <- `names<-`(list(X), deparse(substitute(X)))
      attr(X, "row.names") <- .set_row_names(lx)
      class(X) <- c("data.table","data.frame")
      return(X)
    }
  } else {
    if(inherits(X, "data.table")) return(X)
    if(is.null(attr(X, "names"))) attr(X, "names") <- paste0("V", seq_along(X))
    attr(X, "row.names") <- .set_row_names(length(X[[1L]]))
    class(X) <- c("data.table","data.frame")
    return(X)
  }
}
qM <- function(X) {
  if(is.atomic(X)) {
    d <- dim(X)
    ld <- length(d)
    if(ld > 2L) {
      dn <- dimnames(X)
      dim(X) <- c(d[1L], prod(d[-1L]))
      if(!is.null(dn)) {
        for (i in 2L:ld) if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
        dimnames(X) <- list(dn[[1L]], interaction(expand.grid(dn[-1L]))) # Good??
      }
      return(X)
    } else if(ld == 2L) return(X) else
      return(matrix(X, ncol = 1, dimnames = list(NULL, deparse(substitute(X)))))
  } else {
    rn <- attr(X, "row.names")
    res <- do.call(cbind, X)
    if(!(is.null(rn) || is.integer(rn))) dimnames(res) <- list(rn, names(X))
    return(res)
  }
}


