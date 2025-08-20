

qDF <- function(X, row.names.col = FALSE, keep.attr = FALSE, class = "data.frame") {
  if(is.atomic(X)) {
    d <- dim(X)
    ld <- length(d)
    if(ld > 1L) {
      if(ld > 2L) {
        dn <- dimnames(X)
        dim(X) <- c(d[1L], bprod(d[-1L]))
        if(length(dn)) {
          for (i in 2L:ld) if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
          dimnames(X) <- list(dn[[1L]], interact_names(dn[-1L])) # Good?
        }
      }
      if(!isFALSE(row.names.col) && length(force(dn <- dimnames(X))[[1L]])) {
        res <- c(list(dn[[1L]]), .Call(Cpp_mctl, X, FALSE, 0L))
        names(res) <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", dn[[2L]])
        attr(res, "row.names") <- .set_row_names(length(dn[[1L]]))
      } else res <- .Call(Cpp_mctl, X, TRUE, 1L)
      oldClass(res) <- if(length(class)) class else "data.frame"
      if(!keep.attr) return(res)
      ax <- attributes(X)
      axoth <- names(ax) %!in% c("dim", "dimnames", "class")
      if(any(axoth)) return(addAttributes(res, ax[axoth])) else return(res)
    }
    nam <- names(X)
    if(is.null(nam) || isFALSE(row.names.col)) {
      if(is.null(nam)) {
        res <- `names<-`(list(X), l1orlst(as.character(substitute(X))))
        attr(res, "row.names") <- .set_row_names(length(X))
      } else {
        res <- `names<-`(list(`names<-`(X, NULL)), l1orlst(as.character(substitute(X))))
        attr(res, "row.names") <- nam
      }
    } else {
      res <- list(nam, `names<-`(X, NULL))
      names(res) <- if(length(row.names.col) == 2L) row.names.col else c(
                    if(is.character(row.names.col)) row.names.col[1L] else "row.names",
                    l1orlst(as.character(substitute(X))))
      attr(res, "row.names") <- .set_row_names(length(X))
    }
    return(`oldClass<-`(res, if(length(class)) class else "data.frame"))
  }
  if(keep.attr) {
    # if(all(class(X) == class)) return(X) # better adjust rows ? -> yes, row.names.col should always work !
    if(is.null(attr(X, "names"))) attr(X, "names") <- paste0("V", seq_along(unclass(X)))
    if(is.null(attr(X, "row.names"))) {
      attr(X, "row.names") <- .set_row_names(fnrow(X))
    } else if(!isFALSE(row.names.col)) {
      ax <- attributes(X)
      X <- c(list(ax[["row.names"]]), X)
      ax[["row.names"]] <- .set_row_names(.Call(C_fnrow, X)) # this is ok, X is a list ...
      ax[["names"]] <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", ax[["names"]])
      setattributes(X, ax)
    }
    if(length(class)) return(`oldClass<-`(X, class))
    if(inherits(X, "data.frame")) return(X)
    return(`oldClass<-`(X, "data.frame"))
  }
  nam <- attr(X, "names")
  rn <- attr(X, "row.names")
  attributes(X) <- NULL
  if(is.null(nam)) nam <- paste0("V", seq_along(X))
  if(is.null(rn) || is.numeric(rn)) {
    rn <- .set_row_names(.Call(C_fnrow, X))
  } else if(!isFALSE(row.names.col)) {
    X <- c(list(rn), X)
    rn <- .set_row_names(.Call(C_fnrow, X))
    nam <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", nam)
  }
  # slower: !!
  # setAttributes(X, pairlist(names = nam, row.names = rn, class = if(length(class)) class else "data.frame"))
  names(X) <- nam
  attr(X, "row.names") <- rn # This can be inefficient for large data.frames if character rn !!
  oldClass(X) <- if(length(class)) class else "data.frame"
  X
}

qDT_raw <- function(X, row.names.col, keep.attr, DT_class, X_nam) {
  if(is.atomic(X)) {
    d <- dim(X)
    ld <- length(d)
    if(ld > 1L) {
      if(ld > 2L) {
        dn <- dimnames(X)
        dim(X) <- c(d[1L], bprod(d[-1L]))
        if(length(dn)) {
          for (i in 2L:ld) if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
          dimnames(X) <- list(dn[[1L]], interact_names(dn[-1L])) # Good?
        }
      }
      if(!isFALSE(row.names.col) && length(force(dn <- dimnames(X))[[1L]])) {
        res <- c(list(dn[[1L]]), .Call(Cpp_mctl, X, FALSE, 0L))
        names(res) <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", dn[[2L]])
        attr(res, "row.names") <- .set_row_names(length(dn[[1L]]))
      } else res <- .Call(Cpp_mctl, X, TRUE, 2L)
      oldClass(res) <- DT_class
      if(!keep.attr) return(res)
      ax <- attributes(X)
      axoth <- names(ax) %!in% c("dim", "dimnames", "class")
      return(if(any(axoth)) addAttributes(res, ax[axoth]) else res)
    }
    if(isFALSE(row.names.col) || is.null(nam <- names(X))) {
      res <- `names<-`(list(X), X_nam)
    } else {
      res <- list(nam, `names<-`(X, NULL))
      names(res) <- if(length(row.names.col) == 2L) row.names.col else c(
        if(is.character(row.names.col)) row.names.col[1L] else "row.names", X_nam)
    }
    attr(res, "row.names") <- .set_row_names(length(X))
    return(`oldClass<-`(res, DT_class))
  }
  if(keep.attr) {
    # if(all(class(X) == DT_class)) return(X) # better adjust rows ? -> yes, row.names.col should always work !
    if(is.null(attr(X, "names"))) attr(X, "names") <- paste0("V", seq_along(unclass(X)))
    if(!isFALSE(row.names.col) && length(rn <- attr(X, "row.names"))) {
      ax <- attributes(X)
      X <- c(list(rn), X)
      ax[["names"]] <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", ax[["names"]])
      setattributes(X, ax)
    }
    if(!length(DT_class) && inherits(X, c("data.table", "data.frame"))) return(X)
    attr(X, "row.names") <- .set_row_names(fnrow(X))
  } else {
    nam <- attr(X, "names")
    rncol <- !isFALSE(row.names.col) && length(rn <- attr(X, "row.names"))
    attributes(X) <- NULL
    if(is.null(nam)) nam <- paste0("V", seq_along(X))
    if(rncol) {
      X <- c(list(rn), X)
      nam <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", nam)
    }
    names(X) <- nam
    attr(X, "row.names") <- .set_row_names(.Call(C_fnrow, X))
  }
  return(`oldClass<-`(X, DT_class))
}

qDT <- function(X, row.names.col = FALSE, keep.attr = FALSE, class = c("data.table", "data.frame")) {
   alc(qDT_raw(X, row.names.col, keep.attr,
               if(length(class) || keep.attr) class else c("data.table", "data.frame"),
               if(is.atomic(X) && !is.matrix(X)) l1orlst(as.character(substitute(X))) else NULL))
}

qTBL <- function(X, row.names.col = FALSE, keep.attr = FALSE, class = c("tbl_df", "tbl", "data.frame")) {
       qDT_raw(X, row.names.col, keep.attr,
               if(length(class) || keep.attr) class else c("tbl_df", "tbl", "data.frame"),
               if(is.atomic(X) && !is.matrix(X)) l1orlst(as.character(substitute(X))) else NULL)
}

qM <- function(X, row.names.col = NULL, keep.attr = FALSE, class = NULL, sep = ".") {
  if(keep.attr) {
    if(is.atomic(X)) {
      if(length(class)) oldClass(X) <- class
      if(is.matrix(X)) return(X)
      if(is.array(X)) {
        d <- dim(X)
        dn <- dimnames(X)
        dim(X) <- c(d[1L], bprod(d[-1L]))
        if(length(dn)) {
          for (i in 2L:length(d)) if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
          dimnames(X) <- list(dn[[1L]], interact_names(dn[-1L])) # Good?
        }
      } else {
        nam <- l1orlst(as.character(substitute(X))) # needed before X is changed !!
        dim(X) <- c(length(X), 1L)
        dimnames(X) <- list(names(X), nam)
        names(X) <- NULL
        # if(is.object(X)) oldClass(X) <- NULL Necessary ? Can also have factor or date matrices. Check this !
        # -> qM(wlddev$date, TRUE) is a vector !!
      }
      return(X)
    }
    ax <- attributes(X)
    if(length(row.names.col)) {
      rnc <- cols2int(row.names.col, X, ax[["names"]])
      res <- do.call(cbind, .subset(X, -rnc))
      dimnames(res)[[1L]] <- if(length(rnc) == 1L) .subset2(X, rnc) else
        do.call(paste, c(.subset(X, rnc), list(sep = sep)))
    } else {
      res <- do.call(cbind, X)
      rn <- ax[["row.names"]]
      if(!(is.numeric(rn) || is.null(rn) || rn[1L] == "1"))
         dimnames(res) <- list(rn, ax[["names"]])
    }
    if(length(class)) oldClass(res) <- class
    axoth <- names(ax) %!in% c("names", "row.names", "class")
    if(any(axoth)) return(addAttributes(res, ax[axoth]))
    return(res)
  }
  if(is.atomic(X)) {
    if(!is.array(X)) {
      r <- matrix(X, ncol = 1, dimnames = list(names(X), l1orlst(as.character(substitute(X)))))
      if(is.null(class)) return(r) else return(`oldClass<-`(r, class))
    }
    d <- dim(X)
    dn <- dimnames(X)
    attributes(X) <- NULL
    ld <- length(d)
    if(ld == 2L) {
      # setattributes(X, pairlist(dim = d, dimnames = dn)) # Not faster !
      dim(X) <- d
      dimnames(X) <- dn
    } else {
      dim(X) <- c(d[1L], bprod(d[-1L]))
      if(length(dn)) {
        for (i in 2L:ld) if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
        dimnames(X) <- list(dn[[1L]], interact_names(dn[-1L])) # Good?
      }
    }
    if(length(class)) oldClass(X) <- class
    return(X)
  }
  if(length(row.names.col)) {
    rnc <- cols2int(row.names.col, X, attr(X, "names"))
    res <- do.call(cbind, .subset(X, -rnc))
    if(is.object(res)) attributes(res) <- attributes(res)[c("dim", "dimnames")]
    dimnames(res)[[1L]] <- if(length(rnc) == 1L) .subset2(X, rnc) else
      do.call(paste, c(.subset(X, rnc), list(sep = sep)))
  } else {
    rn <- attr(X, "row.names")
    res <- do.call(cbind, X)
    if(is.object(res)) attributes(res) <- attributes(res)[c("dim", "dimnames")] # if X is list of time-series, do.call(cbind, X) creates ts-matrix.
    if(!(is.numeric(rn) || is.null(rn) || rn[1L] == "1"))
      dimnames(res) <- list(rn, attr(X, "names"))
  }
  if(length(class)) oldClass(res) <- class
  res
}

# Same speed
# tf1 <- function(res)  {
#   res <- do.call(cbind, res)
#   if(is.object(res)) attributes(res) <- attributes(res)[c("dim", "dimnames")]
#   res
# }
#
# tf2 <- function(res)  {
#   res <- do.call(cbind, res)
#   if(is.object(res)) setAttributes(res, attributes(res)[c("dim", "dimnames")])
# }
