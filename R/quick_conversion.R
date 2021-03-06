

qDF <- function(X, row.names.col = FALSE, keep.attr = FALSE, class = "data.frame") {
  if(is.atomic(X)) {
    d <- dim(X)
    ld <- length(d)
    if(ld > 1L) {
      if(ld > 2L) {
        dn <- dimnames(X)
        dim(X) <- c(d[1L], prod(d[-1L]))
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
      names(res) <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", l1orlst(as.character(substitute(X))))
      attr(res, "row.names") <- .set_row_names(length(X))
    }
    return(`oldClass<-`(res, if(length(class)) class else "data.frame"))
  }
  if(keep.attr) {
    # if(all(class(X) == class)) return(X) # better adjust rows ? -> yes, row.names.col should always work !
    if(is.null(attr(X, "names"))) attr(X, "names") <- paste0("V", seq_along(unclass(X)))
    if(is.null(attr(X, "row.names"))) {
      attr(X, "row.names") <- .set_row_names(length(.subset2(X, 1L)))
    } else if(!isFALSE(row.names.col)) {
      ax <- attributes(X)
      X <- c(list(ax[["row.names"]]), X)
      ax[["row.names"]] <- .set_row_names(length(X[[1L]])) # this is ok, X is a list ...
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
    rn <- .set_row_names(length(X[[1L]]))
  } else if(!isFALSE(row.names.col)) {
    X <- c(list(rn), X)
    rn <- .set_row_names(length(X[[1L]]))
    nam <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", nam)
  }
  # slower: !!
  # setAttributes(X, pairlist(names = nam, row.names = rn, class = if(length(class)) class else "data.frame"))
  names(X) <- nam
  attr(X, "row.names") <- rn # This can be inefficient for large data.frames if character rn !!
  oldClass(X) <- if(length(class)) class else "data.frame"
  X
}


qDT <- function(X, row.names.col = FALSE, keep.attr = FALSE, class = c("data.table", "data.frame")) {
  if(is.atomic(X)) {
    d <- dim(X)
    ld <- length(d)
    if(ld > 1L) {
      if(ld > 2L) {
        dn <- dimnames(X)
        dim(X) <- c(d[1L], prod(d[-1L]))
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
      oldClass(res) <- if(length(class)) class else c("data.table", "data.frame")
      if(!keep.attr) return(res)
      ax <- attributes(X)
      axoth <- names(ax) %!in% c("dim", "dimnames", "class")
      if(any(axoth)) return(addAttributes(res, ax[axoth])) else return(res)
    }
    if(isFALSE(row.names.col) || is.null(nam <- names(X))) {
      res <- `names<-`(list(X), l1orlst(as.character(substitute(X))))
    } else {
      res <- list(nam, `names<-`(X, NULL))
      names(res) <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", l1orlst(as.character(substitute(X))))
    }
    attr(res, "row.names") <- .set_row_names(length(X))
    return(`oldClass<-`(res, if(length(class)) class else c("data.table", "data.frame")))
  }
  if(keep.attr) {
    # if(all(class(X) == class)) return(X) # better adjust rows ? -> yes, row.names.col should always work !
    if(is.null(attr(X, "names"))) attr(X, "names") <- paste0("V", seq_along(unclass(X)))
    if(!isFALSE(row.names.col) && length(rn <- attr(X, "row.names"))) {
      ax <- attributes(X)
      X <- c(list(rn), X)
      ax[["names"]] <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", ax[["names"]])
      setattributes(X, ax)
    }
    if(!length(class) && inherits(X, c("data.table", "data.frame"))) return(X)
    attr(X, "row.names") <- .set_row_names(length(.subset2(X, 1L)))
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
    attr(X, "row.names") <- .set_row_names(length(X[[1L]]))
  }
  oldClass(X) <- if(length(class)) class else c("data.table", "data.frame")
  X
}

qTBL <- function(X, row.names.col = FALSE, keep.attr = FALSE)
  qDT(X, row.names.col, keep.attr, c("tbl_df", "tbl", "data.frame"))


qM <- function(X, keep.attr = FALSE, class = NULL) {
  if(keep.attr) {
    if(is.atomic(X)) {
      if(length(class)) oldClass(X) <- class
      if(is.matrix(X)) return(X)
      if(is.array(X)) {
        d <- dim(X)
        dn <- dimnames(X)
        dim(X) <- c(d[1L], prod(d[-1L]))
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
    res <- do.call(cbind, X)
    rn <- ax[["row.names"]]
    if(!(is.numeric(rn) || is.null(rn) || rn[1L] == "1"))
       dimnames(res) <- list(rn, ax[["names"]])
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
      dim(X) <- c(d[1L], prod(d[-1L]))
      if(length(dn)) {
        for (i in 2L:ld) if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
        dimnames(X) <- list(dn[[1L]], interact_names(dn[-1L])) # Good?
      }
    }
    if(length(class)) oldClass(X) <- class
    return(X)
  }
  rn <- attr(X, "row.names")
  res <- do.call(cbind, X)
  if(is.object(res)) attributes(res) <- attributes(res)[c("dim", "dimnames")] # if X is list of time-series, do.call(cbind, X) creates ts-matrix.
  if(!(is.numeric(rn) || is.null(rn) || rn[1L] == "1"))
    dimnames(res) <- list(rn, attr(X, "names"))
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


## Old Versions:

# # before collapse 1.4.0:
# qM_old <- function(X) {
#   if(is.atomic(X)) {
#     d <- dim(X)
#     ld <- length(d)
#     if(ld > 2L) {
#       dn <- dimnames(X)
#       dim(X) <- c(d[1L], prod(d[-1L]))
#       if(length(dn)) {
#         for (i in 2L:ld) if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
#         dimnames(X) <- list(dn[[1L]], interact_names(dn[-1L])) # Good?
#       }
#       return(X)
#     }
#     if(ld == 2L) return(X)
#     return(matrix(X, ncol = 1, dimnames = list(names(X), l1orlst(as.character(substitute(X))))))
#   }
#   rn <- attr(X, "row.names")
#   res <- do.call(cbind, X)
#   if(!(is.numeric(rn) || is.null(rn) || rn[1L] == "1")) dimnames(res) <- list(rn, attr(X, "names"))
#   res
# }


# # before collapse 1.4.0:
# qDF <- function(X, row.names.col = FALSE) {
#   if(is.atomic(X)) {
#     d <- dim(X)
#     ld <- length(d)
#     if(ld >= 2L) {
#       if(ld != 2L) {
#         dn <- dimnames(X)
#         dim(X) <- c(d[1L], prod(d[-1L]))
#         if(length(dn)) {
#           for (i in 2L:ld) if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
#           dimnames(X) <- list(dn[[1L]], interact_names(dn[-1L])) # Good?
#         }
#       }
#       if(!isFALSE(row.names.col) && length(force(dn <- dimnames(X))[[1L]])) {
#         res <- c(list(dn[[1L]]), .Call(Cpp_mctl, X, FALSE, 0L))
#         names(res) <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", dn[[2L]])
#         attr(res, "row.names") <- .set_row_names(length(dn[[1L]]))
#         return(`oldClass<-`(res, "data.frame"))
#       }
#       return(.Call(Cpp_mctl, X, TRUE, 1L))
#     }
#     nam <- names(X)
#     if(isFALSE(row.names.col) || is.null(nam)) {
#       if(is.null(nam)) {
#         res <- `names<-`(list(X), l1orlst(as.character(substitute(X))))
#         attr(res, "row.names") <- .set_row_names(length(X))
#       } else {
#         res <- `names<-`(list(`names<-`(X, NULL)), l1orlst(as.character(substitute(X))))
#         attr(res, "row.names") <- nam
#       }
#     } else {
#       res <- list(nam, `names<-`(X, NULL))
#       names(res) <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", l1orlst(as.character(substitute(X))))
#       attr(res, "row.names") <- .set_row_names(length(X))
#     }
#     return(`oldClass<-`(res, "data.frame"))
#   }
#   # if(inherits(X, "data.frame")) return(X)
#   if(is.null(attr(X, "names"))) attr(X, "names") <- paste0("V", seq_along(unclass(X)))
#   if(is.null(attr(X, "row.names"))) {
#     attr(X, "row.names") <- .set_row_names(fnrow2(X))
#   } else if(!isFALSE(row.names.col)) {
#     ax <- attributes(X)
#     X <- c(list(ax[["row.names"]]), X) # best ??
#     ax[["row.names"]] <- .set_row_names(length(X[[1L]])) # this is ok, X is a list ...
#     ax[["names"]] <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", ax[["names"]])
#     setattributes(X, ax)
#   }
#   oldClass(X) <- "data.frame"
#   X
# }


# # before collapse 1.4.0:
# qDT <- function(X, row.names.col = FALSE) {
#   if(is.atomic(X)) {
#     d <- dim(X)
#     ld <- length(d)
#     if(ld >= 2L) {
#       if(ld != 2L) {
#         dn <- dimnames(X)
#         dim(X) <- c(d[1L], prod(d[-1L]))
#         if(length(dn)) {
#           for (i in 2L:ld) if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
#           dimnames(X) <- list(dn[[1L]], interact_names(dn[-1L])) # Good?
#         }
#       }
#       if(!isFALSE(row.names.col) && length(force(dn <- dimnames(X))[[1L]])) {
#         res <- c(list(dn[[1L]]), .Call(Cpp_mctl, X, FALSE, 0L))
#         names(res) <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", dn[[2L]])
#         attr(res, "row.names") <- .set_row_names(length(dn[[1L]]))
#         return(`oldClass<-`(res, c("data.table","data.frame")))
#       }
#       return(.Call(Cpp_mctl, X, TRUE, 2L))
#     }
#     if(isFALSE(row.names.col) || is.null(nam <- names(X))) {
#       res <- `names<-`(list(X), l1orlst(as.character(substitute(X))))
#     } else {
#       res <- list(nam, `names<-`(X, NULL))
#       names(res) <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", l1orlst(as.character(substitute(X))))
#     }
#     attr(res, "row.names") <- .set_row_names(length(X))
#     return(`oldClass<-`(res, c("data.table","data.frame")))
#   }
#   if(inherits(X, "data.table")) return(X)
#   if(is.null(attr(X, "names"))) attr(X, "names") <- paste0("V", seq_along(unclass(X)))
#   if(!isFALSE(row.names.col) && length(rn <- attr(X, "row.names"))) {
#     ax <- attributes(X)
#     X <- c(list(rn), X)
#     ax[["names"]] <- c(if(is.character(row.names.col)) row.names.col[1L] else "row.names", ax[["names"]])
#     setattributes(X, ax)
#   }
#   attr(X, "row.names") <- .set_row_names(fnrow2(X))
#   oldClass(X) <- c("data.table","data.frame")
#   X
# }


