library(Rcpp)
sourceCpp("R/C++/mrtl_type_dispatch_final.cpp")
sourceCpp("R/C++/qFqG.cpp", rebuild = TRUE) # https://gallery.rcpp.org/articles/fast-factor-generation/
qF <- function(x, ordered = TRUE) {
  if(is.factor(x)) return(x)
  qFCpp(x, ordered)
}
qG <- function(x, ordered = TRUE) {
  if(is.factor(x)) return(x)
  qGCpp(x, ordered)
}
qDF <- function(X) { 
  if(is.atomic(X)) {
    d <- dim(X)
    ld <- length(d)
    if(ld == 2L)
      mctl(X, names = TRUE, ret = 1L) else if (ld > 2L) {
        dn <- dimnames(X)
        dim(X) <- c(d[1L], prod(d[-1L]))
        if (!is.null(dn)) {
          for (i in 2L:ld) if (is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
          dn[[2]] <- interaction(expand.grid(dn[-1L])) # Good??
        }
        mctl(X, names = TRUE, ret = 1L) 
      } else {
        lx <- length(X)
        X <- setNames(list(X), deparse(substitute(X)))
        attr(X,"row.names") <- .set_row_names(lx)
        class(X) <- "data.frame"
        X
      }
  } else {
    if(inherits(X, "data.frame")) return(X)
    if (is.null(attr(X,"names"))) attr(X,"names") <- paste0("V", seq_along(X))
    if (is.null(attr(X,"row.names"))) attr(X,"row.names") <- .set_row_names(length(X[[1]]))
    class(X) <- "data.frame"
    X
  }
}
qDT <- function(X) { # what if already DT ?? return ??  
  if(is.atomic(X)) {
    d <- dim(X)
    ld <- length(d)
    if(ld == 2L) 
      mctl(X, names = TRUE, ret = 2L) else if (ld > 2L) {
        dn <- dimnames(X)
        dim(X) <- c(d[1L], prod(d[-1L]))
        if (!is.null(dn)) {
          for (i in 2L:ld) if (is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
          dn[[2]] <- interaction(expand.grid(dn[-1L])) # Good??
        }
        mctl(X, names = TRUE, ret = 2L)
      } else {
        lx <- length(X)
        X <- setNames(list(X), deparse(substitute(X)))
        attr(X,"row.names") <- .set_row_names(lx)
        class(X) <- c("data.table","data.frame")
        X
      }
  } else {
    if(inherits(X, "data.table")) return(X)
    if (is.null(attr(X, "names"))) attr(X, "names") <- paste0("V",seq_along(X))
    attr(X, "row.names") <- .set_row_names(length(X[[1]]))
    class(X) <- c("data.table","data.frame")
    X
  } 
}
qM <- function(X) {
  if (is.atomic(X)) {
    d <- dim(X)
    ld <- length(d)
    if (ld > 2L) {
      dn <- dimnames(X)
      dim(X) <- c(d[1L], prod(d[-1L]))
      if (!is.null(dn)) {
        if (length(dn[[1L]])) rownames(X) <- dn[[1L]]
        for (i in 2L:ld) if (is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
        colnames(X) <- interaction(expand.grid(dn[-1L])) # Good??
      }
      X
    } else if (ld == 2L) X else 
      matrix(X, ncol = 1, dimnames = list(NULL, deparse(substitute(X))))
  } else {
    rn <- attr(X, "row.names")
    # # if (is.null(cn)) cn <- paste0("C",seq_len(lx)) # Not really necessary
    res <- do.call(cbind, X) # faster !1
    if(!(is.null(rn) || is.integer(rn))) dimnames(res) <- list(rn, names(X))
    res
    # rn <- attr(X, "row.names") # other solution: slightly slower !!
    # res <- unlist(X, use.names = FALSE, recursive = FALSE)
    # dim(res) <- dim.data.frame(X)
    # dimnames(res) <- if(rn[1] == "1") list(NULL, names(X)) else list(rn, names(X))
    # res
  }
}


# # Quick conversions to DF, DT and M
# setdndf <- function(x, value) {
#   if (is.null(value)) {
#     value <- list(.set_row_names(length(x[[1]])), paste0("V",seq_along(x)))
#   } else {
#     if(is.null(value[[1]])) value[[1]] <- .set_row_names(length(x[[1]]))  
#     if(is.null(value[[2]])) value[[2]] <- paste0("V",seq_along(x))
#   }
#   attr(x, "names") <- value[[2]]
#   attr(x, "row.names") <- value[[1]]
#   x
# }
# setdndt <- function(x, value) {
#   attr(x, "names") <- if(is.null(value[[2]])) paste0("V",seq_along(x)) else value[[2]]
#   attr(x, "row.names") <- .set_row_names(length(x[[1]]))
#   x
# }
# qDFOld <- function(X) { 
#   if (is.atomic(X)) {
#     d <- dim(X)
#     ld <- length(d)
#     if (ld == 2L) {
#       X <- setdndf(mctl(X), dimnames(X)) 
#     } else if (ld > 2L) {
#       dn <- dimnames(X)
#       dim(X) <- c(d[1L], prod(d[-1L]))
#       if (!is.null(dn)) {
#         for (i in 2L:ld) if (is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
#         dn[[2]] <- interaction(expand.grid(dn[-1L])) # Good??
#       }
#       X <- setdndf(mctl(X), dn) # dimnames option for mctl??
#     } else {
#       lx <- length(x)
#       X <- setn(list(X), deparse(substitute(X)))
#       attr(X,"row.names") <- .set_row_names(lx)
#     }
#   } else {
#    if (is.null(attr(X,"names"))) attr(X,"names") <- paste0("V", seq_along(X))
#    if (is.null(attr(X,"row.names"))) attr(X,"row.names") <- .set_row_names(length(X[[1]]))
#   }
#   class(X) <- "data.frame" # skip if X is already a data.frame??
#   # data.table::setattr(X, "class", "data.frame")
#   # attr(X, "row.names") <- .set_row_names(length(X[[1]])) # use setattr ?? -> not faster!!
#   # data.table::setattr(X, "row.names", .set_row_names(length(X[[1]])))
#   X
# }
# qDTOld <- function(X) { 
#   if (is.atomic(X)) {
#     d <- dim(X)
#     ld <- length(d)
#     if (ld == 2L) {
#       X <- setdndt(mctl(X), dimnames(X)) 
#       #cn <- dimnames(X)[[2]]
#       #if(is.null(cn)) cn <- paste0("V",seq_len(d[2]))
#       #X <- setn(mctl(X), cn) # build your own faster setNames with attr(x, "names) ??
#     } else if (ld > 2L) {
#       dn <- dimnames(X)
#       dim(X) <- c(d[1L], prod(d[-1L]))
#       if (!is.null(dn)) {
#         for (i in 2L:ld) if (is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
#         dn[[2]] <- interaction(expand.grid(dn[-1L])) # Good??
#       }
#       X <- setdndt(mctl(X), dn) # X <- setn(mctl(X), dn[[2]])
#     } else X <- setn(list(X), deparse(substitute(X)))
#   } else {
#     if (is.null(attr(X, "names"))) attr(X, "names") <- paste0("V",seq_along(X))
#     attr(X, "row.names") <- .set_row_names(length(X[[1]]))
#   } 
#   class(X) <- c("data.table","data.frame") # skip if X is already a data.table??
#   #attr(X, "row.names") <- .set_row_names(length(X[[1]])) # use setattr ?? -> nope, slower
#   X
# }

# Final Versions Were Here !!!!!!!!

# qF2 <- function(x, ordered = FALSE) { # not so great, unless speed up both unique and match !!
#   if(is.factor(x)) return(x)
#   if(is.numeric(x)) {
#     o = .Call(data.table:::Cforder, x, NULL, TRUE, !ordered, 1L, FALSE) # good??
#     s = attr(o, "starts")
#     f = .Call(data.table:::Cfrank, o, s, .Call(data.table:::Cuniqlengths, s, length(o)), "dense")
#     attr(f, "levels") <- as.character(x[o[s]]) # // This conversion takes a lot of time !!!
#   } else {
#     un <- if(ordered) .Internal(unique(x, FALSE, FALSE, length(x))) else sort.default(.Internal(unique(x, FALSE, FALSE, length(x))), method = "quick")
#     f <- if (is.character(x)) .Call(data.table:::Cchmatch, x, un, NA_integer_) else match(x, un)
#     attr(f, "levels") <- as.character(un) # This conversion takes a lot of time !!!
#   }
#   attr(f, "class") <- if(ordered) c("ordered","factor") else "factor"
#   f
# }


# microbenchmark(qF(WDI$countrycode),qF2(WDI$countrycode),factor(WDI$countrycode),as.factor(WDI$countrycode))
# microbenchmark(qF(x),qF2(x),factor(x),as.factor(x), times = 1)

# to do: mctl with dimnames option?? faster??
# qM if numeric row.names discard!!
