
fsplit <- function(x, f, drop, ...) if(drop && is.factor(f))
  split(x, .Call(Cpp_fdroplevels, f, !inherits(f, "na.included")), drop = FALSE, ...) else
    split(x, qF(f), drop = FALSE, ...)

t_list2 <- function(x) .Call(Cpp_mctl, do.call(rbind, x), TRUE, 0L)

# This is for export
t_list <- function(l) {
  lmat <- do.call(rbind, l)
  rn <- dimnames(lmat)[[1L]]
  .Call(C_copyMostAttrib,
        lapply(.Call(Cpp_mctl, lmat, TRUE, 0L), `names<-`, rn), l)
}


# fsplit_DF2 <- function(x, f, drop, ...) {
#
#   if(!length(x)) return(x)
#
#   j <- seq_along(unclass(x))
#
#   ind <- if(drop && is.factor(f))
#     split.default(seq_along(.subset2(x, 1L)), .Call(Cpp_fdroplevels, f, !inherits(f, "na.included")), drop = FALSE, ...) else
#       split.default(seq_along(.subset2(x, 1L)), qF(f), drop = FALSE, ...)
#
#   lapply(ind, function(i) .Call(C_subsetDT, x, i, j)) # tryCatch(, error = function(e) print(list(x = x, ind = ind, j = j)))
# }

# rsplit_default2(mtcars$mpg, slt(mtcars, cyl, vs)) -> Still find source of Bug !! see if it is faster !!

# rsplit_default2 <- function(x, fl, drop = TRUE, flatten = FALSE, ...) { # , check = TRUE
#   if(is.atomic(fl)) return(fsplit(x, fl, drop, ...))
#   if(flatten && length(unclass(fl)) > 1L) return(fsplit(x, finteraction(fl), drop, ...))
#   attributes(fl) <- NULL
#   # if(check) fl <- lapply(fl, qF) # necessary ? -> split.default is actually faster on non-factor variables !
#   rspl <- function(y, fly) {
#     if(length(fly) == 1L) return(fsplit(y, fly[[1L]], drop, ...))
#     mapply(rspl, y = fsplit(y, fly[[1L]], drop, ...),
#            fly = fsplit_DF2(fly[-1L], fly[[1L]], drop, ...), SIMPLIFY = FALSE) # Possibility to avoid transpose ? C_subsetDT ??
#   }
#   rspl(x, fl)
# }


rsplit <- function(x, ...) UseMethod("rsplit")

rsplit.default <- function(x, fl, drop = TRUE, flatten = FALSE, ...) { # , check = TRUE
  if(is.atomic(fl)) return(fsplit(x, fl, drop, ...))
  if(flatten && length(unclass(fl)) > 1L) return(fsplit(x, finteraction(fl), drop, ...))
  attributes(fl) <- NULL
  # if(check) fl <- lapply(fl, qF) # necessary ? -> split.default is actually faster on non-factor variables !
  rspl <- function(y, fly) {
    if(length(fly) == 1L) return(fsplit(y, fly[[1L]], drop, ...))
    mapply(rspl, y = fsplit(y, fly[[1L]], drop, ...),
           fly = t_list2(lapply(fly[-1L], fsplit, fly[[1L]], drop, ...)), SIMPLIFY = FALSE) # Possibility to avoid transpose ? C_subsetDT ??
  }
  rspl(x, fl)
}



# From stackoverflow package:
# rsplit <- function (x, by, drop = FALSE)
# {
#   if (is.atomic(by))
#     return(split(x, by, drop = drop))
#   attributes(by) <- NULL
#   if (length(by) == 1L)
#     return(split(x, by[[1L]], drop = drop))
#   mapply(rsplit, x = split(x, by[[1L]], drop = drop), by = t(lapply(by[-1L], split, by[[1L]], drop = drop)), drop = drop,
#          SIMPLIFY = FALSE)
# }

rsplit.data.frame <- function(x, by, drop = TRUE, flatten = FALSE, # check = TRUE,
                              cols = NULL, keep.by = FALSE, simplify = TRUE, ...) {

  if(is.call(by)) {
    nam <- attr(x, "names")
    if(length(by) == 3L) {
      byn <- ckmatch(all.vars(by[[3L]]), nam)
      cols <- ckmatch(all.vars(by[[2L]]), nam)
    } else { # keep.by always added: Same behavior as L or W !!
      byn <- ckmatch(all.vars(by), nam)
      if(!(is.null(cols) && keep.by))
        cols <- if(is.null(cols)) -byn else cols2int(cols, x, nam, FALSE)
    }
    by <- .subset(x, byn)
    if(length(cols)) x <- fcolsubset(x, if(keep.by) c(byn, cols) else cols)
  } else if(length(cols))
    x <- fcolsubset(x, cols2int(cols, x, attr(x, "names"), FALSE))

  if(simplify && length(unclass(x)) == 1L)
    return(rsplit.default(.subset2(x, 1L), by, drop, flatten, ...))  # , check
  # Note there is a data.table method: split.data.table, which can also do recursive splitting..

  j <- seq_along(unclass(x))
  rn <- attr(x, "row.names")
  if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") {
    fsplit_DF <- function(x, f, ...)
      lapply(fsplit(seq_along(.subset2(x, 1L)), f, drop, ...),
             function(i) .Call(C_subsetDT, x, i, j))
  } else {
    fsplit_DF <- function(x, f, ...)
      lapply(fsplit(seq_along(.subset2(x, 1L)), f, drop, ...),
             function(i) `attr<-`(.Call(C_subsetDT, x, i, j), "row.names", attr(x, "row.names")[i]))
  }

  if(is.atomic(by)) return(fsplit_DF(x, by, ...))
  if(flatten && length(unclass(by)) > 1L) return(fsplit_DF(x, finteraction(by), ...))
  attributes(by) <- NULL
  # if(check) by <- lapply(by, qF) # necessary ?
  rspl_DF <- function(y, fly) {
    if(length(fly) == 1L) return(fsplit_DF(y, fly[[1L]], ...))
    mapply(rspl_DF, y = fsplit_DF(y, fly[[1L]], ...),
           fly = t_list2(lapply(fly[-1L], fsplit, fly[[1L]], drop, ...)), SIMPLIFY = FALSE) # Possibility to avoid transpose ?
  }                # use C_subsetDT here as well ??? what is faster ???
  rspl_DF(x, by)
}


# Misc trial:
#
#     fli <- lapply(fl[-1L], split, fl[[1L]], ...)
#     for(i in seq_len(length(fl)-1L)) {
#       r <- rapply(r, split, fl[[1L]], how = "list")
#       fli <- rapply(fli[-1L], split, fli[[2L]], ..., how = "list")
#     }
#   }
#   r
