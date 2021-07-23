# Row-operations (documented under data transformations...) ... see if any other package has it (i.e. matrixStats etc..)
# or wirhout r ??? look for %+% function on Rducumentation.. rdio.

"%rr%" <- function(X, v) if(is.atomic(X) || is.atomic(v) || inherits(X, "data.frame")) TRA(X, v, "replace_fill") else # outer(rep.int(1L, dim(X)[2L]), v)
  duplAttributes(mapply(function(x, y) TRA(x, y, "replace_fill"), unattrib(X), unattrib(v),
                        USE.NAMES = FALSE, SIMPLIFY = FALSE), X)
"%r+%" <- function(X, v) if(is.atomic(X) || is.atomic(v) || inherits(X, "data.frame")) TRA(X, v, "+") else
  duplAttributes(mapply(function(x, y) TRA(x, y, "+"), unattrib(X), unattrib(v),
                        USE.NAMES = FALSE, SIMPLIFY = FALSE), X)
"%r-%" <- function(X, v) if(is.atomic(X) || is.atomic(v) || inherits(X, "data.frame")) TRA(X, v, "-") else
  duplAttributes(mapply(function(x, y) TRA(x, y, "-"), unattrib(X), unattrib(v),
                        USE.NAMES = FALSE, SIMPLIFY = FALSE), X)
"%r*%" <- function(X, v) if(is.atomic(X) || is.atomic(v) || inherits(X, "data.frame")) TRA(X, v, "*") else
  duplAttributes(mapply(function(x, y) TRA(x, y, "*"), unattrib(X), unattrib(v),
                        USE.NAMES = FALSE, SIMPLIFY = FALSE), X)
"%r/%" <- function(X, v) if(is.atomic(X) || is.atomic(v) || inherits(X, "data.frame")) TRA(X, v, "/") else
  duplAttributes(mapply(function(x, y) TRA(x, y, "/"), unattrib(X), unattrib(v),
                        USE.NAMES = FALSE, SIMPLIFY = FALSE), X)



# othidentity <- function(x, y) y
"%cr%" <- function(X, v) if(is.atomic(X)) return(duplAttributes(rep(v, NCOL(X)), X)) else # outer(rep.int(1L, dim(X)[2L]), v)
  if(is.atomic(v)) return(duplAttributes(lapply(vector("list", length(unclass(X))), function(z) v), X)) else
    copyAttrib(v, X) # copyAttrib first makes a shallow copy of v
"%c+%" <- function(X, v) if(is.atomic(X)) return(X + v) else
  duplAttributes(if(is.atomic(v)) lapply(unattrib(X), `+`, v) else
    mapply(`+`, unattrib(X), unattrib(v), USE.NAMES = FALSE, SIMPLIFY = FALSE), X)
"%c-%" <- function(X, v) if(is.atomic(X)) return(X - v) else
  duplAttributes(if(is.atomic(v)) lapply(unattrib(X), `-`, v) else
    mapply(`-`, unattrib(X), unattrib(v), USE.NAMES = FALSE, SIMPLIFY = FALSE), X)
"%c*%" <- function(X, v) if(is.atomic(X)) return(X * v) else
  duplAttributes(if(is.atomic(v)) lapply(unattrib(X), `*`, v) else
    mapply(`*`, unattrib(X), unattrib(v), USE.NAMES = FALSE, SIMPLIFY = FALSE), X)
"%c/%" <- function(X, v) if(is.atomic(X)) return(X / v) else  # or * 1L/v ??
  duplAttributes(if(is.atomic(v)) lapply(unattrib(X), `/`, v) else
    mapply(`/`, unattrib(X), unattrib(v), USE.NAMES = FALSE, SIMPLIFY = FALSE), X)



getenvFUN <- function(nam, efmt1 = "For this method need to install.packages('%s'), then unload [detach('package:collapse', unload = TRUE)] and reload [library(collapse)].")
  if(is.null(FUN <- .collapse_env[[nam]])) stop(sprintf(efmt1, strsplit(nam, "_", fixed = TRUE)[[1L]][1L])) else FUN


# qM2 <- function(x) if(is.list(x)) do.call(cbind, x) else x

null2NA <- function(x) if(is.null(x)) NA_character_ else x

# flapply <- function(x, FUN, ...) lapply(unattrib(x), FUN, ...) # not really needed ...

vlabels <- function(X, attrn = "label") {
  if(is.atomic(X)) return(null2NA(attr(X, attrn)))
  res <- lapply(X, attr, attrn) # unattrib(X): no names
  res[vapply(res, is.null, TRUE)] <- NA_character_
  unlist(res)
}

# Slower on WDI !!!
# vlabels2 <- function(X, attrn = "label") {
#   if(is.atomic(X)) return(null2NA(attr(X, attrn)))
#   vapply(X, function(x) if(is.null(a <- attr(x, attrn))) NA_character_ else a, character(1L))
# }

"vlabels<-" <- function(X, attrn = "label", value) {
  names(value) <- NULL
  if(is.atomic(X)) return(`attr<-`(X, attrn, value))
  clx <- oldClass(X)
  oldClass(X) <- NULL
  if(is.null(value)) {
    for (i in seq_along(X)) attr(X[[i]], attrn) <- NULL
  } else {
    if(length(X) != length(value)) stop("length(X) must match length(value)")
    for (i in seq_along(value)) attr(X[[i]], attrn) <- value[[i]]
  }
  if(any(clx == "data.table")) return(alc(`oldClass<-`(X, clx)))
  `oldClass<-`(X, clx)
}

setLabels <- function(X, value, attrn = "label") `vlabels<-`(X, attrn, value)

# Also slower on WDI !!
# "vlabels2<-" <- function(X, attrn = "label", value) {
#   names(value) <- NULL
#   if(is.atomic(X)) return(`attr<-`(X, attrn, value))
#   duplAttributes(mapply(function(x, y) `attr<-`(x, attrn, y), `attributes<-`(X, NULL), as.vector(value, "list"),
#                         SIMPLIFY = FALSE, USE.NAMES = FALSE), X)
# }

.c <- function(...) as.character(substitute(c(...))[-1L])


strclp <- function(x) if(length(x) > 1L) paste(x, collapse = " ") else x

pasteclass <- function(x) if(length(cx <- class(x)) > 1L) paste(cx, collapse = " ") else cx

vclasses <- function(X) {
  if(is.atomic(X)) return(pasteclass(X))
  vapply(X, pasteclass, character(1L)) # unattrib(X): no names
}

vtypes <- function(X) {
  if(is.atomic(X)) return(typeof(X))
  vapply(X, typeof, character(1L)) # unattrib(X): no names
}

namlab <- function(X, class = FALSE, attrn = "label") {
  if(!is.list(X)) stop("namlab only works with lists")
  res <- if(class) list(attr(X, "names"), vapply(unattrib(X), pasteclass, character(1)), vlabels(X, attrn)) else
                   list(attr(X, "names"), vlabels(X, attrn)) # could do unattrib(vlabels(X, attrn)) ??
  attributes(res) <- list(names = if(class) c("Variable","Class","Label") else c("Variable","Label"),
                          row.names = .set_row_names(length(unclass(X))),
                          class = "data.frame")
  res
}

add_stub <- function(X, stub, pre = TRUE) {
  if(!is.character(stub)) return(X)
  if(is.atomic(X) && is.array(X)) {
    if(length(dim(X)) > 2L) stop("Can't stub higher dimensional arrays!")
    dn <- dimnames(X)
    cn <- dn[[2L]]
    if(length(cn)) dimnames(X) <- list(dn[[1L]], if(pre) paste0(stub, cn) else paste0(cn, stub))
  } else {
    nam <- attr(X, "names")
    if(length(nam)) {
      attr(X, "names") <- if(pre) paste0(stub, nam) else paste0(nam, stub)
      if(inherits(X, "data.table")) X <- alc(X)
    }
  }
  X
}

rm_stub <- function(X, stub, pre = TRUE, regex = FALSE, ...) {
  if(!is.character(stub)) return(X)
  if(regex)
    rmstubFUN <- function(x) {
      gsub(stub, "", x, ...)
    } else if(pre)
    rmstubFUN <- function(x) { # much faster than using sub!
      v <- startsWith(x, stub)
      x[v] <- substr(x[v], nchar(stub)+1L, 1000000L)
      x
    } else
    rmstubFUN <- function(x) { # much faster than using sub!
      v <- endsWith(x, stub)
      xv <- x[v] # faster ..
      x[v] <- substr(xv, 0L, nchar(xv)-nchar(stub))
      x
    }
  if(is.atomic(X)) {
    d <- dim(X)
    if(is.null(d)) if(is.character(X)) return(rmstubFUN(X)) else stop("Cannot modify a vector that is not character")
    if(length(d) > 2L) stop("Can't remove stub from higher dimensional arrays!")
    dn <- dimnames(X)
    dimnames(X) <- list(dn[[1L]], rmstubFUN(dn[[2L]]))
  } else {
    attr(X, "names") <- rmstubFUN(attr(X, "names"))
    if(inherits(X, "data.table")) X <- alc(X)
  }
  X
}

setRownames <- function(object, nm = if(is.atomic(object)) seq_row(object) else NULL) {
  if(is.list(object)) {
    l <- length(.subset2(object, 1L))
    if(is.null(nm)) nm <- .set_row_names(l) else if(length(nm) != l) stop("supplied row-names must match list extent")
    attr(object, "row.names") <- nm
    if(inherits(object, "data.table")) return(alc(object))
    return(object)
  }
  if(!is.array(object)) stop("Setting row-names only supported on arrays and lists")
  dn <- dimnames(object)
 `dimnames<-`(object, c(list(nm), dn[-1L]))
}

setColnames <- function(object, nm) {
  if(is.atomic(object) && is.array(object))
    dimnames(object)[[2L]] <- nm else {
    attr(object, "names") <- nm
    if(inherits(object, "data.table")) return(alc(object))
  }
  object
}

setDimnames <- function(object, dn, which = NULL) {
  if(is.null(which)) return(`dimnames<-`(object, dn))
  if(is.atomic(dn)) dimnames(object)[[which]] <- dn else
                    dimnames(object)[which] <- dn
  object
}

all_identical <- function(...) {
  if(...length() == 1L && is.list(...)) return(all(vapply(unattrib(...)[-1L], identical, TRUE, .subset2(..., 1L))))
  l <- list(...)
  all(vapply(l[-1L], identical, TRUE, l[[1L]]))
}

all_obj_equal <- function(...) {
  if(...length() == 1L && is.list(...))
    r <- unlist(lapply(unattrib(...)[-1L], all.equal, .subset2(..., 1L)), use.names = FALSE) else {
    l <- list(...)
    r <- unlist(lapply(l[-1L], all.equal, l[[1L]]), use.names = FALSE)
  }
  is.logical(r)
}

cinv <- function(X) chol2inv(chol(X))

interact_names <- function(l) do.call(paste, c(expand.grid(l, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE), list(sep = ".")))

# set over-allocation for data.table's
alc <- function(x) .Call(C_alloccol, x)
condalc <- function(x, DT) if(DT) .Call(C_alloccol, x) else x
alcSA <- function(x, a) .Call(C_alloccol, .Call(C_setAttributes, x, a))
condalcSA <- function(x, a, DT) if(DT) .Call(C_alloccol, .Call(C_setAttributes, x, a)) else .Call(C_setAttributes, x, a)

unattrib <- function(object) `attributes<-`(object, NULL)

# Both equally efficient and therefore redundant !
# setAttr <- function(object, a, v) .Call(C_setAttr, object, a, v)
# setAttrR <- function(object, a, v) `attr<-`(object, a, v)

setAttrib <- function(object, a) .Call(C_setAttrib, object, a)
# setAttribR <- function(object, a) `attributes<-`(object, x)

copyAttrib <- function(to, from) .Call(C_copyAttrib, to, from)
# copyAttribR <- function(to, from) `attributes<-`(to, attributes(from))


copyMostAttrib <- function(to, from) .Call(C_copyMostAttrib, to, from)
# copyMostAttribR <- function(to, from) `mostattributes<-`(to, attributes(from))

addAttributes <- function(x, a) .Call(C_setAttributes, x, c(attributes(x), a))


is_categorical <- function(x) !is.numeric(x)
is.categorical <- is_categorical

is_date <- function(x) inherits(x, c("Date","POSIXlt","POSIXct"))
is.Date <- is_date

"%!in%" <- function(x, table) match(x, table, nomatch = 0L) == 0L

# more consistent with base than na_rm
# na.rm <- function(x) { # cpp version available, but not faster !
#   if(length(attr(x, "names"))) { # gives corruped time-series !
#     ax <- attributes(x)
#     r <- x[!is.na(x)]
#     ax[["names"]] <- names(r)
#     setAttributes(r, ax)
#   } else duplAttributes(x[!is.na(x)], x)
# }

alloc <- function(value, n) .Call(C_alloc, value, n)

allNA <- function(x) .Call(C_allNA, x, TRUE) # True means give error for unsupported vector types, not FALSE.

missing_cases <- function(X, cols = NULL) {
  if(is.list(X)) return(.Call(C_dt_na, X, if(is.null(cols)) seq_along(unclass(X)) else cols2int(cols, X, attr(X, "names"))))
  if(is.matrix(X)) return(if(is.null(cols)) !complete.cases(X) else !complete.cases(X[, cols]))
  is.na(X)
}

na_rm <- function(x) .Call(C_na_rm, x)  # x[!is.na(x)]

na_omit <- function(X, cols = NULL, na.attr = FALSE) {
  if(is.list(X)) {
    iX <- seq_along(unclass(X))
    rl <- if(is.null(cols)) !.Call(C_dt_na, X, iX) else
      !.Call(C_dt_na, X, cols2int(cols, X, attr(X, "names"))) # gives error if X not list
    rkeep <- which(rl)
    if(length(rkeep) == fnrow2(X)) return(condalc(X, inherits(X, "data.table")))
    res <- .Call(C_subsetDT, X, rkeep, iX, FALSE)
    rn <- attr(X, "row.names")
    if(!(is.numeric(rn) || is.null(rn) || rn[1L] == "1")) attr(res, "row.names") <- rn[rkeep]
    if(na.attr) {
      attr(res, "na.action") <- `oldClass<-`(which(!rl), "omit")
      if(inherits(res, "data.table")) return(alc(res))
    }
  } else {
    rl <- if(is.null(cols)) complete.cases(X) else complete.cases(X[, cols])
    rkeep <- which(rl)
    if(length(rkeep) == NROW(X)) return(X)
    res <- if(is.matrix(X)) X[rkeep, , drop = FALSE] else X[rkeep]
    if(na.attr) attr(res, "na.action") <- `oldClass<-`(which(!rl), "omit")
  }
  res
}

na_insert <- function(X, prop = 0.1) {
  if(is.list(X)) {
    n <- fnrow2(X)
    nmiss <- floor(n * prop)
    res <- duplAttributes(lapply(unattrib(X), function(y) `[<-`(y, sample.int(n, nmiss), value = NA)), X)
    return(if(inherits(X, "data.table")) alc(res) else res)
  }
  if(!is.atomic(X)) stop("X must be an atomic vector, array or data.frame")
  l <- length(X)
  X[sample.int(l, floor(l * prop))] <- NA
  X
}

fdapply <- function(X, FUN, ...) duplAttributes(lapply(`attributes<-`(X, NULL), FUN, ...), X)

fnlevels <- function(x) length(attr(x, "levels"))

# flevels <- function(x) attr(x, "levels")

fnrow <- function(X) if(is.list(X)) length(.subset2(X, 1L)) else dim(X)[1L]

fnrow2 <- function(X) length(.subset2(X, 1L))

fncol <- function(X) if(is.list(X)) length(unclass(X)) else dim(X)[2L]

fNCOL <- function(X) if(is.list(X)) length(unclass(X)) else NCOL(X)

fdim <- function(X) {
   if(is.atomic(X)) return(dim(X)) # or if !is.list ?
   oldClass(X) <- NULL
   c(length(X[[1L]]), length(X)) # Faster than c(length(.subset2(X, 1L)), length(unclass(X)))
}

seq_row <- function(X) if(is.list(X)) seq_along(.subset2(X, 1L)) else seq_len(nrow(X))

seq_col <- function(X) if(is.list(X)) seq_along(unclass(X)) else seq_len(ncol(X))

# na.last is false !!
forder.int <- function(x) .Call(C_radixsort, FALSE, FALSE, FALSE, FALSE, TRUE, pairlist(x)) # if(is.unsorted(x)) .Call(C_forder, x, NULL, FALSE, TRUE, 1L, TRUE) else seq_along(x) # since forder gives integer(0) if sorted !
ford <- function(x, g = NULL) {
  if(!is.null(g)) {
    x <- c(if(is.atomic(g)) list(g) else if(is_GRP(g)) g[2L] else g,
           if(is.atomic(x)) list(x) else x, list(method = "radix"))
    return(do.call(order, x))
  }
  if(is.list(x)) return(do.call(order, c(x, list(method = "radix"))))
  if(length(x) < 1000L) .Call(C_radixsort, TRUE, FALSE, FALSE, FALSE, TRUE, pairlist(x)) else order(x, method = "radix")
}

fsetdiff <- function(x, y) x[match(x, y, 0L) == 0L] # not unique !

ffka <- function(x, f) {
   ax <- attributes(x)
  `attributes<-`(f(ax[["levels"]])[x],
   ax[names(ax) %!in% c("levels", "class")])
}


as_numeric_factor <- function(X, keep.attr = TRUE) {
  if(is.atomic(X)) if(keep.attr) return(ffka(X, as.numeric)) else
    return(as.numeric(attr(X, "levels"))[X])
  res <- duplAttributes(lapply(unattrib(X),
    if(keep.attr) (function(y) if(is.factor(y)) ffka(y, as.numeric) else y) else
                  (function(y) if(is.factor(y)) as.numeric(attr(y, "levels"))[y] else y)), X)
  if(inherits(X, "data.table")) return(alc(res))
  res
}

as_character_factor <- function(X, keep.attr = TRUE) {
  if(is.atomic(X)) if(keep.attr) return(ffka(X, tochar)) else
    return(as.character.factor(X))
  res <- duplAttributes(lapply(unattrib(X),
         if(keep.attr) (function(y) if(is.factor(y)) ffka(y, tochar) else y) else
                       (function(y) if(is.factor(y)) as.character.factor(y) else y)), X)
  if(inherits(X, "data.table")) return(alc(res))
  res
}

as.numeric_factor <- as_numeric_factor
as.character_factor <- as_character_factor

setRnDF <- function(df, nm) `attr<-`(df, "row.names", nm)

TtI <- function(x)
  switch(x, replace_fill = 1L, replace = 2L, `-` = 3L, `-+` = 4L, `/` = 5L, `%` = 6L, `+` = 7L, `*` = 8L, `%%` = 9L, `-%%` = 10L,
            stop("Unknown transformation!"))

condsetn <- function(x, value, cond) {
  if(cond) attr(x, "names") <- value
  x
}

setnck <- function(x, value) {
  if(is.null(value)) return(x)
  ren <- nzchar(value)
  if(all(ren)) names(x) <- value else names(x)[ren] <- value[ren]
  x
}

# give_nam <- function(x, gn, stub) {
#   if(!gn) return(x)
#   attr(x, "names") <- paste0(stub, attr(x, "names"))
#   x
# }

ckmatch <- function(x, table, e = "Unknown columns:") if(anyNA(m <- match(x, table))) stop(paste(e, paste(x[is.na(m)], collapse = ", "))) else m

# anyNAerror <- function(x, e) if(anyNA(x)) stop(e) else x

cols2int <- function(cols, x, nam, topos = TRUE) {
 if(is.numeric(cols)) {
   l <- length(unclass(x)) # length(nam) ?
   if(cols[1L] < 0L) { # This is sufficient to check negative indices: No R function allows subsetting mixing positive and negative indices.
     if(-min(cols) > l) stop("Index out of range abs(1:length(x))")
     if(topos) return(seq_len(l)[cols])
     # cols <- seq_len(l)[cols]
     # if(!length(cols) || anyNA(cols)) stop("Index out of range abs(1:length(x))") -> used to put earlier check after if(topos) and use this one instead. But turns out that doesn't always work well.
     # return(cols)
   } else if(max(cols) > l) stop("Index out of range abs(1:length(x))")
   # if(max(abs(cols)) > length(unclass(x))) stop("Index out of range abs(1:length(x))") # Before collapse 1.4.0 !
   return(as.integer(cols)) # as.integer is necessary (for C_subsetCols), and at very little cost..
 }
 if(is.character(cols)) return(ckmatch(cols, nam))
 if(is.function(cols)) return(which(vapply(unattrib(x), cols, TRUE)))
 if(is.logical(cols)) {
  if(length(cols) != length(unclass(x))) stop("Logical subsetting vector must match columns!") # length(nam) ?
  return(which(cols))
 }
 stop("cols must be a function, character vector, numeric indices or logical vector!")
}

# Not needed anymore !!
# cols2log <- function(cols, x, nam) {
#   lx <- length(unclass(x))
#   if(is.logical(cols)) if(length(cols) == lx) return(cols) else stop("Logical subsetting vector must match columns!")
#   if(is.function(cols)) return(vapply(unattrib(x), cols, TRUE))
#   r <- logical(lx)
#   if(is.character(cols)) {
#     r[ckmatch(cols, nam)] <- TRUE
#   } else if(is.numeric(cols)) {
#     if(max(abs(cols)) > lx) stop("Index out of range abs(1:length(x))")
#     r[cols] <- TRUE
#   } else stop("cols must be a function, character vector, numeric indices or logical vector!")
#   r
# }

colsubset <- function(x, ind, checksf = FALSE) {
  if(is.numeric(ind)) return(.Call(C_subsetCols, x, as.integer(ind), checksf))
  if(is.logical(ind)) {
    nc <- length(unclass(x))
    if(length(ind) != nc) stop("Logical subsetting vector must match length(x)")
    ind <- which(ind)
    if(length(ind) == nc) return(x)
    return(.Call(C_subsetCols, x, ind, checksf))
  }
  ind <- if(is.character(ind)) ckmatch(ind, attr(x, "names")) else which(vapply(`attributes<-`(x, NULL), ind, TRUE))
  return(.Call(C_subsetCols, x, ind, checksf))
}

# Previously Fastest! even though it involves code duplication..
# colsubset <- function(x, ind) {
#   ax <- attributes(x)
#   if(is.numeric(ind)) {
#     attributes(x) <- NULL # note: attributes(x) <- NULL is very slightly faster than class(x) <- NULL
#     if(max(abs(ind)) > length(x)) stop("Index out of range abs(1:length(x))")
#     ax[["names"]] <- ax[["names"]][ind]
#     return(.Call(C_setAttributes, x[ind], ax))
#   }
#   if(is.logical(ind)) {
#     attributes(x) <- NULL
#     if(length(ind) != length(x)) stop("Logical subsetting vector must match length(x)")
#     ax[["names"]] <- ax[["names"]][ind]
#     return(.Call(C_setAttributes, x[ind], ax))
#   }
#   ind <- if(is.character(ind)) ckmatch(ind, ax[["names"]]) else vapply(`attributes<-`(x, NULL), ind, TRUE)
#   ax[["names"]] <- ax[["names"]][ind]
#   .Call(C_setAttributes, .subset(x, ind), ax)
# }


fcolsubset <- function(x, ind, checksf = FALSE) { # fastest !
  .Call(C_subsetCols, x, if(is.logical(ind)) which(ind) else as.integer(ind), checksf)
  # Fastet! becore C version:
  # ax <- attributes(x)
  # ax[["names"]] <- ax[["names"]][ind]
  # .Call(C_setAttributes, .subset(x, ind), ax)
}

# Sorted out 1.5.3 -> 1.6.0:
# Fastest because vapply runs faster on a list without any attributes !
# colsubsetFUN <- function(x, FUN) {
#   .Call(C_subsetCols, x, which(vapply(`attributes<-`(x, NULL), FUN, TRUE)))
#   # Fastet! becore C version:
#   # ax <- attributes(x)
#   # attributes(x) <- NULL
#   # ind <- vapply(x, FUN, TRUE)
#   # ax[["names"]] <- ax[["names"]][ind]
#   # .Call(C_setAttributes, x[ind], ax)
# }

at2GRP <- function(x) {
  if(is.nmfactor(x)) return(list(length(attr(x, "levels")), x, NULL))
  res <- list(NULL, NULL, NULL)
  res[[2L]] <- qG(x, sort = FALSE, na.exclude = FALSE)
  res[[1L]] <- attr(res[[2L]], "N.groups")
  res
}

G_t <- function(x, wm = 1L) {
  if(is.null(x)) {
    if(wm > 0L) message(switch(wm, "Panel-lag computed without timevar: Assuming ordered data",
                             "Panel-difference computed without timevar: Assuming ordered data",
                             "Panel-growth rate computed without timevar: Assuming ordered data"))
    return(x)
  } # If integer time variable contains NA, noes not break C++ code..
  if(is.atomic(x)) if(is.integer(unclass(x))) return(x) else return(qG(x, na.exclude = FALSE, sort = TRUE, method = "hash")) # make sure it is sorted ! qG already checks factor !
  if(is_GRP(x)) return(x[[2L]]) else return(GRP.default(x, return.groups = FALSE, sort = TRUE, call = FALSE)[[2L]])
}

# Not currently used !!
# G_t2 <- function(x) {
#   if(is.atomic(x)) if(is.integer(unclass(x))) return(x) else return(qG(x, sort = TRUE, na.exclude = FALSE, method = "hash")) # Hashing seems generally faster for time-variables !!
#   if(is_GRP(x)) return(x[[2L]]) else return(GRP.default(x, return.groups = FALSE, sort = TRUE, call = FALSE)[[2L]])
# }


rgrep <- function(exp, nam, ..., sort = TRUE) if(length(exp) == 1L) grep(exp, nam, ...) else .Call(Cpp_funique, unlist(lapply(exp, grep, nam, ...), use.names = FALSE), sort)
rgrepl <- function(exp, nam, ...) if(length(exp) == 1L) grepl(exp, nam, ...) else Reduce(`|`, lapply(exp, grepl, nam, ...))

fanyDuplicated <- function(x) if(length(x) < 100L) anyDuplicated.default(x) > 0L else .Call(Cpp_fndistinct,x,0L,0L,NULL,FALSE) != length(x)

# NROW2 <- function(x, d) if(length(d)) d[1L] else length(x)
# NCOL2 <- function(d, ilv) if(ilv) d[2L] else 1L

charorNULL <- function(x) if(is.character(x)) x else NULL

tochar <- function(x) if(is.character(x)) x else as.character(x)

# more security here?
# unique_factor <- function(x) {  # Still needed with new collap solution ? -> Nope !
#   res <- seq_along(attr(x, "levels"))
#   .Call(C_duplAttributes, res, x)
# }
# dotstostr <- function(...) {
#   args <- deparse(substitute(c(...)))
#   nc <- nchar(args)
#   substr(args, 2, nc) # 3, nc-1 for no brackets !
# }

unused_arg_action <- function(call, ...) {
  wo <- switch(getOption("collapse_unused_arg_action"), none = 0L, message = 1L, warning = 2L, error = 3L,
               stop("Unused argument encountered. Please instruct collapse what to do about unused arguments by setting
                     options(collapse_unused_arg_action = 'warning'), or 'error', or 'message' or 'none'."))
  if(wo != 0L) {
    args <- deparse(substitute(c(...)))
    nc <- nchar(args)
    args <- substr(args, 2, nc) # 3, nc-1 for no brackets !
    msg <- paste("Unused argument", args, "passed to", as.character(call[[1L]]))
    switch(wo, message(msg), warning(msg), stop(msg))
  }
}

is.nmfactor <- function(x) inherits(x, "factor") && (inherits(x, "na.included") || !anyNA(unclass(x)))

addNA2 <- function(x) {
  if(!anyNA(unclass(x))) return(x)
  clx <- oldClass(x)
  oldClass(x) <- NULL
  if(!anyNA(lev <- attr(x, "levels"))) {
    attr(x, "levels") <- c(lev, NA_character_)
    x[is.na(x)] <- length(lev) + 1L
  } else x[is.na(x)] <- length(lev)
  `oldClass<-`(x, clx)
}

# addNA2 <- function(x) {
#   clx <- c(class(x), "na.included")
#   if(!anyNA(unclass(x))) return(`oldClass<-`(x, clx))
#   ll <- attr(x, "levels")
#   if(!anyNA(ll)) ll <- c(ll, NA)
#   return(`oldClass<-`(factor(x, levels = ll, exclude = NULL), clx))
# }


l1orn <- function(x, nam) if(length(x) == 1L) x else nam
l1orlst <- function(x) if(length(x) == 1L) x else x[length(x)]

fsimplify2array <- function(l) {
  res <- do.call(cbind, l) # lapply(l, `dimnames<-`, NULL) # also faster than unlist..
  dim(res) <- c(dim(l[[1L]]), length(l))
  dimnames(res) <- c(if(length(dn <- dimnames(l[[1L]]))) dn else list(NULL, NULL), list(names(l)))
  res
}


# fss <- function(x, i, j) {
#   rn <- attr(x, "row.names")
#   if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(.Call(C_subsetDT, x, i, j))
#   return(`attr<-`(.Call(C_subsetDT, x, i, j), "row.names", rn[r]))
# }

# fsplit_DF <- function(x, j, f, rnl, ...) {
#   j <- seq_along(unclass(x))
#   rn <- attr(x, "row.names")
#   if(is.numeric(rn) || is.null(rn) || rn[1L] == "1")
#     return(lapply(split.default(seq_along(.subset2(x, 1L)), f, ...),
#                   function(i) .Call(C_subsetDT, x, i, j)))
#   lapply(split.default(seq_along(.subset2(x, 1L)), f, ...),
#          function(i) `attr<-`(.Call(C_subsetDT, x, i, j), "row.names", rn[i]))
# }






