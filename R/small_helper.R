
null2NA <- function(x) if(is.null(x)) NA_character_ else x

# flapply <- function(x, FUN, ...) lapply(unattrib(x), FUN, ...) # not really needed ...

vlabels <- function(X, attrn = "label") {
  if(is.atomic(X)) return(null2NA(attr(X, attrn)))
  res <- lapply(X, attr, attrn) # unattrib(X): no names
  res[vapply(res, is.null, TRUE)] <- NA_character_
  unlist(res)
}

"vlabels<-" <- function(X, attrn = "label", value) {
  if(is.atomic(X)) return(`attr<-`(X, attrn, value))
  clx <- oldClass(X)
  oldClass(X) <- NULL
  if(length(X) != length(value)) stop("length(X) must match length(value)")
  for (i in seq_along(value)) attr(X[[i]], attrn) <- value[i]
  `oldClass<-`(X, clx)
}

# or cns ??
.c <- function(...) as.character(substitute(c(...))[-1L])
# formula = formula(paste0("~ ", paste(as.character(substitute(c(...)))[-1], collapse = " + ")))) # could speed up ??

strclp <- function(x) if(length(x) > 1L) paste(x, collapse = " ") else x

pasteclass <- function(x) if(length(cx <- class(x)) > 1L) paste(cx, collapse = " ") else cx # Faster if length(class(x)) == 1L

vclasses <- function(X) {
  if(is.atomic(X)) return(pasteclass(X))
  vapply(X, pasteclass, character(1)) # unattrib(X): no names !
}

vtypes <- function(X) {
  if(is.atomic(X)) return(typeof(X))
  vapply(X, typeof, character(1)) # unattrib(X): no names !
}

namlab <- function(X, class = FALSE, attrn = "label") {
  if(!is.list(X)) stop("namlab only works with lists")
  res <- if(class) list(attr(X, "names"), vapply(unattrib(X), pasteclass, character(1)), vlabels(X, attrn)) else list(attr(X, "names"), vlabels(X, attrn))
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
    dimnames(X) <- list(dn[[1L]], if(pre) paste0(stub, dn[[2L]]) else paste0(dn[[2L]], stub))
  } else attr(X, "names") <- if(pre) paste0(stub, attr(X, "names")) else paste0(attr(X, "names"), stub)
  X
}

rm_stub <- function(X, stub, pre = TRUE) {
  if(!is.character(stub)) return(X)
  if(pre)
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
  } else attr(X, "names") <- rmstubFUN(attr(X, "names"))
  X
}

setRownames <- function(object, nm = seq_row(object)) {
  if(is.list(object)) {
    l <- length(.subset2(object, 1L))
    if(is.null(nm)) nm <- .set_row_names(l) else if(length(nm) != l) stop("supplied row-names must match list extent")
    return(`attr<-`(object, "row.names", nm))
  }
  if(!is.array(object)) stop("Setting row-names only supported on arrays and lists")
  dn <- dimnames(object)
 `dimnames<-`(object, c(list(nm), dn[-1L]))
}

setColnames <- function(object, nm) {
  if(is.atomic(object) && is.array(object))
    dimnames(object)[[2L]] <- nm else
    attr(object, "names") <- nm
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
  if(...length() == 1L && is.list(...)) return(all(vapply(unattrib(...)[-1L], all.equal, TRUE, .subset2(..., 1L))))
  l <- list(...)
  all(vapply(l[-1L], all.equal, TRUE, l[[1L]]))
}

interact_names <- function(l) do.call(paste, c(expand.grid(l, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE), list(sep = ".")))

unattrib <- function(object) `attributes<-`(object, NULL)

is.categorical <- function(x) !is.numeric(x)

is.Date <- function(x) inherits(x, c("Date","POSIXlt","POSIXct"))

"%!in%" <- function(x, table) match(x, table, nomatch = 0L) == 0L

na_rm <- function(x) x[!is.na(x)]

# more consistent with base than na_rm
# na.rm <- function(x) { # cpp version available, but not faster !
#   if(!is.null(attr(x, "names"))) { # gives corruped time-series !
#     ax <- attributes(x)
#     r <- x[!is.na(x)]
#     ax[["names"]] <- names(r)
#     setAttributes(r, ax)
#   } else duplAttributes(x[!is.na(x)], x)
# }

# fast na.omit !
na_omit <- function(X, cols = NULL, na.attr = FALSE) {
  if(is.list(X)) {
    iX <- seq_along(unclass(X))
    rl <- if(is.null(cols)) !.Call(C_dt_na, X, iX) else !.Call(C_dt_na, X, cols2int(cols, X, attr(X, "names"))) # gives error if X not list
    rn <- attr(X, "row.names")
    if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") res <- .Call(C_subsetDT, X, which(rl), iX) else {
      rkeep <- which(rl)
      res <- .Call(C_subsetDT, X, rkeep, iX)
      attr(res, "row.names") <- rn[rkeep]
    }
  } else {
    rl <- if(is.null(cols)) complete.cases(X) else complete.cases(X[, cols])
    res <- if(is.matrix(X)) X[rl, , drop = FALSE] else X[rl]
  }
  if(na.attr) attr(res, "na.action") <- `oldClass<-`(which(!rl), "omit")
  res
}

na_insert <- function(X, prop = 0.1) {
  if(is.list(X)) {
    n <- fnrow2(X)
    nmiss <- floor(n * prop)
    return(duplAttributes(lapply(unattrib(X), function(x) `[<-`(x, sample.int(n, nmiss), value = NA)), X))
  } else if(!is.null(d <- fdim(X))) {
    n <- d[1L]
    p <- d[2L]
    NAloc <- rep(FALSE, n * p)
    NAloc[sample.int(n * p, floor(n * p * prop))] <- TRUE
    X[matrix(NAloc, nrow = n, ncol = p)] <- NA
  } else if(is.atomic(X)) {
    l <- length(X)
    X[sample.int(l, floor(l * prop))] <- NA
  } else stop("X must be an atomic vector, matrix or data.frame")
  X
}

fdapply <- function(X, FUN, ...) duplAttributes(lapply(`attributes<-`(X, NULL), FUN, ...), X)

fnlevels <- function(x) length(attr(x, "levels"))

# flevels <- function(x) attr(x, "levels")

fnrow <- function(X) if(is.list(X)) length(.subset2(X, 1L)) else dim(X)[1L]

fnrow2 <- function(x) length(.subset2(x, 1L))

fncol <- function(X) if(is.list(X)) length(unclass(X)) else dim(X)[2L]

fNCOL <- function(X) if(is.list(X)) length(unclass(X)) else NCOL(X)

fdim <- function(X) {
   if(is.atomic(X)) return(dim(X)) # or if !is.list ?
   oldClass(X) <- NULL
   c(length(X[[1L]]), length(X))
}

seq_row <- function(X) if(is.list(X)) seq_along(.subset2(X, 1L)) else seq_len(nrow(X))

seq_col <- function(X) if(is.list(X)) seq_along(unclass(X)) else seq_len(ncol(X))

# na.last is false !!
forder.int <- function(x) .Call(C_radixsort, FALSE, FALSE, FALSE, FALSE, TRUE, pairlist(x)) # if(is.unsorted(x)) .Call(C_forder, x, NULL, FALSE, TRUE, 1L, TRUE) else seq_along(x) # since forder gives integer(0) if sorted !
forder.vec <- function(x) if(length(x) < 1000L) .Call(C_radixsort, TRUE, FALSE, FALSE, FALSE, TRUE, pairlist(x)) else order(x, method = "radix")

fsetdiff <- function(x, y) x[match(x, y, 0L) == 0L] # not unique !

ffka <- function(x, f) {
   ax <- attributes(x)
  `attributes<-`(f(ax[["levels"]])[x],
   ax[names(ax) %!in% c("levels", "class")])
}

as.numeric_factor <- function(X, keep.attr = TRUE) {
  if(is.atomic(X)) if(keep.attr) return(ffka(X, as.numeric)) else
    return(as.numeric(attr(X, "levels"))[X])
  fcts <- vapply(unattrib(X), is.factor, TRUE)
  clx <- oldClass(X)
  oldClass(X) <- NULL
  X[fcts] <- if(keep.attr) lapply(X[fcts], ffka, as.numeric) else
    lapply(X[fcts], function(x) as.numeric(attr(x, "levels"))[x])
  return(`oldClass<-`(X, clx))
}

as.character_factor <- function(X, keep.attr = TRUE) {
  if(is.atomic(X)) if(keep.attr) return(ffka(X, tochar)) else
    return(as.character.factor(X))
  fcts <- vapply(unattrib(X), is.factor, TRUE)
  clx <- oldClass(X)
  oldClass(X) <- NULL
  X[fcts] <- if(keep.attr) lapply(X[fcts], ffka, tochar) else
    lapply(X[fcts], as.character.factor)
  return(`oldClass<-`(X, clx))
}

setRnDF <- function(df, nm) `attr<-`(df, "row.names", nm)

addAttributes <- function(x, a) .Call(C_setAttributes, x, c(attributes(x), a))

TtI <- function(x)
  switch(x, replace_fill = 1L, replace = 2L, `-` = 3L, `-+` = 4L, `/` = 5L, `%` = 6L, `+` = 7L, `*` = 8L, `%%` = 9L, `-%%` = 10L,
            stop("Unknown transformation!"))

condsetn <- function(x, value, cond) {
  if(cond) attr(x, "names") <- value
  x
}

# give_nam <- function(x, gn, stub) {
#   if(!gn) return(x)
#   attr(x, "names") <- paste0(stub, attr(x, "names"))
#   x
# }

ckmatch <- function(x, table, e = "Unknown columns:") if(anyNA(m <- match(x, table))) stop(paste(e, paste(x[is.na(m)], collapse = ", "))) else m

# anyNAerror <- function(x, e) if(anyNA(x)) stop(e) else x

cols2int <- function(cols, x, nam) {
 if(is.numeric(cols)) {
  if(max(abs(cols)) > length(unclass(x))) stop("Index out of range abs(1:length(x))")
  return(cols)
 }
 if(is.character(cols)) return(ckmatch(cols, nam))
 if(is.function(cols)) return(which(vapply(unattrib(x), cols, TRUE)))
 if(is.logical(cols)) {
  if(length(cols) != length(unclass(x))) stop("Logical subsetting vector must match columns!")
  return(which(cols))
 }
 stop("cols must be a function, character vector, numeric indices or logical vector!")
}

cols2log <- function(cols, x, nam) {
  lx <- length(unclass(x))
  if(is.logical(cols)) if(length(cols) == lx) return(cols) else stop("Logical subsetting vector must match columns!")
  if(is.function(cols)) return(vapply(unattrib(x), cols, TRUE))
  r <- logical(lx)
  if(is.character(cols)) {
    r[ckmatch(cols, nam)] <- TRUE
  } else if(is.numeric(cols)) {
    if(max(abs(cols)) > lx) stop("Index out of range abs(1:length(x))")
    r[cols] <- TRUE
  } else stop("cols must be a function, character vector, numeric indices or logical vector!")
  r
}

# Fastest! even though it involves code duplication..
colsubset <- function(x, ind) {
  ax <- attributes(x)
  if(is.numeric(ind)) {
    attributes(x) <- NULL # note: attributes(x) <- NULL is very slightly faster than class(x) <- NULL
    if(max(abs(ind)) > length(x)) stop("Index out of range abs(1:length(x))")
    ax[["names"]] <- ax[["names"]][ind]
    return(.Call(C_setAttributes, x[ind], ax))
  }
  if(is.logical(ind)) {
    attributes(x) <- NULL
    if(length(ind) != length(x)) stop("Logical subsetting vector must match length(x)")
    ax[["names"]] <- ax[["names"]][ind]
    return(.Call(C_setAttributes, x[ind], ax))
  }
  ind <- if(is.character(ind)) ckmatch(ind, ax[["names"]]) else vapply(`attributes<-`(x, NULL), ind, TRUE)
  ax[["names"]] <- ax[["names"]][ind]
  .Call(C_setAttributes, .subset(x, ind), ax)
}

# Previous: slower..
# colsubset <- function(x, ind) {
#   ax <- attributes(x)
#   if(is.numeric(ind)) {
#     if(max(abs(ind)) > length(unclass(x))) stop("Index out of range abs(1:length(x))")
#   } else if(is.logical(ind)) {
#     if(length(ind) != length(unclass(x))) stop("Logical subsetting vector must match length(x)")
#   } else ind <- if(is.character(ind)) ckmatch(ind, ax[["names"]]) else vapply(`attributes<-`(x, NULL), ind, TRUE)
#   ax[["names"]] <- ax[["names"]][ind]
#   .Call(C_setAttributes, .subset(x, ind), ax) # return(`attributes<-`(x[ind], ax)) # This is slow on large data -> a lot of checks !
# }


fcolsubset <- function(x, ind) { # fastest !
  ax <- attributes(x)
  ax[["names"]] <- ax[["names"]][ind]
  .Call(C_setAttributes, .subset(x, ind), ax)
}

# Fastest because vapply runs faster on a list without any attributes !
colsubsetFUN <- function(x, FUN) {
  ax <- attributes(x)
  attributes(x) <- NULL
  ind <- vapply(x, FUN, TRUE)
  ax[["names"]] <- ax[["names"]][ind]
  .Call(C_setAttributes, x[ind], ax)
}

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
  if(is.GRP(x)) return(x[[2L]]) else return(GRP.default(x, return.groups = FALSE, sort = TRUE, call = FALSE)[[2L]])
}

# Not currently used !!
# G_t2 <- function(x) {
#   if(is.atomic(x)) if(is.integer(unclass(x))) return(x) else return(qG(x, sort = TRUE, na.exclude = FALSE, method = "hash")) # Hashing seems generally faster for time-variables !!
#   if(is.GRP(x)) return(x[[2L]]) else return(GRP.default(x, return.groups = FALSE, sort = TRUE, call = FALSE)[[2L]])
# }


rgrep <- function(exp, nam, ...) if(length(exp) == 1L) grep(exp, nam, ...) else .Call(Cpp_funique, unlist(lapply(exp, grep, nam, ...), use.names = FALSE), TRUE)
rgrepl <- function(exp, nam, ...) if(length(exp) == 1L) grepl(exp, nam, ...) else Reduce(`|`, lapply(exp, grepl, nam, ...))

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
  dimnames(res) <- c(dimnames(l[[1L]]), list(names(l)))
  res
}


# Experimental:

# faster pipe: more challenging than it seems...
# `%>>%` <- function(lhs, rhs) {
#     rhs_call <- substitute(rhs)
#     eval(rhs_call, envir = lhs, enclos = parent.frame())
# }
#
# > pipeR::`%>>%`
# function (x, expr)
# {
#   x
#   expr <- substitute(expr)
#   envir <- parent.frame()
#   switch(class(expr), `NULL` = NULL, character = {
#     cat(expr, "\n")
#     x
#   }, `{` = pipe_dot(x, expr, envir), `(` = pipe_fun(x,
#                                                     expr[[2L]], envir), pipe_symbol(x, expr, envir, TRUE,
#                                                                                     pipe_first))
# }

