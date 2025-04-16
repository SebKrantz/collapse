# Functions needed for internal use because of option(collapse_mask = "fast-stat-fun")
bsum <- base::sum
bprod <- base::prod
bmin <- base::min
bmax <- base::max

# Row-operations (documented under data transformations...) ...

"%rr%" <- function(X, v) if(is.atomic(X) || is.atomic(v) || inherits(X, "data.frame")) TRA(X, v, "replace_fill") else # outer(rep.int(1L, dim(X)[2L]), v)
  duplAttributes(.mapply(function(x, y) TRA(x, y, "replace_fill"), list(unattrib(X), unattrib(v)), NULL), X)
"%r+%" <- function(X, v) if(is.atomic(X) || is.atomic(v) || inherits(X, "data.frame")) TRA(X, v, "+") else
  duplAttributes(.mapply(function(x, y) TRA(x, y, "+"), list(unattrib(X), unattrib(v)), NULL), X)
"%r-%" <- function(X, v) if(is.atomic(X) || is.atomic(v) || inherits(X, "data.frame")) TRA(X, v, "-") else
  duplAttributes(.mapply(function(x, y) TRA(x, y, "-"), list(unattrib(X), unattrib(v)), NULL), X)
"%r*%" <- function(X, v) if(is.atomic(X) || is.atomic(v) || inherits(X, "data.frame")) TRA(X, v, "*") else
  duplAttributes(.mapply(function(x, y) TRA(x, y, "*"), list(unattrib(X), unattrib(v)), NULL), X)
"%r/%" <- function(X, v) if(is.atomic(X) || is.atomic(v) || inherits(X, "data.frame")) TRA(X, v, "/") else
  duplAttributes(.mapply(function(x, y) TRA(x, y, "/"), list(unattrib(X), unattrib(v)), NULL), X)


"%cr%" <- function(X, V) if(is.atomic(X)) return(duplAttributes(rep(V, NCOL(X)), X)) else # outer(rep.int(1L, dim(X)[2L]), V)
  if(is.atomic(V)) return(duplAttributes(lapply(vector("list", length(unclass(X))), function(z) V), X)) else
    copyAttrib(V, X) # copyAttrib first makes a shallow copy of V
"%c+%" <- function(X, V) if(is.atomic(X)) return(X + V) else
  duplAttributes(if(is.atomic(V)) lapply(unattrib(X), `+`, V) else
    .mapply(`+`, list(unattrib(X), unattrib(V)), NULL), X)
"%c-%" <- function(X, V) if(is.atomic(X)) return(X - V) else
  duplAttributes(if(is.atomic(V)) lapply(unattrib(X), `-`, V) else
    .mapply(`-`, list(unattrib(X), unattrib(V)), NULL), X)
"%c*%" <- function(X, V) if(is.atomic(X)) return(X * V) else
  duplAttributes(if(is.atomic(V)) lapply(unattrib(X), `*`, V) else
    .mapply(`*`, list(unattrib(X), unattrib(V)), NULL), X)
"%c/%" <- function(X, V) if(is.atomic(X)) return(X / V) else  # or * 1L/V ??
  duplAttributes(if(is.atomic(V)) lapply(unattrib(X), `/`, V) else
    .mapply(`/`, list(unattrib(X), unattrib(V)), NULL), X)


# Multiple-assignment
"%=%" <- function(nam, values) invisible(.Call(C_multiassign, nam, values, parent.frame()))
massign <- function(nam, values, envir = parent.frame()) invisible(.Call(C_multiassign, nam, values, envir))

# R implementation:
# "%=%" <- function(lhs, rhs) {
#   if(!is.character(lhs)) stop("lhs needs to be character")
#   if(!is.list(rhs)) rhs <- as.vector(rhs, "list")
#   if(length(lhs) != length(rhs)) stop("length(lhs) not equal to length(rhs)")
#   list2env(`names<-`(rhs, lhs), envir = parent.frame(),
#            parent = NULL, hash = FALSE, size = 0L)
#   invisible()
# }


getenvFUN <- function(nam, efmt1 = "For this method need to install.packages('%s'), then unload [detach('package:collapse', unload = TRUE)] and reload [library(collapse)].")
{
  if(is.null(FUN <- .collapse_env[[nam]])) {
    v <- strsplit(nam, "_", fixed = TRUE)[[1L]]
    .collapse_env[[nam]] <- FUN <- if(requireNamespace(v[1L], quietly = TRUE))
           get0(v[2L], envir = getNamespace(v[1L])) else NULL
    if(is.null(FUN)) stop(sprintf(efmt1, v[1L]))
  }
  FUN
}
# qM2 <- function(x) if(is.list(x)) do.call(cbind, x) else x

null2NA <- function(x) if(is.null(x)) NA_character_ else x

# flapply <- function(x, FUN, ...) lapply(unattrib(x), FUN, ...) # not really needed ...

vlabels <- function(X, attrn = "label", use.names = TRUE) .Call(C_vlabels, X, attrn, use.names) # {
  # if(is.atomic(X)) return(null2NA(attr(X, attrn)))
  # res <- lapply(X, attr, attrn) # unattrib(X): no names
  # res[vapply(res, is.null, TRUE)] <- NA_character_
  # unlist(res)
# }

"vlabels<-" <- function(X, attrn = "label", value) {
  if(is.atomic(X)) return(`attr<-`(X, attrn, value))
  .Call(C_setvlabels, X, attrn, value, NULL)
}

# "vlabels<-" <- function(X, attrn = "label", value) {
#   names(value) <- NULL
#   if(is.atomic(X)) return(`attr<-`(X, attrn, value))
#   clx <- oldClass(X)
#   oldClass(X) <- NULL
#   if(is.null(value)) {
#     for (i in seq_along(X)) attr(X[[i]], attrn) <- NULL
#   } else {
#     if(length(X) != length(value)) stop("length(X) must match length(value)")
#     for (i in seq_along(value)) attr(X[[i]], attrn) <- value[[i]]
#   }
#   if(any(clx == "data.table")) return(alc(`oldClass<-`(X, clx)))
#   `oldClass<-`(X, clx)
# }

# Note: Shallow copy does not work as it only copies the list, but the attribute is a feature of the atomic elements inside...
setLabels <- function(X, value = NULL, attrn = "label", cols = NULL) { # , sc = TRUE
  if(is.atomic(X)) return(`attr<-`(X, attrn, value))
  .Call(C_setvlabels, X, attrn, value, as.integer(cols))
}

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

vclasses <- function(X, use.names = TRUE) {
  if(is.atomic(X)) return(pasteclass(X))
  vapply(X, pasteclass, "", USE.NAMES = use.names) # unattrib(X): no names
}

# https://github.com/wch/r-source/blob/4a409a1a244d842a3098d2783c5b63c9661fc6be/src/main/util.c
R_types <- c("NULL",	      # NILSXP
             "symbol",      # SYMSXP
             "pairlist",	  # LISTSXP
             "closure",	    # CLOSXP
             "environment", # ENVSXP
             "promise",	    # PROMSXP
             "language",	  # LANGSXP
             "special",	    # SPECIALSXP
             "builtin",	    # BUILTINSXP
             "char",		    # CHARSXP
             "logical",	    # LGLSXP
             "",
             "",
             "integer",	    # INTSXP
             "double",	    # REALSXP
             "complex",	    # CPLXSXP
             "character",	  # STRSXP
             "...",		      # DOTSXP
             "any",		      # ANYSXP
             "list",		    # VECSXP
             "expression",	# EXPRSXP
             "bytecode",	  # BCODESXP
             "externalptr",	# EXTPTRSXP
             "weakref",	    # WEAKREFSXP
             "raw",		      # RAWSXP
             "S4")		      # S4SXP
# /* aliases : */
# { "numeric",	REALSXP	   },
# { "name",		SYMSXP	   },


vtypes <- function(X, use.names = TRUE) {
  if(is.atomic(X)) return(typeof(X))
  res <- R_types[.Call(C_vtypes, X, 0L)]
  if(use.names) names(res) <- attr(X, "names")
  res
  # vapply(X, typeof, "") # unattrib(X): no names
}

vlengths <- function(X, use.names = TRUE) .Call(C_vlengths, X, use.names)

namlab <- function(X, class = FALSE, attrn = "label", N = FALSE, Ndistinct = FALSE) {
  if(!is.list(X)) stop("namlab only works with lists")
  res <- list(Variable = attr(X, "names"))
  attributes(X) <- NULL
  if(class) res$Class <- vapply(X, pasteclass, "", USE.NAMES = FALSE)
  if(N) res$N <- fnobs.data.frame(X)
  if(Ndistinct) res$Ndist <- fndistinct.data.frame(X, na.rm = TRUE)
  res$Label <- vlabels(X, attrn, FALSE)
  attr(res, "row.names") <- c(NA_integer_, -length(X))
  oldClass(res) <- "data.frame"
  res
}

add_stub <- function(X, stub, pre = TRUE, cols = NULL) {
  if(!is.character(stub)) return(X)
  if(is.atomic(X) && is.array(X)) {
    if(length(dim(X)) > 2L) stop("Can't stub higher dimensional arrays!")
    dn <- dimnames(X)
    cn <- dn[[2L]]
    if(length(cn)) {
      if(length(cols)) cn[cols] <- if(pre) paste0(stub, cn[cols]) else paste0(cn[cols], stub)
      else cn <- if(pre) paste0(stub, cn) else paste0(cn, stub)
      dimnames(X) <- list(dn[[1L]], cn)
    }
  } else {
    nam <- attr(X, "names")
    if(length(nam)) {
      if(length(cols)) attr(X, "names")[cols] <- if(pre) paste0(stub, nam[cols]) else paste0(nam[cols], stub)
      else attr(X, "names") <- if(pre) paste0(stub, nam) else paste0(nam, stub)
      if(inherits(X, "data.table")) X <- alc(X)
    }
  }
  X
}

rm_stub <- function(X, stub, pre = TRUE, regex = FALSE, cols = NULL, ...) {
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
    if(is.null(d)) if(is.character(X)) return(if(length(cols)) replace(X, cols, rmstubFUN(X[cols])) else rmstubFUN(X)) else stop("Cannot modify a vector that is not character")
    if(length(d) > 2L) stop("Can't remove stub from higher dimensional arrays!")
    dn <- dimnames(X)
    cn <- dn[[2L]]
    dimnames(X) <- list(dn[[1L]], if(length(cols)) replace(cn, cols, rmstubFUN(cn[cols])) else rmstubFUN(cn))
  } else {
    nam <- attr(X, "names")
    attr(X, "names") <- if(length(cols)) replace(nam, cols, rmstubFUN(nam[cols])) else rmstubFUN(nam)
    if(inherits(X, "data.table")) X <- alc(X)
  }
  X
}

setRownames <- function(object, nm = if(is.atomic(object)) seq_row(object) else NULL) {
  if(is.list(object)) {
    l <- .Call(C_fnrow, object)
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

all_funs <- function(expr) .Call(C_all_funs, expr)

cinv <- function(x) chol2inv(chol(x))

vec <- function(X) {
  if(is.atomic(X)) return(`attributes<-`(X, NULL))
  .Call(C_pivot_long, X, NULL, FALSE)
}

interact_names <- function(l) {
  oldClass(l) <- NULL
  if(length(l) == 2L) return(`dim<-`(outer(l[[1L]], l[[2L]], paste, sep = "."), NULL))
  do.call(paste, c(expand.grid(l, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE), list(sep = ".")))
}

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
setattrib <- function(object, a) {
  .Call(C_setattributes, object, a)
  return(invisible(object))
}

# setAttribR <- function(object, a) `attributes<-`(object, x)

copyAttrib <- function(to, from) .Call(C_copyAttrib, to, from)
# copyAttribR <- function(to, from) `attributes<-`(to, attributes(from))


copyMostAttrib <- function(to, from) .Call(C_copyMostAttrib, to, from)
# copyMostAttribR <- function(to, from) `mostattributes<-`(to, attributes(from))

addAttributes <- function(x, a) .Call(C_setAttributes, x, c(attributes(x), a))


is_categorical <- function(x) !is.numeric(x)
# is.categorical <- function(x) {
#   .Deprecated(msg = "'is.categorical' was renamed to 'is_categorical'. It will be removed end of 2023, see help('collapse-renamed').")
#   !is.numeric(x)
# }

is_date <- function(x) inherits(x, c("Date","POSIXlt","POSIXct"))
# is.Date <- function(x) {
#   .Deprecated(msg = "'is.Date' was renamed to 'is_date'. It will be removed end of 2023, see help('collapse-renamed').")
#   inherits(x, c("Date","POSIXlt","POSIXct"))
# }

# more consistent with base than na_rm
# na.rm <- function(x) { # cpp version available, but not faster !
#   if(length(attr(x, "names"))) { # gives corruped time-series !
#     ax <- attributes(x)
#     r <- x[!is.na(x)]
#     ax[["names"]] <- names(r)
#     setAttributes(r, ax)
#   } else duplAttributes(x[!is.na(x)], x)
# }

whichv <- function(x, value, invert = FALSE) .Call(C_whichv, x, value, invert)
"%==%" <- function(x, value) .Call(C_whichv, x, value, FALSE)
"%!=%" <- function(x, value) .Call(C_whichv, x, value, TRUE)
whichNA <- function(x, invert = FALSE) .Call(C_whichv, x, NA, invert)

frange <- function(x, na.rm = .op[["na.rm"]], finite = FALSE) .Call(C_frange, x, na.rm, finite)
.range <- function(x, na.rm = TRUE, finite = FALSE) .Call(C_frange, x, na.rm, finite)
alloc <- function(value, n, simplify = TRUE) .Call(C_alloc, value, n, simplify)
vgcd <- function(x) .Call(C_vecgcd, x)
fdist <- function(x, v = NULL, ..., method = "euclidean", nthreads = .op[["nthreads"]]) .Call(C_fdist, if(is.atomic(x)) x else qM(x), v, method, nthreads)

allNA <- function(x) .Call(C_allNA, x, TRUE) # True means give error for unsupported vector types, not FALSE.
anyv <- function(x, value) .Call(C_anyallv, x, value, FALSE)
allv <- function(x, value) .Call(C_anyallv, x, value, TRUE)


copyv <- function(X, v, R, ..., invert = FALSE, vind1 = FALSE, xlist = FALSE) {
  if(is.list(X, ...) && !xlist) { # Making sure some error is produced if dots are used
    if(is.list(R)) {
      res <- .mapply(function(x, r) .Call(C_setcopyv, x, v, r, invert, FALSE, vind1),
                     list(unattrib(X), unattrib(R)), NULL)
    } else {
      res <- lapply(unattrib(X), function(x) .Call(C_setcopyv, x, v, R, invert, FALSE, vind1))
    }
    return(condalc(duplAttributes(res, X), inherits(X, "data.table")))
  }
  .Call(C_setcopyv, X, v, R, invert, FALSE, vind1)
}
setv  <- function(X, v, R, ..., invert = FALSE, vind1 = FALSE, xlist = FALSE) {
  if(is.list(X, ...) && !xlist) { # Making sure some error is produced if dots are used
    if(is.list(R)) {
      .mapply(function(x, r) .Call(C_setcopyv, x, v, r, invert, TRUE, vind1),
              list(unattrib(X), unattrib(R)), NULL)
    } else {
      lapply(unattrib(X), function(x) .Call(C_setcopyv, x, v, R, invert, TRUE, vind1))
    }
    return(invisible(X))
  }
  invisible(.Call(C_setcopyv, X, v, R, invert, TRUE, vind1))
}

setop <- function(X, op, V, ..., rowwise = FALSE) # Making sure some error is produced if dots are used
  invisible(.Call(C_setop, X, V, switch(op, "+" = 1L, "-" = 2L, "*" = 3L, "/" = 4L, stop("Unsupported operation:", op)), rowwise), ...)

"%+=%" <- function(X, V) invisible(.Call(C_setop, X, V, 1L, FALSE))
"%-=%" <- function(X, V) invisible(.Call(C_setop, X, V, 2L, FALSE))
"%*=%" <- function(X, V) invisible(.Call(C_setop, X, V, 3L, FALSE))
"%/=%" <- function(X, V) invisible(.Call(C_setop, X, V, 4L, FALSE))

# Internal functions
missDF <- function(x, cols = seq_along(unclass(x))) .Call(C_dt_na, x, cols, 0, FALSE)
frowSums <- function(x) {
  nr <- dim(x)[1L]
  .rowSums(x, nr, length(x)/nr)
}
fcolSums <- function(x) {
  nr <- dim(x)[1L]
  .colSums(x, nr, length(x)/nr)
}

missing_cases <- function(X, cols = NULL, prop = 0, count = FALSE) {
  if(is.list(X)) return(.Call(C_dt_na, X, if(is.null(cols)) seq_along(unclass(X)) else cols2int(cols, X, attr(X, "names")), prop, count))
  if(is.matrix(X)) {
    if(length(cols)) X <- X[, cols]
    if(is.matrix(X)) return(if(count) as.integer(frowSums(is.na(X))) else if(prop > 0) # as.integer() needed to establish consistency (integer output)
      frowSums(is.na(X)) >= bmax(as.integer(prop * NCOL(X)), 1L) else !complete.cases(X))
  }
  if(count) as.integer(is.na(X)) else is.na(X) # Note: as.integer() here is inefficient, but storage.mode() <- "integer" is also. Would have to export a R wrapper to C function SET_TYPEOF()... but this is probably never invoked anyway.
}

na_rm <- function(x) .Call(C_na_rm, x)  # x[!is.na(x)]
# Also takes names along, whereas na_rm does not preserve names of list
null_rm <- function(l) if(!all(ind <- vlengths(l, FALSE) > 0L)) .subset(l, ind) else l

all_eq <- function(x) .Call(C_anyallv, x, x[1L], TRUE)

na_omit <- function(X, cols = NULL, na.attr = FALSE, prop = 0, ...) {
  if(is.list(X)) {
    iX <- seq_along(unclass(X))
    rl <- .Call(C_dt_na, X, if(is.null(cols)) iX else cols2int(cols, X, attr(X, "names")), prop, FALSE)
    rkeep <- whichv(rl, FALSE)
    if(length(rkeep) == fnrow(X)) return(condalc(X, inherits(X, "data.table")))
    res <- .Call(C_subsetDT, X, rkeep, iX, FALSE) # This allocates data.tables...
    rn <- attr(X, "row.names")
    if(!(is.numeric(rn) || is.null(rn) || rn[1L] == "1")) attr(res, "row.names") <- Csv(rn, rkeep)
    if(na.attr) {
      attr(res, "na.action") <- `oldClass<-`(which(rl), "omit")
      if(inherits(res, "data.table") && !inherits(X, "pdata.frame")) return(alc(res))
    }
    if(inherits(X, "pdata.frame")) {
      index <- findex(X)
      index_omit <- droplevels_index(.Call(C_subsetDT, index, rkeep, seq_along(unclass(index)), FALSE), ...)
      if(inherits(X, "indexed_frame")) return(reindex(res, index_omit)) # data.table handled here
      attr(res, "index") <- index_omit
    }
  } else {
    Xcols <- if(is.null(cols)) X else X[, cols]
    rl <- if(prop > 0 && is.matrix(Xcols)) frowSums(is.na(Xcols)) < bmax(as.integer(prop * ncol(Xcols)), 1L) else complete.cases(Xcols)
    rkeep <- which(rl)
    if(length(rkeep) == NROW(X)) return(X)
    res <- if(is.matrix(X)) X[rkeep, , drop = FALSE, ...] else X[rkeep, ...]
    if(na.attr) attr(res, "na.action") <- `oldClass<-`(whichv(rl, FALSE), "omit")
  }
  res
}

na_insert <- function(X, prop = 0.1, value = NA, set = FALSE) {
  if(is.list(X)) {
    n <- fnrow(X)
    nmiss <- floor(n * prop)
    if(set) {
      lapply(unattrib(X), function(y) scv(y, sample.int(n, nmiss), value, TRUE))
      return(invisible(X))
    }
    res <- duplAttributes(lapply(unattrib(X), function(y) scv(y, sample.int(n, nmiss), value)), X)
    return(if(inherits(X, "data.table")) alc(res) else res)
  }
  if(!is.atomic(X)) stop("X must be an atomic vector, array or data.frame")
  l <- length(X)
  if(set) {
    scv(X, sample.int(l, floor(l * prop)), value, TRUE)
    return(invisible(X))
  }
  return(scv(X, sample.int(l, floor(l * prop)), value))
}

fdapply <- function(X, FUN, ...) duplAttributes(lapply(`attributes<-`(X, NULL), FUN, ...), X)

fnlevels <- function(x) length(attr(x, "levels"))

# flevels <- function(x) attr(x, "levels")

fnrow <- function(X) .Call(C_fnrow, X)  # if(is.list(X)) length(.subset2(X, 1L)) else dim(X)[1L]

fncol <- function(X) if(is.list(X)) length(unclass(X)) else dim(X)[2L]

fNCOL <- function(X) if(is.list(X)) length(unclass(X)) else NCOL(X)

fdim <- function(X) {
   if(is.atomic(X)) return(dim(X)) # or if !is.list ?
   c(.Call(C_fnrow, X), length(unclass(X)))
}

seq_row <- function(X) seq_len(.Call(C_fnrow, X))

seq_col <- function(X) if(is.list(X)) seq_along(unclass(X)) else seq_len(dim(X)[2L])

# na.last = TRUE, same default as order():
forder.int <- function(x, na.last = TRUE, decreasing = FALSE) .Call(C_radixsort, na.last, decreasing, FALSE, FALSE, TRUE, pairlist(x)) # if(is.unsorted(x)) .Call(C_forder, x, NULL, FALSE, TRUE, 1L, TRUE) else seq_along(x) # since forder gives integer(0) if sorted !

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

as_integer_factor <- function(X, keep.attr = TRUE) {
  if(is.atomic(X)) if(keep.attr) return(ffka(X, as.integer)) else
    return(as.integer(attr(X, "levels"))[X])
  res <- duplAttributes(lapply(unattrib(X),
    if(keep.attr) (function(y) if(is.factor(y)) ffka(y, as.integer) else y) else
                  (function(y) if(is.factor(y)) as.integer(attr(y, "levels"))[y] else y)), X)
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

# as.numeric_factor <- function(X, keep.attr = TRUE) {
#   .Deprecated(msg = "'as.numeric_factor' was renamed to 'as_numeric_factor'. It will be removed end of 2023, see help('collapse-renamed').")
#   as_numeric_factor(X, keep.attr)
# }
#
# as.character_factor <- function(X, keep.attr = TRUE) {
#   .Deprecated(msg = "'as.character_factor' was renamed to 'as_character_factor'. It will be removed end of 2023, see help('collapse-renamed').")
#   as_character_factor(X, keep.attr)
# }


setRnDF <- function(df, nm) `attr<-`(df, "row.names", nm)

# TtI <- function(x)
#   switch(x, replace_fill = 1L, replace = 2L, `-` = 3L, `-+` = 4L, `/` = 5L, `%` = 6L, `+` = 7L, `*` = 8L, `%%` = 9L, `-%%` = 10L,
#             stop("Unknown transformation!"))

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

do_stub <- function(stub, nam, default) {
  if(is.character(stub)) return(paste0(stub, nam))
  if(isTRUE(stub)) paste0(default, nam) else nam
}
# give_nam <- function(x, gn, stub) {
#   if(!gn) return(x)
#   attr(x, "names") <- paste0(stub, attr(x, "names"))
#   x
# }

fmatch <- function(x, table, nomatch = NA_integer_, count = FALSE, overid = 1L) .Call(C_fmatch, x, table, nomatch, count, overid)
ckmatch <- function(x, table, e = "Unknown columns:", ...) if(anyNA(m <- fmatch(x, table, NA_integer_, ...))) stop(paste(e, if(is.list(x)) paste(c("\n", capture.output(ss(x, is.na(m)))), collapse = "\n") else paste(x[is.na(m)], collapse = ", "))) else m
"%fin%" <- function(x, table) as.logical(fmatch(x, table, 0L, overid = 2L)) # export through set_collapse(mask = "%in%")
"%!in%" <- function(x, table) is.na(fmatch(x, table, overid = 2L))
"%!iin%" <- function(x, table) whichNA(fmatch(x, table, overid = 2L))
"%iin%" <- function(x, table) whichNA(fmatch(x, table, overid = 2L), invert = TRUE)
# anyNAerror <- function(x, e) if(anyNA(x)) stop(e) else x

cols2int <- function(cols, x, nam, topos = TRUE) {
 if(is.numeric(cols)) {
   if(length(cols) == 0L) return(integer(0L))
   l <- length(unclass(x)) # length(nam) ?
   if(cols[1L] < 0L) { # This is sufficient to check negative indices: No R function allows subsetting mixing positive and negative indices.
     if(-bmin(cols) > l) stop("Index out of range abs(1:length(x))")
     if(topos) return(seq_len(l)[cols])
     # cols <- seq_len(l)[cols]
     # if(!length(cols) || anyNA(cols)) stop("Index out of range abs(1:length(x))") -> used to put earlier check after if(topos) and use this one instead. But turns out that doesn't always work well.
     # return(cols)
   } else if(bmax(cols) > l) stop("Index out of range abs(1:length(x))")
   # if(bmax(abs(cols)) > length(unclass(x))) stop("Index out of range abs(1:length(x))") # Before collapse 1.4.0 !
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

# Needed for fmutate
cols2char <- function(cols, x, nam) {
  if(is.character(cols)) return(cols)
  if(!length(cols)) return("") # Needed if NULL is passed
  if(is.numeric(cols)) {
    l <- length(nam)
    if(cols[1L] < 0L) {
      if(-bmin(cols) > l) stop("Index out of range abs(1:length(x))")
    } else if(bmax(cols) > l) stop("Index out of range abs(1:length(x))")
    return(nam[cols])
  }
  if(is.function(cols)) return(nam[vapply(unattrib(x), cols, TRUE)])
  if(is.logical(cols)) {
    if(length(cols) != length(nam)) stop("Logical subsetting vector must match columns!")
    return(nam[cols])
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
#     if(bmax(abs(cols)) > lx) stop("Index out of range abs(1:length(x))")
#     r[cols] <- TRUE
#   } else stop("cols must be a function, character vector, numeric indices or logical vector!")
#   r
# }

# Helper for operator functions...
cols2intrmgn <- function(gn, cols, x) {
  if(is.function(cols)) {
    cols <- if(identical(cols, is.numeric)) .Call(C_vtypes, x, 1L) else vapply(unattrib(x), cols, TRUE)
    cols[gn] <- FALSE
    return(which(cols))
  }
  if(is.null(cols)) return(seq_along(unclass(x))[-gn])
  if(is.numeric(cols) && length(cols) && cols[1L] < 0) {
    res <- logical(length(unclass(x)))
    res[cols] <- TRUE
    res[gn] <- FALSE
    return(which(res))
  }
  cols2int(cols, x, attr(x, "names"), FALSE)
}

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
#     if(bmax(abs(ind)) > length(x)) stop("Index out of range abs(1:length(x))")
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

rgrep <- function(exp, nam, ..., sort = TRUE) if(length(exp) == 1L) grep(exp, nam, ...) else funique.default(unlist(lapply(exp, grep, nam, ...), use.names = FALSE), sort)
rgrepl <- function(exp, nam, ...) if(length(exp) == 1L) grepl(exp, nam, ...) else Reduce(`|`, lapply(exp, grepl, nam, ...))

fanyDuplicated <- function(x) if(length(x) < 100L) anyDuplicated.default(x) > 0L else .Call(C_fndistinct,x,NULL,FALSE,1L) != length(x)

# NROW2 <- function(x, d) if(length(d)) d[1L] else length(x)
# NCOL2 <- function(d, ilv) if(ilv) d[2L] else 1L

issorted <- function(x, strictly = FALSE) .Call(C_issorted, x, strictly)

charorNULL <- function(x) if(is.character(x)) x else NULL

tochar <- function(x) if(is.character(x)) x else as.character(x)  # if(is.object(x)) as.character(x) else .Call(C_aschar, x)

# dotstostr <- function(...) {
#   args <- deparse(substitute(c(...)))
#   nc <- nchar(args)
#   substr(args, 2, nc) # 3, nc-1 for no brackets !
# }

switch_msg <- function(msg, which = NULL) {
  if(is.null(which)) stop(msg)
  switch(which, error = stop(msg), message = message(msg), warning = warning(msg))
}

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
    .Call(C_setcopyv, x, NA_integer_, length(lev) + 1L, FALSE, TRUE, FALSE) # x[is.na(x)] <- length(lev) + 1L
  } else .Call(C_setcopyv, x, NA_integer_, length(lev), FALSE, TRUE, FALSE) # x[is.na(x)] <- length(lev)
  oldClass(x) <- clx
  x
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
