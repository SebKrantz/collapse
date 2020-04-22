# library(Rcpp)
# sourceCpp("C++/small_helper.cpp")
null2NA <- function(x) if(is.null(x)) NA_character_ else x
# flapply <- function(x, FUN, ...) lapply(unattrib(x), FUN, ...) # not really needed ...
vlabels <- function(X, attrn = "label") {
  if(is.atomic(X)) return(null2NA(attr(X, attrn)))
  res <- lapply(unattrib(X), attr, attrn)
  res[vapply(res, is.null, TRUE)] <- NA_character_
  unlist(res)
}
"vlabels<-" <- function(X, attrn = "label", value) {
  if(is.atomic(X)) return(`attr<-`(X, attrn, value))
  clx <- class(X)
  class(X) <- NULL
  for (i in seq_along(value)) attr(X[[i]], attrn) <- value[i]
  `oldClass<-`(X, clx)
}
strclp <- function(x) if(length(x) > 1L) paste(x, collapse = " ") else x
pasteclass <- function(x) if(length(cx <- class(x)) > 1L) paste(cx, collapse = " ") else cx # Faster if length(class(x)) == 1L
vclasses <- function(X) {
  if(is.atomic(X)) return(pasteclass(X))
  vapply(unattrib(X), pasteclass, character(1))
}
vtypes <- function(X) {
  if(is.atomic(X)) return(typeof(X))
  vapply(unattrib(X), typeof, character(1))
}
namlab <- function(X, class = FALSE, attrn = "label") {
  if(!is.list(X)) stop("namlab only works with lists")
  res <- if(class) list(attr(X, "names"), vapply(unattrib(X), pasteclass, character(1)), vlabels(X, attrn)) else list(attr(X, "names"), vlabels(X, attrn))
  attributes(res) <- list(names = if(class) c("Variable","Class","Label") else c("Variable","Label"),
                          row.names = .set_row_names(length(unclass(X))),
                          class = "data.frame")
  return(res)
}
add_stub <- function(X, stub, pre = TRUE) {
  if(!is.character(stub)) return(X)
  if(is.array(X)) {
    if(length(dim(X)) > 2L) stop("Can't stub higher dimensional arrays!")
    dn <- dimnames(X)
    dimnames(X) <- list(dn[[1L]], if(pre) paste0(stub, dn[[2L]]) else paste0(dn[[2L]], stub))
  } else attr(X, "names") <- if(pre) paste0(stub, attr(X, "names")) else paste0(attr(X, "names"), stub)
  X
}
rm_stub <- function(X, stub, ...) {
  if(!is.character(stub)) return(X)
  if(is.array(X)) {
    if(length(dim(X)) > 2L) stop("Can't remove stub from higher dimensional arrays!")
    dn <- dimnames(X)
    dimnames(X) <- list(dn[[1L]], sub(stub, "", dn[[2L]], ...)) # use ^ or $ to restrict to pre or port matching.
  } else attr(X, "names") <- sub(stub, "", attr(X, "names"), ...)
  X
}

setRownames <- function(object = nm, nm = seq_row(object)) `rownames<-`(object, nm)
setColnames <- function(object = nm, nm) `colnames<-`(object, nm)
setDimnames <- function(object = dn, dn) `dimnames<-`(object, dn)
all_identical <- function(...) {
  if(length(list(...)) == 1L && is.list(...)) { # if(length(match.call())-1L == 1L && is.list(...)) # https://stackoverflow.com/questions/44011918/count-number-of-arguments-passed-to-function
    all(unlist(lapply(...[-1L], identical, ...[[1L]]), use.names = FALSE)) # use vapply ??
  } else {
    l <- list(...)
    all(unlist(lapply(l[-1L], identical, l[[1L]]), use.names = FALSE)) # use vapply ??
  }
}
all_obj_equal <- function(...) {
  if(length(list(...)) == 1L && is.list(...)) { # if(length(match.call())-1L == 1L && is.list(...)) # https://stackoverflow.com/questions/44011918/count-number-of-arguments-passed-to-function
    all(unlist(lapply(...[-1L], all.equal, ...[[1L]]), use.names = FALSE)) # use vapply ??
  } else {
    l <- list(...)
    all(unlist(lapply(l[-1L], all.equal, l[[1L]]), use.names = FALSE)) # use vapply ??
  }
}

finteraction <- function(...) { # does it drop levels ?? -> Yes !!!
  ll <- length(list(...))
  if(ll == 1L && is.list(...)) return(as.factor.GRP(GRP.default(...)))
  as.factor.GRP(GRP.default(list(...)))
}

interact_names <- function(l) do.call(paste, c(expand.grid(l, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE), list(sep = ".")))
unattrib <- function(object) `attributes<-`(object, NULL)
is.categorical <- function(x) !is.numeric(x)
is.Date <- function(x) inherits(x, c("Date","POSIXlt","POSIXct"))
"%!in%" <- function(x, table) match(x, table, nomatch = 0L) == 0L
na_rm <- function(x) x[!is.na(x)] # more consistent with base than na_rm !!! if not Cpp version that's fine !!
# na.rm <- function(x) { # cpp version available, but not faster !!
#   if(!is.null(attr(x, "names"))) { # gives corruped time-series !!
#     ax <- attributes(x)
#     r <- x[!is.na(x)]
#     ax[["names"]] <- names(r)
#     setAttributes(r, ax)
#   } else duplAttributes(x[!is.na(x)], x)
# }

# fast na.omit !!
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
  if(na.attr) attr(res, "na.action") <- `class<-`(which(!rl), "omit")
  return(res)
}
na_insert <- function(X, prop = 0.1) {
  if(!is.null(d <- dim(X))) {
    n <- d[1L]
    p <- d[2L]
    NAloc <- rep(FALSE, n * p)
    NAloc[sample.int(n * p, floor(n * p * prop))] <- TRUE
    X[matrix(NAloc, nrow = n, ncol = p)] <- NA
  } else if(is.atomic(X)) {
    l <- length(X)
    X[sample.int(l, floor(l * prop))] <- NA
  } else stop("X must be an atomic vector, matrix or data.frame")
  return(X)
}
fnlevels <- function(x) length(attr(x, "levels")) # make cpp version ?? -> nope, slower !!
# flevels <- function(x) attr(x, "levels")
fnrow <- function(x) if(is.list(x)) length(unclass(x)[[1L]]) else dim(x)[1L]
fnrow2 <- function(x) length(unclass(x)[[1L]])
fncol <- function(x) if(is.list(x)) length(unclass(x)) else dim(x)[2L]
fdim <- function(x) {
   if(is.atomic(x)) return(dim(x)) # or if !is.list ???
   class(x) <- NULL
   c(length(x[[1L]]), length(x))
}
seq_row <- function(X) if(is.list(X)) seq_along(unclass(X)[[1L]]) else seq_len(nrow(X))
seq_col <- function(X) if(is.list(X)) seq_along(unclass(X)) else seq_len(ncol(X))

forder.int <- function(x) if(is.unsorted(x)) .Call(C_forder, x, NULL, FALSE, TRUE, 1L, TRUE) else seq_along(x) # since forder gives integer(0) if sorted !!
fsetdiff <- function(x, y) x[match(x, y, 0L) == 0L] # not unique !!
as.numeric_factor <- function(X) {
  if(is.atomic(X)) return(as.numeric(attr(X, "levels"))[X])
    fcts <- vapply(unattrib(X), is.factor, TRUE)
    # if(all(fcts)) return(dapply(X, function(x) as.numeric(attr(x, "levels"))[x]))
    clx <- class(X)
    class(X) <- NULL
    X[fcts] <- lapply(X[fcts], function(x) as.numeric(attr(x, "levels"))[x])
    return(`oldClass<-`(X, clx))
}
as.character_factor <- function(X) {
  if(is.atomic(X)) return(as.character(attr(X, "levels"))[X])
  fcts <- vapply(unattrib(X), is.factor, TRUE)
  clx <- class(X)
  class(X) <- NULL
  X[fcts] <- lapply(X[fcts], function(x) as.character(attr(x, "levels"))[x])
  return(`oldClass<-`(X, clx))
}

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


setRow.names <- function(df, nm) `attr<-`(df, "row.names", nm)
addAttributes <- function(x, a) .Call(Cpp_setAttributes, x, c(attributes(x), a))

TRAtoInt <- function(x) # A lot faster than match based verion !!!
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
anyNAerror <- function(x, e) if(anyNA(x)) stop(e) else x
ckmatch <- function(x, table, e = "Unknown columns:") if(anyNA(m <- match(x, table))) stop(paste(e, paste(x[is.na(m)], collapse = ", "))) else m
cols2int <- function(cols, x, nam) {
 if(is.numeric(cols)) {
  if(max(abs(cols)) > length(unclass(x))) stop("Index out of range abs(1:length(x))")
  return(cols)
 } else if(is.function(cols))
  return(which(vapply(unattrib(x), cols, TRUE))) else if(is.character(cols))
  return(ckmatch(cols, nam)) else if(is.logical(cols)) {
    if(length(cols) != length(unclass(x))) stop("Logical subsetting vector must match columns!")
    return(which(cols))
  } else stop("cols must be a function, character vector, numeric indices or logical vector!")
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
  return(r)
}
colsubset <- function(x, ind) { # also works for grouped tibbles !!
  ax <- attributes(x)
  attributes(x) <- NULL #  and good here since vapply without attributes is faster...
  if(is.numeric(ind)) { # faster using switch(mode(ind), ...) ??
    if(max(abs(ind)) > length(x)) stop("Index out of range abs(1:length(x))")
  } else if(is.logical(ind)) {
    if(length(ind) != length(x)) stop("Logical subsetting vector must match length(x)")
  } else ind <- if(is.character(ind)) ckmatch(ind, ax[["names"]]) else
                vapply(x, ind, TRUE)
  ax[["names"]] <- ax[["names"]][ind]
  .Call(Cpp_setAttributes, x[ind], ax) # return(`attributes<-`(x[ind], ax)) # This is slow on large data -> a lot of checks !!!
}

fcolsubset <- function(x, ind) {
  ax <- attributes(x)
  attributes(x) <- NULL
  ax[["names"]] <- ax[["names"]][ind]
  .Call(Cpp_setAttributes, x[ind], ax)
}

# fcolsubset.int <- function(x, ind) { # same speed,, perhaps a bit slower... (the above is faster !!)
#   ax <- attributes(x)
#   class(x) <- NULL
#   ax[["names"]] <- names(x)[ind]
#   return(.Call(Cpp_setAttributes, x[ind], ax))
# }
# fcolsubset.int2 <- function(x, ind) { # same speed, all really neglible differences...
#   ax <- attributes(x)
#   ax[["names"]] <- ax[["names"]][ind]
#   return(.Call(Cpp_setAttributes, unclass(x)[ind], ax))
# }

# Fastest because vapply runs faster on a list without any attributes !!!
colsubsetFUN <- function(x, FUN) {
  ax <- attributes(x)
  attributes(x) <- NULL
  ind <- vapply(x, FUN, TRUE)
  ax[["names"]] <- ax[["names"]][ind]
  .Call(Cpp_setAttributes, x[ind], ax)
}


at2GRP <- function(x) {
  if(is.nmfactor(x))
  return(list(length(attr(x, "levels")), x, NULL)) else {
    res <- list(NULL, NULL, NULL)
    res[[2L]] <- qG(x, ordered = FALSE, na.exclude = FALSE)
    res[[1L]] <- attr(res[[2L]], "N.groups")
    return(res)
  }
}
# ret2int <- function(x) match(ret, c("cols","all","add")) # return = c("cols","all","add")
G_t <- function(x, m = TRUE, wm = 1L) {
  if(is.null(x)) {
    if(m) message(switch(wm, "Panel-lag computed without timevar: Assuming ordered data",
                             "Panel-difference computed without timevar: Assuming ordered data",
                             "Panel-growth rate computed without timevar: Assuming ordered data"))
    return(x)
  } else if(is.atomic(x)) {
    if(is.integer(x)) return(x) else return(qG(x, na.exclude = FALSE)) # make sure it is ordered !!! qG already ckecks factor !!
  } else if(is.GRP(x)) return(x[[2L]]) else return(GRP.default(x, return.groups = FALSE)[[2L]])
}
rgrep <- function(exp, nam, ...) if(length(exp) == 1L) grep(exp, nam, ...) else funique(unlist(lapply(exp, grep, nam, ...), use.names = FALSE))
# NROW2 <- function(x, d) if(length(d)) d[1L] else length(x)
# NCOL2 <- function(d, ilv) if(ilv) d[2L] else 1L
charorNULL <- function(x) if(is.character(x)) x else NULL
# more security here??
# unique_factor <- function(x) {  # Still needed with new collap solution ?? -> Nope !!
#   res <- seq_along(attr(x, "levels"))
#   .Call(Cpp_duplAttributes, res, x)
# }
# dotstostr <- function(...) {
#   args <- deparse(substitute(c(...)))
#   nc <- nchar(args)
#   substr(args, 2, nc) # 3, nc-1 for no brackets !!
# }

unused_arg_warning <- function(call, ...) {
  args <- deparse(substitute(c(...)))
  nc <- nchar(args)
  args <- substr(args, 2, nc) # 3, nc-1 for no brackets !
  warning(paste("Unused arguments",args,"passed to", as.character(call[[1L]])))
}

is.nmfactor <- function(x) inherits(x, "factor") && (inherits(x, "na.included") || !anyNA(unclass(x)))
addNA2 <- function(x) {
  if(!anyNA(unclass(x))) return(x)
  clx <- class(x)
  class(x) <- NULL
  if(!anyNA(lev <- attr(x, "levels"))) {
    attr(x, "levels") <- c(lev, NA_character_)
    x[is.na(x)] <- length(lev) + 1L
  } else x[is.na(x)] <- length(lev)
  return(`oldClass<-`(x, clx))
}

l1orn <- function(x, nam) if(length(x) == 1L) x else nam

fsimplify2array <- function(l) {
  res <- do.call(cbind, l) # lapply(l, `dimnames<-`, NULL) # also faster than unlist..
  dim(res) <- c(dim(l[[1L]]), length(l))
  dimnames(res) <- c(dimnames(l[[1L]]), list(names(l)))
  res
}


# addNA2 <- function(x) {
#   clx <- c(class(x), "na.included")
#   if(!anyNA(unclass(x))) return(`oldClass<-`(x, clx))
#   ll <- attr(x, "levels")
#   if(!anyNA(ll)) ll <- c(ll, NA)
#   return(`oldClass<-`(factor(x, levels = ll, exclude = NULL), clx))
# }
