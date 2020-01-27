# library(Rcpp)
# sourceCpp("C++/small_helper.cpp")

# Export --------------------------------------
vlabels <- function(X, attrn = "label") {
  if(is.atomic(X)) {
    res <- attr(X, attrn)
    if(is.null(res)) NA else res
  } else {
    res <- lapply(X, attr, attrn)
    res[vapply(res, is.null, TRUE)] <- NA
    unlist(res)
  }
}
"vlabels<-" <- function(X, attrn = "label", value) {
  if(is.atomic(X)) {
    attr(X, attrn) <- value
  } else {
    for (i in seq_along(value)) attr(X[[i]], attrn) <- value[i]
  }
  X
}
namlab <- function(X, class = FALSE, attrn = "label") {
  pasteclass <- function(x) paste(class(x), collapse = " ")
  res <- if(class) list(names(X), vapply(X, pasteclass, character(1)), vlabels(X, attrn)) else list(names(X), vlabels(X, attrn))
  attributes(res) <- list(names = if(class) c("Variable","Class","Label") else c("Variable","Label"),
                          row.names = .set_row_names(length(X)),
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
seq_row <- function(X) seq_len(nrow(X))
seq_col <- function(X) seq_len(ncol(X))
setRownames <- function(object = nm, nm = seq_row(object)) {
  rownames(object) <- nm
  object
}
setColnames <- function(object = nm, nm) {
  colnames(object) <- nm
  object
}
setDimnames <- function(object = dn, dn) {
  dimnames(object) <- dn
  object
}
pwcor <- function(X, ...) `class<-`(cor(X, ..., use = "pairwise.complete.obs"), c("pwcor","matrix"))
pwcov <- function(X, ...) `class<-`(cov(X, ..., use = "pairwise.complete.obs"), c("pwcov","matrix"))
print.pwcov <- print.qsu
print.pwcor <- function(x, digits = 2) print.table(dapply(round(x, digits), function(j) sub("^(-?)0.", "\\1.", j)), right = TRUE) # print.table(, right = TRUE)

all_identical <- function(...) {
  if(length(match.call())-1L == 1L && is.list(...)) { # https://stackoverflow.com/questions/44011918/count-number-of-arguments-passed-to-function
    all(unlist(lapply(...[-1L], identical, ...[[1L]]), use.names = FALSE)) # use vapply ??
  } else {
    l <- list(...)
    all(unlist(lapply(l[-1L], identical, l[[1L]]), use.names = FALSE)) # use vapply ??
  }
}
all.identical <- all_identical
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
as.numeric_factor <- function(X) {
  if(is.list(X)) {
    fcts <- vapply(X, is.factor, TRUE, USE.NAMES = FALSE)
    get_vars(X, fcts) <- lapply(colsubset(X, fcts), function(x) as.numeric(attr(x, "levels"))[x])
    return(X)
  } else return(as.numeric(attr(X, "levels"))[X])
}




setRow.names <- function(df, nm) {
  attr(df, "row.names") <- nm
  df
}
remove_attributes <- function(x) {
  attributes(x) <- NULL
  x
}
addAttributes <- function(x, a) {
  ax <- attributes(x)
  .Call(Cpp_setAttributes, x, c(ax, a))
}
TRAtoInt <- function(x) # A lot faster than match based verion !!!
  switch(x, replace_fill = 1L, replace = 2L, `-` = 3L, `-+` = 4L, `/` = 5L, `%` = 6L, `+` = 7L, `*` = 8L,
            stop("Unknown transformation!"))
condsetn <- function(x, value, cond) {
  if(cond) attr(x, "names") <- value
  x
}
give_nam <- function(x, gn, stub) {
  if(!gn) return(x)
  attr(x, "names") <- paste0(stub, attr(x, "names"))
  x
}
anyNAerror <- function(x, e) if(anyNA(x)) stop(e) else x
cols2int <- function(cols, x, nam) {
 if(is.numeric(cols)) {
  if(max(abs(cols)) > length(x)) stop("Index out of range abs(1:length(x))")
  return(cols)
 } else if(is.function(cols))
  return(which(vapply(x, cols, TRUE))) else if(is.character(cols))
  return(anyNAerror(match(cols, nam), "Unknown column names!")) else if(is.logical(cols)) {
    if(length(cols) != length(x)) stop("Logical subsetting vector must match columns!")
    return(which(cols))
  } else stop("cols must be a function, character vector, numeric indices or logical vector!")
}
cols2log <- function(cols, x, nam) {
  if(is.logical(cols)) if(length(cols) == length(x)) return(cols) else stop("Logical subsetting vector must match columns!")
  if(is.function(cols)) return(vapply(x, cols, TRUE))
  r <- logical(length(x))
  if(is.character(cols)) {
    r[anyNAerror(match(cols, nam), "Unknown column names!")] <- TRUE
  } else if(is.numeric(cols)) {
    if(max(abs(cols)) > length(r)) stop("Index out of range abs(1:length(x))")
    r[cols] <- TRUE
  } else stop("cols must be a function, character vector, numeric indices or logical vector!")
  return(r)
}
colsubset <- function(x, ind) { # also works for grouped tibbles !!
  ax <- attributes(x)
  attributes(x) <- NULL # faster than unclass !! and good here since vapply on unclassed is faster !!
  if(is.numeric(ind)) {
    if(max(abs(ind)) > length(x)) stop("Index out of range abs(1:length(x))")
  } else if(is.logical(ind)) {
    if(length(ind) != length(x)) stop("Logical subsetting vector must match length(x)")
  } else {
    ind <- if(is.function(ind)) vapply(x, ind, TRUE, USE.NAMES = FALSE) else
      anyNAerror(match(ind, ax[["names"]]), "Unknown column names!")
  }
  ax[["names"]] <- ax[["names"]][ind]
  return(.Call(Cpp_setAttributes, x[ind], ax)) # return(`attributes<-`(x[ind], ax)) # This is slow on large data -> a lot of checks !!!
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
G_t <- function(x, m = TRUE) {
  if(is.null(x)) {
    if(m) message("Panel-lag computed without timevar: Assuming ordered data")
    return(x)
  } else if(is.atomic(x)) {
    if(is.integer(x)) return(x) else return(qG(x, na.exclude = FALSE)) # make sure it is ordered !!! qG already ckecks factor !!
  } else if(is.GRP(x)) return(x[[2L]]) else return(GRP(x, return.groups = FALSE)[[2L]])
}
rgrep <- function(exp, nam, ...) if(length(exp) > 1L) funique(unlist(lapply(exp, grep, nam, ...), use.names = FALSE)) else grep(exp, nam, ...)
NROW2 <- function(x, d) if(length(d)) d[1L] else length(x)
NCOL2 <- function(d, ilv) if(ilv) d[2L] else 1L
charorNULL <- function(x) if(is.character(x)) x else NULL
# more security here??
unique_factor <- function(x) {
  res <- seq_along(attr(x, "levels"))
  .Call(Cpp_duplAttributes, res, x)
}
dotstostr <- function(...) {
  args <- deparse(substitute(c(...)))
  nc <- nchar(args)
  substr(args, 2, nc) # 3, nc-1 for no brackets !!
}
is.nmfactor <- function(x) inherits(x, "factor") && (inherits(x, "na.included") || !anyNA(unclass(x)))
addNA2 <- function(x) {
  if(!anyNA(unclass(x))) return(x)
  ax <- attributes(x)
  if(!anyNA(ax[["levels"]])) ax[["levels"]] <- c(ax[["levels"]], NA)
  attributes(x) <- NULL
  x[is.na(x)] <- length(ax[["levels"]])
  return(setAttributes(x, ax))
}

# addNA2 <- function(x) {
#   clx <- c(class(x), "na.included")
#   if(!anyNA(unclass(x))) return(`oldClass<-`(x, clx))
#   ll <- attr(x, "levels")
#   if(!anyNA(ll)) ll <- c(ll, NA)
#   return(`oldClass<-`(factor(x, levels = ll, exclude = NULL), clx))
# }
