library(Rcpp)
sourceCpp("C++/small_helper.cpp")
vlabels <- function(X) {
  if(is.atomic(X)) {
    res = attr(X, "label")
    if(is.null(res)) NA else res
  } else {
    res = lapply(X, attr, "label")
    res[vapply(res, is.null, TRUE)] = NA
    unlist(res)
  }
}
"vlabels<-" <- function(X, value) {
  if(is.atomic(X)) {
    attr(X, "label") = value
  } else {
    for (i in seq_along(value)) attr(X[[i]],"label") = value[i]
  }
  X
}
namlab <- function(X, class = FALSE) {
  res <- if(class) list(names(X), vapply(X, class, character(1)), vlabels(X)) else list(names(X), vlabels(X))
  attributes(res) <- list(names = if(class) c("Variable","Class","Label") else c("Variable","Label"),
                          row.names = .set_row_names(length(X)),
                          class = "data.frame")
  return(res)
}
seq_row <- function(x) seq_len(nrow(x))
seq_col <- function(x) seq_len(ncol(x))
setRow.names <- function(df, nm) {
  attr(df, "row.names") <- nm
  df
}
setRownames <- function(object = nm, nm = NULL) {
  if (is.null(nm)) nm <- seq_row(object)
  rownames(object) <- nm
  object
}
setColnames <- function(object = nm, nm) {
  colnames(object) <- nm
  object
}
setDimnames <- function(object = nm, nm) {
  dimnames(object) <- nm
  object
}
remove_attributes <- function(x) {
  attributes(x) <- NULL
  x
}
pwcor <- function(x, ...) cor(x, use = "pairwise.complete.obs", ...)
pwcov <- function(x, ...) cov(x, use = "pairwise.complete.obs", ...)
na.rm <- function(x) x[!is.na(x)] # cpp version available !!
fnlevels <- function(x) length(attr(x, "levels")) # make cpp version ??
TRAtoInt <- function(x) # A lot faster than match based verion !!!
  switch(x, replace_fill = 1L, replace = 2L, `-` = 3L, `-+` = 4L, `/` = 5L, `%` = 6L, `+` = 7L, `*` = 8L,
            stop("Unknown transformation!"))
all.identical <- function(...) {
  if(length(match.call())-1 == 1 && is.list(...)) { # https://stackoverflow.com/questions/44011918/count-number-of-arguments-passed-to-function
    all(unlist(lapply(...[-1],identical,...[[1]]), use.names = FALSE)) # use vapply ??
  } else {
    l <- list(...)
    all(unlist(lapply(l[-1],identical,l[[1]]), use.names = FALSE)) # use vapply ??
  }
}
give_nam <- function(x, gn, stub) {
  if(!gn) return(x)
  attr(x, "names") <- paste0(stub, attr(x, "names"))
  x
}
cols2int <- function(cols, x, nam) {
  if(is.function(cols)) which(vapply(x, cols, TRUE)) else if(is.character(cols))
    match(cols, nam) else cols
}
at2GRP <- function(x) {
  if(is.factor(x))
  return(list(length(attr(x, "levels")), x, NULL)) else {
    res <- list(NULL, NULL, NULL)
    res[[2]] <- qG(x, ordered = FALSE)
    res[[1]] <- attr(res[[2]], "N.groups")
    return(res)
  }
}
# ret2int <- function(x) match(ret, c("cols","all","add")) # return = c("cols","all","add")
G_t <- function(x, m = TRUE) {
  if(is.null(x)) {
    if(m) message("Panel-lag computed without timevar: Assuming ordered data")
    return(x)
  } else if(is.atomic(x)) {
    if(is.integer(x)) return(x) else return(qG(x)) # make sure it is ordered !!! qG already ckecks factor !!
  } else if(all(class(l) == "GRP")) return(l[[2]]) else return(GRP(l, return.groups = FALSE)[[2]])
}
anyNAerror <- function(x, e) if(anyNA(x)) stop(e) else x
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
  return(setAttributes(x[ind], ax)) # return(`attributes<-`(x[ind], ax)) # This is slow on large data -> a lot of checks !!!
}
rgrep <- function(exp, nam, ...) if(length(exp) > 1L) sort.int(unique.default(vapply(exp, grep, 1L, nam, ...))) else grep(exp, nam, ...)
is.categorical <- function(x) !is.numeric(x)
is.Date <- function(x) inherits(x, c("Date","POSIXlt","POSIXct"))
NROW2 <- function(x, d) if(length(d)) d[1L] else length(x)
NCOL2 <- function(d, ilv) if(ilv) d[2L] else 1L
charorNULL <- function(x) if(is.character(x)) x else NULL
