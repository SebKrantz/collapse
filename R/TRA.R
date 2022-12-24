

TRA <- function(x, STATS, FUN = "-", ...) UseMethod("TRA") # , x

setTRA <- function(x, STATS, FUN = "-", ...) invisible(TRA(x, STATS, FUN, ..., set = TRUE))

TRA.default <- function(x, STATS, FUN = "-", g = NULL, set = FALSE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(TRA.matrix(x, STATS, FUN, g, set, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(C_TRA,x,STATS,0L,FUN,set))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) {
      if(fnlevels(g) != length(STATS)) stop("number of groups must match length(STATS)")
    } else {
      g <- qG(g, na.exclude = FALSE) # needs to be ordered to be compatible with fast functions !!
      if(attr(g, "N.groups") != length(STATS)) stop("number of groups must match length(STATS)")
    }
    return(.Call(C_TRA,x,STATS,g,FUN,set))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  if(g[[1L]] != length(STATS)) stop("number of groups must match length(STATS)")
  .Call(C_TRA,x,STATS,g[[2L]],FUN,set)
}

TRA.matrix <- function(x, STATS, FUN = "-", g = NULL, set = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(C_TRAm,x,STATS,0L,FUN,set))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) {
      if(fnlevels(g) != nrow(STATS)) stop("number of groups must match nrow(STATS)")
    } else {
      g <- qG(g, na.exclude = FALSE) # needs to be ordered to be compatible with fast functions !!
      if(attr(g, "N.groups") != nrow(STATS)) stop("number of groups must match nrow(STATS)")
    }
    return(.Call(C_TRAm,x,STATS,g,FUN,set))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  if(g[[1L]] != nrow(STATS)) stop("number of groups must match nrow(STATS)")
  .Call(C_TRAm,x,STATS,g[[2L]],FUN,set)
}

TRA.data.frame <- function(x, STATS, FUN = "-", g = NULL, set = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(C_TRAl,x,STATS,0L,FUN,set))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) {
      if(fnlevels(g) != fnrow(STATS)) stop("number of groups must match nrow(STATS)")
    } else {
      g <- qG(g, na.exclude = FALSE) # needs to be ordered to be compatible with fast functions !!
      if(attr(g, "N.groups") != fnrow(STATS)) stop("number of groups must match nrow(STATS)")
    }
    return(.Call(C_TRAl,x,STATS,g,FUN,set))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  if(g[[1L]] != fnrow(STATS)) stop("number of groups must match nrow(STATS)")
  .Call(C_TRAl,x,STATS,g[[2L]],FUN,set)
}

TRA.list <- function(x, ...) TRA.data.frame(x, ...)

TRA.grouped_df <- function(x, STATS, FUN = "-", keep.group_vars = TRUE, set = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  clx <- oldClass(x)
  oldClass(x) <- NULL
  oldClass(STATS) <- NULL
  if(g[[1L]] != length(STATS[[1L]])) stop("number of groups must match nrow(STATS)")
  nognst <- names(STATS) %!in% g[[5L]]
  mt <- ckmatch(names(STATS), names(x), "Variables in STATS not found in x:")
  mt <- mt[nognst]
  x[mt] <- .Call(C_TRAl,x[mt],STATS[nognst],g[[2L]],FUN,set)
  if(!keep.group_vars) x[names(x) %in% g[[5L]]] <- NULL
  oldClass(x) <- clx
  x
}
