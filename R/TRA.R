

TRA <- function(x, STATS, FUN = "-", ...) UseMethod("TRA") # , x

TRA.default <- function(x, STATS, FUN = "-", g = NULL, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(TRA.matrix(x, STATS, FUN, g, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_TRA,x,STATS,0L,TtI(FUN)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) {
      if(fnlevels(g) != length(STATS)) stop("number of groups must match length(STATS)")
    } else {
      g <- qG(g, na.exclude = FALSE) # needs to be ordered to be compatible with fast functions !!
      if(attr(g, "N.groups") != length(STATS)) stop("number of groups must match length(STATS)")
    }
    return(.Call(Cpp_TRA,x,STATS,g,TtI(FUN)))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  if(g[[1L]] != length(STATS)) stop("number of groups must match length(STATS)")
  .Call(Cpp_TRA,x,STATS,g[[2L]],TtI(FUN))
}

TRA.matrix <- function(x, STATS, FUN = "-", g = NULL, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_TRAm,x,STATS,0L,TtI(FUN)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) {
      if(fnlevels(g) != nrow(STATS)) stop("number of groups must match nrow(STATS)")
    } else {
      g <- qG(g, na.exclude = FALSE) # needs to be ordered to be compatible with fast functions !!
      if(attr(g, "N.groups") != nrow(STATS)) stop("number of groups must match nrow(STATS)")
    }
    return(.Call(Cpp_TRAm,x,STATS,g,TtI(FUN)))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  if(g[[1L]] != nrow(STATS)) stop("number of groups must match nrow(STATS)")
  .Call(Cpp_TRAm,x,STATS,g[[2L]],TtI(FUN))
}

TRA.data.frame <- function(x, STATS, FUN = "-", g = NULL, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_TRAl,x,STATS,0L,TtI(FUN)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) {
      if(fnlevels(g) != fnrow2(STATS)) stop("number of groups must match nrow(STATS)")
    } else {
      g <- qG(g, na.exclude = FALSE) # needs to be ordered to be compatible with fast functions !!
      if(attr(g, "N.groups") != fnrow2(STATS)) stop("number of groups must match nrow(STATS)")
    }
    return(.Call(Cpp_TRAl,x,STATS,g,TtI(FUN)))
  }
  if(!is_GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  if(g[[1L]] != fnrow2(STATS)) stop("number of groups must match nrow(STATS)")
  .Call(Cpp_TRAl,x,STATS,g[[2L]],TtI(FUN))
}

TRA.list <- function(x, STATS, FUN = "-", g = NULL, ...) TRA.data.frame(x, STATS, FUN, g, ...)

TRA.grouped_df <- function(x, STATS, FUN = "-", keep.group_vars = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  clx <- oldClass(x)
  oldClass(x) <- NULL
  oldClass(STATS) <- NULL
  if(g[[1L]] != length(STATS[[1L]])) stop("number of groups must match nrow(STATS)")
  nognst <- names(STATS) %!in% g[[5L]]
  mt <- ckmatch(names(STATS), names(x), "Variables in STATS not found in x:")
  mt <- mt[nognst]
  x[mt] <- .Call(Cpp_TRAl,x[mt],STATS[nognst],g[[2L]],TtI(FUN))
  if(!keep.group_vars) x[names(x) %in% g[[5L]]] <- NULL
  oldClass(x) <- clx
  x
}




# sourceCpp('R/C++/TRAset.cpp')
# sourceCpp('R/C++/TRAsetl.cpp')
# sourceCpp('R/C++/TRAseta.cpp')
#
# setTRA <- function(X, X_ag, g = 0L, trans = "replace", ...) {
#   UseMethod("setTRA", X)
# }
# setTRA.default <- function(X, X_ag, g = 0L, trans = "replace", ...) {
#   # if(!(is.atomic(X_ag) && is.null(dim(X_ag)))) stop("X_ag must be a vector") # Cpp already gives error !! matrix takes first element..
#   if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#   if(is.list(g)) setTRACpp(X, X_ag, g[[2]], trans) else setTRACpp(X, X_ag, g, trans)
# }
# setTRA.matrix <- function(X, X_ag, g = 0L, trans = "replace", ...) {
#   if(!is.atomic(X_ag)) stop("X_ag must be a vector or matrix")
#   if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#   if(is.list(g)) setTRAmCpp(X, X_ag, g[[2]], trans) else setTRAmCpp(X, X_ag, g, trans)
# }
# setTRA.data.frame <- function(X, X_ag, g = 0L, trans = "replace", ...) {
#   if(is.array(X_ag)) stop("X_ag must be a vetor or list / data.frame")
#   if(is.list(X_ag) && length(X_ag[[1]]) == 1) X_ag <- unlist(X_ag, use.names = FALSE)
#   if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#   if(is.list(g)) setTRAlCpp(X, X_ag, g[[2]], trans) else setTRAlCpp(X, X_ag, g, trans)
# }
# setTRA.list <- function(X, X_ag, g = 0L, trans = "replace", ...) {
#   if(is.array(X_ag)) stop("X_ag must be a vetor or list / data.frame")
#   if(is.list(X_ag) && length(X_ag[[1]]) == 1) X_ag <- unlist(X_ag, use.names = FALSE)
#   if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#   if(is.list(g)) setTRAlCpp(X, X_ag, g[[2]], trans) else setTRAlCpp(X, X_ag, g, trans)
# }
