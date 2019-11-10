library(Rcpp)
sourceCpp('R/C++/TRA.cpp')
sourceCpp('R/C++/TRAl.cpp')
sourceCpp('R/C++/TRAa.cpp')

TRA <- function(X, X_ag, g = 0L, trans = "replace", ...) {
  UseMethod("TRA", X)
}
TRA.default <- function(X, X_ag, g = 0L, trans = "replace", ...) {
  # if(!(is.atomic(X_ag) && is.null(dim(X_ag)))) stop("X_ag must be a vector") # Cpp already gives error !! matrix takes first element.. 
  if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
  if(is.list(g)) TRACpp(X, X_ag, g[[2]], trans) else TRACpp(X, X_ag, g, trans)
}
TRA.matrix <- function(X, X_ag, g = 0L, trans = "replace", ...) {
  if(!is.atomic(X_ag)) stop("X_ag must be a vector or matrix")
  if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
  if(is.list(g)) TRAmCpp(X, X_ag, g[[2]], trans) else TRAmCpp(X, X_ag, g, trans)
}
TRA.data.frame <- function(X, X_ag, g = 0L, trans = "replace", ...) {
  if(is.array(X_ag)) stop("X_ag must be a vetor or list / data.frame")
  if(is.list(X_ag) && length(X_ag[[1]]) == 1) X_ag <- unlist(X_ag, use.names = FALSE)
  if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
  res <- if(is.list(g)) TRAlCpp(X, X_ag, g[[2]], trans) else TRAlCpp(X, X_ag, g, trans)
  attributes(res) <- attributes(X)
  res
}
TRA.list <- function(X, X_ag, g = 0L, trans = "replace", ...) {
  if(is.array(X_ag)) stop("X_ag must be a vetor or list / data.frame")
  if(is.list(X_ag) && length(X_ag[[1]]) == 1) X_ag <- unlist(X_ag, use.names = FALSE)
  if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
  res <- if(is.list(g)) TRAlCpp(X, X_ag, g[[2]], trans) else TRAlCpp(X, X_ag, g, trans)
  attributes(res) <- attributes(X)
  res
}


sourceCpp('R/C++/TRAset.cpp')
sourceCpp('R/C++/TRAsetl.cpp')
sourceCpp('R/C++/TRAseta.cpp')

setTRA <- function(X, X_ag, g = 0L, trans = "replace", ...) {
  UseMethod("setTRA", X)
}
setTRA.default <- function(X, X_ag, g = 0L, trans = "replace", ...) {
  # if(!(is.atomic(X_ag) && is.null(dim(X_ag)))) stop("X_ag must be a vector") # Cpp already gives error !! matrix takes first element.. 
  if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
  if(is.list(g)) setTRACpp(X, X_ag, g[[2]], trans) else setTRACpp(X, X_ag, g, trans)
}
setTRA.matrix <- function(X, X_ag, g = 0L, trans = "replace", ...) {
  if(!is.atomic(X_ag)) stop("X_ag must be a vector or matrix")
  if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
  if(is.list(g)) setTRAmCpp(X, X_ag, g[[2]], trans) else setTRAmCpp(X, X_ag, g, trans)
}
setTRA.data.frame <- function(X, X_ag, g = 0L, trans = "replace", ...) {
  if(is.array(X_ag)) stop("X_ag must be a vetor or list / data.frame")
  if(is.list(X_ag) && length(X_ag[[1]]) == 1) X_ag <- unlist(X_ag, use.names = FALSE)
  if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
  if(is.list(g)) setTRAlCpp(X, X_ag, g[[2]], trans) else setTRAlCpp(X, X_ag, g, trans)
}
setTRA.list <- function(X, X_ag, g = 0L, trans = "replace", ...) {
  if(is.array(X_ag)) stop("X_ag must be a vetor or list / data.frame")
  if(is.list(X_ag) && length(X_ag[[1]]) == 1) X_ag <- unlist(X_ag, use.names = FALSE)
  if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
  if(is.list(g)) setTRAlCpp(X, X_ag, g[[2]], trans) else setTRAlCpp(X, X_ag, g, trans)
}