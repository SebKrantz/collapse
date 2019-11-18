# library(Rcpp)
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# sdfsdfsd
# Note: ng is not supplied to TRACpp for easier stat functions, thus you need to check in R !!

TRA <- function(x, STATS, FUN = "-", g = NULL, keep.group_keys = TRUE, ...) {
  UseMethod("TRA", x)
}
TRA.default <- function(x, STATS, FUN = "-", g = NULL, ...) {
  if(is.null(g)) return(TRACpp(x,STATS,0L,TRAtoInt(FUN))) else if(is.atomic(g)) {
    if(is.factor(g)) {
      if(fnlevels(g) != length(STATS)) stop("number of groups must match length(STATS)")
    } else {
      g <- qG(g) # needs to be ordered to be compatible with fast functions !!
      if(attr(g, "N.groups") != length(STATS)) stop("number of groups must match length(STATS)")
    }
    return(TRACpp(x,STATS,g,TRAtoInt(FUN)))
  } else {
    if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
    if(g[[1L]] != length(STATS)) stop("number of groups must match length(STATS)")
    return(TRACpp(x,STATS,g[[2L]],TRAtoInt(FUN)))
  }
}
TRA.matrix <- function(x, STATS, FUN = "-", g = NULL, ...) {
  if(is.null(g)) return(TRAmCpp(x,STATS,0L,TRAtoInt(FUN))) else if(is.atomic(g)) {
    if(is.factor(g)) {
      if(fnlevels(g) != nrow(STATS)) stop("number of groups must match nrow(STATS)")
    } else {
      g <- qG(g) # needs to be ordered to be compatible with fast functions !!
      if(attr(g, "N.groups") != nrow(STATS)) stop("number of groups must match nrow(STATS)")
    }
    return(TRAmCpp(x,STATS,g,TRAtoInt(FUN)))
  } else {
    if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
    if(g[[1L]] != nrow(STATS)) stop("number of groups must match nrow(STATS)")
    return(TRAmCpp(x,STATS,g[[2L]],TRAtoInt(FUN)))
  }
}
TRA.data.frame <- function(x, STATS, FUN = "-", g = NULL, ...) {
  if(is.null(g)) return(TRAlCpp(x,STATS,0L,TRAtoInt(FUN))) else if(is.atomic(g)) {
    if(is.factor(g)) {
      if(fnlevels(g) != nrow(STATS)) stop("number of groups must match nrow(STATS)")
    } else {
      g <- qG(g) # needs to be ordered to be compatible with fast functions !!
      if(attr(g, "N.groups") != nrow(STATS)) stop("number of groups must match nrow(STATS)")
    }
    return(TRAlCpp(x,STATS,g,TRAtoInt(FUN)))
  } else {
    if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
    if(g[[1L]] != nrow(STATS)) stop("number of groups must match nrow(STATS)")
    return(TRAlCpp(x,STATS,g[[2L]],TRAtoInt(FUN)))
  }
}
# TRA.grouped_df <- function(x, STATS, FUN = "-", keep.group_keys = TRUE, ...) {
#   g <- GRP.grouped_df(x)
#   gn <- which(names(x) %in% g[[5L]])
#   if(g[[1L]] != nrow(STATS)) stop("number of groups must match nrow(STATS)")
#   if(length(gn) > 0L && length(STATS) != length(x)) {
#     ax <- attributes(x)
#     attributes(x) <- NULL
#     namst <- names(STATS) # Finish !! shoudl sweep out only the indicated column / and setcolorder afterwards ??
#     ginst <- namst %in% g[[5L]]
#     if(anyNA(mt <- match(namst, ax[["names"]]))) stop("the variable names of x and STATS must match")
#     nomt <- which(!(ax[["names"]] %in% namst))
#     if(!all(nomg <- nomatch %in% gn)) {
#       gn <- c(gn, nomatch[!nomg])
#       STATS <- unclass(STATS)[-match(g[[5L]], namst)]
#     }
#
#     if(any(ginst)) STATS <- unclass(STATS)[!ginst]
#
#     if(keep.group_keys) {
#       ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
#       return(setAttributes(c(x[gn],TRAlCpp(x[-gn],STATS,g[[2L]],TRAtoInt(FUN))), ax))
#     } else {
#       ax[["names"]] <- ax[["names"]][-gn]
#       return(setAttributes(TRAlCpp(x[-gn],STATS,g[[2L]],TRAtoInt(FUN)), ax))
#     }
#   } else return(TRAlCpp(x,STATS,g[[2L]],TRAtoInt(FUN)))
# }


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
