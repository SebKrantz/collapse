# library(Rcpp)
# sourceCpp('src/fdiff.cpp', rebuild = TRUE)
# sourceCpp('src/fdiffa.cpp', rebuild = TRUE)
# sourceCpp('src/fdiffl.cpp', rebuild = TRUE) # On large test data (10 mio obs), unordered panel-difference still gave error !!!!!
# source("R/GRP.R")
# source("R/small_helper.R")
# source("R/quick_conversion.R")

# For principle innovations of this code see flag.R and flag.cpp # stubs instead of give.names !!!

fdiff <- function(x, n = 1, diff = 1, ...) { # , g = NULL, t = NULL, fill = NA, stubs = TRUE
  UseMethod("fdiff", x)
}
fdiff.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(is.null(g))
    return(fdiffCpp(x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),stubs)) else if(is.atomic(g)) {
      if(is.factor(g)) nl <- fnlevels(g) else {
        g <- qG(g)
        nl <- attr(g, "N.groups")
      }
      return(fdiffCpp(x,n,diff,fill,nl,g,NULL,G_t(t),stubs))
    } else {
      if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
      return(fdiffCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs))
    }
}
fdiff.pseries <- function(x, n = 1, diff = 1, fill = NA, stubs = TRUE, ...) {
  index <- attr(x, "index")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  if(is.matrix(x))
    fdiffmCpp(x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs) else
      fdiffCpp(x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs)
}
fdiff.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(is.null(g))
    return(fdiffmCpp(x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),stubs)) else if(is.atomic(g)) {
      if(is.factor(g)) nl <- fnlevels(g) else {
        g <- qG(g)
        nl <- attr(g, "N.groups")
      }
      fdiffmCpp(x,n,diff,fill,nl,g,NULL,G_t(t),stubs)
    } else {
      if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
      fdiffmCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs)
    }
}
fdiff.grouped_df <- function(x, n = 1, diff = 1, t = NULL, fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {
  g <- GRP.grouped_df(x)
  tsym <- deparse(substitute(t))
  nam <- names(x)
  gn <- which(nam %in% g[[5L]])
  if(!(tsym == "NULL" || is.na(tn <- match(tsym, nam)))) {
    if(any(gn == tn)) stop("timevar coincides with grouping variables!")
    t <- x[[tn]]
    gn <- c(gn, tn)
  }
  if(length(gn)) {
    if(!keep.ids)
      return(fdifflCpp(x[-gn],n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs)) else {
        ax <- attributes(x)
        class(x) <- NULL # Works for multiple lags !!
        res <- c(x[gn],fdifflCpp(x[-gn],n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs))
        ax[["names"]] <- names(res)
        return(setAttributes(res, ax))
      }
  } else return(fdifflCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs))
}
fdiff.data.frame <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(is.null(g))
    return(fdifflCpp(x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),stubs)) else if(is.atomic(g)) {
      if(is.factor(g)) nl <- fnlevels(g) else {
        g <- qG(g)
        nl <- attr(g, "N.groups")
      }
      fdifflCpp(x,n,diff,fill,nl,g,NULL,G_t(t),stubs)
    } else {
      if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
      fdifflCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs)
    }
}
fdiff.pdata.frame <- function(x, n = 1, diff = 1, fill = NA, stubs = TRUE, ...) {
  index <- attr(x, "index")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  fdifflCpp(x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs)
}

# use xt instead of by ???
D <- function(x, n = 1, diff = 1, ...) { # g = NULL, t = NULL, fill = NA, stubs = TRUE
  UseMethod("D", x)
}
D.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  fdiff.default(x, n, diff, g, t, fill, stubs, ...)
}
D.pseries <- function(x, n = 1, diff = 1, fill = NA, stubs = TRUE, ...) {
  fdiff.pseries(x, n, diff, fill, stubs, ...)
}
D.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  fdiff.matrix(x, n, diff, g, t, fill, stubs, ...)
}
D.grouped_df <- fdiff.grouped_df
D.data.frame <- function(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric,
                         fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {

  if(is.call(by) || is.call(t)) {
    ax <- attributes(x)
    class(x) <- NULL
    nam <- ax[["names"]]

    if(is.call(by)) {
      if(length(by) == 3L) {
        cols <- anyNAerror(match(all.vars(by[[2L]]), nam), "Unknown variables passed to by!")
        gn <- anyNAerror(match(all.vars(by[[3L]]), nam), "Unknown variables passed to by!")
      } else {
        gn <- anyNAerror(match(all.vars(by), nam), "Unknown variables passed to by!")
        cols <- if(is.null(cols)) seq_along(x)[-gn] else cols2int(cols, x, nam)
      }
      by <- if(length(gn) == 1L) at2GRP(x[[gn]]) else GRP(x, gn, return.groups = FALSE)
      if(!keep.ids) gn <- NULL
    } else {
      gn <- NULL
      if(!is.null(cols)) cols <- cols2int(cols, x, nam)
      if(!is.GRP(by)) by <- if(is.null(by)) list(0L, 0L, NULL) else if(is.atomic(by)) # Necessary for if by is passed externally !!
        at2GRP(by) else GRP(by, return.groups = FALSE)
    }

    if(is.call(t)) {
      t <- all.vars(t)
      tn <- anyNAerror(match(t, nam), "Unknown variables passed to t!")
      t <- x[[tn]]
      cols <- if(is.null(cols)) seq_along(x)[-tn] else cols[cols != tn]
      if(keep.ids) gn <- c(gn, tn)
    }

    res <- if(length(gn))
      c(x[gn], fdifflCpp(x[cols],n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),stubs)) else
        fdifflCpp(x[cols],n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),stubs)
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(!is.null(cols)) {
    ax <- attributes(x)
    x <- if(is.function(cols)) unclass(x)[vapply(x, cols, TRUE)] else unclass(x)[cols]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }

  if(is.null(by))
    return(fdifflCpp(x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),stubs)) else if(is.atomic(by)) {
      if(is.factor(by)) nl <- fnlevels(by) else {
        by <- qG(by)
        nl <- attr(by, "N.groups")
      }
      fdifflCpp(x,n,diff,fill,nl,by,NULL,G_t(t),stubs)
    } else {
      if(!is.GRP(by)) by <- GRP(by, return.groups = FALSE)
      fdifflCpp(x,n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),stubs)
    }
}
D.pdata.frame <- function(x, n = 1, diff = 1, cols = is.numeric, fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {

  ax <- attributes(x)
  nam <- ax[["names"]]
  index <- ax[["index"]]

  if(keep.ids) {
    gn <- which(nam %in% names(index))
    if(length(gn) && is.null(cols)) cols <- seq_along(x)[-gn]
  } else gn <- NULL

  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])

  if(!is.null(cols)) cols <- cols2int(cols, x, nam)

  if(length(gn) && !is.null(cols)) {
    class(x) <- NULL # Works for multiple lags !!
    res <- c(x[gn], fdifflCpp(x[cols],n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs))
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(!length(gn)) # could speed up ??
    return(fdifflCpp(x[cols],n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs)) else
      return(fdifflCpp(x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs))
}


# Previous versions were deleted -> see previous verstions of L.data.frame and L.pdata.frame in flag.R !!

# Old versions: several functions for t, and no Xcols...
# is_factor_G <- function(x) {
#   if(is.factor(x) || is.integer(x)) x else qG(x)
# }
# is_GRP_G <- function(l) {
#   if(all(class(l) == "GRP")) l[[2L]] else GRP(l, return.groups = FALSE)[[2L]]
# }
# fdiff <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
#   UseMethod("fdiff", x)
# }
# fdiff.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
#   if(is.null(g)) {
#     if(is.null(t)) fdiffCpp(x,n,diff,fill, names = stubs) else if(is.atomic(t))
#       fdiffCpp(x,n,diff,fill, t = is_factor_G(t), names = stubs) else
#         fdiffCpp(x,n,diff,fill, t = is_GRP_G(t), names = stubs)
#   } else if(is.atomic(g)) {
#     if(!is.factor(g)) g = qF(g)
#     if(is.null(t)) {
#       message("Panel-difference computed without timevar: Assuming ordered data")
#       fdiffCpp(x,n,diff,fill,nlevels(g),g, names = stubs)
#     } else if(is.atomic(t)) fdiffCpp(x,n,diff,fill,nlevels(g),g,NULL,is_factor_G(t),stubs) else
#       fdiffCpp(x,n,diff,fill,nlevels(g),g,NULL,is_GRP_G(t),stubs)
#   } else {
#     if(!all(class(g) == "GRP")) g = GRP(g, return.groups = FALSE)
#     if(is.null(t)) {
#       message("Panel-difference computed without timevar: Assuming ordered data")
#       fdiffCpp(x,n,diff,fill,g[[1L]],g[[2L]], names = stubs)
#     } else if(is.atomic(t)) fdiffCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],is_factor_G(t),stubs) else
#       fdiffCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],is_GRP_G(t),stubs)
#   }
# }
# fdiff.pseries <- function(x, n = 1, diff = 1, fill = NA, stubs = TRUE, ...) {
#   index <- attr(x, "index")
#   fdiffCpp(x,n,diff,fill,nlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs)
# }
# fdiff.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
#   if(is.null(g)) {
#     if(is.null(t)) fdiffmCpp(x,n,diff,fill, names = stubs) else if(is.atomic(t))
#       fdiffmCpp(x,n,diff,fill, t = is_factor_G(t), names = stubs) else
#         fdiffmCpp(x,n,diff,fill, t = is_GRP_G(t), names = stubs)
#   } else if(is.atomic(g)) {
#     if(!is.factor(g)) g = qF(g)
#     if(is.null(t)) {
#       message("Panel-difference computed without timevar: Assuming ordered data")
#       fdiffmCpp(x,n,diff,fill,nlevels(g),g,NULL,NULL,stubs)
#     } else if(is.atomic(t)) fdiffmCpp(x,n,diff,fill,nlevels(g),g,NULL,is_factor_G(t),stubs) else
#       fdiffmCpp(x,n,diff,fill,nlevels(g),g,NULL,is_GRP_G(t),stubs)
#   } else {
#     if(!all(class(g) == "GRP")) g = GRP(g, return.groups = FALSE)
#     if(is.null(t)) {
#       message("Panel-difference computed without timevar: Assuming ordered data")
#       fdiffmCpp(x,n,diff,fill,g[[1L]],g[[2L]],NULL,NULL,stubs)
#     } else if(is.atomic(t)) fdiffmCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],is_factor_G(t),stubs) else
#       fdiffmCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],is_GRP_G(t),stubs)
#   }
# }
# fdiff.data.frame <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
#   if(is.null(g)) {
#     if(is.null(t)) fdifflCpp(x,n,diff,fill, names = stubs) else if(is.atomic(t))
#       fdifflCpp(x,n,diff,fill, t = is_factor_G(t), names = stubs) else
#         fdifflCpp(x,n,diff,fill, t = is_GRP_G(t), names = stubs)
#   } else if(is.atomic(g)) {
#     if(!is.factor(g)) g = qF(g)
#     if(is.null(t)) {
#       message("Panel-difference computed without timevar: Assuming ordered data")
#       fdifflCpp(x,n,diff,fill,nlevels(g),g,NULL,NULL,stubs)
#     } else if(is.atomic(t)) fdifflCpp(x,n,diff,fill,nlevels(g),g,NULL,is_factor_G(t),stubs) else
#       fdifflCpp(x,n,diff,fill,nlevels(g),g,NULL,is_GRP_G(t),stubs)
#   } else {
#     if(!all(class(g) == "GRP")) g = GRP(g, return.groups = FALSE)
#     if(is.null(t)) {
#       message("Panel-difference computed without timevar: Assuming ordered data")
#       fdifflCpp(x,n,diff,fill,g[[1L]],g[[2L]],NULL,NULL,stubs)
#     } else if(is.atomic(t)) fdifflCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],is_factor_G(t),stubs) else
#       fdifflCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],is_GRP_G(t),stubs)
#   }
# }
# fdiff.pdata.frame <- function(x, n = 1, diff = 1, fill = NA, stubs = TRUE, ...) {
#   index <- attr(x, "index")
#   fdifflCpp(x,n,diff,fill,nlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs)
# }





# Previous attempts: Generating names in R !!!
# diffnames <- function(n, diff, nam = NULL) {
#   pos <- abs(n) > 1
#   neg <- n < -1
#
#   lapply(n[!pos], switch, `-1` = "FD", `0` = "", `1` = "D")
#   # as.vector(t(outer(n, diff, paste, sep="."))) # https://stackoverflow.com/questions/16143700/pasting-two-vectors-with-combinations-of-all-vectors-elements/23533365
#   paste(rep(n, each = length(diff)), diff, sep = ".") # fastest !!
#
#   res <- character(length(n))
#   res[pos] <- paste0("D",n[pos],".")
#   res[neg] <- paste0("FD",abs(n[neg]),".")
#   paste0(rep(res, length(nam)), rep(nam, each = length(n)))
# }
#
# diffn <- function(n, diff) {
#   pos <- abs(n) > 1
#   neg <- n < -1
#
#   lapply(n[!pos], switch, `-1` = "FD", `0` = "", `1` = "D")
#   # as.vector(t(outer(n, diff, paste, sep="."))) # https://stackoverflow.com/questions/16143700/pasting-two-vectors-with-combinations-of-all-vectors-elements/23533365
#   paste(rep(n, each = length(diff)), diff, sep = ".") # fastest !!
#
#   res <- character(length(n))
#   res[pos] <- paste0("D",n[pos],".")
#   res[neg] <- paste0("FD",abs(n[neg]),".")
#   paste0(rep(res, length(nam)), rep(nam, each = length(n)))
# }
#
# # Also make compatible with all types of data !!!
# # Also if fill = NULL -> Delete !!
#
# fdiff <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, ...) {
#   UseMethod("fdiff", x)
# }
# fdiff.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, ...) {
#    if(is.null(g)) {
#     if(is.null(t)) {
#        res <- fdiffCpp(x,n,fill)
#        if(diff > 1) { # this will be too much code !! -> do in C++ !!
#          i <- 1
#          repeat {
#            i = i + 1
#            if(i == diff) break
#            res <- fdiffCpp(res,n,fill)
#          }
#        }
#       } else { # best way to organize this code ??
#       if(is.factor(t)) fdiffCpp(x,n,fill,0L,0L,NULL,t) else if (inherits(t,"GRP"))
#         fdiffCpp(x,n,fill,0L,0L,NULL,t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#     }
#   } else { # best way to organize this code ??
#     if(is.factor(g)) {
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         fdiffCpp(x,n,fill,nlevels(g),g)
#       } else {
#         if(is.factor(t)) fdiffCpp(x,n,fill,nlevels(g),g,NULL,t) else if (inherits(t,"GRP"))
#           fdiffCpp(x,n,fill,nlevels(g),g,NULL,t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     } else {
#       if(!inherits(g,"GRP")) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         fdiffCpp(x,n,fill,g[[1L]],g[[2L]])
#       } else {
#         if(is.factor(t)) fdiffCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t) else if (inherits(t,"GRP"))
#           fdiffCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     }
#   }
#   if(length(n) == 1) attributes(res) <- attributes(x) else dimnames(res) <- list(names(x), n) # lagn(n)) # lagn(n)) # # paste0(ifelse(n>=0,"L","F"),abs(n)) # fastest way ?? or ax["dim"] <- attr(res,"dim") # list(as.character(ax[["names"]]), paste0(ifelse(n>=0,"L","F"),abs(n)))
#   res
# }
# fdiff.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, ...) {
#   res <- if(is.null(g)) {
#     if(is.null(t)) fdiffmCpp(x,n,fill) else { # best way to organize this code ??
#       if(is.factor(t)) fdiffmCpp(x,n,fill,0L,0L,NULL,t) else if (inherits(t,"GRP"))
#         fdiffmCpp(x,n,fill,0L,0L,NULL,t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#     }
#   } else { # best way to organize this code ??
#     if (is.factor(g)) {
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         fdiffmCpp(x,n,fill,nlevels(g),g)
#       } else {
#         if(is.factor(t)) fdiffmCpp(x,n,fill,nlevels(g),g,NULL,t) else if (inherits(t,"GRP"))
#           fdiffmCpp(x,n,fill,nlevels(g),g,NULL,t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     } else {
#       if(!inherits(g,"GRP")) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         fdiffmCpp(x,n,fill,g[[1L]],g[[2L]])
#       } else {
#         if(is.factor(t)) fdiffmCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t) else if (inherits(t,"GRP"))
#           fdiffmCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     }
#   }
#   if(length(n)>1) {
#     dx <- dimnames(x)
#     dx[[2L]] <- lagnames(dx[[2L]], n)
#     dimnames(res) <- dx
#   } else dimnames(res) <- dimnames(x)
#   res
# }
# fdiff.data.frame <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, ...) {
#   ax <- attributes(x)
#   res <- if(is.null(g)) {
#     if(is.null(t)) fdifflCpp(x,n,fill) else { # best way to organize this code ??
#       if(is.factor(t)) fdifflCpp(x,n,fill,0L,0L,NULL,t) else if (inherits(t,"GRP"))
#         fdifflCpp(x,n,fill,0L,0L,NULL,t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#     }
#   } else { # best way to organize this code ??
#     if (is.factor(g)) {
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         fdifflCpp(x,n,fill,nlevels(g),g)
#       } else {
#         if(is.factor(t)) fdifflCpp(x,n,fill,nlevels(g),g,NULL,t) else if (inherits(t,"GRP"))
#           fdifflCpp(x,n,fill,nlevels(g),g,NULL,t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     } else {
#       if(!inherits(g,"GRP")) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         fdifflCpp(x,n,fill,g[[1L]],g[[2L]])
#       } else {
#         if(is.factor(t)) fdifflCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t) else if (inherits(t,"GRP"))
#           fdifflCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     }
#   }
#   if(length(n)>1) ax[["names"]] <- lagnames(ax[["names"]], n)
#   attributes(res) <- ax
#   res
# }
