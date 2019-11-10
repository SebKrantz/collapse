library(Rcpp)
sourceCpp('R/C++/fdiff.cpp', rebuild = TRUE) 
sourceCpp('R/C++/fdiffa.cpp', rebuild = TRUE)
sourceCpp('R/C++/fdiffl.cpp', rebuild = TRUE) # On large test data (10 mio obs), unordered panel-difference still gave error !!!!!
source("R/collapse R/GRP.R")
source("R/collapse R/small_helper.R")
source("R/collapse R/quick_conversion.R")

# For principle innovations of this code see flag.R and flag.cpp

fdiff <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...) {
  UseMethod("fdiff", x)
}
fdiff.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...) { 
  if(is.null(g)) 
    return(fdiffCpp(x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),give.names)) else if(is.atomic(g)) {
      if(is.factor(g)) nl <- fnlevels(g) else {
        g <- qG(g) 
        nl <- attr(g, "N.groups")
      }
      return(fdiffCpp(x,n,diff,fill,nl,g,NULL,G_t(t),give.names))
    } else {
      if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
      return(fdiffCpp(x,n,diff,fill,g[[1]],g[[2]],g[[3]],G_t(t),give.names))
    }
}
fdiff.pseries <- function(x, n = 1, diff = 1, fill = NA, give.names = TRUE, ...) { 
  index <- attr(x, "index")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  if(is.matrix(x))
    fdiffmCpp(x,n,diff,fill,fnlevels(index[[1]]),index[[1]],NULL,index[[2]],give.names) else
      fdiffCpp(x,n,diff,fill,fnlevels(index[[1]]),index[[1]],NULL,index[[2]],give.names) 
}
fdiff.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...) { 
  if(is.null(g)) 
    return(fdiffmCpp(x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),give.names)) else if(is.atomic(g)) {
      if(is.factor(g)) nl <- fnlevels(g) else {
        g <- qG(g) 
        nl <- attr(g, "N.groups")
      }
      fdiffmCpp(x,n,diff,fill,nl,g,NULL,G_t(t),give.names)
    } else {
      if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
      fdiffmCpp(x,n,diff,fill,g[[1]],g[[2]],g[[3]],G_t(t),give.names) 
    }
}
fdiff.grouped_df <- function(x, n = 1, diff = 1, t = NULL, fill = NA, give.names = TRUE, keep.ids = TRUE, ...) { 
  g <- GRP.grouped_df(x)
  tsym <- deparse(substitute(t))
  nam <- names(x)
  gn <- na.rm(match(names(g[[4]]), nam))
  if(!(tsym == "NULL" || is.na(tn <- match(tsym, nam)))) {
    if(any(gn == tn)) stop("timevar coincides with grouping variables!")
    t <- x[[tn]]
    gn <- c(gn, tn)
  }
  if(length(gn)) {
    if(!keep.ids) 
      return(fdifflCpp(x[-gn],n,diff,fill,g[[1]],g[[2]],g[[3]],G_t(t),give.names)) else {
        ax <- attributes(x)
        class(x) <- NULL # Works for multiple lags !!
        res <- c(x[gn],fdifflCpp(x[-gn],n,diff,fill,g[[1]],g[[2]],g[[3]],G_t(t),give.names))
        ax[["names"]] <- names(res)
        return(`attributes<-`(res, ax))
      }
  } else return(fdifflCpp(x,n,diff,fill,g[[1]],g[[2]],g[[3]],G_t(t),give.names))
}
fdiff.data.frame <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...) { 
  if(is.null(g)) 
    return(fdifflCpp(x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),give.names)) else if(is.atomic(g)) {
      if(is.factor(g)) nl <- fnlevels(g) else {
        g <- qG(g) 
        nl <- attr(g, "N.groups")
      }
      fdifflCpp(x,n,diff,fill,nl,g,NULL,G_t(t),give.names)
    } else {
      if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
      fdifflCpp(x,n,diff,fill,g[[1]],g[[2]],g[[3]],G_t(t),give.names) 
    }
}
fdiff.pdata.frame <- function(x, n = 1, diff = 1, fill = NA, give.names = TRUE, ...) { 
  index <- attr(x, "index")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  fdifflCpp(x,n,diff,fill,fnlevels(index[[1]]),index[[1]],NULL,index[[2]],give.names) 
}

# use xt instead of by ??? 
D <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...) {
  UseMethod("D", x)
}
D.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...) {
  fdiff.default(x, n, diff, g, t, fill, give.names, ...)
}
D.pseries <- function(x, n = 1, diff = 1, fill = NA, give.names = TRUE, ...) { 
  fdiff.pseries(x, n, diff, fill, give.names, ...)
}
D.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...) {
  fdiff.matrix(x, n, diff, g, t, fill, give.names, ...)
}
D.grouped_df <- fdiff.grouped_df
D.data.frame <- function(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric, 
                         fill = NA, give.names = TRUE, keep.ids = TRUE, ...) { 
  
  if(is.call(by) || is.call(t)) {
    ax <- attributes(x)
    class(x) <- NULL
    nam <- ax[["names"]]
    
    if(is.call(by)) {
      if(length(by) == 3L) {
        cols <- match(all.vars(by[[2L]]), nam)
        gn <- match(all.vars(by[[3L]]), nam)
      } else {
        gn <- match(all.vars(by), nam)
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
      tn <- match(t, nam)
      t <- x[[tn]]
      cols <- if(is.null(cols)) seq_along(x)[-tn] else cols[cols != tn]
      if(keep.ids) gn <- c(gn, tn)
    } 
    
    res <- if(length(gn)) 
      c(x[gn], fdifflCpp(x[cols],n,diff,fill,by[[1]],by[[2]],by[[3]],G_t(t),give.names)) else
        fdifflCpp(x[cols],n,diff,fill,by[[1]],by[[2]],by[[3]],G_t(t),give.names)
    ax[["names"]] <- names(res)
    return(`attributes<-`(res, ax))
  } else if(!is.null(cols)) {
    ax <- attributes(x)
    x <- if(is.function(cols)) unclass(x)[vapply(x, cols, TRUE)] else unclass(x)[cols]
    ax[["names"]] <- names(x)
    attributes(x) <- ax
  } 
  
  if(is.null(by)) 
    return(fdifflCpp(x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),give.names)) else if(is.atomic(by)) {
      if(is.factor(by)) nl <- fnlevels(by) else {
        by <- qG(by) 
        nl <- attr(by, "N.groups")
      }
      fdifflCpp(x,n,diff,fill,nl,by,NULL,G_t(t),give.names)
    } else {
      if(!is.GRP(by)) by <- GRP(by, return.groups = FALSE)
      fdifflCpp(x,n,diff,fill,by[[1]],by[[2]],by[[3]],G_t(t),give.names) 
    }
}
D.pdata.frame <- function(x, n = 1, diff = 1, cols = is.numeric, fill = NA, give.names = TRUE, keep.ids = TRUE, ...) { 
  
  ax <- attributes(x)
  nam <- ax[["names"]]
  index <- ax[["index"]]
  
  if(keep.ids) {
    gn <- match(names(index), nam)
    gn <- gn[!is.na(gn)]
    if(length(gn) && is.null(cols)) cols <- seq_along(x)[-gn] 
  } else gn <- NULL
  
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  
  if(!is.null(cols)) cols <- cols2int(cols, x, nam)
  
  if(length(gn) && !is.null(cols)) {
    class(x) <- NULL # Works for multiple lags !!
    res <- c(x[gn], fdifflCpp(x[cols],n,diff,fill,fnlevels(index[[1]]),index[[1]],NULL,index[[2]],give.names))
    ax[["names"]] <- names(res)
    return(`attributes<-`(res, ax))
  } else if(!length(gn)) # could speed up ?? 
    return(fdifflCpp(x[cols],n,diff,fill,fnlevels(index[[1]]),index[[1]],NULL,index[[2]],give.names)) else
      return(fdifflCpp(x,n,diff,fill,fnlevels(index[[1]]),index[[1]],NULL,index[[2]],give.names))
}


# Previous versions were deleted -> see previous verstions of L.data.frame and L.pdata.frame in flag.R !!

# Old versions: several functions for t, and no Xcols...
# is_factor_G <- function(x) {
#   if(is.factor(x) || is.integer(x)) x else qG(x)
# }
# is_GRP_G <- function(l) {
#   if(all(class(l) == "GRP")) l[[2]] else GRP(l, return.groups = FALSE)[[2]]
# }
# fdiff <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...) {
#   UseMethod("fdiff", x)
# }
# fdiff.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...) { 
#   if(is.null(g)) {
#     if(is.null(t)) fdiffCpp(x,n,diff,fill, names = give.names) else if(is.atomic(t)) 
#       fdiffCpp(x,n,diff,fill, t = is_factor_G(t), names = give.names) else 
#         fdiffCpp(x,n,diff,fill, t = is_GRP_G(t), names = give.names)
#   } else if(is.atomic(g)) {
#     if(!is.factor(g)) g = qF(g)
#     if(is.null(t)) {
#       message("Panel-difference computed without timevar: Assuming ordered data")
#       fdiffCpp(x,n,diff,fill,nlevels(g),g, names = give.names) 
#     } else if(is.atomic(t)) fdiffCpp(x,n,diff,fill,nlevels(g),g,NULL,is_factor_G(t),give.names) else 
#       fdiffCpp(x,n,diff,fill,nlevels(g),g,NULL,is_GRP_G(t),give.names)
#   } else {
#     if(!all(class(g) == "GRP")) g = GRP(g, return.groups = FALSE)
#     if(is.null(t)) {
#       message("Panel-difference computed without timevar: Assuming ordered data")
#       fdiffCpp(x,n,diff,fill,g[[1]],g[[2]], names = give.names) 
#     } else if(is.atomic(t)) fdiffCpp(x,n,diff,fill,g[[1]],g[[2]],g[[3]],is_factor_G(t),give.names) else 
#       fdiffCpp(x,n,diff,fill,g[[1]],g[[2]],g[[3]],is_GRP_G(t),give.names) 
#   }
# }
# fdiff.pseries <- function(x, n = 1, diff = 1, fill = NA, give.names = TRUE, ...) { 
#   index <- attr(x, "index")
#   fdiffCpp(x,n,diff,fill,nlevels(index[[1]]),index[[1]],NULL,index[[2]],give.names) 
# }
# fdiff.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...) { 
#   if(is.null(g)) {
#     if(is.null(t)) fdiffmCpp(x,n,diff,fill, names = give.names) else if(is.atomic(t)) 
#       fdiffmCpp(x,n,diff,fill, t = is_factor_G(t), names = give.names) else 
#         fdiffmCpp(x,n,diff,fill, t = is_GRP_G(t), names = give.names)
#   } else if(is.atomic(g)) {
#     if(!is.factor(g)) g = qF(g)
#     if(is.null(t)) {
#       message("Panel-difference computed without timevar: Assuming ordered data")
#       fdiffmCpp(x,n,diff,fill,nlevels(g),g,NULL,NULL,give.names) 
#     } else if(is.atomic(t)) fdiffmCpp(x,n,diff,fill,nlevels(g),g,NULL,is_factor_G(t),give.names) else 
#       fdiffmCpp(x,n,diff,fill,nlevels(g),g,NULL,is_GRP_G(t),give.names)
#   } else {
#     if(!all(class(g) == "GRP")) g = GRP(g, return.groups = FALSE)
#     if(is.null(t)) {
#       message("Panel-difference computed without timevar: Assuming ordered data")
#       fdiffmCpp(x,n,diff,fill,g[[1]],g[[2]],NULL,NULL,give.names) 
#     } else if(is.atomic(t)) fdiffmCpp(x,n,diff,fill,g[[1]],g[[2]],g[[3]],is_factor_G(t),give.names) else 
#       fdiffmCpp(x,n,diff,fill,g[[1]],g[[2]],g[[3]],is_GRP_G(t),give.names) 
#   }
# }
# fdiff.data.frame <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...) { 
#   if(is.null(g)) {
#     if(is.null(t)) fdifflCpp(x,n,diff,fill, names = give.names) else if(is.atomic(t)) 
#       fdifflCpp(x,n,diff,fill, t = is_factor_G(t), names = give.names) else 
#         fdifflCpp(x,n,diff,fill, t = is_GRP_G(t), names = give.names)
#   } else if(is.atomic(g)) {
#     if(!is.factor(g)) g = qF(g)
#     if(is.null(t)) {
#       message("Panel-difference computed without timevar: Assuming ordered data")
#       fdifflCpp(x,n,diff,fill,nlevels(g),g,NULL,NULL,give.names) 
#     } else if(is.atomic(t)) fdifflCpp(x,n,diff,fill,nlevels(g),g,NULL,is_factor_G(t),give.names) else 
#       fdifflCpp(x,n,diff,fill,nlevels(g),g,NULL,is_GRP_G(t),give.names)
#   } else {
#     if(!all(class(g) == "GRP")) g = GRP(g, return.groups = FALSE)
#     if(is.null(t)) {
#       message("Panel-difference computed without timevar: Assuming ordered data")
#       fdifflCpp(x,n,diff,fill,g[[1]],g[[2]],NULL,NULL,give.names) 
#     } else if(is.atomic(t)) fdifflCpp(x,n,diff,fill,g[[1]],g[[2]],g[[3]],is_factor_G(t),give.names) else 
#       fdifflCpp(x,n,diff,fill,g[[1]],g[[2]],g[[3]],is_GRP_G(t),give.names) 
#   }
# }
# fdiff.pdata.frame <- function(x, n = 1, diff = 1, fill = NA, give.names = TRUE, ...) { 
#   index <- attr(x, "index")
#   fdifflCpp(x,n,diff,fill,nlevels(index[[1]]),index[[1]],NULL,index[[2]],give.names) 
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
#         fdiffCpp(x,n,fill,0L,0L,NULL,t[[2]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#     }
#   } else { # best way to organize this code ??
#     if(is.factor(g)) {
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         fdiffCpp(x,n,fill,nlevels(g),g) 
#       } else {
#         if(is.factor(t)) fdiffCpp(x,n,fill,nlevels(g),g,NULL,t) else if (inherits(t,"GRP"))
#           fdiffCpp(x,n,fill,nlevels(g),g,NULL,t[[2]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     } else {
#       if(!inherits(g,"GRP")) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         fdiffCpp(x,n,fill,g[[1]],g[[2]]) 
#       } else {
#         if(is.factor(t)) fdiffCpp(x,n,fill,g[[1]],g[[2]],g[[3]],t) else if (inherits(t,"GRP"))
#           fdiffCpp(x,n,fill,g[[1]],g[[2]],g[[3]],t[[2]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
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
#         fdiffmCpp(x,n,fill,0L,0L,NULL,t[[2]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#     }
#   } else { # best way to organize this code ??
#     if (is.factor(g)) {
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         fdiffmCpp(x,n,fill,nlevels(g),g) 
#       } else {
#         if(is.factor(t)) fdiffmCpp(x,n,fill,nlevels(g),g,NULL,t) else if (inherits(t,"GRP"))
#           fdiffmCpp(x,n,fill,nlevels(g),g,NULL,t[[2]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     } else {
#       if(!inherits(g,"GRP")) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         fdiffmCpp(x,n,fill,g[[1]],g[[2]]) 
#       } else {
#         if(is.factor(t)) fdiffmCpp(x,n,fill,g[[1]],g[[2]],g[[3]],t) else if (inherits(t,"GRP"))
#           fdiffmCpp(x,n,fill,g[[1]],g[[2]],g[[3]],t[[2]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     }
#   }
#   if(length(n)>1) {
#     dx <- dimnames(x)
#     dx[[2]] <- lagnames(dx[[2]], n)
#     dimnames(res) <- dx
#   } else dimnames(res) <- dimnames(x) 
#   res
# }
# fdiff.data.frame <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, ...) { 
#   ax <- attributes(x)
#   res <- if(is.null(g)) {
#     if(is.null(t)) fdifflCpp(x,n,fill) else { # best way to organize this code ??
#       if(is.factor(t)) fdifflCpp(x,n,fill,0L,0L,NULL,t) else if (inherits(t,"GRP"))
#         fdifflCpp(x,n,fill,0L,0L,NULL,t[[2]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#     }
#   } else { # best way to organize this code ??
#     if (is.factor(g)) {
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         fdifflCpp(x,n,fill,nlevels(g),g) 
#       } else {
#         if(is.factor(t)) fdifflCpp(x,n,fill,nlevels(g),g,NULL,t) else if (inherits(t,"GRP"))
#           fdifflCpp(x,n,fill,nlevels(g),g,NULL,t[[2]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     } else {
#       if(!inherits(g,"GRP")) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         fdifflCpp(x,n,fill,g[[1]],g[[2]]) 
#       } else {
#         if(is.factor(t)) fdifflCpp(x,n,fill,g[[1]],g[[2]],g[[3]],t) else if (inherits(t,"GRP"))
#           fdifflCpp(x,n,fill,g[[1]],g[[2]],g[[3]],t[[2]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     }
#   }
#   if(length(n)>1) ax[["names"]] <- lagnames(ax[["names"]], n)
#   attributes(res) <- ax
#   res
# }
