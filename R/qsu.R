# library(Rcpp)
# sourceCpp('R/C++/fbstats.cpp')

# Note: for principal innovations of this code see fsum.R !!
# nv <- function(x) unclass(x)[vapply(x, is.numeric, TRUE, USE.NAMES = FALSE)]
# cols2int <- function(x, cols) if(is.function(cols)) which(vapply(x, cols, TRUE)) else if(is.character(cols)) match(cols, names(x)) else cols

qsu <- function(x, ...) { # g = NULL, pid = NULL, w = NULL, higher = FALSE, array = TRUE,
  UseMethod("qsu", x)
}
qsu.default <- function(x, g = NULL, pid = NULL, w = NULL, higher = FALSE, array = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    if(is.null(g)) {
      if(is.null(pid)) return(fbstatsCpp(x,higher, w = w)) else if(is.atomic(pid)) {
        if(is.nmfactor(pid)) nid <- fnlevels(pid) else {
          pid <- qG(pid, na.exclude = FALSE)
          nid <- attr(pid, "N.groups")
        }
        return(fbstatsCpp(x,higher,0L,0L,nid,pid,w))
      } else if(!is.GRP(pid)) pid <- GRP(pid, return.groups = FALSE)
        return(fbstatsCpp(x,higher,0L,0L,pid[[1L]],pid[[2L]],w))
    } else if (is.atomic(g)) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        if(is.null(pid)) return(fbstatsCpp(x,higher,length(lev),g,0L,0L,w,TRUE,TRUE,lev)) else if(is.atomic(pid)) {
          if(is.nmfactor(pid)) nid <- fnlevels(pid) else {
            pid <- qG(pid, na.exclude = FALSE)
            nid <- attr(pid, "N.groups")
          }
          return(fbstatsCpp(x,higher,length(lev),g,nid,pid,w,array,TRUE,lev))
        } else if(!is.GRP(pid)) pid <- GRP(pid, return.groups = FALSE)
        return(fbstatsCpp(x,higher,length(lev),g,pid[[1L]],pid[[2L]],w,array,TRUE,lev))
    } else {
      if(!is.GRP(g)) g <- GRP(g)
      if(is.null(pid)) return(fbstatsCpp(x,higher,g[[1L]],g[[2L]],0L,0L,w,TRUE,TRUE,group_names.GRP(g))) else if(is.atomic(pid)) {
        if(is.nmfactor(pid)) nid <- fnlevels(pid) else {
          pid <- qG(pid, na.exclude = FALSE)
          nid <- attr(pid, "N.groups")
        }
        return(fbstatsCpp(x,higher,g[[1L]],g[[2L]],nid,pid,w,array,TRUE,group_names.GRP(g)))
      } else if(!is.GRP(pid)) pid <- GRP(pid, return.groups = FALSE)
      return(fbstatsCpp(x,higher,g[[1L]],g[[2L]],pid[[1L]],pid[[2L]],w,array,TRUE,group_names.GRP(g)))
    }
}
qsu.pseries <- function(x, g = NULL, w = NULL, higher = FALSE, array = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  index <- attr(x, "index")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  if(is.null(g))
  return(fbstatsCpp(x,higher,0L,0L,fnlevels(index[[1L]]),index[[1L]],w)) else if (is.atomic(g)) {
    if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
    lev <- attr(g, "levels")
    return(fbstatsCpp(x,higher,length(lev),g,fnlevels(index[[1L]]),index[[1L]],w,array,TRUE,lev))
  } else if(!is.GRP(g)) g <- GRP(g)
    return(fbstatsCpp(x,higher,g[[1L]],g[[2L]],fnlevels(index[[1L]]),index[[1L]],w,array,TRUE,group_names.GRP(g)))
}
qsu.matrix <- function(x, g = NULL, pid = NULL, w = NULL, higher = FALSE, array = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(is.null(g)) {
    if(is.null(pid)) return(fbstatsmCpp(x,higher, w = w)) else if(is.atomic(pid)) {
      if(is.nmfactor(pid)) nid <- fnlevels(pid) else {
        pid <- qG(pid, na.exclude = FALSE)
        nid <- attr(pid, "N.groups")
      }
      return(fbstatsmCpp(x,higher,0L,0L,nid,pid,w,array))
    } else if(!is.GRP(pid)) pid <- GRP(pid, return.groups = FALSE)
    return(fbstatsmCpp(x,higher,0L,0L,pid[[1L]],pid[[2L]],w,array))
  } else if (is.atomic(g)) {
    if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
    lev <- attr(g, "levels")
    if(is.null(pid)) return(fbstatsmCpp(x,higher,length(lev),g,0L,0L,w,array,lev)) else if(is.atomic(pid)) {
      if(is.nmfactor(pid)) nid <- fnlevels(pid) else {
        pid <- qG(pid, na.exclude = FALSE)
        nid <- attr(pid, "N.groups")
      }
      return(fbstatsmCpp(x,higher,length(lev),g,nid,pid,w,array,lev))
    } else if(!is.GRP(pid)) pid <- GRP(pid, return.groups = FALSE)
    return(fbstatsmCpp(x,higher,length(lev),g,pid[[1L]],pid[[2L]],w,array,lev))
  } else {
    if(!is.GRP(g)) g <- GRP(g)
    if(is.null(pid)) return(fbstatsmCpp(x,higher,g[[1L]],g[[2L]],0L,0L,w,array,group_names.GRP(g))) else if(is.atomic(pid)) {
      if(is.nmfactor(pid)) nid <- fnlevels(pid) else {
        pid <- qG(pid, na.exclude = FALSE)
        nid <- attr(pid, "N.groups")
      }
      return(fbstatsmCpp(x,higher,g[[1L]],g[[2L]],nid,pid,w,array,group_names.GRP(g)))
    } else if(!is.GRP(pid)) pid <- GRP(pid, return.groups = FALSE)
    return(fbstatsmCpp(x,higher,g[[1L]],g[[2L]],pid[[1L]],pid[[2L]],w,array,group_names.GRP(g)))
  }
}
qsu.data.frame <- function(x, by = NULL, pid = NULL, w = NULL, cols = NULL, higher = FALSE, array = TRUE, vlabels = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  formby <- is.call(by)
  formpid <- is.call(pid)

  # fastest solution!! (see checks below !!)
  if(formby || formpid) {
    v <- NULL
    class(x) <- NULL
    if(formby) {
      if(length(by) == 3L) {
        v <- all.vars(by[[2L]])
        namby <- all.vars(by[[3L]])
      } else namby <- all.vars(by)
      by <- if(length(namby) == 1L) x[[namby]] else GRP(x, namby)
    } else namby <- NULL
    if(formpid) {
      if(length(pid) == 3L) {
        v <- all.vars(pid[[2L]])
        nampid <- all.vars(pid[[3L]])
      } else nampid <- all.vars(pid)
      pid <- if(length(nampid) == 1L) x[[nampid]] else GRP(x, nampid)
    } else nampid <- NULL
    if(is.null(v)) {
      x <- if(is.null(cols)) x[-match(c(namby,nampid), names(x))] else x[cols2int(cols, x, names(x))]
    } else x <- x[v]
  } else if(!is.null(cols)) x <- unclass(x)[cols2int(cols, x, names(x))]


  # Get vlabels !!
  if(vlabels) names(x) <- paste(names(x), vlabels(x), sep = ": ")

  # original code:
  if(is.null(by)) {
    if(is.null(pid)) return(fbstatslCpp(x,higher, w = w)) else if(is.atomic(pid)) {
      if(is.nmfactor(pid)) nid <- fnlevels(pid) else {
        pid <- qG(pid, na.exclude = FALSE)
        nid <- attr(pid, "N.groups")
      }
      return(drop(fbstatslCpp(x,higher,0L,0L,nid,pid,w,array)))
    } else if(!is.GRP(pid)) pid <- GRP(pid, return.groups = FALSE)
    return(drop(fbstatslCpp(x,higher,0L,0L,pid[[1L]],pid[[2L]],w,array)))
  } else if (is.atomic(by)) {
    if(!is.nmfactor(by)) by <- qF(by, na.exclude = FALSE)
    lev <- attr(by, "levels")
    if(is.null(pid)) return(drop(fbstatslCpp(x,higher,length(lev),by,0L,0L,w,array,lev))) else if(is.atomic(pid)) {
      if(is.nmfactor(pid)) nid <- fnlevels(pid) else {
        pid <- qG(pid, na.exclude = FALSE)
        nid <- attr(pid, "N.groups")
      }
      return(drop(fbstatslCpp(x,higher,length(lev),by,nid,pid,w,array,lev)))
    } else if(!is.GRP(pid)) pid <- GRP(pid, return.groups = FALSE)
    return(drop(fbstatslCpp(x,higher,length(lev),by,pid[[1L]],pid[[2L]],w,array,lev)))
  } else {
    if(!is.GRP(by)) by <- GRP(by)
    if(is.null(pid)) return(drop(fbstatslCpp(x,higher,by[[1L]],by[[2L]],0L,0L,w,array,group_names.GRP(by)))) else if(is.atomic(pid)) {
      if(is.nmfactor(pid)) nid <- fnlevels(pid) else {
        pid <- qG(pid, na.exclude = FALSE)
        nid <- attr(pid, "N.groups")
      }
      return(drop(fbstatslCpp(x,higher,by[[1L]],by[[2L]],nid,pid,w,array,group_names.GRP(by))))
    } else if(!is.GRP(pid)) pid <- GRP(pid, return.groups = FALSE)
    return(drop(fbstatslCpp(x,higher,by[[1L]],by[[2L]],pid[[1L]],pid[[2L]],w,array,group_names.GRP(by))))
  }
}
qsu.pdata.frame <- function(x, by = NULL, w = NULL, cols = NULL, higher = FALSE, array = TRUE, vlabels = FALSE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  index <- attr(x, "index")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])

  # fastest solution!
  if(is.call(by)) {
    class(x) <- NULL
    if(length(by) == 3L) {
      v <- all.vars(by[[2L]])
      namby <- all.vars(by[[3L]])
    } else {
      namby <- all.vars(by)
      v <- if(is.null(cols)) -match(namby, names(x)) else cols2int(cols, x, names(x))
    }
    by <- if(length(namby) == 1L) x[[namby]] else GRP(x, namby)
    x <- x[v]
  } else if(!is.null(cols)) x <- unclass(x)[cols2int(cols, x, names(x))]

  if(vlabels) names(x) <- paste(names(x), vlabels(x), sep = ": ")

  if(is.null(by))
    return(drop(fbstatslCpp(x,higher,0L,0L,fnlevels(index[[1L]]),index[[1L]],w,array))) else if (is.atomic(by)) {
      if(!is.nmfactor(by)) by <- qF(by, na.exclude = FALSE)
      lev <- attr(by, "levels")
      return(drop(fbstatslCpp(x,higher,length(lev),by,fnlevels(index[[1L]]),index[[1L]],w,array,lev)))
    } else if(!is.GRP(by)) by <- GRP(by)
    return(drop(fbstatslCpp(x,higher,by[[1L]],by[[2L]],fnlevels(index[[1L]]),index[[1L]],w,array,group_names.GRP(by))))
}

# give class "qsu" !! Note: round gives error of list output!
# Also thing about aperm options !!
# print.qsu <- function(x) {
#   print.default(round(x, 2), na.print = "-", digits = 16)
# }

print.qsu <- function(x, digits = 2, nonsci.digits = 9, na.print = "-", return = FALSE, ...) {
  formatfun <- function(x) { # , drop0trailing = FALSE redundat ??
    xx <- formatC(unclass(round(x, digits)), format = "g", digits = nonsci.digits, big.mark = ",", big.interval = 6, ...) # format(unclass(round(x,2)), digits = digits, drop0trailing = TRUE, big.mark = ",", big.interval = 6, scientific = FALSE)
    if(any(ina <- is.na(x))) xx[ina] <- na.print
    return(xx)
  }
  xx <- if(is.atomic(x)) formatfun(x) else rapply(x, formatfun, how = "list") # No longer necessary, but keep, maybe you want to print lists using print.qsu.
  if(return) return(xx) else print.default(xx, quote = FALSE, right = TRUE)
  invisible(x)
}
# View.qsu <- function(x) View(unclass(x))


# testing formula inputs:

# best!!! (subsetting unclassed objects is better), and reassigning x does not take more memry than deleting columns, it is often faster !!, also GRP does not mind unclassed objects !!
# formtest2 <- function(x, by = NULL, xt = NULL, cols = NULL) {
#   formby <- is.call(by)
#   formxt <- is.call(xt)
#
#   # fastest solution: (check: is reassigning x memory efficient ?? should not rather delete columns ??)
#   if(formby || formxt) {
#     v <- NULL
#     class(x) <- NULL # this is faster !!
#     if(formby) {
#       if(length(by) == 3L) {
#         v <- all.vars(by[[2L]])
#         namby <- all.vars(by[[3L]])
#       } else namby <- all.vars(by)
#       by <- if(length(namby) == 1L) x[[namby]] else GRP(x, namby)
#     } else namby <- NULL
#     if(formxt) {
#       if(length(xt) == 3L) {
#         v <- all.vars(xt[[2L]])
#         namxt <- all.vars(xt[[3L]])
#       } else namxt <- all.vars(xt)
#       xt <- if(length(namxt) == 1L) x[[namxt]] else GRP(x, namxt)
#     } else namxt <- NULL
#     if(is.null(v)) { # reassign ?? or set NULL ??? what is more memory efficient ??
#       x <- if(is.null(cols)) x[-match(c(namby,namxt), names(x))] else if(is.function(cols))
#             x[vapply(x, cols, TRUE)] else x[cols]
#     } else x <- x[v]
#   } else if(!is.null(cols)) {
#     # class(x) <- NULL # but unclass is faster !!
#     x <- if(is.function(cols)) unclass(x)[vapply(x, cols, TRUE)] else unclass(x)[cols]
#   }
#   return(list(x, by, xt))
# }
#
# formtest1 <- function(x, by = NULL, xt = NULL, cols = NULL) {
#   formby <- is.call(by)
#   formxt <- is.call(xt)
#
#   # fastest solution: (check: is reassigning x memory efficient ?? should not rather delete columns ??)
#   if(formby || formxt) {
#     v <- NULL
#     if(formby) {
#       if(length(by) == 3L) {
#         v <- all.vars(by[[2L]])
#         namby <- all.vars(by[[3L]])
#       } else namby <- all.vars(by)
#       by <- if(length(namby) == 1L) x[[namby]] else GRP(x, namby)
#     } else namby <- NULL
#     if(formxt) {
#       if(length(xt) == 3L) {
#         v <- all.vars(xt[[2L]])
#         namxt <- all.vars(xt[[3L]])
#       } else namxt <- all.vars(xt)
#       xt <- if(length(namxt) == 1L) x[[namxt]] else GRP(x, namxt)
#     } else namxt <- NULL
#     if(is.null(v)) { # reassign ?? or set NULL ??? what is more memory efficient ??
#       x <- if(is.null(cols)) unclass(x)[-match(c(namby,namxt), names(x))] else if(is.function(cols))
#         unclass(x)[vapply(x, cols, TRUE)] else unclass(x)[cols]
#     } else x <- unclass(x)[v]
#   } else if(!is.null(cols)) {
#     x <- if(is.function(cols)) unclass(x)[vapply(x, cols, TRUE)] else unclass(x)[cols]
#   }
#   return(list(x, by, xt))
# }
#
# formtest3 <- function(x, by = NULL, xt = NULL, cols = NULL) {
#   formby <- is.call(by)
#   formxt <- is.call(xt)
#
#   # fastest solution: (check: is reassigning x memory efficient ?? should not rather delete columns ??)
#   if(formby || formxt) {
#     v <- NULL
#     class(x) <- NULL
#     if(formby) {
#       if(length(by) == 3L) {
#         v <- all.vars(by[[2L]])
#         namby <- all.vars(by[[3L]])
#       } else namby <- all.vars(by)
#       by <- if(length(namby) == 1L) x[[namby]] else GRP(x, namby)
#     } else namby <- NULL
#     if(formxt) {
#       if(length(xt) == 3L) {
#         v <- all.vars(xt[[2L]])
#         namxt <- all.vars(xt[[3L]])
#       } else namxt <- all.vars(xt)
#       xt <- if(length(namxt) == 1L) x[[namxt]] else GRP(x, namxt)
#     } else namxt <- NULL
#     if(is.null(v)) { # reassign ?? or set NULL ??? what is more memory efficient ??
#       if(is.null(cols)) x[match(c(namby,namxt), names(x))] <- NULL else if(is.function(cols))
#         x[!vapply(x, cols, TRUE)] <- NULL else x[-cols] <- NULL
#     } else x[-v] <- NULL
#   } else if(!is.null(cols)) {
#     if(is.function(cols)) x[!vapply(x, cols, TRUE)] <- NULL else x[-cols] <- NULL
#   }
#   return(list(x, by, xt))
# }
#


# fmean.data.frame <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE) {
#     if(is.null(g)) return(fmeanlCpp(x,0L,0L,NULL,w,na.rm,drop)) else if (is.atomic(g)) {
#       if(use.g.names && !inherits(x, "data.table")) {
#         if(!is.factor(g)) g <- qF(g)
#         lev <- attr(g, "levels")
#         return(setRow.names(fmeanlCpp(x,length(lev),g,NULL,w,na.rm), lev))
#       } else {
#         if(is.factor(g)) return(fmeanlCpp(x,fnlevels(g),g,NULL,w,na.rm)) else {
#           g <- qG(g)
#           return(fmeanlCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm))
#         }
#       }
#     } else {
#       if(!all(class(g) == "GRP")) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
#       if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- names.GRP(g)))
#         return(setRow.names(fmeanlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm), if(is.double(groups)) paste0(groups) else groups)) else
#           return(fmeanlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm))
#     }
#   } else {
#     if(is.null(g)) return(TRAlCpp(x,fmeanlCpp(x,0L,0L,NULL,w,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
#       if(is.factor(g)) return(TRAlCpp(x,fmeanlCpp(x,fnlevels(g),g,NULL,w,na.rm),g,TRAtoInt(TRA))) else {
#         g <- qG(g)
#         return(TRAlCpp(x,fmeanlCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm),g,TRAtoInt(TRA)))
#       }
#     } else {
#       if(!all(class(g) == "GRP")) g <- GRP(g, return.groups = FALSE)
#       return(TRAlCpp(x,fmeanlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm),g[[2]],TRAtoInt(TRA)))
#     }
#   }
# }
# fmean.grouped_df <- function(x, w = NULL, TRA = FALSE, na.rm = TRUE, drop.groups = FALSE, drop.w = TRUE, ...) {
#   g <- GRP.grouped_df(x)
#   wsym <- deparse(substitute(w))
#   nam <- names(x)
#   gn <- match(names(g[[4]]), nam)
#   gn2 = gn <- gn[!is.na(gn)]
#   if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
#     w <- x[[wn]]
#     if(any(gn == wn)) stop("Weights coincide with grouping variables!")
#     if(drop.w) if(drop.groups) gn <- c(gn,wn) else gn2 <- c(gn2,wn)
#   }
#   if(length(gn)) {
#     if(drop.groups) {
#       if(TRA == FALSE) return(fmeanlCpp(x[-gn],g[[1]],g[[2]],g[[3]],w,na.rm)) else {
#         x <- x[-gn]
#         return(TRAlCpp(x,fmeanlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm),g[[2]],TRAtoInt(TRA)))
#       }
#     } else {
#       ax <- attributes(x)
#       attributes(x) <- NULL
#       ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn2])
#       if(TRA == FALSE) {
#         ax[["row.names"]] <- .set_row_names(g[[1]])
#         return(`attributes<-`(c(g[[4]],fmeanlCpp(x[-gn2],g[[1]],g[[2]],g[[3]],w,na.rm)), ax))
#       } else
#         return(`attributes<-`(c(x[gn],TRAlCpp(x[-gn2],fmeanlCpp(x[-gn2],g[[1]],g[[2]],g[[3]],w,na.rm),g[[2]],TRAtoInt(TRA))), ax))
#     }
#   } else {
#     if(TRA == FALSE)
#       return(fmeanlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm)) else
#         return(TRAlCpp(x,fmeanlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm),g[[2]],TRAtoInt(TRA)))
#   }
# }

