# library(Rcpp)
# sourceCpp('src/fsum.cpp')
# sourceCpp('src/fmean.cpp')
# sourceCpp('src/fmeana.cpp')
# sourceCpp('src/fmeanl.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# Note: for principal innovations of this code see fsum.R !!

fmean <- function(x, ...) { # g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, keep.group_keys = TRUE, keep.w = TRUE,
  UseMethod("fmean", x)
}
fmean.default <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fmeanCpp(x,0L,0L,NULL,w,na.rm)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fmeanCpp(x,length(lev),g,NULL,w,na.rm), lev))
      } else {
        if(is.factor(g)) return(fmeanCpp(x,fnlevels(g),g,NULL,w,na.rm)) else {
          g <- qG(g)
          return(fmeanCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fmeanCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), group.names.GRP(g))) else
        return(fmeanCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fmeanCpp(x,0L,0L,NULL,w,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fmeanCpp(x,fnlevels(g),g,NULL,w,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fmeanCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fmeanCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmean.matrix <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fmeanmCpp(x,0L,0L,NULL,w,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fmeanmCpp(x,length(lev),g,NULL,w,na.rm), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(fmeanmCpp(x,fnlevels(g),g,NULL,w,na.rm)) else {
          g <- qG(g)
          return(fmeanmCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fmeanmCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), list(group.names.GRP(g), dimnames(x)[[2L]]))) else
        return(fmeanmCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fmeanmCpp(x,0L,0L,NULL,w,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fmeanmCpp(x,fnlevels(g),g,NULL,w,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fmeanmCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fmeanmCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmean.data.frame <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fmeanlCpp(x,0L,0L,NULL,w,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fmeanlCpp(x,length(lev),g,NULL,w,na.rm), lev))
      } else {
        if(is.factor(g)) return(fmeanlCpp(x,fnlevels(g),g,NULL,w,na.rm)) else {
          g <- qG(g)
          return(fmeanlCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group.names.GRP(g)))
        return(setRow.names(fmeanlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), groups)) else
          return(fmeanlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm))
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fmeanlCpp(x,0L,0L,NULL,w,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fmeanlCpp(x,fnlevels(g),g,NULL,w,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fmeanlCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fmeanlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmean.grouped_df <- function(x, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = FALSE,
                             keep.group_keys = TRUE, keep.w = TRUE, ...) {
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- names(x)
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- TRA == FALSE
  sumw <- NULL

  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- x[[wn]]
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    gn <- c(gn, wn)
    if(keep.w) {
      if(nTRAl) sumw <- `names<-`(list(fsumCpp(w,g[[1L]],g[[2L]],na.rm)), paste0("sum.", wsym)) else if(keep.group_keys)
        gn2 <- gn else sumw <- gn2 <- wn
    }
  }

  gl <- length(gn) > 0L # necessary here, not before !!!

  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      ax[["groups"]] <- NULL
      ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
      ax[["row.names"]] <- if(use.g.names) group.names.GRP(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_keys) {
          ax[["names"]] <- c(g[[5L]], names(sumw), ax[["names"]][-gn])
          return(setAttributes(c(g[[4L]], sumw, fmeanlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm)), ax))
        } else {
          ax[["names"]] <- c(names(sumw), ax[["names"]][-gn])
          return(setAttributes(c(sumw, fmeanlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm)), ax))
        }
      } else return(setAttributes(fmeanlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), ax))
    } else if(keep.group_keys || (keep.w && length(sumw))) {
      ax[["names"]] <- c(ax[["names"]][gn2], ax[["names"]][-gn])
      return(setAttributes(c(x[gn2],TRAlCpp(x[-gn],fmeanlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(TRAlCpp(x[-gn],fmeanlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(TRAlCpp(x,fmeanlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)))
}



# Previous Versions: Only factor and GRP objects allowed !!
# fmean <- function(x, g = NULL, w = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   UseMethod("fmean", x)
# }
# fmean.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE) {
#     if(is.null(g)) fmeanCpp(x,0L,0L,0L,w,na.rm) else if (is.factor(g)) {
#       if(use.g.names) {
#         nam <- levels(g)
#         res <- if(na.rm) fmeanCpp(x,length(nam),g,0L,w,TRUE) else
#           fmeanCpp(x,length(nam),g,.Internal(tabulate(g,length(nam))),w,FALSE)
#         names(res) <- nam
#         res
#       } else if(na.rm) fmeanCpp(x,nlevels(g),g,0L,w,TRUE) else fmeanCpp(x,nlevels(g),g,.Internal(tabulate(g,nlevels(g))),w,FALSE)
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(use.g.names) {
#         if(is.atomic(g[["groups"]])) `names<-`(fmeanCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[["groups"]]) else `names<-`(fmeanCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), do.call(paste,c(g[["groups"]],list(sep = "."))))
#       } else fmeanCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm) # speed loss ??
#     }
#   } else {
#     if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#     if(is.null(g)) TRACpp(x,fmeanCpp(x,0L,0L,0L,w,na.rm),0L,TRA) else if (is.factor(g)) {
#       if(na.rm) TRACpp(x,fmeanCpp(x,nlevels(g),g,0L,w,TRUE),g,TRA) else TRACpp(x,fmeanCpp(x,nlevels(g),g,.Internal(tabulate(g,nlevels(g))),w,FALSE),g,TRA)
#     } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRACpp(x,fmeanCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRA)
#       }
#   }
# }
# fmean.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE) {
#     if(is.null(g)) fmeanmCpp(x,0L,0L,0L,w,na.rm,drop) else if (is.factor(g)) {
#       if(use.g.names) {
#         dn <- list(levels(g),dimnames(x)[[2L]])
#         res <- if(na.rm) fmeanmCpp(x,length(dn[[1L]]),g,0L,w,TRUE,drop) else
#           fmeanmCpp(x,length(dn[[1L]]),g,.Internal(tabulate(g,length(dn[[1L]]))),w,FALSE,drop)
#         dimnames(res) <- dn
#         res
#       } else {
#         res <- if(na.rm) fmeanmCpp(x,nlevels(g),g,0L,w,TRUE,drop) else fmeanmCpp(x,nlevels(g),g,.Internal(tabulate(g,nlevels(g))),w,FALSE,drop)
#         dimnames(res) <- list(NULL,dimnames(x)[[2L]])
#         res
#       }
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       res <- fmeanmCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop)
#       dimnames(res) <- if(use.g.names) {
#         if(is.atomic(g[["groups"]])) list(g[["groups"]], dimnames(x)[[2L]]) else list(do.call(paste,c(g[["groups"]],list(sep = "."))), dimnames(x)[[2L]])
#       } else list(NULL,dimnames(x)[[2L]])
#       res
#     }
#   } else {
#     if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#     if(is.null(g)) TRAmCpp(x,fmeanmCpp(x,0L,0L,0L,w,na.rm,TRUE),0L,TRA) else if (is.factor(g)) {
#       if(na.rm) TRAmCpp(x,fmeanmCpp(x,nlevels(g),g,0L,w,TRUE,TRUE),g,TRA) else TRAmCpp(x,fmeanmCpp(x,nlevels(g),g,.Internal(tabulate(g,nlevels(g))),w,FALSE,TRUE),g,TRA)
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       TRAmCpp(x,fmeanmCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,TRUE),g[[2L]],TRA)
#     }
#   }
# }
# fmean.data.frame <- function(x, g = NULL, w = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE && is.null(g) && drop) fmeanlCpp(x,0L,0L,0L,w,na.rm,TRUE) else {
#     if(TRA == FALSE) {
#       ax <- attributes(x)
#       if(is.null(g) && !drop) {
#         res <- fmeanlCpp(x,0L,0L,0L,w,na.rm,FALSE)
#         ax[["row.names"]] <- 1L
#       } else if (is.factor(g)) {
#         lev <- levels(g)
#         res <- if(na.rm) fmeanlCpp(x,length(lev),g,0L,w,TRUE,drop) else fmeanlCpp(x,length(lev),g,.Internal(tabulate(g,length(lev))),w,FALSE,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table")) lev else .set_row_names(length(lev))
#       } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         res <- fmeanlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table") && !is.null(g[["groups"]])) { # necessary, else corrupted df
#           if(is.atomic(g[["groups"]])) paste0(g[["groups"]]) else do.call(paste,c(g[["groups"]],list(sep = ".")))
#         } else .set_row_names(g[[1L]])
#       }
#       attributes(res) <- ax
#       res
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#       if(is.null(g)) TRAlCpp(x,fmeanlCpp(x,0L,0L,0L,w,na.rm,TRUE),0L,TRA) else if (is.factor(g)) {
#         if(na.rm) TRAlCpp(x,fmeanlCpp(x,nlevels(g),g,0L,w,TRUE,TRUE),g,TRA) else TRAlCpp(x,fmeanlCpp(x,nlevels(g),g,.Internal(tabulate(g,nlevels(g))),w,FALSE,TRUE),g,TRA)
#       } else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fmeanlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,TRUE),g[[2L]],TRA)
#       }
#     }
#   }
# }
# fmean.list <- function(x, g = NULL, w = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, ...) {
#   if(TRA == FALSE && is.null(g) && drop) fmeanlCpp(x,0L,0L,0L,w,na.rm,TRUE) else {
#     ax <- attributes(x)
#     if(TRA == FALSE) {
#       if(is.null(g) && !drop) res <- fmeanlCpp(x,0L,0L,0L,w,na.rm,FALSE) else if (is.factor(g)) {
#         res <- if(na.rm) fmeanlCpp(x,length(lev),g,0L,w,TRUE,drop) else fmeanlCpp(x,length(lev),g,.Internal(tabulate(g,length(lev))),w,FALSE,drop)
#         } else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           res <- fmeanlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop)
#         }
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#       res <- if(is.null(g)) TRAlCpp(x,fmeanlCpp(x,0L,0L,0L,w,na.rm,TRUE),0L,TRA) else if (is.factor(g)) {
#         if(na.rm) TRAlCpp(x,fmeanlCpp(x,nlevels(g),g,0L,w,TRUE,TRUE),g,TRA) else TRAlCpp(x,fmeanlCpp(x,nlevels(g),g,.Internal(tabulate(g,nlevels(g))),w,FALSE,TRUE),g,TRA)
#       } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAlCpp(x,fmeanlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,TRUE),g[[2L]],TRA)
#       }
#     }
#     attributes(res) <- ax
#     res
#   }
# }

# Old:
#     if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))



# Previous Version:
# library(Rcpp)
# sourceCpp('src/fgwmean.cpp')
# sourceCpp('src/fgwmeana.cpp')
# sourceCpp('src/fgwmeanl.cpp')
#
#
# fmean <- function(x, g = 0L, w = NULL, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L, gs = NULL) {
#   UseMethod("fmean", x)
# }
# fmean.default <- function(x, g = 0L, w = NULL, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L, gs = NULL) {
#   # if (ng == 0L && is.null(w)) .External(data.table:::Cfastmean, x, na.rm) else
#   if(is.atomic(g)) fgwmeanCpp(x,ng,g,gs,w,na.rm,return,fill) else fgwmeanCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,return,fill)
# }
# fmean.matrix <- function(x, g = 0L, w = NULL, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L, gs = NULL) {
#   if(is.atomic(g)) {
#     if(ng == 0L)
#       `names<-`(fgwmeanmCpp(x,ng,g,gs,w,na.rm,return,fill), dimnames(x)[[2L]]) else {
#     ax <- attributes(x)
#     res <- fgwmeanmCpp(x,ng,g,gs,w,na.rm,return,fill)
#     if(return == 0L) {
#       ax[["dimnames"]][[1L]] <- character(0)
#       ax[["dim"]][[1L]] <- ng
#     }
#     attributes(res) <- ax
#     res
#     }
#   } else {
#     ax <- attributes(x)
#     res <- fgwmeanmCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,return,fill)
#     if(return == 0L) {
#       ax[["dimnames"]][[1L]] <- character(0)
#       ax[["dim"]][[1L]] <- g[[1L]]
#     }
#     attributes(res) <- ax
#     res
#   }
# }
# fmean.data.frame <- function(x, g = 0L, w = NULL, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L, gs = NULL) {
#   ax <- attributes(x)
#   if(is.atomic(g)) {
#     res <- fgwmeanlCpp(x,ng,g,gs,w,na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(ifelse(ng == 0L,1L,ng))
#   } else {
#     res <- fgwmeanlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(g[[1L]])
#   }
#   attributes(res) <- ax
#   res
# }
# fmean.list <- function(x, g = 0L, w = NULL, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L, gs = NULL) {
#   ax <- attributes(x)
#   res <- if(is.atomic(g)) fgwmeanlCpp(x,ng,g,gs,w,na.rm,return,fill) else fgwmeanlCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,return,fill)
#   attributes(res) <- ax
#   res
# } # make it work for lists with other classes??


# setwd("C:/Users/Sebastian Krantz/Documents/R")
# library(Rcpp)
# sourceCpp('fgmean.cpp')
# sourceCpp('fgmeana.cpp')
# sourceCpp('fgmeanl.cpp')
#
# fmeanOld <- function(x, g = 0L, gs = 0L, na.rm = TRUE, return = 0L) { # parallel argument!!
#   if (is.atomic(x)) { # method dispatch???
#     ld <- length(dim(x))
#     if (ld == 2L) {
#       if (gs == 0L) colMeans(x, na.rm) else fgmeanacpp(x,g,gs,na.rm,return)
#     } else if (ld > 2L) stop("fmean does not support higher dimensional arrays") else {
#      if (gs == 0L) .External(data.table:::Cfastmean, x, na.rm) else fgmeancpp(x,g,gs,na.rm,return)
#     }
#   } else fgmeanlcpp(x,g,gs,na.rm,return)
# }
#
# # method dispatch here slower...
# fmean <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L, gs = 0L) {
#   UseMethod("fmean", x)
# }
# fmean.default <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L, gs = 0L) {
#   if(is.atomic(g)) {
#   if (gs == 0L) .External(data.table:::Cfastmean, x, na.rm) else fgmeancpp(x,g,ng,gs,na.rm,return,fill)
#   } else fgmeancpp(x,g[[2L]],g[[1L]],g[[3L]],na.rm,return,fill)
# }
# fmean.matrix <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L, gs = 0L) {
#   if(is.atomic(g)) {
#   if (ng == 0L && return == 0L) colMeans(x, na.rm) else fgmeanacpp(x,g,ng,gs,na.rm,return,fill)
#   } else fgmeanacpp(x,g[[2L]],g[[1L]],g[[3L]],na.rm,return,fill)
# }
# fmean.data.frame <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L, gs = 0L) {
#   ax <- attributes(x)
#   if(is.atomic(g)) {
#     res <- fgmeanlcpp(x,g,ng,gs,na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(ifelse(ng == 0L,1L,ng))
#   } else { # if (gs == 0L) colMeans(x, na.rm) # slower!!
#     res <- fgmeanlcpp(x,g[[2L]],g[[1L]],g[[3L]],na.rm,return,fill)
#    if(return == 0L) ax[["row.names"]] <- .set_row_names(g[[1L]])
#   }
#   attributes(res) <- ax
#   res
# }
# fmean.list <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L, gs = 0L) {
#   ax <- attributes(x)
#   res <- if(is.atomic(g)) fgmeanlcpp(x,g,ng,gs,na.rm,return,fill) else fgmeanlcpp(x,g[[2L]],g[[1L]],g[[3L]],na.rm,return,fill)
#   attributes(res) <- ax
#   res
# } # make it work for lists with other classes??
#
# fmean <- function(x, g = 0L, gs = 0L, na.rm = TRUE, return = 0L) { # parallel argument!!
#   if (is.atomic(x)) {
#     ld <- length(dim(x))
#     if (ld == 2L) {
#       if (gs == 0L) 1 else 1
#     } else if (ld > 2L) stop("fmean does not support higher dimensional arrays") else {
#       if (gs == 0L) 1 else 1
#     }
#   } # else fgmeanlcpp(x,g,gs,na.rm,return)
# }
#
# m = as.matrix(mtcars)
