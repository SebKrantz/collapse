# library(Rcpp)
# sourceCpp('src/fmax.cpp')
# sourceCpp('src/fmaxa.cpp')
# sourceCpp('src/fmaxl.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# For foundational changes to this code see fsum.R !!

fmax <- function(x, ...) { # g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, keep.group_vars = TRUE,
  UseMethod("fmax", x)
}
fmax.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(TRA == FALSE) {
    if(is.null(g)) return(.Call(Cpp_fmax,x,0L,0L,na.rm)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fmax,x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fmax,x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(.Call(Cpp_fmax,x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fmax,x,g[[1L]],g[[2L]],na.rm), group_names.GRP(g))) else
        return(.Call(Cpp_fmax,x,g[[1L]],g[[2L]],na.rm))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fmax,x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fmax,x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fmax,x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fmax,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmax.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(TRA == FALSE) {
    if(is.null(g)) return(.Call(Cpp_fmaxm,x,0L,0L,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fmaxm,x,length(lev),g,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fmaxm,x,fnlevels(g),g,na.rm,FALSE)) else {
          g <- qG(g)
          return(.Call(Cpp_fmaxm,x,attr(g,"N.groups"),g,na.rm,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fmaxm,x,g[[1L]],g[[2L]],na.rm,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fmaxm,x,g[[1L]],g[[2L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fmaxm,x,0L,0L,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fmaxm,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fmaxm,x,attr(g,"N.groups"),g,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fmaxm,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmax.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(TRA == FALSE) {
    if(is.null(g)) return(.Call(Cpp_fmaxl,x,0L,0L,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fmaxl,x,length(lev),g,na.rm,FALSE), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fmaxl,x,fnlevels(g),g,na.rm,FALSE)) else {
          g <- qG(g)
          return(.Call(Cpp_fmaxl,x,attr(g,"N.groups"),g,na.rm,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE), groups)) else
          return(.Call(Cpp_fmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fmaxl,x,0L,0L,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fmaxl,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fmaxl,x,attr(g,"N.groups"),g,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmax.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- GRP.grouped_df(x)
  gn <- which(names(x) %in% g[[5L]])
  nTRAl <- TRA == FALSE
  gl <- length(gn) > 0L
  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      ax[["groups"]] <- NULL
      ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
      ax[["row.names"]] <- if(use.g.names) group_names.GRP(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], ax[["names"]][-gn])
          return(setAttributes(c(g[[4L]],.Call(Cpp_fmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE)), ax))
        } else {
          ax[["names"]] <- ax[["names"]][-gn]
          return(setAttributes(.Call(Cpp_fmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE), ax))
        }
      } else return(setAttributes(.Call(Cpp_fmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fmaxl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fmaxl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
}




# Previous Version:
# fmax <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) {
#   UseMethod("fmax", x)
# }
# fmax.default <- function(x, g = NULL, na.rm = TRUE, trans = FALSE, use.g.names = TRUE, ...) {
#   if(trans == FALSE) {
#     if(is.null(g)) fmaxCpp(x,0L,0L,na.rm) else if (is.factor(g)) {
#       if(use.g.names) {
#         nam <- levels(g)
#         res <- fmaxCpp(x,length(nam),g,na.rm)
#         names(res) <- nam
#         res
#       } else fmaxCpp(x,nlevels(g),g,na.rm)
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(use.g.names) {
#         if(is.atomic(g[["groups"]])) `names<-`(fmaxCpp(x,g[[1L]],g[[2L]],na.rm),g[["groups"]]) else `names<-`(fmaxCpp(x,g[[1L]],g[[2L]],na.rm), do.call(paste,c(g[["groups"]],list(sep = "."))))
#       } else fmaxCpp(x,g[[1L]],g[[2L]],na.rm) # speed loss ??
#     }
#   } else {
#     if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRACpp(x,fmaxCpp(x,0L,0L,na.rm),0L,trans) else if (is.factor(g))
#       TRACpp(x,fmaxCpp(x,nlevels(g),g,na.rm),g,trans) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRACpp(x,fmaxmCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],trans)
#       }
#   }
# }
# fmax.matrix <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) {
#   if(trans == FALSE) {
#     if(is.null(g)) fmaxmCpp(x,0L,0L,na.rm,drop) else if (is.factor(g)) {
#       if(use.g.names) {
#         dn <- list(levels(g),dimnames(x)[[2L]])
#         res <- fmaxmCpp(x,length(dn[[1L]]),g,na.rm,drop)
#         dimnames(res) <- dn
#         res
#       } else {
#         res <- fmaxmCpp(x,nlevels(g),g,na.rm,drop)
#         dimnames(res) <- list(NULL,dimnames(x)[[2L]])
#         res
#       }
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       res <- fmaxmCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#       dimnames(res) <- if(use.g.names) {
#         if(is.atomic(g[["groups"]])) list(g[["groups"]], dimnames(x)[[2L]]) else list(do.call(paste,c(g[["groups"]],list(sep = "."))), dimnames(x)[[2L]])
#       } else list(NULL,dimnames(x)[[2L]])
#       res
#     }
#   } else {
#     if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRAmCpp(x,fmaxmCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#       TRAmCpp(x,fmaxmCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAmCpp(x,fmaxmCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],trans)
#       }
#   }
# }
# fmax.data.frame <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) {
#   if(trans == FALSE && is.null(g) && drop) fmaxlCpp(x,0L,0L,na.rm,drop) else {
#     ax <- attributes(x)
#     if(trans == FALSE) {
#       if(is.null(g) && !drop) {
#         res <- fmaxlCpp(x,0L,0L,na.rm,drop)
#         ax[["row.names"]] <- 1L
#       } else if (is.factor(g)) {
#         lev <- levels(g)
#         res <- fmaxlCpp(x,length(lev),g,na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table")) lev else .set_row_names(length(lev))
#       } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         res <- fmaxlCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table") && !is.null(g[["groups"]])) { # necessary, else corrupted df
#           if(is.atomic(g[["groups"]])) paste0(g[["groups"]]) else do.call(paste,c(g[["groups"]],list(sep = ".")))
#         } else .set_row_names(g[[1L]])
#       }
#     } else {
#       if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fmaxlCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#         TRAlCpp(x,fmaxlCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fmaxlCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],trans)
#         }
#     }
#     attributes(res) <- ax
#     res
#   }
# }
# fmax.list <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, ...) {
#   if(trans == FALSE && is.null(g) && drop) fmaxlCpp(x,0L,0L,na.rm,drop) else {
#     ax <- attributes(x)
#     if(trans == FALSE) {
#       if(is.null(g) && !drop) res <- fmaxlCpp(x,0L,0L,na.rm,drop) else if (is.factor(g))
#         res <- fmaxlCpp(x,nlevels(g),g,na.rm,drop) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           res <- fmaxlCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#         }
#     } else {
#       if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fmaxlCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#         TRAlCpp(x,fmaxlCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fmaxlCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],trans)
#         }
#     }
#     attributes(res) <- ax
#     res
#   }
# }



# Previous version (Old):
# library(Rcpp)
# sourceCpp('src/fgmax.cpp')
# sourceCpp('src/fgmaxa.cpp')
# sourceCpp('src/fgmaxl.cpp')
#
#
# fmax <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   UseMethod("fmax", x)
# }
# fmax.default <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   if(is.atomic(g)) fgmaxCpp(x,ng,g,na.rm,return,fill) else fgmaxCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
# }
# fmax.matrix <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   if(is.atomic(g)) {
#     if(ng == 0L)
#       `names<-`(fgmaxmCpp(x,ng,g,na.rm,return,fill), dimnames(x)[[2L]]) else {
#         ax <- attributes(x)
#         res <- fgmaxmCpp(x,ng,g,na.rm,return,fill)
#         if(return == 0L) {
#           ax[["dimnames"]][[1L]] <- character(0)
#           ax[["dim"]][[1L]] <- ng
#         }
#         attributes(res) <- ax
#         res
#       }
#   } else {
#     ax <- attributes(x)
#     res <- fgmaxmCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
#     if(return == 0L) {
#       ax[["dimnames"]][[1L]] <- character(0)
#       ax[["dim"]][[1L]] <- g[[1L]]
#     }
#     attributes(res) <- ax
#     res
#   }
# }
# fmax.data.frame <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   ax <- attributes(x)
#   if(is.atomic(g)) {
#     res <- fgmaxlCpp(x,ng,g,na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(ifelse(ng == 0L,1L,ng))
#   } else {
#     res <- fgmaxlCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(g[[1L]])
#   }
#   attributes(res) <- ax
#   res
# }
# fmax.list <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   ax <- attributes(x)
#   res <- if(is.atomic(g)) fgmaxlCpp(x,ng,g,na.rm,return,fill) else fgmaxlCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
#   attributes(res) <- ax
#   res
# }
