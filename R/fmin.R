# library(Rcpp)
# sourceCpp('src/fmin.cpp')
# sourceCpp('src/fmina.cpp')
# sourceCpp('src/fminl.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# For foundational changes to this code see fsum.R !!

fmin <- function(x, ...) { # g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, keep.group_vars = TRUE,
  UseMethod("fmin", x)
}
fmin.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(TRA == FALSE) {
    if(is.null(g)) return(.Call(Cpp_fmin,x,0L,0L,na.rm)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fmin,x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fmin,x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(.Call(Cpp_fmin,x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fmin,x,g[[1L]],g[[2L]],na.rm), group_names.GRP(g))) else
        return(.Call(Cpp_fmin,x,g[[1L]],g[[2L]],na.rm))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fmin,x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fmin,x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fmin,x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fmin,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmin.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(TRA == FALSE) {
    if(is.null(g)) return(.Call(Cpp_fminm,x,0L,0L,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fminm,x,length(lev),g,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fminm,x,fnlevels(g),g,na.rm,FALSE)) else {
          g <- qG(g)
          return(.Call(Cpp_fminm,x,attr(g,"N.groups"),g,na.rm,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fminm,x,g[[1L]],g[[2L]],na.rm,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fminm,x,g[[1L]],g[[2L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fminm,x,0L,0L,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fminm,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fminm,x,attr(g,"N.groups"),g,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fminm,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmin.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(TRA == FALSE) {
    if(is.null(g)) return(.Call(Cpp_fminl,x,0L,0L,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fminl,x,length(lev),g,na.rm,FALSE), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fminl,x,fnlevels(g),g,na.rm,FALSE)) else {
          g <- qG(g)
          return(.Call(Cpp_fminl,x,attr(g,"N.groups"),g,na.rm,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE), groups)) else
          return(.Call(Cpp_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fminl,x,0L,0L,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fminl,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fminl,x,attr(g,"N.groups"),g,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmin.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
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
          return(setAttributes(c(g[[4L]],.Call(Cpp_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE)), ax))
        } else {
          ax[["names"]] <- ax[["names"]][-gn]
          return(setAttributes(.Call(Cpp_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE), ax))
        }
      } else return(setAttributes(.Call(Cpp_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fminl,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fminl,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
}




# Previous Version
# fmin <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) {
#   UseMethod("fmin", x)
# }
# fmin.default <- function(x, g = NULL, na.rm = TRUE, trans = FALSE, use.g.names = TRUE, ...) {
#   if(trans == FALSE) {
#     if(is.null(g)) fminCpp(x,0L,0L,na.rm) else if (is.factor(g)) {
#       if(use.g.names) {
#         nam <- levels(g)
#         res <- fminCpp(x,length(nam),g,na.rm)
#         names(res) <- nam
#         res
#       } else fminCpp(x,nlevels(g),g,na.rm)
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(use.g.names) {
#         if(is.atomic(g[["groups"]])) `names<-`(fminCpp(x,g[[1L]],g[[2L]],na.rm),g[["groups"]]) else `names<-`(fminCpp(x,g[[1L]],g[[2L]],na.rm), do.call(paste,c(g[["groups"]],list(sep = "."))))
#       } else fminCpp(x,g[[1L]],g[[2L]],na.rm) # speed loss ??
#     }
#   } else {
#     if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRACpp(x,fminCpp(x,0L,0L,na.rm),0L,trans) else if (is.factor(g))
#       TRACpp(x,fminCpp(x,nlevels(g),g,na.rm),g,trans) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRACpp(x,fminmCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],trans)
#       }
#   }
# }
# fmin.matrix <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) {
#   if(trans == FALSE) {
#     if(is.null(g)) fminmCpp(x,0L,0L,na.rm,drop) else if (is.factor(g)) {
#       if(use.g.names) {
#         dn <- list(levels(g),dimnames(x)[[2L]])
#         res <- fminmCpp(x,length(dn[[1L]]),g,na.rm,drop)
#         dimnames(res) <- dn
#         res
#       } else {
#         res <- fminmCpp(x,nlevels(g),g,na.rm,drop)
#         dimnames(res) <- list(NULL,dimnames(x)[[2L]])
#         res
#       }
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       res <- fminmCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#       dimnames(res) <- if(use.g.names) {
#         if(is.atomic(g[["groups"]])) list(g[["groups"]], dimnames(x)[[2L]]) else list(do.call(paste,c(g[["groups"]],list(sep = "."))), dimnames(x)[[2L]])
#       } else list(NULL,dimnames(x)[[2L]])
#       res
#     }
#   } else {
#     if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRAmCpp(x,fminmCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#       TRAmCpp(x,fminmCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAmCpp(x,fminmCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],trans)
#       }
#   }
# }
# fmin.data.frame <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) {
#   if(trans == FALSE && is.null(g) && drop) fminlCpp(x,0L,0L,na.rm,drop) else {
#     ax <- attributes(x)
#     if(trans == FALSE) {
#       if(is.null(g) && !drop) {
#         res <- fminlCpp(x,0L,0L,na.rm,drop)
#         ax[["row.names"]] <- 1L
#       } else if (is.factor(g)) {
#         lev <- levels(g)
#         res <- fminlCpp(x,length(lev),g,na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table")) lev else .set_row_names(length(lev))
#       } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         res <- fminlCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table") && !is.null(g[["groups"]])) { # necessary, else corrupted df
#           if(is.atomic(g[["groups"]])) paste0(g[["groups"]]) else do.call(paste,c(g[["groups"]],list(sep = ".")))
#         } else .set_row_names(g[[1L]])
#       }
#     } else {
#       if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fminlCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#         TRAlCpp(x,fminlCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fminlCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],trans)
#         }
#     }
#     attributes(res) <- ax
#     res
#   }
# }
# fmin.list <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, ...) {
#   if(trans == FALSE && is.null(g) && drop) fminlCpp(x,0L,0L,na.rm,drop) else {
#     ax <- attributes(x)
#     if(trans == FALSE) {
#       if(is.null(g) && !drop) res <- fminlCpp(x,0L,0L,na.rm,drop) else if (is.factor(g))
#         res <- fminlCpp(x,nlevels(g),g,na.rm,drop) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           res <- fminlCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#         }
#     } else {
#       if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fminlCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#         TRAlCpp(x,fminlCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fminlCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],trans)
#         }
#     }
#     attributes(res) <- ax
#     res
#   }
# }



# Previous Version (old)
# library(Rcpp)
# sourceCpp('src/fgmin.cpp')
# sourceCpp('src/fgmina.cpp')
# sourceCpp('src/fgminl.cpp')
#
#
# fmin <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   UseMethod("fmin", x)
# }
# fmin.default <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   if(is.atomic(g)) fgminCpp(x,ng,g,na.rm,return,fill) else fgminCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
# }
# fmin.matrix <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   if(is.atomic(g)) {
#     if(ng == 0L)
#       `names<-`(fgminmCpp(x,ng,g,na.rm,return,fill), dimnames(x)[[2L]]) else {
#         ax <- attributes(x)
#         res <- fgminmCpp(x,ng,g,na.rm,return,fill)
#         if(return == 0L) {
#           ax[["dimnames"]][[1L]] <- character(0)
#           ax[["dim"]][[1L]] <- ng
#         }
#         attributes(res) <- ax
#         res
#       }
#   } else {
#     ax <- attributes(x)
#     res <- fgminmCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
#     if(return == 0L) {
#       ax[["dimnames"]][[1L]] <- character(0)
#       ax[["dim"]][[1L]] <- g[[1L]]
#     }
#     attributes(res) <- ax
#     res
#   }
# }
# fmin.data.frame <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   ax <- attributes(x)
#   if(is.atomic(g)) {
#     res <- fgminlCpp(x,ng,g,na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(ifelse(ng == 0L,1L,ng))
#   } else {
#     res <- fgminlCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(g[[1L]])
#   }
#   attributes(res) <- ax
#   res
# }
# fmin.list <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   ax <- attributes(x)
#   res <- if(is.atomic(g)) fgminlCpp(x,ng,g,na.rm,return,fill) else fgminlCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
#   attributes(res) <- ax
#   res
# }
