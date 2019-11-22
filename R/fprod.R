# library(Rcpp)
# sourceCpp('src/fprod.cpp')
# sourceCpp('src/fproda.cpp')
# sourceCpp('src/fprodl.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# For foundational changes to this code see fsum.R !!

fprod <- function(x, ...) { # g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, keep.group_vars = TRUE,
  UseMethod("fprod", x)
}
fprod.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fprodCpp(x,0L,0L,na.rm)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fprodCpp(x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(fprodCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(fprodCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fprodCpp(x,g[[1L]],g[[2L]],na.rm), group_names.GRP(g))) else
        return(fprodCpp(x,g[[1L]],g[[2L]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fprodCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fprodCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fprodCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fprodCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fprod.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fprodmCpp(x,0L,0L,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fprodmCpp(x,length(lev),g,na.rm), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(fprodmCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(fprodmCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fprodmCpp(x,g[[1L]],g[[2L]],na.rm), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(fprodmCpp(x,g[[1L]],g[[2L]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fprodmCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fprodmCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fprodmCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fprodmCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fprod.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fprodlCpp(x,0L,0L,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fprodlCpp(x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(fprodlCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(fprodlCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(fprodlCpp(x,g[[1L]],g[[2L]],na.rm), groups)) else
          return(fprodlCpp(x,g[[1L]],g[[2L]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fprodlCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fprodlCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fprodlCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fprodlCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fprod.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
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
          return(setAttributes(c(g[[4L]],fprodlCpp(x[-gn],g[[1L]],g[[2L]],na.rm)), ax))
        } else {
          ax[["names"]] <- ax[["names"]][-gn]
          return(setAttributes(fprodlCpp(x[-gn],g[[1L]],g[[2L]],na.rm), ax))
        }
      } else return(setAttributes(fprodlCpp(x,g[[1L]],g[[2L]],na.rm), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      return(setAttributes(c(x[gn],TRAlCpp(x[-gn],fprodlCpp(x[-gn],g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(TRAlCpp(x[-gn],fprodlCpp(x[-gn],g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(TRAlCpp(x,fprodlCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
}




# Previous Version
# fprod <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) {
#   UseMethod("fprod", x)
# }
# fprod.default <- function(x, g = NULL, na.rm = TRUE, trans = FALSE, use.g.names = TRUE, ...) {
#   if(trans == FALSE) {
#     if(is.null(g)) fprodCpp(x,0L,0L,na.rm) else if (is.factor(g)) {
#       if(use.g.names) {
#         nam <- levels(g)
#         res <- fprodCpp(x,length(nam),g,na.rm)
#         names(res) <- nam
#         res
#       } else fprodCpp(x,nlevels(g),g,na.rm)
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(use.g.names) {
#         if(is.atomic(g[["groups"]])) `names<-`(fprodCpp(x,g[[1L]],g[[2L]],na.rm),g[["groups"]]) else `names<-`(fprodCpp(x,g[[1L]],g[[2L]],na.rm), do.call(paste,c(g[["groups"]],list(sep = "."))))
#       } else fprodCpp(x,g[[1L]],g[[2L]],na.rm) # speed loss ??
#     }
#   } else {
#     if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRACpp(x,fprodCpp(x,0L,0L,na.rm),0L,trans) else if (is.factor(g))
#       TRACpp(x,fprodCpp(x,nlevels(g),g,na.rm),g,trans) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRACpp(x,fprodmCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],trans)
#       }
#   }
# }
# fprod.matrix <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) {
#   if(trans == FALSE) {
#     if(is.null(g)) fprodmCpp(x,0L,0L,na.rm,drop) else if (is.factor(g)) {
#       if(use.g.names) {
#         dn <- list(levels(g),dimnames(x)[[2L]])
#         res <- fprodmCpp(x,length(dn[[1L]]),g,na.rm,drop)
#         dimnames(res) <- dn
#         res
#       } else {
#         res <- fprodmCpp(x,nlevels(g),g,na.rm,drop)
#         dimnames(res) <- list(NULL,dimnames(x)[[2L]])
#         res
#       }
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       res <- fprodmCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#       dimnames(res) <- if(use.g.names) {
#         if(is.atomic(g[["groups"]])) list(g[["groups"]], dimnames(x)[[2L]]) else list(do.call(paste,c(g[["groups"]],list(sep = "."))), dimnames(x)[[2L]])
#       } else list(NULL,dimnames(x)[[2L]])
#       res
#     }
#   } else {
#     if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRAmCpp(x,fprodmCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#       TRAmCpp(x,fprodmCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAmCpp(x,fprodmCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],trans)
#       }
#   }
# }
# fprod.data.frame <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) {
#   if(trans == FALSE && is.null(g) && drop) fprodlCpp(x,0L,0L,na.rm,drop) else {
#     ax <- attributes(x)
#     if(trans == FALSE) {
#       if(is.null(g) && !drop) {
#         res <- fprodlCpp(x,0L,0L,na.rm,drop)
#         ax[["row.names"]] <- 1L
#       } else if (is.factor(g)) {
#         lev <- levels(g)
#         res <- fprodlCpp(x,length(lev),g,na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table")) lev else .set_row_names(length(lev))
#       } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         res <- fprodlCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table") && !is.null(g[["groups"]])) { # necessary, else corrupted df
#           if(is.atomic(g[["groups"]])) paste0(g[["groups"]]) else do.call(paste,c(g[["groups"]],list(sep = ".")))
#         } else .set_row_names(g[[1L]])
#       }
#     } else {
#       if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fprodlCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#         TRAlCpp(x,fprodlCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fprodlCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],trans)
#         }
#     }
#     attributes(res) <- ax
#     res
#   }
# }
# fprod.list <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, ...) {
#   if(trans == FALSE && is.null(g) && drop) fprodlCpp(x,0L,0L,na.rm,drop) else {
#     ax <- attributes(x)
#     if(trans == FALSE) {
#       if(is.null(g) && !drop) res <- fprodlCpp(x,0L,0L,na.rm,drop) else if (is.factor(g))
#         res <- fprodlCpp(x,nlevels(g),g,na.rm,drop) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           res <- fprodlCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#         }
#     } else {
#       if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fprodlCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#         TRAlCpp(x,fprodlCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fprodlCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],trans)
#         }
#     }
#     attributes(res) <- ax
#     res
#   }
# }






# Previous Version (old)
# library(Rcpp)
# sourceCpp('R/C++/fgprod.cpp')
# sourceCpp('R/C++/fgproda.cpp')
# sourceCpp('R/C++/fgprodl.cpp')
#
#
# fprod <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   UseMethod("fprod", x)
# }
# fprod.default <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   if(is.atomic(g)) fgprodCpp(x,ng,g,na.rm,return,fill) else fgprodCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
# }
# fprod.matrix <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   if(is.atomic(g)) {
#     if(ng == 0L)
#       `names<-`(fgprodmCpp(x,ng,g,na.rm,return,fill), dimnames(x)[[2L]]) else {
#         ax <- attributes(x)
#         res <- fgprodmCpp(x,ng,g,na.rm,return,fill)
#         if(return == 0L) {
#           ax[["dimnames"]][[1L]] <- character(0)
#           ax[["dim"]][[1L]] <- ng
#         }
#         attributes(res) <- ax
#         res
#       }
#   } else {
#     ax <- attributes(x)
#     res <- fgprodmCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
#     if(return == 0L) {
#       ax[["dimnames"]][[1L]] <- character(0)
#       ax[["dim"]][[1L]] <- g[[1L]]
#     }
#     attributes(res) <- ax
#     res
#   }
# }
# fprod.data.frame <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   ax <- attributes(x)
#   if(is.atomic(g)) {
#     res <- fgprodlCpp(x,ng,g,na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(ifelse(ng == 0L,1L,ng))
#   } else {
#     res <- fgprodlCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(g[[1L]])
#   }
#   attributes(res) <- ax
#   res
# }
# fprod.list <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   ax <- attributes(x)
#   res <- if(is.atomic(g)) fgprodlCpp(x,ng,g,na.rm,return,fill) else fgprodlCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
#   attributes(res) <- ax
#   res
# }
