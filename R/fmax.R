
# For foundational changes to this code see fsum.R !!

fmax <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  UseMethod("fmax", x)
}
fmax.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fmaxCpp(x,0L,0L,na.rm)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fmaxCpp(x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(fmaxCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(fmaxCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fmaxCpp(x,g[[1]],g[[2]],na.rm), group.names.GRP(g))) else 
        return(fmaxCpp(x,g[[1]],g[[2]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fmaxCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fmaxCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fmaxCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fmaxCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fmax.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fmaxmCpp(x,0L,0L,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fmaxmCpp(x,length(lev),g,na.rm), list(lev, dimnames(x)[[2]])))
      } else {
        if(is.factor(g)) return(fmaxmCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(fmaxmCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fmaxmCpp(x,g[[1]],g[[2]],na.rm), list(group.names.GRP(g), dimnames(x)[[2]]))) else
        return(fmaxmCpp(x,g[[1]],g[[2]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fmaxmCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fmaxmCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fmaxmCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fmaxmCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fmax.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fmaxlCpp(x,0L,0L,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fmaxlCpp(x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(fmaxlCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(fmaxlCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group.names.GRP(g))) 
        return(setRow.names(fmaxlCpp(x,g[[1]],g[[2]],na.rm), if(is.double(groups)) paste0(groups) else groups)) else 
          return(fmaxlCpp(x,g[[1]],g[[2]],na.rm)) 
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fmaxlCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fmaxlCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fmaxlCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fmaxlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fmax.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, drop.groups = FALSE, ...) { 
  g <- GRP.grouped_df(x)
  gn <- match(names(g[[4]]), names(x))
  gn <- gn[!is.na(gn)]
  if(length(gn)) {
    if(drop.groups) {
      if(TRA == FALSE) return(fmaxlCpp(x[-gn],g[[1]],g[[2]],na.rm)) else {
        x <- x[-gn] 
        return(TRAlCpp(x,fmaxlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
      }
    } else {
      ax <- attributes(x)
      attributes(x) <- NULL 
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      if(TRA == FALSE) {
        ax[["row.names"]] <- .set_row_names(g[[1]])
        return(`attributes<-`(c(g[[4]],fmaxlCpp(x[-gn],g[[1]],g[[2]],na.rm)), ax))
      } else 
        return(`attributes<-`(c(x[gn],TRAlCpp(x[-gn],fmaxlCpp(x[-gn],g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA))), ax))
    }
  } else {
    if(TRA == FALSE)
      return(fmaxlCpp(x,g[[1]],g[[2]],na.rm)) else 
        return(TRAlCpp(x,fmaxlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
  }
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
#         if(is.atomic(g[["groups"]])) `names<-`(fmaxCpp(x,g[[1]],g[[2]],na.rm),g[["groups"]]) else `names<-`(fmaxCpp(x,g[[1]],g[[2]],na.rm), do.call(paste,c(g[["groups"]],list(sep = ".")))) 
#       } else fmaxCpp(x,g[[1]],g[[2]],na.rm) # speed loss ??
#     }
#   } else {
#     if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRACpp(x,fmaxCpp(x,0L,0L,na.rm),0L,trans) else if (is.factor(g))
#       TRACpp(x,fmaxCpp(x,nlevels(g),g,na.rm),g,trans) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRACpp(x,fmaxmCpp(x,g[[1]],g[[2]],na.rm),g[[2]],trans)
#       }
#   }
# }
# fmax.matrix <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) { 
#   if(trans == FALSE) {
#     if(is.null(g)) fmaxmCpp(x,0L,0L,na.rm,drop) else if (is.factor(g)) {
#       if(use.g.names) {
#         dn <- list(levels(g),dimnames(x)[[2]]) 
#         res <- fmaxmCpp(x,length(dn[[1]]),g,na.rm,drop)
#         dimnames(res) <- dn
#         res
#       } else {
#         res <- fmaxmCpp(x,nlevels(g),g,na.rm,drop)
#         dimnames(res) <- list(NULL,dimnames(x)[[2]])
#         res
#       }
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       res <- fmaxmCpp(x,g[[1]],g[[2]],na.rm,drop)
#       dimnames(res) <- if(use.g.names) {
#         if(is.atomic(g[["groups"]])) list(g[["groups"]], dimnames(x)[[2]]) else list(do.call(paste,c(g[["groups"]],list(sep = "."))), dimnames(x)[[2]])
#       } else list(NULL,dimnames(x)[[2]])
#       res
#     }
#   } else {
#     if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRAmCpp(x,fmaxmCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#       TRAmCpp(x,fmaxmCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAmCpp(x,fmaxmCpp(x,g[[1]],g[[2]],na.rm,TRUE),g[[2]],trans)
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
#         res <- fmaxlCpp(x,g[[1]],g[[2]],na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table") && !is.null(g[["groups"]])) { # necessary, else corrupted df 
#           if(is.atomic(g[["groups"]])) paste0(g[["groups"]]) else do.call(paste,c(g[["groups"]],list(sep = "."))) 
#         } else .set_row_names(g[[1]])
#       }
#     } else {
#       if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fmaxlCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#         TRAlCpp(x,fmaxlCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fmaxlCpp(x,g[[1]],g[[2]],na.rm,TRUE),g[[2]],trans)
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
#           res <- fmaxlCpp(x,g[[1]],g[[2]],na.rm,drop)
#         }
#     } else {
#       if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fmaxlCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#         TRAlCpp(x,fmaxlCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fmaxlCpp(x,g[[1]],g[[2]],na.rm,TRUE),g[[2]],trans)
#         }
#     }
#     attributes(res) <- ax
#     res
#   }
# }



# Previous version (Old):
# library(Rcpp)
# sourceCpp('R/C++/fgmax.cpp')
# sourceCpp('R/C++/fgmaxa.cpp')
# sourceCpp('R/C++/fgmaxl.cpp')
# 
# 
# fmax <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) { 
#   UseMethod("fmax", x)
# }
# fmax.default <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) { 
#   if(is.atomic(g)) fgmaxCpp(x,ng,g,na.rm,return,fill) else fgmaxCpp(x,g[[1]],g[[2]],na.rm,return,fill)
# }
# fmax.matrix <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) { 
#   if(is.atomic(g)) { 
#     if(ng == 0L) 
#       `names<-`(fgmaxmCpp(x,ng,g,na.rm,return,fill), dimnames(x)[[2]]) else {
#         ax <- attributes(x)  
#         res <- fgmaxmCpp(x,ng,g,na.rm,return,fill)
#         if(return == 0L) {
#           ax[["dimnames"]][[1]] <- character(0)
#           ax[["dim"]][[1]] <- ng
#         }
#         attributes(res) <- ax
#         res
#       }
#   } else { 
#     ax <- attributes(x)
#     res <- fgmaxmCpp(x,g[[1]],g[[2]],na.rm,return,fill)
#     if(return == 0L) {
#       ax[["dimnames"]][[1]] <- character(0)
#       ax[["dim"]][[1]] <- g[[1]]
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
#     res <- fgmaxlCpp(x,g[[1]],g[[2]],na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(g[[1]])
#   }
#   attributes(res) <- ax
#   res
# }
# fmax.list <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) { 
#   ax <- attributes(x)
#   res <- if(is.atomic(g)) fgmaxlCpp(x,ng,g,na.rm,return,fill) else fgmaxlCpp(x,g[[1]],g[[2]],na.rm,return,fill)
#   attributes(res) <- ax
#   res
# }
