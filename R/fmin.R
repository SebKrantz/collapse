library(Rcpp)
sourceCpp('R/C++/fmin.cpp')
sourceCpp('R/C++/fmina.cpp')
sourceCpp('R/C++/fminl.cpp')
sourceCpp('R/C++/TRA.cpp')
sourceCpp('R/C++/TRAl.cpp')
sourceCpp('R/C++/TRAa.cpp')

# For foundational changes to this code see fsum.R !!

fmin <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  UseMethod("fmin", x)
}
fmin.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fminCpp(x,0L,0L,na.rm)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fminCpp(x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(fminCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(fminCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fminCpp(x,g[[1]],g[[2]],na.rm), group.names.GRP(g))) else 
        return(fminCpp(x,g[[1]],g[[2]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fminCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fminCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fminCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fminCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fmin.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fminmCpp(x,0L,0L,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fminmCpp(x,length(lev),g,na.rm), list(lev, dimnames(x)[[2]])))
      } else {
        if(is.factor(g)) return(fminmCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(fminmCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fminmCpp(x,g[[1]],g[[2]],na.rm), list(group.names.GRP(g), dimnames(x)[[2]]))) else
        return(fminmCpp(x,g[[1]],g[[2]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fminmCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fminmCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fminmCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fminmCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fmin.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fminlCpp(x,0L,0L,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fminlCpp(x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(fminlCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(fminlCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group.names.GRP(g))) 
        return(setRow.names(fminlCpp(x,g[[1]],g[[2]],na.rm), if(is.double(groups)) paste0(groups) else groups)) else 
          return(fminlCpp(x,g[[1]],g[[2]],na.rm)) 
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fminlCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fminlCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fminlCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fminlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fmin.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, drop.groups = FALSE, ...) { 
  g <- GRP.grouped_df(x)
  gn <- match(names(g[[4]]), names(x))
  gn <- gn[!is.na(gn)]
  if(length(gn)) {
    if(drop.groups) {
      if(TRA == FALSE) return(fminlCpp(x[-gn],g[[1]],g[[2]],na.rm)) else {
        x <- x[-gn] 
        return(TRAlCpp(x,fminlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
      }
    } else {
      ax <- attributes(x)
      attributes(x) <- NULL 
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      if(TRA == FALSE) {
        ax[["row.names"]] <- .set_row_names(g[[1]])
        return(`attributes<-`(c(g[[4]],fminlCpp(x[-gn],g[[1]],g[[2]],na.rm)), ax))
      } else 
        return(`attributes<-`(c(x[gn],TRAlCpp(x[-gn],fminlCpp(x[-gn],g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA))), ax))
    }
  } else {
    if(TRA == FALSE)
      return(fminlCpp(x,g[[1]],g[[2]],na.rm)) else 
        return(TRAlCpp(x,fminlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
  }
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
#         if(is.atomic(g[["groups"]])) `names<-`(fminCpp(x,g[[1]],g[[2]],na.rm),g[["groups"]]) else `names<-`(fminCpp(x,g[[1]],g[[2]],na.rm), do.call(paste,c(g[["groups"]],list(sep = ".")))) 
#       } else fminCpp(x,g[[1]],g[[2]],na.rm) # speed loss ??
#     }
#   } else {
#     if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRACpp(x,fminCpp(x,0L,0L,na.rm),0L,trans) else if (is.factor(g))
#       TRACpp(x,fminCpp(x,nlevels(g),g,na.rm),g,trans) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRACpp(x,fminmCpp(x,g[[1]],g[[2]],na.rm),g[[2]],trans)
#       }
#   }
# }
# fmin.matrix <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) { 
#   if(trans == FALSE) {
#     if(is.null(g)) fminmCpp(x,0L,0L,na.rm,drop) else if (is.factor(g)) {
#       if(use.g.names) {
#         dn <- list(levels(g),dimnames(x)[[2]]) 
#         res <- fminmCpp(x,length(dn[[1]]),g,na.rm,drop)
#         dimnames(res) <- dn
#         res
#       } else {
#         res <- fminmCpp(x,nlevels(g),g,na.rm,drop)
#         dimnames(res) <- list(NULL,dimnames(x)[[2]])
#         res
#       }
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       res <- fminmCpp(x,g[[1]],g[[2]],na.rm,drop)
#       dimnames(res) <- if(use.g.names) {
#         if(is.atomic(g[["groups"]])) list(g[["groups"]], dimnames(x)[[2]]) else list(do.call(paste,c(g[["groups"]],list(sep = "."))), dimnames(x)[[2]])
#       } else list(NULL,dimnames(x)[[2]])
#       res
#     }
#   } else {
#     if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRAmCpp(x,fminmCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#       TRAmCpp(x,fminmCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAmCpp(x,fminmCpp(x,g[[1]],g[[2]],na.rm,TRUE),g[[2]],trans)
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
#         res <- fminlCpp(x,g[[1]],g[[2]],na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table") && !is.null(g[["groups"]])) { # necessary, else corrupted df 
#           if(is.atomic(g[["groups"]])) paste0(g[["groups"]]) else do.call(paste,c(g[["groups"]],list(sep = "."))) 
#         } else .set_row_names(g[[1]])
#       }
#     } else {
#       if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fminlCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#         TRAlCpp(x,fminlCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fminlCpp(x,g[[1]],g[[2]],na.rm,TRUE),g[[2]],trans)
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
#           res <- fminlCpp(x,g[[1]],g[[2]],na.rm,drop)
#         }
#     } else {
#       if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fminlCpp(x,narm = na.rm),0L,trans) else if (is.factor(g))
#         TRAlCpp(x,fminlCpp(x,nlevels(g),g,na.rm,TRUE),g,trans) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fminlCpp(x,g[[1]],g[[2]],na.rm,TRUE),g[[2]],trans)
#         }
#     }
#     attributes(res) <- ax
#     res
#   }
# }



# Previous Version (old)
# library(Rcpp)
# sourceCpp('R/C++/fgmin.cpp')
# sourceCpp('R/C++/fgmina.cpp')
# sourceCpp('R/C++/fgminl.cpp')
# 
# 
# fmin <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) { 
#   UseMethod("fmin", x)
# }
# fmin.default <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) { 
#   if(is.atomic(g)) fgminCpp(x,ng,g,na.rm,return,fill) else fgminCpp(x,g[[1]],g[[2]],na.rm,return,fill)
# }
# fmin.matrix <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) { 
#   if(is.atomic(g)) { 
#     if(ng == 0L) 
#       `names<-`(fgminmCpp(x,ng,g,na.rm,return,fill), dimnames(x)[[2]]) else {
#         ax <- attributes(x)  
#         res <- fgminmCpp(x,ng,g,na.rm,return,fill)
#         if(return == 0L) {
#           ax[["dimnames"]][[1]] <- character(0)
#           ax[["dim"]][[1]] <- ng
#         }
#         attributes(res) <- ax
#         res
#       }
#   } else { 
#     ax <- attributes(x)
#     res <- fgminmCpp(x,g[[1]],g[[2]],na.rm,return,fill)
#     if(return == 0L) {
#       ax[["dimnames"]][[1]] <- character(0)
#       ax[["dim"]][[1]] <- g[[1]]
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
#     res <- fgminlCpp(x,g[[1]],g[[2]],na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(g[[1]])
#   }
#   attributes(res) <- ax
#   res
# }
# fmin.list <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) { 
#   ax <- attributes(x)
#   res <- if(is.atomic(g)) fgminlCpp(x,ng,g,na.rm,return,fill) else fgminlCpp(x,g[[1]],g[[2]],na.rm,return,fill)
#   attributes(res) <- ax
#   res
# }
