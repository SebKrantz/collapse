library(Rcpp)
sourceCpp('R/C++/fsum.cpp')
sourceCpp('R/C++/fsuma.cpp')
sourceCpp('R/C++/fsuml.cpp')
sourceCpp('R/C++/TRA.cpp')
sourceCpp('R/C++/TRAl.cpp')
sourceCpp('R/C++/TRAa.cpp')


fsum <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  UseMethod("fsum", x)
}
fsum.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fsumCpp(x,0L,0L,na.rm)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fsumCpp(x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(fsumCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(fsumCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
      # } else if (.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   if(use.g.names) { 
      #     lev <- attr(g, "levels")
      #     return(`names<-`(fsumCpp(x,length(lev),g,na.rm), lev))
      #   } else return(fsumCpp(x,fnlevels(g),g,na.rm))
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fsumCpp(x,g[[1]],g[[2]],na.rm), group.names.GRP(g))) else 
        return(fsumCpp(x,g[[1]],g[[2]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fsumCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fsumCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fsumCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
    #   g <- interaction(g)
    #   return(TRACpp(x,fsumCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fsumCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fsum.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fsummCpp(x,0L,0L,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fsummCpp(x,length(lev),g,na.rm), list(lev, dimnames(x)[[2]])))
      } else {
        if(is.factor(g)) return(fsummCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(fsummCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
      # } else if (.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   if(use.g.names) { 
      #     lev <- attr(g, "levels")
      #     return(`dimnames<-`(fsummCpp(x,length(lev),g,na.rm), list(lev, dimnames(x)[[2]])))
      #   } else return(fsummCpp(x,fnlevels(g),g,na.rm))
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fsummCpp(x,g[[1]],g[[2]],na.rm), list(group.names.GRP(g), dimnames(x)[[2]]))) else
        return(fsummCpp(x,g[[1]],g[[2]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fsummCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fsummCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fsummCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
      # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   return(TRAmCpp(x,fsummCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fsummCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fsum.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fsumlCpp(x,0L,0L,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fsumlCpp(x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(fsumlCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(fsumlCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
      # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   if(use.g.names && !inherits(x, "data.table")) { 
      #     lev <- attr(g, "levels")
      #     return(setRow.names(fsumlCpp(x,length(lev),g,na.rm), lev))
      #   } else return(fsumlCpp(x,fnlevels(g),g,na.rm))
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group.names.GRP(g))) 
        return(setRow.names(fsumlCpp(x,g[[1]],g[[2]],na.rm), if(is.double(groups)) paste0(groups) else groups)) else 
          return(fsumlCpp(x,g[[1]],g[[2]],na.rm)) 
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fsumlCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fsumlCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fsumlCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
      # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   return(TRAlCpp(x,fsumlCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fsumlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fsum.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, drop.groups = FALSE, ...) { # drop grouping columns argument ?? 
  g <- GRP.grouped_df(x)
  gn <- match(names(g[[4]]), names(x))
  gn <- gn[!is.na(gn)]
  if(length(gn)) {
    if(drop.groups) { # No improvements possibly, else corrupted data.frame !!
      if(TRA == FALSE) return(fsumlCpp(x[-gn],g[[1]],g[[2]],na.rm)) else {
        x <- x[-gn] # A lot faster !! (not doing it two times !!)
        return(TRAlCpp(x,fsumlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
      }
    } else {
      ax <- attributes(x)
      attributes(x) <- NULL # remove attributes is more memory efficient than unclass !!
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      if(TRA == FALSE) {
        ax[["row.names"]] <- .set_row_names(g[[1]])
        return(`attributes<-`(c(g[[4]],fsumlCpp(x[-gn],g[[1]],g[[2]],na.rm)), ax))
      } else 
        return(`attributes<-`(c(x[gn],TRAlCpp(x[-gn],fsumlCpp(x[-gn],g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA))), ax))
    }
  } else {
    if(TRA == FALSE)
      return(fsumlCpp(x,g[[1]],g[[2]],na.rm)) else 
        return(TRAlCpp(x,fsumlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
  }
}



# Older Verions:

# fsum.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) { 
#   if(TRA == FALSE) {
#     if(is.null(g)) fsumCpp(x,0L,0L,na.rm) else if (is.factor(g)) {
#       if(use.g.names) {
#         nam <- levels(g)
#         res <- fsumCpp(x,length(nam),g,na.rm)
#         names(res) <- nam
#         res
#       } else fsumCpp(x,nlevels(g),g,na.rm)
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(use.g.names) { 
#         if(is.atomic(g[["groups"]])) `names<-`(fsumCpp(x,g[[1]],g[[2]],na.rm),g[["groups"]]) else `names<-`(fsumCpp(x,g[[1]],g[[2]],na.rm), do.call(paste,c(g[["groups"]],list(sep = ".")))) 
#       } else fsumCpp(x,g[[1]],g[[2]],na.rm) # speed loss ??
#     }
#   } else {
#     if(is.null(g)) TRACpp(x,fsumCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA)) else if (is.factor(g))
#       TRACpp(x,fsumCpp(x,nlevels(g),g,na.rm),g,TRA) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRACpp(x,fsumCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRA)
#       }
#   }
#   # if(is.factor(g)) g = list(nlevels(g),g) else if (!is.list(g)) stop("g must be a a factor, or a list in which the first element indicates the number of groups, and the second element an integer group-id")
#   # if(TRA == FALSE) fsumCpp(x,g[[1]],g[[2]],na.rm) else {
#   #   if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#   #   TRACpp(x,fsumCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRA)
#   # }
# }

# Another approach to fsum.default specification, but not really 
#   if(is.null(g)) { 
#    if(TRA == FALSE) fsumCpp(x,0L,0L,na.rm) else TRACpp(x,fsumCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))
#   } else if(is.atomic(g)) {
#     if(is.factor(g)) ng <- nlevels(g) else {
#       g <- qG(g)
#       ng <- attr(g, "N.groups")
#     }
#    if(TRA == FALSE) fsumCpp(x,ng,g,na.rm) else 
#    TRACpp(x,fsumCpp(x,ng,g,na.rm),g,TRAtoInt(TRA))
#   } else {
#     if(is.GRP(g)) {
#       if(TRA == FALSE) fsumCpp(x,g[[1]],g[[2]],na.rm) else 
#       TRACpp(x,fsumCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA))
#     } else if(.Internal(islistfactor(g, FALSE))) {
#         g <- interaction(g)
#         ng <- nlevels(g)
#       if(TRA == FALSE) fsumCpp(x,ng,g,na.rm) else 
#       TRACpp(x,fsumCpp(x,ng,g,na.rm),g,TRAtoInt(TRA))
#     } else {
#         g <- GRP(g, return.groups = FALSE)
#       if(TRA == FALSE) fsumCpp(x,g[[1]],g[[2]],na.rm) else 
#       TRACpp(x,fsumCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA))
#     } 
#   }
# }


# fsum.matrix <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, use.g.names = TRUE, ...) { 
#   #  if(return == 0L) fsummCpp(x,g[[1]],g[[2]],na.rm,drop) else TRAmCpp(x,fsummCpp(x,g[[1]],g[[2]],na.rm,TRUE),g[[2]],return) 
#   if(TRA == FALSE) {
#     if(is.null(g)) fsummCpp(x,0L,0L,na.rm,drop) else if (is.factor(g)) {
#       if(use.g.names) {
#         dn <- list(levels(g),dimnames(x)[[2]]) 
#         res <- fsummCpp(x,length(dn[[1]]),g,na.rm,drop)
#         dimnames(res) <- dn
#         res
#       } else {
#         res <- fsummCpp(x,nlevels(g),g,na.rm,drop)
#         dimnames(res) <- list(NULL,dimnames(x)[[2]])
#         res
#       }
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       res <- fsummCpp(x,g[[1]],g[[2]],na.rm,drop)
#       dimnames(res) <- if(use.g.names) {
#         if(is.atomic(g[["groups"]])) list(g[["groups"]], dimnames(x)[[2]]) else list(do.call(paste,c(g[["groups"]],list(sep = "."))), dimnames(x)[[2]])
#       } else list(NULL,dimnames(x)[[2]])
#       res
#     }
#   } else {
#     if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRAmCpp(x,fsummCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#       TRAmCpp(x,fsummCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAmCpp(x,fsummCpp(x,g[[1]],g[[2]],na.rm,TRUE),g[[2]],TRA)
#       }
#   }
#   # mostattributes(res) <- ax # attributes(x) # This is very costly !!
#   # dimnames(res) <- dimnames(x) # Fastest ?? -> do in C++?? -> yes!!
# }
# fsum.data.frameold <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE && is.null(g) && drop) fsumlCpp(x,0L,0L,na.rm,drop) else {
#     ax <- attributes(x)
#     if(TRA == FALSE) {
#       if(is.null(g) && !drop) {
#         res <- fsumlCpp(x,0L,0L,na.rm,drop)
#         ax[["row.names"]] <- 1L
#       } else if (is.factor(g)) {
#         lev <- levels(g)
#         res <- fsumlCpp(x,length(lev),g,na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table")) lev else # tryCatch??
#           .set_row_names(length(lev))
#       } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         res <- fsumlCpp(x,g[[1]],g[[2]],na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table") && !is.null(g[["groups"]])) { # necessary, else corrupted df
#           if(is.atomic(g[["groups"]])) paste0(g[["groups"]]) else do.call(paste,c(g[["groups"]],list(sep = ".")))
#         } else .set_row_names(g[[1]])
#       }
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fsumlCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#         TRAlCpp(x,fsumlCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fsumlCpp(x,g[[1]],g[[2]],na.rm,TRUE),g[[2]],TRA)
#         }
#     }
#     attributes(res) <- ax
#     # if(is.factor(g)) g = list(nlevels(g),g) else if (!is.list(g)) stop("g must be a a factor, or a list in which the first element indicates the number of groups, and the second element an integer group-id")
#     # ng = g[[1]]
#     # if(TRA == FALSE) {
#     #   res <- fsumlCpp(x,ng,g[[2]],na.rm,drop)
#     #   if(!drop || ng != 0L) {
#     #     ax <- attributes(x)
#     #     ax[["row.names"]] <- .set_row_names(ifelse(ng == 0L,1L,ng))
#     #     attributes(res) <- ax
#     #   }
#     # } else {
#     #   if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     #   res <- TRAlCpp(x,fsumlCpp(x,ng,g[[2]],na.rm,TRUE),g[[2]],TRA)
#     #   attributes(res) <- attributes(x)
#     # }
#     res
#   }
# }


# fsum.list <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, ...) { 
#   if(TRA == FALSE && is.null(g) && drop) fsumlCpp(x,0L,0L,na.rm,drop) else { 
#     ax <- attributes(x)
#     if(TRA == FALSE) {
#       if(is.null(g) && !drop) res <- fsumlCpp(x,0L,0L,na.rm,drop) else if (is.factor(g)) 
#         res <- fsumlCpp(x,nlevels(g),g,na.rm,drop) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           res <- fsumlCpp(x,g[[1]],g[[2]],na.rm,drop)
#         }
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fsumlCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#         TRAlCpp(x,fsumlCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fsumlCpp(x,g[[1]],g[[2]],na.rm,TRUE),g[[2]],TRA)
#         }
#     }
#     attributes(res) <- ax
#     res
#   }
#   # if(is.factor(g)) g = list(nlevels(g),g) else if (!is.list(g)) stop("g must be a a factor, or a list in which the first element indicates the number of groups, and the second element an integer group-id")
#   # ng = g[[1]]
#   # if(TRA == FALSE) {
#   #   res <- fsumlCpp(x,ng,g[[2]],na.rm,drop)
#   #   if(!drop || ng != 0L) attributes(res) <- attributes(x)
#   # } else {
#   #   if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#   #   res <- TRAlCpp(x,fsumlCpp(x,ng,g[[2]],na.rm,TRUE),g[[2]],TRA)
#   #   attributes(res) <- attributes(x)
#   # }
#   # res
# }


# OLD:
# fsum <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) { 
#   UseMethod("fsum", x)
# }
# fsum.default <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) { 
#   if(is.atomic(g)) fsumCpp(x,ng,g,na.rm,return,fill) else fsumCpp(x,g[[1]],g[[2]],na.rm,return,fill)
# }
# fsum.matrix <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) { 
#   if(is.atomic(g)) { 
#     if(ng == 0L) 
#       `names<-`(fsummCpp(x,ng,g,na.rm,return,fill), dimnames(x)[[2]]) else {
#         ax <- attributes(x)  
#         res <- fsummCpp(x,ng,g,na.rm,return,fill)
#         if(return == 0L) {
#           ax[["dimnames"]][[1]] <- character(0)
#           ax[["dim"]][[1]] <- ng
#         }
#         attributes(res) <- ax
#         res
#       }
#   } else { 
#     ax <- attributes(x)
#     res <- fsummCpp(x,g[[1]],g[[2]],na.rm,return,fill)
#     if(return == 0L) {
#       ax[["dimnames"]][[1]] <- character(0)
#       ax[["dim"]][[1]] <- g[[1]]
#     }
#     attributes(res) <- ax
#     res
#   }
# }
# fsum.data.frame <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) { 
#   ax <- attributes(x)
#   if(is.atomic(g)) {
#     res <- fsumlCpp(x,ng,g,na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(ifelse(ng == 0L,1L,ng))
#   } else { 
#     res <- fsumlCpp(x,g[[1]],g[[2]],na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(g[[1]])
#   }
#   attributes(res) <- ax
#   res
# }
# fsum.list <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) { 
#   ax <- attributes(x)
#   res <- if(is.atomic(g)) fsumlCpp(x,ng,g,na.rm,return,fill) else fsumlCpp(x,g[[1]],g[[2]],na.rm,return,fill)
#   attributes(res) <- ax
#   res
# }
