# library(Rcpp)
# sourceCpp('src/fsum.cpp')
# sourceCpp('src/fsuma.cpp')
# sourceCpp('src/fsuml.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')


fsum <- function(x, ...) { # g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, keep.group_vars = TRUE
  UseMethod("fsum", x)
}
fsum.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(TRA == FALSE) {
    if(is.null(g)) return(.Call(Cpp_fsum,x,0L,0L,na.rm)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fsum,x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fsum,x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(.Call(Cpp_fsum,x,attr(g,"N.groups"),g,na.rm))
        }
      }
      # } else if (.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   if(use.g.names) {
      #     lev <- attr(g, "levels")
      #     return(`names<-`(.Call(Cpp_fsum,x,length(lev),g,na.rm), lev))
      #   } else return(.Call(Cpp_fsum,x,fnlevels(g),g,na.rm))
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_fsum,x,g[[1L]],g[[2L]],na.rm), group_names.GRP(g))) else
        return(.Call(Cpp_fsum,x,g[[1L]],g[[2L]],na.rm))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
    #   g <- interaction(g)
    #   return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_fsum,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsum.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(TRA == FALSE) {
    if(is.null(g)) return(.Call(Cpp_fsumm,x,0L,0L,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fsumm,x,length(lev),g,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fsumm,x,fnlevels(g),g,na.rm,FALSE)) else {
          g <- qG(g)
          return(.Call(Cpp_fsumm,x,attr(g,"N.groups"),g,na.rm,FALSE))
        }
      }
      # } else if (.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   if(use.g.names) {
      #     lev <- attr(g, "levels")
      #     return(`dimnames<-`(.Call(Cpp_fsumm,x,length(lev),g,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      #   } else return(.Call(Cpp_fsumm,x,fnlevels(g),g,na.rm,FALSE))
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_fsumm,x,g[[1L]],g[[2L]],na.rm,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_fsumm,x,g[[1L]],g[[2L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,0L,0L,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,attr(g,"N.groups"),g,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
      # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_fsumm,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsum.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(TRA == FALSE) {
    if(is.null(g)) return(.Call(Cpp_fsuml,x,0L,0L,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_fsuml,x,length(lev),g,na.rm,FALSE), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_fsuml,x,fnlevels(g),g,na.rm,FALSE)) else {
          g <- qG(g)
          return(.Call(Cpp_fsuml,x,attr(g,"N.groups"),g,na.rm,FALSE))
        }
      }
      # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   if(use.g.names && !inherits(x, "data.table")) {
      #     lev <- attr(g, "levels")
      #     return(setRow.names(.Call(Cpp_fsuml,x,length(lev),g,na.rm,FALSE), lev))
      #   } else return(.Call(Cpp_fsuml,x,fnlevels(g),g,na.rm,FALSE))
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],na.rm,FALSE), groups)) else
          return(.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,0L,0L,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,attr(g,"N.groups"),g,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
      # } else if(.Internal(islistfactor(g, FALSE))) { -> slower than GRP for large data !!
      #   g <- interaction(g)
      #   return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA)))
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fsum.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) { # drop grouping columns argument ??
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  g <- GRP.grouped_df(x)
  gn <- which(names(x) %in% g[[5L]]) # faster than na.rm(match(names(g[[4L]]), names(x)))
  nTRAl <- TRA == FALSE
  gl <- length(gn) > 0L
  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL # remove attributes is more memory efficient than unclass !!
    # if(gl) ax[["names"]] <- if(keep.group_vars) c(g[[5L]], ax[["names"]][-gn]) else ax[["names"]][-gn]
    if(nTRAl) {
      ax[["groups"]] <- NULL # best way ??? -> ye, not slower than following (tested GGDC) ngp <- name(ax) != "groups"
      ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
      ax[["row.names"]] <- if(use.g.names) group_names.GRP(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], ax[["names"]][-gn])
          return(setAttributes(c(g[[4L]],.Call(Cpp_fsuml,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE)), ax))
        } else {
          ax[["names"]] <- ax[["names"]][-gn]
          return(setAttributes(.Call(Cpp_fsuml,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE), ax))
        }
      } else return(setAttributes(.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],na.rm,FALSE), ax))
    } else if(keep.group_vars) {
         ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
         return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fsuml,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
         ax[["names"]] <- ax[["names"]][-gn]
         return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fsuml,x[-gn],g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fsuml,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
}


# Previous grouped_df version.
# fsum.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, drop.groups = FALSE, ...) {
#   g <- GRP.grouped_df(x)
#   gn <- match(names(g[[4]]), names(x))
#   gn <- gn[!is.na(gn)]
#   if(length(gn)) {
#     if(drop.groups) {
#       if(TRA == FALSE) return(fsumlCpp(x[-gn],g[[1]],g[[2]],na.rm)) else {
#         x <- x[-gn]
#         return(TRAlCpp(x,fsumlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
#       }
#     } else {
#       ax <- attributes(x)
#       attributes(x) <- NULL
#       ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
#       if(TRA == FALSE) {
#         ax[["row.names"]] <- .set_row_names(g[[1]])
#         return(`attributes<-`(c(g[[4]],fsumlCpp(x[-gn],g[[1]],g[[2]],na.rm)), ax))
#       } else
#         return(`attributes<-`(c(x[gn],TRAlCpp(x[-gn],fsumlCpp(x[-gn],g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA))), ax))
#     }
#   } else {
#     if(TRA == FALSE)
#       return(fsumlCpp(x,g[[1]],g[[2]],na.rm)) else
#         return(TRAlCpp(x,fsumlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
#   }
# }



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
#         if(is.atomic(g[["groups"]])) `names<-`(fsumCpp(x,g[[1L]],g[[2L]],na.rm),g[["groups"]]) else `names<-`(fsumCpp(x,g[[1L]],g[[2L]],na.rm), do.call(paste,c(g[["groups"]],list(sep = "."))))
#       } else fsumCpp(x,g[[1L]],g[[2L]],na.rm) # speed loss ??
#     }
#   } else {
#     if(is.null(g)) TRACpp(x,fsumCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA)) else if (is.factor(g))
#       TRACpp(x,fsumCpp(x,nlevels(g),g,na.rm),g,TRA) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRACpp(x,fsumCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRA)
#       }
#   }
#   # if(is.factor(g)) g = list(nlevels(g),g) else if (!is.list(g)) stop("g must be a a factor, or a list in which the first element indicates the number of groups, and the second element an integer group-id")
#   # if(TRA == FALSE) fsumCpp(x,g[[1L]],g[[2L]],na.rm) else {
#   #   if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#   #   TRACpp(x,fsumCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRA)
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
#       if(TRA == FALSE) fsumCpp(x,g[[1L]],g[[2L]],na.rm) else
#       TRACpp(x,fsumCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA))
#     } else if(.Internal(islistfactor(g, FALSE))) {
#         g <- interaction(g)
#         ng <- nlevels(g)
#       if(TRA == FALSE) fsumCpp(x,ng,g,na.rm) else
#       TRACpp(x,fsumCpp(x,ng,g,na.rm),g,TRAtoInt(TRA))
#     } else {
#         g <- GRP(g, return.groups = FALSE)
#       if(TRA == FALSE) fsumCpp(x,g[[1L]],g[[2L]],na.rm) else
#       TRACpp(x,fsumCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA))
#     }
#   }
# }


# fsum.matrix <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   #  if(return == 0L) fsummCpp(x,g[[1L]],g[[2L]],na.rm,drop) else TRAmCpp(x,fsummCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],return)
#   if(TRA == FALSE) {
#     if(is.null(g)) fsummCpp(x,0L,0L,na.rm,drop) else if (is.factor(g)) {
#       if(use.g.names) {
#         dn <- list(levels(g),dimnames(x)[[2L]])
#         res <- fsummCpp(x,length(dn[[1L]]),g,na.rm,drop)
#         dimnames(res) <- dn
#         res
#       } else {
#         res <- fsummCpp(x,nlevels(g),g,na.rm,drop)
#         dimnames(res) <- list(NULL,dimnames(x)[[2L]])
#         res
#       }
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       res <- fsummCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#       dimnames(res) <- if(use.g.names) {
#         if(is.atomic(g[["groups"]])) list(g[["groups"]], dimnames(x)[[2L]]) else list(do.call(paste,c(g[["groups"]],list(sep = "."))), dimnames(x)[[2L]])
#       } else list(NULL,dimnames(x)[[2L]])
#       res
#     }
#   } else {
#     if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRAmCpp(x,fsummCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#       TRAmCpp(x,fsummCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAmCpp(x,fsummCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],TRA)
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
#         res <- fsumlCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table") && !is.null(g[["groups"]])) { # necessary, else corrupted df
#           if(is.atomic(g[["groups"]])) paste0(g[["groups"]]) else do.call(paste,c(g[["groups"]],list(sep = ".")))
#         } else .set_row_names(g[[1L]])
#       }
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fsumlCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#         TRAlCpp(x,fsumlCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fsumlCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],TRA)
#         }
#     }
#     attributes(res) <- ax
#     # if(is.factor(g)) g = list(nlevels(g),g) else if (!is.list(g)) stop("g must be a a factor, or a list in which the first element indicates the number of groups, and the second element an integer group-id")
#     # ng = g[[1L]]
#     # if(TRA == FALSE) {
#     #   res <- fsumlCpp(x,ng,g[[2L]],na.rm,drop)
#     #   if(!drop || ng != 0L) {
#     #     ax <- attributes(x)
#     #     ax[["row.names"]] <- .set_row_names(ifelse(ng == 0L,1L,ng))
#     #     attributes(res) <- ax
#     #   }
#     # } else {
#     #   if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     #   res <- TRAlCpp(x,fsumlCpp(x,ng,g[[2L]],na.rm,TRUE),g[[2L]],TRA)
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
#           res <- fsumlCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#         }
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fsumlCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#         TRAlCpp(x,fsumlCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fsumlCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],TRA)
#         }
#     }
#     attributes(res) <- ax
#     res
#   }
#   # if(is.factor(g)) g = list(nlevels(g),g) else if (!is.list(g)) stop("g must be a a factor, or a list in which the first element indicates the number of groups, and the second element an integer group-id")
#   # ng = g[[1L]]
#   # if(TRA == FALSE) {
#   #   res <- fsumlCpp(x,ng,g[[2L]],na.rm,drop)
#   #   if(!drop || ng != 0L) attributes(res) <- attributes(x)
#   # } else {
#   #   if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#   #   res <- TRAlCpp(x,fsumlCpp(x,ng,g[[2L]],na.rm,TRUE),g[[2L]],TRA)
#   #   attributes(res) <- attributes(x)
#   # }
#   # res
# }


# OLD:
# fsum <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   UseMethod("fsum", x)
# }
# fsum.default <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   if(is.atomic(g)) fsumCpp(x,ng,g,na.rm,return,fill) else fsumCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
# }
# fsum.matrix <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   if(is.atomic(g)) {
#     if(ng == 0L)
#       `names<-`(fsummCpp(x,ng,g,na.rm,return,fill), dimnames(x)[[2L]]) else {
#         ax <- attributes(x)
#         res <- fsummCpp(x,ng,g,na.rm,return,fill)
#         if(return == 0L) {
#           ax[["dimnames"]][[1L]] <- character(0)
#           ax[["dim"]][[1L]] <- ng
#         }
#         attributes(res) <- ax
#         res
#       }
#   } else {
#     ax <- attributes(x)
#     res <- fsummCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
#     if(return == 0L) {
#       ax[["dimnames"]][[1L]] <- character(0)
#       ax[["dim"]][[1L]] <- g[[1L]]
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
#     res <- fsumlCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
#     if(return == 0L) ax[["row.names"]] <- .set_row_names(g[[1L]])
#   }
#   attributes(res) <- ax
#   res
# }
# fsum.list <- function(x, g = 0L, na.rm = TRUE, return = 0L, fill = return == 0L, ng = 0L) {
#   ax <- attributes(x)
#   res <- if(is.atomic(g)) fsumlCpp(x,ng,g,na.rm,return,fill) else fsumlCpp(x,g[[1L]],g[[2L]],na.rm,return,fill)
#   attributes(res) <- ax
#   res
# }
