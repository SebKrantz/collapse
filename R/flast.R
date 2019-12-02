# library(Rcpp)
# sourceCpp('src/flast.cpp')
# sourceCpp('src/flasta.cpp')
# sourceCpp('src/flastl.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# For foundational changes to this code see fsum.R !!

flast <- function(x, ...) { # g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, keep.group_vars = TRUE,
  UseMethod("flast", x)
}
flast.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(TRA == FALSE) {
    if(is.null(g)) return(.Call(Cpp_flast,x,0L,0L,na.rm)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_flast,x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_flast,x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(.Call(Cpp_flast,x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(.Call(Cpp_flast,x,g[[1L]],g[[2L]],na.rm), group_names.GRP(g))) else
        return(.Call(Cpp_flast,x,g[[1L]],g[[2L]],na.rm))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_flast,x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_flast,x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRA,x,.Call(Cpp_flast,x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRA,x,.Call(Cpp_flast,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
flast.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(TRA == FALSE) {
    if(is.null(g)) return(.Call(Cpp_flastm,x,0L,0L,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_flastm,x,length(lev),g,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(.Call(Cpp_flastm,x,fnlevels(g),g,na.rm,FALSE)) else {
          g <- qG(g)
          return(.Call(Cpp_flastm,x,attr(g,"N.groups"),g,na.rm,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(.Call(Cpp_flastm,x,g[[1L]],g[[2L]],na.rm,FALSE), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(.Call(Cpp_flastm,x,g[[1L]],g[[2L]],na.rm,FALSE))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_flastm,x,0L,0L,na.rm,TRUE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_flastm,x,fnlevels(g),g,na.rm,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRAm,x,.Call(Cpp_flastm,x,attr(g,"N.groups"),g,na.rm,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAm,x,.Call(Cpp_flastm,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
flast.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) stop("Unknown argument ", dotstostr(...))
  if(TRA == FALSE) {
    if(is.null(g)) {
      if(drop) return(unlist(.Call(Cpp_flastl,x,0L,0L,na.rm))) else return(.Call(Cpp_flastl,x,0L,0L,na.rm))
    } else if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(.Call(Cpp_flastl,x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(.Call(Cpp_flastl,x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(.Call(Cpp_flastl,x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(.Call(Cpp_flastl,x,g[[1L]],g[[2L]],na.rm), groups)) else
          return(.Call(Cpp_flastl,x,g[[1L]],g[[2L]],na.rm))
    }
  } else {
    if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_flastl,x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_flastl,x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(.Call(Cpp_TRAl,x,.Call(Cpp_flastl,x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(.Call(Cpp_TRAl,x,.Call(Cpp_flastl,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
flast.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
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
          return(setAttributes(c(g[[4L]],.Call(Cpp_flastl,x[-gn],g[[1L]],g[[2L]],na.rm)), ax))
        } else {
          ax[["names"]] <- ax[["names"]][-gn]
          return(setAttributes(.Call(Cpp_flastl,x[-gn],g[[1L]],g[[2L]],na.rm), ax))
        }
      } else return(setAttributes(.Call(Cpp_flastl,x,g[[1L]],g[[2L]],na.rm), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_flastl,x[-gn],g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_flastl,x[-gn],g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_flastl,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
}




# Previous Version:
# # Note: `names<-` is slow !!, but no faster option !!
# flast <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   UseMethod("flast", x)
# }
# flast.default <- function(x, g = NULL, na.rm = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE) {
#     if(is.null(g)) flastCpp(x,0L,0L,na.rm) else if (is.factor(g)) {
#       if(use.g.names) {
#         nam <- attr(g, "levels")
#         res <- flastCpp(x,length(nam),g,na.rm)
#         names(res) <- nam
#         res
#       } else flastCpp(x,nlevels(g),g,na.rm)
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(use.g.names) {
#         if(is.atomic(g[["groups"]])) `names<-`(flastCpp(x,g[[1L]],g[[2L]],na.rm),g[["groups"]]) else `names<-`(flastCpp(x,g[[1L]],g[[2L]],na.rm), do.call(paste,c(g[["groups"]],list(sep = "."))))
#       } else flastCpp(x,g[[1L]],g[[2L]],na.rm) # speed loss ??
#     }
#   } else {
#     if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#     if(is.null(g)) TRACpp(x,flastCpp(x,0L,0L,na.rm),0L,TRA) else if (is.factor(g))
#       TRACpp(x,flastCpp(x,nlevels(g),g,na.rm),g,TRA) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRACpp(x,flastCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRA)
#       }
#   }
# }
# flast.matrix <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE) {
#     if(is.null(g)) flastmCpp(x,0L,0L,na.rm,drop) else if (is.factor(g)) {
#       if(use.g.names) {
#         dn <- list(levels(g),dimnames(x)[[2L]])
#         res <- flastmCpp(x,length(dn[[1L]]),g,na.rm,drop)
#         dimnames(res) <- dn
#         res
#       } else {
#         res <- flastmCpp(x,nlevels(g),g,na.rm,drop)
#         dimnames(res) <- list(NULL,dimnames(x)[[2L]])
#         res
#       }
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       res <- flastmCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#       dimnames(res) <- if(use.g.names) {
#         if(is.atomic(g[["groups"]])) list(g[["groups"]], dimnames(x)[[2L]]) else list(do.call(paste,c(g[["groups"]],list(sep = "."))), dimnames(x)[[2L]])
#       } else list(NULL,dimnames(x)[[2L]])
#       res
#     }
#   } else {
#     if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#     if(is.null(g)) TRAmCpp(x,flastmCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#       TRAmCpp(x,flastmCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAmCpp(x,flastmCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],TRA)
#       }
#   }
# }
# flast.data.frame <- function(x, g = NULL, na.rm = TRUE, drop = FALSE, TRA = FALSE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE && is.null(g) && drop) unlist(flastlCpp(x,0L,0L,na.rm,TRUE)) else {
#     if(TRA == FALSE) {
#       ax <- attributes(x)
#       if(is.null(g) && !drop) {
#         res <- flastlCpp(x,0L,0L,na.rm,drop)
#         ax[["row.names"]] <- 1L
#       } else if (is.factor(g)) {
#         lev <- levels(g)
#         res <- flastlCpp(x,length(lev),g,na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table")) lev else .set_row_names(length(lev))
#       } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         res <- flastlCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table") && !is.null(g[["groups"]])) { # necessary, else corrupted df
#           if(is.atomic(g[["groups"]])) paste0(g[["groups"]]) else do.call(paste,c(g[["groups"]],list(sep = ".")))
#         } else .set_row_names(g[[1L]])
#       }
#       attributes(res) <- ax
#       res
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#       if(is.null(g)) TRAlCpp(x,flastlCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#         TRAlCpp(x,flastlCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,flastlCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],TRA)
#         }
#     }
#   }
# }
# flast.list <- function(x, g = NULL, na.rm = TRUE, drop = FALSE, TRA = FALSE, ...) {
#   if(TRA == FALSE && is.null(g) && drop) unlist(flastlCpp(x,0L,0L,na.rm,TRUE)) else {
#     if(TRA == FALSE) {
#       ax <- attributes(x)
#       if(is.null(g) && !drop) res <- flastlCpp(x,0L,0L,na.rm,drop) else if (is.factor(g))
#         res <- flastlCpp(x,nlevels(g),g,na.rm,drop) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           res <- flastlCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#         }
#       attributes(res) <- ax
#       res
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#       if(is.null(g)) TRAlCpp(x,flastlCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#         TRAlCpp(x,flastlCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,flastlCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],TRA)
#         }
#     }
#   }
# }

# Old:
#     if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
