# library(Rcpp)
# sourceCpp('src/ffirst.cpp')
# sourceCpp('src/ffirsta.cpp')
# sourceCpp('src/ffirstl.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')


# For foundational changes to this code see fsum.R !!

ffirst <- function(x, ...) { # g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, keep.group_vars = TRUE,
  UseMethod("ffirst", x)
}
ffirst.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(ffirstCpp(x,0L,0L,na.rm)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(ffirstCpp(x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(ffirstCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(ffirstCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(ffirstCpp(x,g[[1L]],g[[2L]],na.rm), group_names.GRP(g))) else
        return(ffirstCpp(x,g[[1L]],g[[2L]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,ffirstCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,ffirstCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,ffirstCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,ffirstCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
ffirst.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(ffirstmCpp(x,0L,0L,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(ffirstmCpp(x,length(lev),g,na.rm), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(ffirstmCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(ffirstmCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(ffirstmCpp(x,g[[1L]],g[[2L]],na.rm), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(ffirstmCpp(x,g[[1L]],g[[2L]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,ffirstmCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,ffirstmCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,ffirstmCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,ffirstmCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
ffirst.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) {
      if(drop) return(unlist(ffirstlCpp(x,0L,0L,na.rm))) else return(ffirstlCpp(x,0L,0L,na.rm))
      } else if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(ffirstlCpp(x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(ffirstlCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(ffirstlCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(ffirstlCpp(x,g[[1L]],g[[2L]],na.rm), groups)) else
          return(ffirstlCpp(x,g[[1L]],g[[2L]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,ffirstlCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,ffirstlCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,ffirstlCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,ffirstlCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
ffirst.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
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
          return(setAttributes(c(g[[4L]],ffirstlCpp(x[-gn],g[[1L]],g[[2L]],na.rm)), ax))
        } else {
          ax[["names"]] <- ax[["names"]][-gn]
          return(setAttributes(ffirstlCpp(x[-gn],g[[1L]],g[[2L]],na.rm), ax))
        }
      } else return(setAttributes(ffirstlCpp(x,g[[1L]],g[[2L]],na.rm), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      return(setAttributes(c(x[gn],TRAlCpp(x[-gn],ffirstlCpp(x[-gn],g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(TRAlCpp(x[-gn],ffirstlCpp(x[-gn],g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(TRAlCpp(x,ffirstlCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRAtoInt(TRA)))
}



# Previous Version:
# ffirst <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   UseMethod("ffirst", x)
# }
# ffirst.default <- function(x, g = NULL, na.rm = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE) {
#     if(is.null(g)) ffirstCpp(x,0L,0L,na.rm) else if (is.factor(g)) {
#       if(use.g.names) {
#         nam <- levels(g)
#         res <- ffirstCpp(x,length(nam),g,na.rm)
#         names(res) <- nam
#         res
#       } else ffirstCpp(x,nlevels(g),g,na.rm)
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(use.g.names) {
#         if(is.atomic(g[["groups"]])) `names<-`(ffirstCpp(x,g[[1L]],g[[2L]],na.rm),g[["groups"]]) else `names<-`(ffirstCpp(x,g[[1L]],g[[2L]],na.rm), do.call(paste,c(g[["groups"]],list(sep = "."))))
#       } else ffirstCpp(x,g[[1L]],g[[2L]],na.rm) # speed loss ??
#     }
#   } else {
#     if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#     if(is.null(g)) TRACpp(x,ffirstCpp(x,0L,0L,na.rm),0L,TRA) else if (is.factor(g))
#       TRACpp(x,ffirstCpp(x,nlevels(g),g,na.rm),g,TRA) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRACpp(x,ffirstCpp(x,g[[1L]],g[[2L]],na.rm),g[[2L]],TRA)
#       }
#   }
# }
# ffirst.matrix <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE) {
#     if(is.null(g)) ffirstmCpp(x,0L,0L,na.rm,drop) else if (is.factor(g)) {
#       if(use.g.names) {
#         dn <- list(levels(g),dimnames(x)[[2L]])
#         res <- ffirstmCpp(x,length(dn[[1L]]),g,na.rm,drop)
#         dimnames(res) <- dn
#         res
#       } else {
#         res <- ffirstmCpp(x,nlevels(g),g,na.rm,drop)
#         dimnames(res) <- list(NULL,dimnames(x)[[2L]])
#         res
#       }
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       res <- ffirstmCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#       dimnames(res) <- if(use.g.names) {
#         if(is.atomic(g[["groups"]])) list(g[["groups"]], dimnames(x)[[2L]]) else list(do.call(paste,c(g[["groups"]],list(sep = "."))), dimnames(x)[[2L]])
#       } else list(NULL,dimnames(x)[[2L]])
#       res
#     }
#   } else {
#     if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#     if(is.null(g)) TRAmCpp(x,ffirstmCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#       TRAmCpp(x,ffirstmCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAmCpp(x,ffirstmCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],TRA)
#       }
#   }
# }
# ffirst.data.frame <- function(x, g = NULL, na.rm = TRUE, drop = FALSE, TRA = FALSE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE && is.null(g) && drop) unlist(ffirstlCpp(x,0L,0L,na.rm,TRUE)) else {
#     if(TRA == FALSE) {
#       ax <- attributes(x)
#       if(is.null(g) && !drop) {
#         res <- ffirstlCpp(x,0L,0L,na.rm,drop)
#         ax[["row.names"]] <- 1L
#       } else if (is.factor(g)) {
#         lev <- levels(g)
#         res <- ffirstlCpp(x,length(lev),g,na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table")) lev else .set_row_names(length(lev))
#       } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         res <- ffirstlCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table") && !is.null(g[["groups"]])) { # necessary, else corrupted df
#           if(is.atomic(g[["groups"]])) paste0(g[["groups"]]) else do.call(paste,c(g[["groups"]],list(sep = ".")))
#         } else .set_row_names(g[[1L]])
#       }
#       attributes(res) <- ax
#       res
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#       if(is.null(g)) TRAlCpp(x,ffirstlCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#         TRAlCpp(x,ffirstlCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,ffirstlCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],TRA)
#         }
#     }
#   }
# }
# ffirst.list <- function(x, g = NULL, na.rm = TRUE, drop = FALSE, TRA = FALSE, ...) {
#   if(TRA == FALSE && is.null(g) && drop) unlist(ffirstlCpp(x,0L,0L,na.rm,TRUE)) else {
#     if(TRA == FALSE) {
#       ax <- attributes(x)
#       if(is.null(g) && !drop) res <- ffirstlCpp(x,0L,0L,na.rm,drop) else if (is.factor(g))
#         res <- ffirstlCpp(x,nlevels(g),g,na.rm,drop) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           res <- ffirstlCpp(x,g[[1L]],g[[2L]],na.rm,drop)
#         }
#       attributes(res) <- ax
#       res
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#       if(is.null(g)) TRAlCpp(x,ffirstlCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#         TRAlCpp(x,ffirstlCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,ffirstlCpp(x,g[[1L]],g[[2L]],na.rm,TRUE),g[[2L]],TRA)
#         }
#     }
#   }
# }

# Old:
#     if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
