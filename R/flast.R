library(Rcpp)
sourceCpp('R/C++/flast.cpp')
sourceCpp('R/C++/flasta.cpp')
sourceCpp('R/C++/flastl.cpp')
sourceCpp('R/C++/TRA.cpp')
sourceCpp('R/C++/TRAl.cpp')
sourceCpp('R/C++/TRAa.cpp')

# For foundational changes to this code see fsum.R !!

flast <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  UseMethod("flast", x)
}
flast.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(flastCpp(x,0L,0L,na.rm)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(flastCpp(x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(flastCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(flastCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(flastCpp(x,g[[1]],g[[2]],na.rm), group.names.GRP(g))) else 
        return(flastCpp(x,g[[1]],g[[2]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,flastCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,flastCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,flastCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,flastCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
flast.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(flastmCpp(x,0L,0L,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(flastmCpp(x,length(lev),g,na.rm), list(lev, dimnames(x)[[2]])))
      } else {
        if(is.factor(g)) return(flastmCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(flastmCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(flastmCpp(x,g[[1]],g[[2]],na.rm), list(group.names.GRP(g), dimnames(x)[[2]]))) else
        return(flastmCpp(x,g[[1]],g[[2]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,flastmCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,flastmCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,flastmCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,flastmCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
flast.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) {
      if(drop) return(unlist(flastlCpp(x,0L,0L,na.rm))) else return(flastlCpp(x,0L,0L,na.rm))
    } else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(flastlCpp(x,length(lev),g,na.rm), lev))
      } else {
        if(is.factor(g)) return(flastlCpp(x,fnlevels(g),g,na.rm)) else {
          g <- qG(g)
          return(flastlCpp(x,attr(g,"N.groups"),g,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group.names.GRP(g))) 
        return(setRow.names(flastlCpp(x,g[[1]],g[[2]],na.rm), if(is.double(groups)) paste0(groups) else groups)) else 
          return(flastlCpp(x,g[[1]],g[[2]],na.rm)) 
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,flastlCpp(x,0L,0L,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,flastlCpp(x,fnlevels(g),g,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,flastlCpp(x,attr(g,"N.groups"),g,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,flastlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
flast.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, drop.groups = FALSE, ...) { 
  g <- GRP.grouped_df(x)
  gn <- match(names(g[[4]]), names(x))
  gn <- gn[!is.na(gn)]
  if(length(gn)) {
    if(drop.groups) {
      if(TRA == FALSE) return(flastlCpp(x[-gn],g[[1]],g[[2]],na.rm)) else {
        x <- x[-gn] 
        return(TRAlCpp(x,flastlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
      }
    } else {
      ax <- attributes(x)
      attributes(x) <- NULL 
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      if(TRA == FALSE) {
        ax[["row.names"]] <- .set_row_names(g[[1]])
        return(`attributes<-`(c(g[[4]],flastlCpp(x[-gn],g[[1]],g[[2]],na.rm)), ax))
      } else 
        return(`attributes<-`(c(x[gn],TRAlCpp(x[-gn],flastlCpp(x[-gn],g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA))), ax))
    }
  } else {
    if(TRA == FALSE)
      return(flastlCpp(x,g[[1]],g[[2]],na.rm)) else 
        return(TRAlCpp(x,flastlCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRAtoInt(TRA)))
  }
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
#         if(is.atomic(g[["groups"]])) `names<-`(flastCpp(x,g[[1]],g[[2]],na.rm),g[["groups"]]) else `names<-`(flastCpp(x,g[[1]],g[[2]],na.rm), do.call(paste,c(g[["groups"]],list(sep = ".")))) 
#       } else flastCpp(x,g[[1]],g[[2]],na.rm) # speed loss ??
#     }
#   } else {
#     if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#     if(is.null(g)) TRACpp(x,flastCpp(x,0L,0L,na.rm),0L,TRA) else if (is.factor(g))
#       TRACpp(x,flastCpp(x,nlevels(g),g,na.rm),g,TRA) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRACpp(x,flastCpp(x,g[[1]],g[[2]],na.rm),g[[2]],TRA)
#       }
#   }
# }
# flast.matrix <- function(x, g = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, use.g.names = TRUE, ...) { 
#   if(TRA == FALSE) {
#     if(is.null(g)) flastmCpp(x,0L,0L,na.rm,drop) else if (is.factor(g)) {
#       if(use.g.names) {
#         dn <- list(levels(g),dimnames(x)[[2]]) 
#         res <- flastmCpp(x,length(dn[[1]]),g,na.rm,drop)
#         dimnames(res) <- dn
#         res
#       } else {
#         res <- flastmCpp(x,nlevels(g),g,na.rm,drop)
#         dimnames(res) <- list(NULL,dimnames(x)[[2]])
#         res
#       }
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       res <- flastmCpp(x,g[[1]],g[[2]],na.rm,drop)
#       dimnames(res) <- if(use.g.names) {
#         if(is.atomic(g[["groups"]])) list(g[["groups"]], dimnames(x)[[2]]) else list(do.call(paste,c(g[["groups"]],list(sep = "."))), dimnames(x)[[2]])
#       } else list(NULL,dimnames(x)[[2]])
#       res
#     }
#   } else {
#     if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#     if(is.null(g)) TRAmCpp(x,flastmCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#       TRAmCpp(x,flastmCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAmCpp(x,flastmCpp(x,g[[1]],g[[2]],na.rm,TRUE),g[[2]],TRA)
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
#         res <- flastlCpp(x,g[[1]],g[[2]],na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table") && !is.null(g[["groups"]])) { # necessary, else corrupted df 
#           if(is.atomic(g[["groups"]])) paste0(g[["groups"]]) else do.call(paste,c(g[["groups"]],list(sep = "."))) 
#         } else .set_row_names(g[[1]])
#       }
#       attributes(res) <- ax
#       res
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#       if(is.null(g)) TRAlCpp(x,flastlCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#         TRAlCpp(x,flastlCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,flastlCpp(x,g[[1]],g[[2]],na.rm,TRUE),g[[2]],TRA)
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
#           res <- flastlCpp(x,g[[1]],g[[2]],na.rm,drop)
#         }
#       attributes(res) <- ax
#       res
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#       if(is.null(g)) TRAlCpp(x,flastlCpp(x,narm = na.rm),0L,TRA) else if (is.factor(g))
#         TRAlCpp(x,flastlCpp(x,nlevels(g),g,na.rm,TRUE),g,TRA) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,flastlCpp(x,g[[1]],g[[2]],na.rm,TRUE),g[[2]],TRA)
#         }
#     }
#   }
# }

# Old: 
#     if(is.character(TRA)) TRA <- match(TRA,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))