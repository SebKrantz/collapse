library(Rcpp)
sourceCpp('R/C++/fnobs.cpp')
sourceCpp('R/C++/fnobsa.cpp')
sourceCpp('R/C++/fnobsl.cpp')
sourceCpp('R/C++/TRA.cpp')
sourceCpp('R/C++/TRAl.cpp')
sourceCpp('R/C++/TRAa.cpp')

# For foundational changes to this code see fsum.R !!

# Still check: Some weird things going on with some functions. in default there is still a metric cpp version !!
# Rename to fN instead of fnobs ?? or fNobs -> best of both worlds!!

fNobs <- function(x, g = NULL, TRA = FALSE, use.g.names = TRUE, ...) {
  UseMethod("fNobs", x)
}
fNobs.default <- function(x, g = NULL, TRA = FALSE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fnobsCpp(x,0L,0L)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fnobsCpp(x,length(lev),g), lev))
      } else {
        if(is.factor(g)) return(fnobsCpp(x,fnlevels(g),g)) else {
          g <- qG(g)
          return(fnobsCpp(x,attr(g,"N.groups"),g))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fnobsCpp(x,g[[1]],g[[2]]), group.names.GRP(g))) else 
        return(fnobsCpp(x,g[[1]],g[[2]]))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fnobsCpp(x,0L,0L),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fnobsCpp(x,fnlevels(g),g),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fnobsCpp(x,attr(g,"N.groups"),g),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fnobsCpp(x,g[[1]],g[[2]]),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fNobs.matrix <- function(x, g = NULL, TRA = FALSE, drop = TRUE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fnobsmCpp(x,0L,0L,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fnobsmCpp(x,length(lev),g), list(lev, dimnames(x)[[2]])))
      } else {
        if(is.factor(g)) return(fnobsmCpp(x,fnlevels(g),g)) else {
          g <- qG(g)
          return(fnobsmCpp(x,attr(g,"N.groups"),g))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fnobsmCpp(x,g[[1]],g[[2]]), list(group.names.GRP(g), dimnames(x)[[2]]))) else
        return(fnobsmCpp(x,g[[1]],g[[2]]))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fnobsmCpp(x,0L,0L),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fnobsmCpp(x,fnlevels(g),g),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fnobsmCpp(x,attr(g,"N.groups"),g),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fnobsmCpp(x,g[[1]],g[[2]]),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fNobs.data.frame <- function(x, g = NULL, TRA = FALSE, drop = TRUE, use.g.names = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fnobslCpp(x,0L,0L,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fnobslCpp(x,length(lev),g), lev))
      } else {
        if(is.factor(g)) return(fnobslCpp(x,fnlevels(g),g)) else {
          g <- qG(g)
          return(fnobslCpp(x,attr(g,"N.groups"),g))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group.names.GRP(g))) 
        return(setRow.names(fnobslCpp(x,g[[1]],g[[2]]), if(is.double(groups)) paste0(groups) else groups)) else 
          return(fnobslCpp(x,g[[1]],g[[2]])) 
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fnobslCpp(x,0L,0L),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fnobslCpp(x,fnlevels(g),g),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fnobslCpp(x,attr(g,"N.groups"),g),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fnobslCpp(x,g[[1]],g[[2]]),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fNobs.grouped_df <- function(x, TRA = FALSE, drop.groups = FALSE, ...) { 
  g <- GRP.grouped_df(x)
  gn <- match(names(g[[4]]), names(x))
  gn <- gn[!is.na(gn)]
  if(length(gn)) {
    if(drop.groups) {
      if(TRA == FALSE) return(fnobslCpp(x[-gn],g[[1]],g[[2]])) else {
        x <- x[-gn] 
        return(TRAlCpp(x,fnobslCpp(x,g[[1]],g[[2]]),g[[2]],TRAtoInt(TRA)))
      }
    } else {
      ax <- attributes(x)
      attributes(x) <- NULL 
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      if(TRA == FALSE) {
        ax[["row.names"]] <- .set_row_names(g[[1]])
        return(`attributes<-`(c(g[[4]],fnobslCpp(x[-gn],g[[1]],g[[2]])), ax))
      } else 
        return(`attributes<-`(c(x[gn],TRAlCpp(x[-gn],fnobslCpp(x[-gn],g[[1]],g[[2]]),g[[2]],TRAtoInt(TRA))), ax))
    }
  } else {
    if(TRA == FALSE)
      return(fnobslCpp(x,g[[1]],g[[2]])) else 
        return(TRAlCpp(x,fnobslCpp(x,g[[1]],g[[2]]),g[[2]],TRAtoInt(TRA)))
  }
}



# Previous Versions:
# fNobs <- function(x, g = NULL, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) {
#   UseMethod("fnobs", x)
# }
# fNobs.default <- function(x, g = NULL, trans = FALSE, use.g.names = TRUE, ...) { 
#   if(trans == FALSE) {
#     if(is.null(g)) fnobsCpp(x,0L,0L) else if (is.factor(g)) {
#       if(use.g.names) {
#         nam <- levels(g)
#         res <- fnobsCpp(x,length(nam),g)
#         names(res) <- nam
#         res
#       } else fnobsCpp(x,nlevels(g),g)
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(use.g.names) { 
#         if(is.atomic(g[["groups"]])) `names<-`(fnobsCpp(x,g[[1]],g[[2]]),g[["groups"]]) else `names<-`(fnobsCpp(x,g[[1]],g[[2]]), do.call(paste,c(g[["groups"]],list(sep = ".")))) 
#       } else fnobsCpp(x,g[[1]],g[[2]]) # speed loss ??
#     }
#   } else {
#     if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRACpp(x,fnobsCpp(x,0L,0L),0L,trans) else if (is.factor(g))
#       TRACpp(x,fnobsCpp(x,nlevels(g),g),g,trans) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRACpp(x,fnobsCpp(x,g[[1]],g[[2]]),g[[2]],trans)
#       }
#   }
# }
# fNobs.matrix <- function(x, g = NULL, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) { 
#   if(trans == FALSE) {
#     if(is.null(g)) fnobsmCpp(x,0L,0L,drop) else if (is.factor(g)) {
#       if(use.g.names) {
#         dn <- list(levels(g),dimnames(x)[[2]]) 
#         res <- fnobsmCpp(x,length(dn[[1]]),g,drop)
#         dimnames(res) <- dn
#         res
#       } else {
#         res <- fnobsmCpp(x,nlevels(g),g,drop)
#         dimnames(res) <- list(NULL,dimnames(x)[[2]])
#         res
#       }
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       res <- fnobsmCpp(x,g[[1]],g[[2]],drop)
#       dimnames(res) <- if(use.g.names) {
#         if(is.atomic(g[["groups"]])) list(g[["groups"]], dimnames(x)[[2]]) else list(do.call(paste,c(g[["groups"]],list(sep = "."))), dimnames(x)[[2]])
#       } else list(NULL,dimnames(x)[[2]])
#       res
#     }
#   } else {
#     if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#     if(is.null(g)) TRAmCpp(x,fnobsmCpp(x),0L,trans) else if (is.factor(g))
#       TRAmCpp(x,fnobsmCpp(x,nlevels(g),g,TRUE),g,trans) else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAmCpp(x,fnobsmCpp(x,g[[1]],g[[2]],TRUE),g[[2]],trans)
#       }
#   }
# }
# fNobs.data.frame <- function(x, g = NULL, drop = TRUE, trans = FALSE, use.g.names = TRUE, ...) { 
#   if(trans == FALSE && is.null(g) && drop) fnobslCpp(x) else { 
#     ax <- attributes(x)
#     if(trans == FALSE) {
#       if(is.null(g) && !drop) { 
#         res <- fnobslCpp(x,0L,0L,FALSE) 
#         ax[["row.names"]] <- 1L
#       } else if (is.factor(g)) {
#         lev <- levels(g)
#         res <- fnobslCpp(x,length(lev),g,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table")) lev else .set_row_names(length(lev))
#       } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         res <- fnobslCpp(x,g[[1]],g[[2]],drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table") && !is.null(g[["groups"]])) { # necessary, else corrupted df 
#           if(is.atomic(g[["groups"]])) paste0(g[["groups"]]) else do.call(paste,c(g[["groups"]],list(sep = "."))) 
#         } else .set_row_names(g[[1]])
#       }
#     } else {
#       if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fnobslCpp(x),0L,trans) else if (is.factor(g))
#         TRAlCpp(x,fnobslCpp(x,nlevels(g),g,TRUE),g,trans) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fnobslCpp(x,g[[1]],g[[2]],TRUE),g[[2]],trans)
#         }
#     }
#     attributes(res) <- ax
#     res
#   }
# }
# fNobs.list <- function(x, g = NULL, drop = TRUE, trans = FALSE, ...) { 
#   if(trans == FALSE && is.null(g) && drop) fnobslCpp(x) else { 
#     ax <- attributes(x)
#     if(trans == FALSE) {
#       if(is.null(g) && !drop) res <- fnobslCpp(x,0L,0L,FALSE) else if (is.factor(g)) 
#         res <- fnobslCpp(x,nlevels(g),g,drop) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           res <- fnobslCpp(x,g[[1]],g[[2]],drop)
#         }
#     } else {
#       if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
#       res <- if(is.null(g)) TRAlCpp(x,fnobslCpp(x),0L,trans) else if (is.factor(g))
#         TRAlCpp(x,fnobslCpp(x,nlevels(g),g,TRUE),g,trans) else {
#           if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#           TRAlCpp(x,fnobslCpp(x,g[[1]],g[[2]],TRUE),g[[2]],trans)
#         }
#     }
#     attributes(res) <- ax
#     res
#   }
# }
