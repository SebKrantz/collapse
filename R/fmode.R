# library(Rcpp)
# sourceCpp('src/fmax.cpp')
# sourceCpp('src/fmode.cpp') # Still test thoroughly !!, make better fmodelCpp !!!
# sourceCpp('src/fmodea.cpp') # Still test thoroughly !!
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# Note: for principal innovations of this code see fsum.R !!

fmode <- function(x, ...) { # g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, keep.group_keys = TRUE, keep.w = TRUE,
  UseMethod("fmode", x)
}
fmode.default <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fmodeCpp(x,0L,0L,NULL,w,na.rm)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fmodeCpp(x,length(lev),g,NULL,w,na.rm), lev))
      } else {
        if(is.factor(g)) return(fmodeCpp(x,fnlevels(g),g,NULL,w,na.rm)) else {
          g <- qG(g)
          return(fmodeCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fmodeCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), group.names.GRP(g))) else
        return(fmodeCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fmodeCpp(x,0L,0L,NULL,w,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fmodeCpp(x,fnlevels(g),g,NULL,w,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fmodeCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fmodeCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmode.matrix <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fmodemCpp(x,0L,0L,NULL,w,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fmodemCpp(x,length(lev),g,NULL,w,na.rm), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(fmodemCpp(x,fnlevels(g),g,NULL,w,na.rm)) else {
          g <- qG(g)
          return(fmodemCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fmodemCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), list(group.names.GRP(g), dimnames(x)[[2L]]))) else
        return(fmodemCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fmodemCpp(x,0L,0L,NULL,w,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fmodemCpp(x,fnlevels(g),g,NULL,w,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fmodemCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fmodemCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmode.data.frame <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) {
      if(drop) return(unlist(fmodelCpp(x,0L,0L,NULL,w,na.rm))) else return(fmodelCpp(x,0L,0L,NULL,w,na.rm))
    } else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fmodelCpp(x,length(lev),g,NULL,w,na.rm), lev))
      } else {
        if(is.factor(g)) return(fmodelCpp(x,fnlevels(g),g,NULL,w,na.rm)) else {
          g <- qG(g)
          return(fmodelCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group.names.GRP(g)))
        return(setRow.names(fmodelCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), groups)) else
          return(fmodelCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm))
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fmodelCpp(x,0L,0L,NULL,w,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fmodelCpp(x,fnlevels(g),g,NULL,w,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fmodelCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fmodelCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmode.grouped_df <- function(x, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = FALSE,
                             keep.group_keys = TRUE, keep.w = TRUE, ...) {
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- names(x)
  gn2 <- gn <- which(nam %in% g[[5L]])
  nTRAl <- TRA == FALSE
  sumw <- NULL

  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- x[[wn]]
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    gn <- c(gn, wn)
    if(keep.w) {
      if(nTRAl) sumw <- `names<-`(list(fmaxCpp(w,g[[1L]],g[[2L]],na.rm)), paste0("max.", wsym)) else if(keep.group_keys)
        gn2 <- gn else sumw <- gn2 <- wn
    }
  }

  gl <- length(gn) > 0L # necessary here, not before !!!

  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      ax[["groups"]] <- NULL
      ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
      ax[["row.names"]] <- if(use.g.names) group.names.GRP(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_keys) {
          ax[["names"]] <- c(g[[5L]], names(sumw), ax[["names"]][-gn])
          return(setAttributes(c(g[[4L]], sumw, fmodelCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm)), ax))
        } else {
          ax[["names"]] <- c(names(sumw), ax[["names"]][-gn])
          return(setAttributes(c(sumw, fmodelCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm)), ax))
        }
      } else return(setAttributes(fmodelCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), ax))
    } else if(keep.group_keys || (keep.w && length(sumw))) {
      ax[["names"]] <- c(ax[["names"]][gn2], ax[["names"]][-gn])
      return(setAttributes(c(x[gn2],TRAlCpp(x[-gn],fmodelCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(TRAlCpp(x[-gn],fmodelCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(TRAlCpp(x,fmodelCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRAtoInt(TRA)))
}



# Previous Versions: (also with old unordered_map fMode)

# fMode <- function(x, g = NULL, w = NULL, na.rm = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   UseMethod("fMode", x)
# }
# fMode.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE) {
#     if(is.null(g)) fModeCpp(x,0L,0L,0L,w,na.rm) else if (is.factor(g)) {
#       if(use.g.names) {
#         nam <- levels(g)
#         res <- if(na.rm) fModeCpp(x,length(nam),g,0L,w,TRUE) else
#           fModeCpp(x,length(nam),g,.Internal(tabulate(g,length(nam))),w,FALSE)
#         names(res) <- nam
#         res
#       } else if(na.rm) fModeCpp(x,nlevels(g),g,0L,w,TRUE) else fModeCpp(x,nlevels(g),g,.Internal(tabulate(g,nlevels(g))),w,FALSE)
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(use.g.names) {
#         if(is.atomic(g[["groups"]])) `names<-`(fModeCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[["groups"]]) else `names<-`(fModeCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm), do.call(paste,c(g[["groups"]],list(sep = "."))))
#       } else fModeCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm) # speed loss ??
#     }
#   } else {
#     if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#     if(is.null(g)) TRACpp(x,fModeCpp(x,0L,0L,0L,w,na.rm),0L,TRA) else if (is.factor(g)) {
#       if(na.rm) TRACpp(x,fModeCpp(x,nlevels(g),g,0L,w,TRUE),g,TRA) else TRACpp(x,fModeCpp(x,nlevels(g),g,.Internal(tabulate(g,nlevels(g))),w,FALSE),g,TRA)
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       TRACpp(x,fModeCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm),g[[2L]],TRA)
#     }
#   }
# }
# fMode.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE) {
#     if(is.null(g)) fModemCpp(x,0L,0L,0L,w,na.rm,drop) else if (is.factor(g)) {
#       if(use.g.names) {
#         dn <- list(levels(g),dimnames(x)[[2L]])
#         res <- if(na.rm) fModemCpp(x,length(dn[[1L]]),g,0L,w,TRUE,drop) else
#           fModemCpp(x,length(dn[[1L]]),g,.Internal(tabulate(g,length(dn[[1L]]))),w,FALSE,drop)
#         dimnames(res) <- dn
#         res
#       } else {
#         res <- if(na.rm) fModemCpp(x,nlevels(g),g,0L,w,TRUE,drop) else fModemCpp(x,nlevels(g),g,.Internal(tabulate(g,nlevels(g))),w,FALSE,drop)
#         dimnames(res) <- list(NULL,dimnames(x)[[2L]])
#         res
#       }
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       res <- fModemCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop)
#       dimnames(res) <- if(use.g.names) {
#         if(is.atomic(g[["groups"]])) list(g[["groups"]], dimnames(x)[[2L]]) else list(do.call(paste,c(g[["groups"]],list(sep = "."))), dimnames(x)[[2L]])
#       } else list(NULL,dimnames(x)[[2L]])
#       res
#     }
#   } else {
#     if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#     if(is.null(g)) TRAmCpp(x,fModemCpp(x,0L,0L,0L,w,na.rm,TRUE),0L,TRA) else if (is.factor(g)) {
#       if(na.rm) TRAmCpp(x,fModemCpp(x,nlevels(g),g,0L,w,TRUE,TRUE),g,TRA) else TRAmCpp(x,fModemCpp(x,nlevels(g),g,.Internal(tabulate(g,nlevels(g))),w,FALSE,TRUE),g,TRA)
#     } else {
#       if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       TRAmCpp(x,fModemCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,TRUE),g[[2L]],TRA)
#     }
#   }
# }
# fMode.data.frame <- function(x, g = NULL, w = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, use.g.names = TRUE, ...) {
#   if(TRA == FALSE && is.null(g) && drop) unlist(fModelCpp(x,0L,0L,0L,w,na.rm,TRUE)) else {
#     if(TRA == FALSE) {
#       ax <- attributes(x)
#       if(is.null(g) && !drop) {
#         res <- fModelCpp(x,0L,0L,0L,w,na.rm,FALSE)
#         ax[["row.names"]] <- 1L
#       } else if (is.factor(g)) {
#         lev <- levels(g)
#         res <- if(na.rm) fModelCpp(x,length(lev),g,0L,w,TRUE,drop) else fModelCpp(x,length(lev),g,.Internal(tabulate(g,length(lev))),w,FALSE,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table")) lev else .set_row_names(length(lev))
#       } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         res <- fModelCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop)
#         ax[["row.names"]] <- if(use.g.names && !inherits(x,"data.table") && !is.null(g[["groups"]])) { # necessary, else corrupted df
#           if(is.atomic(g[["groups"]])) paste0(g[["groups"]]) else do.call(paste,c(g[["groups"]],list(sep = ".")))
#         } else .set_row_names(g[[1L]])
#       }
#       attributes(res) <- ax
#       res
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#       if(is.null(g)) TRAlCpp(x,fModelCpp(x,0L,0L,0L,w,na.rm,TRUE),0L,TRA) else if (is.factor(g)) {
#         if(na.rm) TRAlCpp(x,fModelCpp(x,nlevels(g),g,0L,w,TRUE,TRUE),g,TRA) else TRAlCpp(x,fModelCpp(x,nlevels(g),g,.Internal(tabulate(g,nlevels(g))),w,FALSE,TRUE),g,TRA)
#       } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAlCpp(x,fModelCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,TRUE),g[[2L]],TRA)
#       }
#     }
#   }
# }
# fMode.list <- function(x, g = NULL, w = NULL, na.rm = TRUE, drop = TRUE, TRA = FALSE, ...) {
#   if(TRA == FALSE && is.null(g) && drop) unlist(fModelCpp(x,0L,0L,0L,w,na.rm,TRUE)) else {
#     if(TRA == FALSE) {
#       ax <- attributes(x)
#       if(is.null(g) && !drop) res <- fModelCpp(x,0L,0L,0L,w,na.rm,FALSE) else if (is.factor(g)) {
#         res <- if(na.rm) fModelCpp(x,length(lev),g,0L,w,TRUE,drop) else fModelCpp(x,length(lev),g,.Internal(tabulate(g,length(lev))),w,FALSE,drop)
#       } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         res <- fModelCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,drop)
#       }
#       attributes(res) <- ax
#       res
#     } else {
#       if(is.character(TRA)) TRA <- match(TRA,c("replace_fill","replace","-","-+","/","%","+","*"))
#       if(is.null(g)) TRAlCpp(x,fModelCpp(x,0L,0L,0L,w,na.rm,TRUE),0L,TRA) else if (is.factor(g)) {
#         if(na.rm) TRAlCpp(x,fModelCpp(x,nlevels(g),g,0L,w,TRUE,TRUE),g,TRA) else TRAlCpp(x,fModelCpp(x,nlevels(g),g,.Internal(tabulate(g,nlevels(g))),w,FALSE,TRUE),g,TRA)
#       } else {
#         if(!is.list(g)) stop("g must be a a factor, or a GRP() object, see ?GRP")
#         TRAlCpp(x,fModelCpp(x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,TRUE),g[[2L]],TRA)
#       }
#     }
#   }
# }


# Old:
#     if(is.character(trans)) trans <- match(trans,c("replace.na.fill","replace","subtract","subtract.add.avg","divide","percentage","add","multiply"))
