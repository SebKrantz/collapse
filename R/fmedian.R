
# For foundational changes to this code see fsum.R !!

fmedian <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  UseMethod("fmedian", x)
}
fmedian.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fmedianCpp(x,0L,0L,NULL,na.rm)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fmedianCpp(x,length(lev),g,NULL,na.rm), lev))
      } else {
        if(is.factor(g)) return(fmedianCpp(x,fnlevels(g),g,NULL,na.rm)) else {
          g <- qG(g)
          return(fmedianCpp(x,attr(g,"N.groups"),g,NULL,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fmedianCpp(x,g[[1]],g[[2]],g[[3]],na.rm), group.names.GRP(g))) else 
        return(fmedianCpp(x,g[[1]],g[[2]],g[[3]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fmedianCpp(x,0L,0L,NULL,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fmedianCpp(x,fnlevels(g),g,NULL,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fmedianCpp(x,attr(g,"N.groups"),g,NULL,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fmedianCpp(x,g[[1]],g[[2]],g[[3]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fmedian.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fmedianmCpp(x,0L,0L,NULL,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fmedianmCpp(x,length(lev),g,NULL,na.rm), list(lev, dimnames(x)[[2]])))
      } else {
        if(is.factor(g)) return(fmedianmCpp(x,fnlevels(g),g,NULL,na.rm)) else {
          g <- qG(g)
          return(fmedianmCpp(x,attr(g,"N.groups"),g,NULL,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fmedianmCpp(x,g[[1]],g[[2]],g[[3]],na.rm), list(group.names.GRP(g), dimnames(x)[[2]]))) else
        return(fmedianmCpp(x,g[[1]],g[[2]],g[[3]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fmedianmCpp(x,0L,0L,NULL,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fmedianmCpp(x,fnlevels(g),g,NULL,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fmedianmCpp(x,attr(g,"N.groups"),g,NULL,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fmedianmCpp(x,g[[1]],g[[2]],g[[3]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fmedian.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fmedianlCpp(x,0L,0L,NULL,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fmedianlCpp(x,length(lev),g,NULL,na.rm), lev))
      } else {
        if(is.factor(g)) return(fmedianlCpp(x,fnlevels(g),g,NULL,na.rm)) else {
          g <- qG(g)
          return(fmedianlCpp(x,attr(g,"N.groups"),g,NULL,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group.names.GRP(g))) 
        return(setRow.names(fmedianlCpp(x,g[[1]],g[[2]],g[[3]],na.rm), if(is.double(groups)) paste0(groups) else groups)) else 
          return(fmedianlCpp(x,g[[1]],g[[2]],g[[3]],na.rm)) 
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fmedianlCpp(x,0L,0L,NULL,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fmedianlCpp(x,fnlevels(g),g,NULL,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fmedianlCpp(x,attr(g,"N.groups"),g,NULL,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fmedianlCpp(x,g[[1]],g[[2]],g[[3]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fmedian.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, drop.groups = FALSE, ...) { 
  g <- GRP.grouped_df(x)
  gn <- match(names(g[[4]]), names(x))
  gn <- gn[!is.na(gn)]
  if(length(gn)) {
    if(drop.groups) {
      if(TRA == FALSE) return(fmedianlCpp(x[-gn],g[[1]],g[[2]],g[[3]],na.rm)) else {
        x <- x[-gn] 
        return(TRAlCpp(x,fmedianlCpp(x,g[[1]],g[[2]],g[[3]],na.rm),g[[2]],TRAtoInt(TRA)))
      }
    } else {
      ax <- attributes(x)
      attributes(x) <- NULL 
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      if(TRA == FALSE) {
        ax[["row.names"]] <- .set_row_names(g[[1]])
        return(`attributes<-`(c(g[[4]],fmedianlCpp(x[-gn],g[[1]],g[[2]],g[[3]],na.rm)), ax))
      } else 
        return(`attributes<-`(c(x[gn],TRAlCpp(x[-gn],fmedianlCpp(x[-gn],g[[1]],g[[2]],g[[3]],na.rm),g[[2]],TRAtoInt(TRA))), ax))
    }
  } else {
    if(TRA == FALSE)
      return(fmedianlCpp(x,g[[1]],g[[2]],g[[3]],na.rm)) else 
        return(TRAlCpp(x,fmedianlCpp(x,g[[1]],g[[2]],g[[3]],na.rm),g[[2]],TRAtoInt(TRA)))
  }
}
