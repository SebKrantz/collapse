library(Rcpp)
sourceCpp('src/fmedian.cpp')
sourceCpp('src/fmediana.cpp')
sourceCpp('src/fmedianl.cpp')
sourceCpp('src/TRA.cpp')
sourceCpp('src/TRAl.cpp')
sourceCpp('src/TRAa.cpp')

# For foundational changes to this code see fsum.R !!

fmedian <- function(x, ...) { # g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, keep.group_keys = TRUE,
  UseMethod("fmedian", x)
}
fmedian.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fmedianCpp(x,0L,0L,NULL,na.rm)) else if(is.atomic(g)) {
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
      if(use.g.names) return(`names<-`(fmedianCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm), group.names.GRP(g))) else
        return(fmedianCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fmedianCpp(x,0L,0L,NULL,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fmedianCpp(x,fnlevels(g),g,NULL,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fmedianCpp(x,attr(g,"N.groups"),g,NULL,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fmedianCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmedian.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fmedianmCpp(x,0L,0L,NULL,na.rm,drop)) else if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fmedianmCpp(x,length(lev),g,NULL,na.rm), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(fmedianmCpp(x,fnlevels(g),g,NULL,na.rm)) else {
          g <- qG(g)
          return(fmedianmCpp(x,attr(g,"N.groups"),g,NULL,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fmedianmCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm), list(group.names.GRP(g), dimnames(x)[[2L]]))) else
        return(fmedianmCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fmedianmCpp(x,0L,0L,NULL,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fmedianmCpp(x,fnlevels(g),g,NULL,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fmedianmCpp(x,attr(g,"N.groups"),g,NULL,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fmedianmCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmedian.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fmedianlCpp(x,0L,0L,NULL,na.rm,drop)) else if(is.atomic(g)) {
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
        return(setRow.names(fmedianlCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm), groups)) else
          return(fmedianlCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fmedianlCpp(x,0L,0L,NULL,na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fmedianlCpp(x,fnlevels(g),g,NULL,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fmedianlCpp(x,attr(g,"N.groups"),g,NULL,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fmedianlCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fmedian.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, use.g.names = FALSE, keep.group_keys = TRUE, ...) {
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
      ax[["row.names"]] <- if(use.g.names) group.names.GRP(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_keys) {
          ax[["names"]] <- c(g[[5L]], ax[["names"]][-gn])
          return(setAttributes(c(g[[4L]],fmedianlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm)), ax))
        } else {
          ax[["names"]] <- ax[["names"]][-gn]
          return(setAttributes(fmedianlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm), ax))
        }
      } else return(setAttributes(fmedianlCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm), ax))
    } else if(keep.group_keys) {
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      return(setAttributes(c(x[gn],TRAlCpp(x[-gn],fmedianlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(TRAlCpp(x[-gn],fmedianlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(TRAlCpp(x,fmedianlCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TRAtoInt(TRA)))
}
