# library(Rcpp)
# sourceCpp('src/fNdistinct.cpp')
# sourceCpp('src/fNdistincta.cpp')
# sourceCpp('src/TRA.cpp')
# sourceCpp('src/TRAl.cpp')
# sourceCpp('src/TRAa.cpp')

# Better data.frame funtcion ??
# For foundational changes to this code see fsum.R !!
# also: matrix method needs memory equal to size of the object, while data.frame method does not need any memory !!

fNdistinct <- function(x, ...) { # g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, keep.group_vars = TRUE,
  UseMethod("fNdistinct", x)
}
fNdistinct.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fNdistinctCpp(x, narm = na.rm)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fNdistinctCpp(x,length(lev),g,NULL,na.rm), lev))
      } else {
        if(is.factor(g)) return(fNdistinctCpp(x,fnlevels(g),g,NULL,na.rm)) else {
          g <- qG(g)
          return(fNdistinctCpp(x,attr(g,"N.groups"),g,NULL,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fNdistinctCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm), group_names.GRP(g))) else
        return(fNdistinctCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fNdistinctCpp(x, narm = na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fNdistinctCpp(x,fnlevels(g),g,NULL,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fNdistinctCpp(x,attr(g,"N.groups"),g,NULL,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fNdistinctCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fNdistinct.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fNdistinctmCpp(x,0L,0L,NULL,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fNdistinctmCpp(x,length(lev),g,NULL,na.rm), list(lev, dimnames(x)[[2L]])))
      } else {
        if(is.factor(g)) return(fNdistinctmCpp(x,fnlevels(g),g,NULL,na.rm)) else {
          g <- qG(g)
          return(fNdistinctmCpp(x,attr(g,"N.groups"),g,NULL,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fNdistinctmCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm), list(group_names.GRP(g), dimnames(x)[[2L]]))) else
        return(fNdistinctmCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fNdistinctmCpp(x, narm = na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fNdistinctmCpp(x,fnlevels(g),g,NULL,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fNdistinctmCpp(x,attr(g,"N.groups"),g,NULL,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fNdistinctmCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fNdistinct.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fNdistinctlCpp(x,0L,0L,NULL,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fNdistinctlCpp(x,length(lev),g,NULL,na.rm), lev))
      } else {
        if(is.factor(g)) return(fNdistinctlCpp(x,fnlevels(g),g,NULL,na.rm)) else {
          g <- qG(g)
          return(fNdistinctlCpp(x,attr(g,"N.groups"),g,NULL,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group_names.GRP(g)))
        return(setRow.names(fNdistinctlCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm), groups)) else
          return(fNdistinctlCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fNdistinctlCpp(x, narm = na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fNdistinctlCpp(x,fnlevels(g),g,NULL,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fNdistinctlCpp(x,attr(g,"N.groups"),g,NULL,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fNdistinctlCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TRAtoInt(TRA)))
    }
  }
}
fNdistinct.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
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
          return(setAttributes(c(g[[4L]],fNdistinctlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm)), ax))
        } else {
          ax[["names"]] <- ax[["names"]][-gn]
          return(setAttributes(fNdistinctlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm), ax))
        }
      } else return(setAttributes(fNdistinctlCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      return(setAttributes(c(x[gn],TRAlCpp(x[-gn],fNdistinctlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TRAtoInt(TRA))), ax))
    } else {
      ax[["names"]] <- ax[["names"]][-gn]
      return(setAttributes(TRAlCpp(x[-gn],fNdistinctlCpp(x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TRAtoInt(TRA)), ax))
    }
  } else return(TRAlCpp(x,fNdistinctlCpp(x,g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TRAtoInt(TRA)))
}
