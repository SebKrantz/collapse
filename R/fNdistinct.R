

# For foundational changes to this code see fsum.R
# Note: matrix method needs memory equal to size of the object, while data.frame method does not need any memory ?!

fNdistinct <- function(x, ...) UseMethod("fNdistinct") # , x

fNdistinct.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fNdistinct.matrix(x, g, TRA, na.rm, use.g.names, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fNdistinct,x,0L,0L,NULL,na.rm))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fNdistinct,x,length(lev),g,NULL,na.rm), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fNdistinct,x,fnlevels(g),g,NULL,na.rm))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fNdistinct,x,attr(g,"N.groups"),g,NULL,na.rm))
    }
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(Cpp_fNdistinct,x,g[[1L]],g[[2L]],g[[3L]],na.rm), GRPnames(g)))
    return(.Call(Cpp_fNdistinct,x,g[[1L]],g[[2L]],g[[3L]],na.rm))
  }
  if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fNdistinct,x,0L,0L,NULL,na.rm),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fNdistinct,x,fnlevels(g),g,NULL,na.rm),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRA,x,.Call(Cpp_fNdistinct,x,attr(g,"N.groups"),g,NULL,na.rm),g,TtI(TRA)))
  }
  if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRA,x,.Call(Cpp_fNdistinct,x,g[[1L]],g[[2L]],g[[3L]],na.rm),g[[2L]],TtI(TRA))
}

fNdistinct.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fNdistinctm,x,0L,0L,NULL,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fNdistinctm,x,length(lev),g,NULL,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fNdistinctm,x,fnlevels(g),g,NULL,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fNdistinctm,x,attr(g,"N.groups"),g,NULL,na.rm,FALSE))
    }
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(Cpp_fNdistinctm,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(Cpp_fNdistinctm,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE))
  }
  if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fNdistinctm,x,0L,0L,NULL,na.rm,TRUE),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fNdistinctm,x,fnlevels(g),g,NULL,na.rm,FALSE),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRAm,x,.Call(Cpp_fNdistinctm,x,attr(g,"N.groups"),g,NULL,na.rm,FALSE),g,TtI(TRA)))
  }
  if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRAm,x,.Call(Cpp_fNdistinctm,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TtI(TRA))
}

fNdistinct.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fNdistinctl,x,0L,0L,NULL,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(Cpp_fNdistinctl,x,length(lev),g,NULL,na.rm,FALSE), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fNdistinctl,x,fnlevels(g),g,NULL,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fNdistinctl,x,attr(g,"N.groups"),g,NULL,na.rm,FALSE))
    }
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(Cpp_fNdistinctl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE), groups))
    return(.Call(Cpp_fNdistinctl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE))
  }
  if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fNdistinctl,x,0L,0L,NULL,na.rm,TRUE),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fNdistinctl,x,fnlevels(g),g,NULL,na.rm,FALSE),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRAl,x,.Call(Cpp_fNdistinctl,x,attr(g,"N.groups"),g,NULL,na.rm,FALSE),g,TtI(TRA)))
  }
  if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRAl,x,.Call(Cpp_fNdistinctl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TtI(TRA))
}

fNdistinct.list <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)
  fNdistinct.data.frame(x, g, TRA, na.rm, use.g.names, drop, ...)

fNdistinct.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  nam <- attr(x, "names")
  gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  gl <- length(gn) > 0L
  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      ax[["groups"]] <- NULL
      ax[["class"]] <- fsetdiff(ax[["class"]], c("GRP_df", "grouped_df"))
      ax[["row.names"]] <- if(use.g.names) GRPnames(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], nam[-gn])
          return(setAttributes(c(g[[4L]],.Call(Cpp_fNdistinctl,x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE)), ax))
        }
        ax[["names"]] <- nam[-gn]
        return(setAttributes(.Call(Cpp_fNdistinctl,x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(Cpp_fNdistinctl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE)), ax))
      } else return(setAttributes(.Call(Cpp_fNdistinctl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fNdistinctl,x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TtI(TRA))), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fNdistinctl,x[-gn],g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TtI(TRA)), ax))
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fNdistinctl,x,g[[1L]],g[[2L]],g[[3L]],na.rm,FALSE),g[[2L]],TtI(TRA)))
}

