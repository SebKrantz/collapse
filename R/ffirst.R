
# Note: for foundational changes to this code see fsum.R

ffirst <- function(x, ...) UseMethod("ffirst") # , x

ffirst.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_ffirst,x,0L,0L,na.rm))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_ffirst,x,length(lev),g,na.rm), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_ffirst,x,fnlevels(g),g,na.rm))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_ffirst,x,attr(g,"N.groups"),g,na.rm))
    }
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(Cpp_ffirst,x,g[[1L]],g[[2L]],na.rm), GRPnames(g)))
    return(.Call(Cpp_ffirst,x,g[[1L]],g[[2L]],na.rm))
  }
  if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_ffirst,x,0L,0L,na.rm),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_ffirst,x,fnlevels(g),g,na.rm),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRA,x,.Call(Cpp_ffirst,x,attr(g,"N.groups"),g,na.rm),g,TtI(TRA)))
  }
  if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRA,x,.Call(Cpp_ffirst,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TtI(TRA))
}

ffirst.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_ffirstm,x,0L,0L,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_ffirstm,x,length(lev),g,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_ffirstm,x,fnlevels(g),g,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_ffirstm,x,attr(g,"N.groups"),g,na.rm,FALSE))
    }
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(Cpp_ffirstm,x,g[[1L]],g[[2L]],na.rm,FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(Cpp_ffirstm,x,g[[1L]],g[[2L]],na.rm,FALSE))
  }
  if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_ffirstm,x,0L,0L,na.rm,TRUE),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_ffirstm,x,fnlevels(g),g,na.rm,FALSE),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRAm,x,.Call(Cpp_ffirstm,x,attr(g,"N.groups"),g,na.rm,FALSE),g,TtI(TRA)))
  }
  if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRAm,x,.Call(Cpp_ffirstm,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TtI(TRA))
}

ffirst.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) if(drop) return(unlist(.Call(Cpp_ffirstl,x,0L,0L,na.rm))) else return(.Call(Cpp_ffirstl,x,0L,0L,na.rm))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(Cpp_ffirstl,x,length(lev),g,na.rm), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_ffirstl,x,fnlevels(g),g,na.rm))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_ffirstl,x,attr(g,"N.groups"),g,na.rm))
    }
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- GRPnames(g)))
      return(setRnDF(.Call(Cpp_ffirstl,x,g[[1L]],g[[2L]],na.rm), groups))
    return(.Call(Cpp_ffirstl,x,g[[1L]],g[[2L]],na.rm))
  }
  if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_ffirstl,x,0L,0L,na.rm),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_ffirstl,x,fnlevels(g),g,na.rm),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRAl,x,.Call(Cpp_ffirstl,x,attr(g,"N.groups"),g,na.rm),g,TtI(TRA)))
  }
  if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRAl,x,.Call(Cpp_ffirstl,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TtI(TRA))
}

ffirst.list <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)
  ffirst.data.frame(x, g, TRA, na.rm, use.g.names, drop, ...)

ffirst.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
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
      ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
      ax[["row.names"]] <- if(use.g.names) GRPnames(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], nam[-gn])
          return(setAttributes(c(g[[4L]],.Call(Cpp_ffirstl,x[-gn],g[[1L]],g[[2L]],na.rm)), ax))
        }
        ax[["names"]] <- nam[-gn]
        return(setAttributes(.Call(Cpp_ffirstl,x[-gn],g[[1L]],g[[2L]],na.rm), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(Cpp_ffirstl,x,g[[1L]],g[[2L]],na.rm)), ax))
      } else return(setAttributes(.Call(Cpp_ffirstl,x,g[[1L]],g[[2L]],na.rm), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_ffirstl,x[-gn],g[[1L]],g[[2L]],na.rm),g[[2L]],TtI(TRA))), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_ffirstl,x[-gn],g[[1L]],g[[2L]],na.rm),g[[2L]],TtI(TRA)), ax))
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_ffirstl,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TtI(TRA)))
}
