
# Note: for foundational changes to this code see fsum.R

flast <- function(x, ...) UseMethod("flast") # , x

flast.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(flast.matrix(x, g, TRA, na.rm, use.g.names, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(C_flast,x,0L,0L,na.rm))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(C_flast,x,length(lev),g,na.rm), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_flast,x,fnlevels(g),g,na.rm))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_flast,x,attr(g,"N.groups"),g,na.rm))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(C_flast,x,g[[1L]],g[[2L]],na.rm), GRPnames(g)))
    return(.Call(C_flast,x,g[[1L]],g[[2L]],na.rm))
  }
  if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(C_flast,x,0L,0L,na.rm),0L,TtI(TRA)))
  g <- G_guo(g)
  .Call(Cpp_TRA,x,.Call(C_flast,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TtI(TRA))
}

flast.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(C_flastm,x,0L,0L,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(C_flastm,x,length(lev),g,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(C_flastm,x,fnlevels(g),g,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_flastm,x,attr(g,"N.groups"),g,na.rm,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(C_flastm,x,g[[1L]],g[[2L]],na.rm,FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(C_flastm,x,g[[1L]],g[[2L]],na.rm,FALSE))
  }
  if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(C_flastm,x,0L,0L,na.rm,TRUE),0L,TtI(TRA)))
  g <- G_guo(g)
  .Call(Cpp_TRAm,x,.Call(C_flastm,x,g[[1L]],g[[2L]],na.rm,FALSE),g[[2L]],TtI(TRA))
}

flast.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) if(drop) return(unlist(.Call(C_flastl,x,0L,0L,na.rm))) else return(.Call(C_flastl,x,0L,0L,na.rm))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(C_flastl,x,length(lev),g,na.rm), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_flastl,x,fnlevels(g),g,na.rm))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_flastl,x,attr(g,"N.groups"),g,na.rm))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(C_flastl,x,g[[1L]],g[[2L]],na.rm), groups))
    return(.Call(C_flastl,x,g[[1L]],g[[2L]],na.rm))
  }
  if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(C_flastl,x,0L,0L,na.rm),0L,TtI(TRA)))
  g <- G_guo(g)
  .Call(Cpp_TRAl,x,.Call(C_flastl,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TtI(TRA))
}

flast.list <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)
  flast.data.frame(x, g, TRA, na.rm, use.g.names, drop, ...)

flast.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
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
          return(setAttributes(c(g[[4L]],.Call(C_flastl,x[-gn],g[[1L]],g[[2L]],na.rm)), ax))
        }
        ax[["names"]] <- nam[-gn]
        return(setAttributes(.Call(C_flastl,x[-gn],g[[1L]],g[[2L]],na.rm), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(C_flastl,x,g[[1L]],g[[2L]],na.rm)), ax))
      } else return(setAttributes(.Call(C_flastl,x,g[[1L]],g[[2L]],na.rm), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(C_flastl,x[-gn],g[[1L]],g[[2L]],na.rm),g[[2L]],TtI(TRA))), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(C_flastl,x[-gn],g[[1L]],g[[2L]],na.rm),g[[2L]],TtI(TRA)), ax))
  } else return(.Call(Cpp_TRAl,x,.Call(C_flastl,x,g[[1L]],g[[2L]],na.rm),g[[2L]],TtI(TRA)))
}
