
# For foundational changes to this code see fsum.R

fNobs <- function(x, ...) UseMethod("fNobs") # , x

fNobs.default <- function(x, g = NULL, TRA = NULL, use.g.names = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fNobs.matrix(x, g, TRA, na.rm, use.g.names, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fNobs,x,0L,0L))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(Cpp_fNobs,x,length(lev),g), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fNobs,x,fnlevels(g),g))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fNobs,x,attr(g,"N.groups"),g))
    }
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(Cpp_fNobs,x,g[[1L]],g[[2L]]), GRPnames(g)))
    return(.Call(Cpp_fNobs,x,g[[1L]],g[[2L]]))
  }
  if(is.null(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fNobs,x,0L,0L),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRA,x,.Call(Cpp_fNobs,x,fnlevels(g),g),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRA,x,.Call(Cpp_fNobs,x,attr(g,"N.groups"),g),g,TtI(TRA)))
  }
  if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRA,x,.Call(Cpp_fNobs,x,g[[1L]],g[[2L]]),g[[2L]],TtI(TRA))
}

fNobs.matrix <- function(x, g = NULL, TRA = NULL, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fNobsm,x,0L,0L,drop))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(Cpp_fNobsm,x,length(lev),g,FALSE), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fNobsm,x,fnlevels(g),g,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fNobsm,x,attr(g,"N.groups"),g,FALSE))
    }
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(Cpp_fNobsm,x,g[[1L]],g[[2L]],FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(Cpp_fNobsm,x,g[[1L]],g[[2L]],FALSE))
  }
  if(is.null(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fNobsm,x,0L,0L,TRUE),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRAm,x,.Call(Cpp_fNobsm,x,fnlevels(g),g,FALSE),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRAm,x,.Call(Cpp_fNobsm,x,attr(g,"N.groups"),g,FALSE),g,TtI(TRA)))
  }
  if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRAm,x,.Call(Cpp_fNobsm,x,g[[1L]],g[[2L]],FALSE),g[[2L]],TtI(TRA))
}

fNobs.data.frame <- function(x, g = NULL, TRA = NULL, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(TRA)) {
    if(is.null(g)) return(.Call(Cpp_fNobsl,x,0L,0L,drop))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(Cpp_fNobsl,x,length(lev),g,FALSE), lev))
      }
      if(is.nmfactor(g)) return(.Call(Cpp_fNobsl,x,fnlevels(g),g,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(Cpp_fNobsl,x,attr(g,"N.groups"),g,FALSE))
    }
    if(!is.GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(Cpp_fNobsl,x,g[[1L]],g[[2L]],FALSE), groups))
    return(.Call(Cpp_fNobsl,x,g[[1L]],g[[2L]],FALSE))
  }
  if(is.null(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fNobsl,x,0L,0L,TRUE),0L,TtI(TRA)))
  if(is.atomic(g)) {
    if(is.nmfactor(g)) return(.Call(Cpp_TRAl,x,.Call(Cpp_fNobsl,x,fnlevels(g),g,FALSE),g,TtI(TRA)))
    g <- qG(g, na.exclude = FALSE)
    return(.Call(Cpp_TRAl,x,.Call(Cpp_fNobsl,x,attr(g,"N.groups"),g,FALSE),g,TtI(TRA)))
  }
  if(!is.GRP(g)) g <- GRP.default(g, return.groups = FALSE, call = FALSE)
  .Call(Cpp_TRAl,x,.Call(Cpp_fNobsl,x,g[[1L]],g[[2L]],FALSE),g[[2L]],TtI(TRA))
}

fNobs.list <- function(x, g = NULL, TRA = NULL, use.g.names = TRUE, drop = TRUE, ...)
  fNobs.data.frame(x, g, TRA, use.g.names, drop, ...)

fNobs.grouped_df <- function(x, TRA = NULL, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
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
          return(setAttributes(c(g[[4L]],.Call(Cpp_fNobsl,x[-gn],g[[1L]],g[[2L]],FALSE)), ax))
        }
        ax[["names"]] <- nam[-gn]
        return(setAttributes(.Call(Cpp_fNobsl,x[-gn],g[[1L]],g[[2L]],FALSE), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(Cpp_fNobsl,x,g[[1L]],g[[2L]],FALSE)), ax))
      } else return(setAttributes(.Call(Cpp_fNobsl,x,g[[1L]],g[[2L]],FALSE), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fNobsl,x[-gn],g[[1L]],g[[2L]],FALSE),g[[2L]],TtI(TRA))), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(.Call(Cpp_TRAl,x[-gn],.Call(Cpp_fNobsl,x[-gn],g[[1L]],g[[2L]],FALSE),g[[2L]],TtI(TRA)), ax))
  } else return(.Call(Cpp_TRAl,x,.Call(Cpp_fNobsl,x,g[[1L]],g[[2L]],FALSE),g[[2L]],TtI(TRA)))
}
