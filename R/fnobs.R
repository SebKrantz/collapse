
# For foundational changes to this code see fsum.R

fnobs <- function(x, ...) UseMethod("fnobs") # , x

fnobs.default <- function(x, g = NULL, TRA = NULL, use.g.names = TRUE, ...) {
  # if(is.matrix(x) && !inherits(x, "matrix")) return(fnobs.matrix(x, g, TRA, use.g.names, ...))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fnobs,x,0L,0L))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(C_fnobs,x,length(lev),g), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fnobs,x,fnlevels(g),g))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fnobs,x,attr(g,"N.groups"),g))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(C_fnobs,x,g[[1L]],g[[2L]]), GRPnames(g)))
    return(.Call(C_fnobs,x,g[[1L]],g[[2L]]))
  }
  if(is.null(g)) return(TRAC(x,.Call(C_fnobs,x,0L,0L),0L,TRA, ...))
  g <- G_guo(g)
  TRAC(x,.Call(C_fnobs,x,g[[1L]],g[[2L]]),g[[2L]],TRA, ...)
}

fnobs.matrix <- function(x, g = NULL, TRA = NULL, use.g.names = TRUE, drop = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fnobsm,x,0L,0L,drop))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(C_fnobsm,x,length(lev),g,FALSE), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(C_fnobsm,x,fnlevels(g),g,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fnobsm,x,attr(g,"N.groups"),g,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(C_fnobsm,x,g[[1L]],g[[2L]],FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(C_fnobsm,x,g[[1L]],g[[2L]],FALSE))
  }
  if(is.null(g)) return(TRAmC(x,.Call(C_fnobsm,x,0L,0L,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAmC(x,.Call(C_fnobsm,x,g[[1L]],g[[2L]],FALSE),g[[2L]],TRA, ...)
}

fnobs.zoo <- function(x, ...) if(is.matrix(x)) fnobs.matrix(x, ...) else fnobs.default(x, ...)
fnobs.units <- fnobs.zoo

fnobs.data.frame <- function(x, g = NULL, TRA = NULL, use.g.names = TRUE, drop = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_fnobsl,x,0L,0L,drop))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(C_fnobsl,x,length(lev),g,FALSE), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_fnobsl,x,fnlevels(g),g,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_fnobsl,x,attr(g,"N.groups"),g,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(C_fnobsl,x,g[[1L]],g[[2L]],FALSE), groups))
    return(.Call(C_fnobsl,x,g[[1L]],g[[2L]],FALSE))
  }
  if(is.null(g)) return(TRAlC(x,.Call(C_fnobsl,x,0L,0L,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAlC(x,.Call(C_fnobsl,x,g[[1L]],g[[2L]],FALSE),g[[2L]],TRA, ...)
}

fnobs.list <- function(x, ...) fnobs.data.frame(x, ...)

fnobs.grouped_df <- function(x, TRA = NULL, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
  g <- GRP.grouped_df(x, call = FALSE)
  if(is.null(g[[4L]])) keep.group_vars <- FALSE
  nam <- attr(x, "names")
  gn <- which(nam %in% g[[5L]])
  nTRAl <- is.null(TRA)
  gl <- length(gn) > 0L
  if(gl || nTRAl) {
    ax <- attributes(x)
    attributes(x) <- NULL
    if(nTRAl) {
      if(!missing(...)) unused_arg_action(match.call(), ...)
      ax[["groups"]] <- NULL
      ax[["class"]] <- fsetdiff(ax[["class"]], c("GRP_df", "grouped_df"))
      ax[["row.names"]] <- if(use.g.names) GRPnames(g) else .set_row_names(g[[1L]])
      if(gl) {
        if(keep.group_vars) {
          ax[["names"]] <- c(g[[5L]], nam[-gn])
          return(setAttributes(c(g[[4L]],.Call(C_fnobsl,x[-gn],g[[1L]],g[[2L]],FALSE)), ax))
        }
        ax[["names"]] <- nam[-gn]
        return(setAttributes(.Call(C_fnobsl,x[-gn],g[[1L]],g[[2L]],FALSE), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(C_fnobsl,x,g[[1L]],g[[2L]],FALSE)), ax))
      } else return(setAttributes(.Call(C_fnobsl,x,g[[1L]],g[[2L]],FALSE), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],TRAlC(x[-gn],.Call(C_fnobsl,x[-gn],g[[1L]],g[[2L]],FALSE),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(C_fnobsl,x[-gn],g[[1L]],g[[2L]],FALSE),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(C_fnobsl,x,g[[1L]],g[[2L]],FALSE),g[[2L]],TRA, ...))
}

fNobs <- function(x, ...) {
  message("Note that 'fNobs' was renamed to 'fnobs'. The S3 generic will not be removed anytime soon, but please use updated function names in new code, see help('collapse-renamed')")
  UseMethod("fnobs")
}
fNobs.default <- function(x, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fnobs.matrix(x, ...))
  # .Deprecated(msg = "This method belongs to a renamed function and will be removed end of 2022, see help('collapse-renamed')")
  fnobs.default(x, ...)
}
fNobs.matrix <- function(x, ...) {
  # .Deprecated(msg = "This method belongs to a renamed function and will be removed end of 2022, see help('collapse-renamed')")
  fnobs.matrix(x, ...)
}
fNobs.data.frame <- function(x, ...) {
  # .Deprecated(msg = "This method belongs to a renamed function and will be removed end of 2022, see help('collapse-renamed')")
  fnobs.data.frame(x, ...)
}
