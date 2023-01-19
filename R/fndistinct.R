
fndistinct <- function(x, ...) UseMethod("fndistinct") # , x

fndistinct.default <- function(x, g = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = TRUE, nthreads = .op[["nthreads"]], ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fndistinct.matrix(x, g, TRA, na.rm, use.g.names, nthreads = nthreads, ...))
  if(!is.null(g)) g <- GRP(g, return.groups = use.g.names && is.null(TRA), call = FALSE) # sort = FALSE for TRA: not faster here...
  res <- .Call(C_fndistinct,x,g,na.rm,nthreads)
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(res)
    if(use.g.names) names(res) <- GRPnames(g, FALSE)
    return(res)
  }
  TRAC(x,res,g[[2L]],TRA, ...)
}

fndistinct.matrix <- function(x, g = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = TRUE, drop = TRUE, nthreads = .op[["nthreads"]], ...) {
  if(!is.null(g)) g <- GRP(g, return.groups = use.g.names && is.null(TRA), call = FALSE) # sort = FALSE for TRA: not faster here...
  res <- .Call(C_fndistinctm,x,g,na.rm,drop,nthreads)
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(res)
    if(use.g.names) dimnames(res)[[1L]] <- GRPnames(g)
    return(res)
  }
  TRAmC(x,res,g[[2L]],TRA, ...)
}

fndistinct.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = TRUE, drop = TRUE, nthreads = .op[["nthreads"]], ...) {
  if(!is.null(g)) g <- GRP(g, return.groups = use.g.names && is.null(TRA), call = FALSE) # sort = FALSE for TRA: not faster here...
  res <- .Call(C_fndistinctl,x,g,na.rm,drop,nthreads)
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(res)
    if(use.g.names && !inherits(x, "data.table") && length(gn <- GRPnames(g)))
      attr(res, "row.names") <- gn
    return(res)
  }
  TRAlC(x,res,g[[2L]],TRA, ...)
}

fndistinct.list <- function(x, ...) fndistinct.data.frame(x, ...)

fndistinct.grouped_df <- function(x, TRA = NULL, na.rm = .op[["na.rm"]], use.g.names = FALSE, keep.group_vars = TRUE, nthreads = .op[["nthreads"]], ...) {
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
          return(setAttributes(c(g[[4L]],.Call(C_fndistinctl,x[-gn],g,na.rm,FALSE,nthreads)), ax))
        }
        ax[["names"]] <- nam[-gn]
        return(setAttributes(.Call(C_fndistinctl,x[-gn],g,na.rm,FALSE,nthreads), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(C_fndistinctl,x,g,na.rm,FALSE,nthreads)), ax))
      } else return(setAttributes(.Call(C_fndistinctl,x,g,na.rm,FALSE,nthreads), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],TRAlC(x[-gn],.Call(C_fndistinctl,x[-gn],g,na.rm,FALSE,nthreads),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(C_fndistinctl,x[-gn],g,na.rm,FALSE,nthreads),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(C_fndistinctl,x,g,na.rm,FALSE,nthreads),g[[2L]],TRA, ...))
}


fNdistinct <- function(x, ...) {
  message("Note that 'fNdistinct' was renamed to 'fndistinct'. The S3 generic will not be removed anytime soon, but please use updated function names in new code, see help('collapse-renamed')")
  UseMethod("fndistinct")
}
fNdistinct.default <- function(x, ...) {
  # .Deprecated(msg = "This method belongs to a renamed function and will be removed end of 2022, see help('collapse-renamed')")
  fndistinct.default(x, ...)
}
fNdistinct.matrix <- function(x, ...) {
  # .Deprecated(msg = "This method belongs to a renamed function and will be removed end of 2022, see help('collapse-renamed')")
  fndistinct.matrix(x, ...)
}
fNdistinct.data.frame <- function(x, ...) {
  # .Deprecated(msg = "This method belongs to a renamed function and will be removed end of 2022, see help('collapse-renamed')")
  fndistinct.data.frame(x, ...)
}
