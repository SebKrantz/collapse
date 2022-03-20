
# Note: for foundational changes to this code see fsum.R

ffirst <- function(x, ...) UseMethod("ffirst") # , x

ffirst.default <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(ffirst.matrix(x, g, TRA, na.rm, use.g.names, ...))
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_ffirst,x,0L,0L,NULL,na.rm))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`names<-`(.Call(C_ffirst,x,length(lev),g,NULL,na.rm), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_ffirst,x,fnlevels(g),g,NULL,na.rm))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_ffirst,x,attr(g,"N.groups"),g,NULL,na.rm))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`names<-`(.Call(C_ffirst,x,g[[1L]],g[[2L]],g[[8L]],na.rm), GRPnames(g)))
    return(.Call(C_ffirst,x,g[[1L]],g[[2L]],g[[8L]],na.rm))
  }
  if(is.null(g)) return(TRAC(x,.Call(C_ffirst,x,0L,0L,NULL,na.rm),0L,TRA, ...))
  g <- G_guo(g)
  TRAC(x,.Call(C_ffirst,x,g[[1L]],g[[2L]],g$group.starts,na.rm),g[[2L]],TRA, ...)
}

ffirst.matrix <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) return(.Call(C_ffirstm,x,0L,0L,NULL,na.rm,drop))
    if(is.atomic(g)) {
      if(use.g.names) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(`dimnames<-`(.Call(C_ffirstm,x,length(lev),g,NULL,na.rm,FALSE), list(lev, dimnames(x)[[2L]])))
      }
      if(is.nmfactor(g)) return(.Call(C_ffirstm,x,fnlevels(g),g,NULL,na.rm,FALSE))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_ffirstm,x,attr(g,"N.groups"),g,NULL,na.rm,FALSE))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names) return(`dimnames<-`(.Call(C_ffirstm,x,g[[1L]],g[[2L]],g[[8L]],na.rm,FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
    return(.Call(C_ffirstm,x,g[[1L]],g[[2L]],g[[8L]],na.rm,FALSE))
  }
  if(is.null(g)) return(TRAmC(x,.Call(C_ffirstm,x,0L,0L,NULL,na.rm,TRUE),0L,TRA, ...))
  g <- G_guo(g)
  TRAmC(x,.Call(C_ffirstm,x,g[[1L]],g[[2L]],g$group.starts,na.rm,FALSE),g[[2L]],TRA, ...)
}

ffirst.data.frame <- function(x, g = NULL, TRA = NULL, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(is.null(TRA)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(is.null(g)) if(drop) return(unlist(.Call(C_ffirstl,x,0L,0L,NULL,na.rm))) else return(.Call(C_ffirstl,x,0L,0L,NULL,na.rm))
    if(is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
        lev <- attr(g, "levels")
        return(setRnDF(.Call(C_ffirstl,x,length(lev),g,NULL,na.rm), lev))
      }
      if(is.nmfactor(g)) return(.Call(C_ffirstl,x,fnlevels(g),g,NULL,na.rm))
      g <- qG(g, na.exclude = FALSE)
      return(.Call(C_ffirstl,x,attr(g,"N.groups"),g,NULL,na.rm))
    }
    if(!is_GRP(g)) g <- GRP.default(g, return.groups = use.g.names, call = FALSE)
    if(use.g.names && !inherits(x, "data.table") && length(groups <- GRPnames(g)))
      return(setRnDF(.Call(C_ffirstl,x,g[[1L]],g[[2L]],g[[8L]],na.rm), groups))
    return(.Call(C_ffirstl,x,g[[1L]],g[[2L]],g[[8L]],na.rm))
  }
  if(is.null(g)) return(TRAlC(x,.Call(C_ffirstl,x,0L,0L,NULL,na.rm),0L,TRA, ...))
  g <- G_guo(g)
  TRAlC(x,.Call(C_ffirstl,x,g[[1L]],g[[2L]],g$group.starts,na.rm),g[[2L]],TRA, ...)
}

ffirst.list <- function(x, ...) ffirst.data.frame(x, ...)

ffirst.grouped_df <- function(x, TRA = NULL, na.rm = TRUE, use.g.names = FALSE, keep.group_vars = TRUE, ...) {
  g <- GRP.grouped_df(x, call = FALSE)
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
          return(setAttributes(c(g[[4L]],.Call(C_ffirstl,x[-gn],g[[1L]],g[[2L]],g[[8L]],na.rm)), ax))
        }
        ax[["names"]] <- nam[-gn]
        return(setAttributes(.Call(C_ffirstl,x[-gn],g[[1L]],g[[2L]],g[[8L]],na.rm), ax))
      } else if(keep.group_vars) {
        ax[["names"]] <- c(g[[5L]], nam)
        return(setAttributes(c(g[[4L]],.Call(C_ffirstl,x,g[[1L]],g[[2L]],g[[8L]],na.rm)), ax))
      } else return(setAttributes(.Call(C_ffirstl,x,g[[1L]],g[[2L]],g[[8L]],na.rm), ax))
    } else if(keep.group_vars) {
      ax[["names"]] <- c(nam[gn], nam[-gn])
      return(setAttributes(c(x[gn],TRAlC(x[-gn],.Call(C_ffirstl,x[-gn],g[[1L]],g[[2L]],g[[8L]],na.rm),g[[2L]],TRA, ...)), ax))
    }
    ax[["names"]] <- nam[-gn]
    return(setAttributes(TRAlC(x[-gn],.Call(C_ffirstl,x[-gn],g[[1L]],g[[2L]],g[[8L]],na.rm),g[[2L]],TRA, ...), ax))
  } else return(TRAlC(x,.Call(C_ffirstl,x,g[[1L]],g[[2L]],g[[8L]],na.rm),g[[2L]],TRA, ...))
}
