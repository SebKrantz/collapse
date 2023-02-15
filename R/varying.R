
varying <- function(x, ...) UseMethod("varying") # , x

varying.default <- function(x, g = NULL, any_group = TRUE, use.g.names = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(varying.matrix(x, g, any_group, use.g.names, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_varying,x,0L,0L,any_group))
  if(is.atomic(g)) {
    if(use.g.names && !any_group) {
      if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE, sort = !any_group && .op[["sort"]])
      lev <- attr(g, "levels")
      return(`names<-`(.Call(Cpp_varying,x,length(lev),g,any_group), lev))
    }
    if(is.nmfactor(g)) return(.Call(Cpp_varying,x,fnlevels(g),g,any_group))
    g <- qG(g, na.exclude = FALSE, sort = !any_group && .op[["sort"]])
    return(.Call(Cpp_varying,x,attr(g,"N.groups"),g,any_group))
  }
  if(!is_GRP(g)) g <- GRP.default(g, sort = !any_group && .op[["sort"]], return.groups = use.g.names && !any_group, call = FALSE)
  if(use.g.names && !any_group) return(`names<-`(.Call(Cpp_varying,x,g[[1L]],g[[2L]],any_group), GRPnames(g)))
  .Call(Cpp_varying,x,g[[1L]],g[[2L]],any_group)
}

varying.pseries <- function(x, effect = 1L, any_group = TRUE, use.g.names = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- if(length(effect) == 1L) .subset2(findex(x), effect) else finteraction(.subset(findex(x), effect), sort = !any_group && .op[["sort"]])
  if(!any_group && use.g.names) {
    lev <- attr(g, "levels")
    return(`names<-`(.Call(Cpp_varying,x,length(lev),g,any_group), lev))
  }
  .Call(Cpp_varying,x,fnlevels(g),g,any_group)
}

varying.matrix <- function(x, g = NULL, any_group = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_varyingm,x,0L,0L,any_group,drop))
  if(is.atomic(g)) {
    if(use.g.names && !any_group) {
      if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE, sort = !any_group && .op[["sort"]])
      lev <- attr(g, "levels")
      return(`dimnames<-`(.Call(Cpp_varyingm,x,length(lev),g,any_group,FALSE), list(lev, dimnames(x)[[2L]])))
    }
    if(is.nmfactor(g)) return(.Call(Cpp_varyingm,x,fnlevels(g),g,any_group,drop))
    g <- qG(g, na.exclude = FALSE, sort = !any_group && .op[["sort"]])
    return(.Call(Cpp_varyingm,x,attr(g,"N.groups"),g,any_group,drop))
  }
  if(!is_GRP(g)) g <- GRP.default(g, sort = !any_group && .op[["sort"]], return.groups = use.g.names && !any_group, call = FALSE)
  if(use.g.names && !any_group) return(`dimnames<-`(.Call(Cpp_varyingm,x,g[[1L]],g[[2L]],any_group,FALSE), list(GRPnames(g), dimnames(x)[[2L]])))
  .Call(Cpp_varyingm,x,g[[1L]],g[[2L]],any_group,drop)
}

varying.data.frame <- function(x, by = NULL, cols = NULL, any_group = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)

  if(is.call(by)) {
    nam <- attr(x, "names")
    if(length(by) == 3L) {
      cols <- ckmatch(all.vars(by[[2L]]), nam)
      gn <- ckmatch(all.vars(by[[3L]]), nam)
    } else {
      gn <- ckmatch(all.vars(by), nam)
      cols <- if(is.null(cols)) seq_along(unclass(x))[-gn] else cols2int(cols, x, nam, FALSE)
    }
    by <- if(length(gn) == 1L) .subset2(x, gn) else GRP.default(x, gn, sort = !any_group && .op[["sort"]], return.groups = use.g.names && !any_group, call = FALSE)
    x <- fcolsubset(x, cols)
  } else if(length(cols)) x <- colsubset(x, cols)

  if(is.null(by)) return(.Call(Cpp_varyingl,x,0L,0L,any_group,drop))
  if(is.atomic(by)) {
    if(use.g.names && !any_group && !inherits(x, "data.table")) {
      if(!is.nmfactor(by)) by <- qF(by, na.exclude = FALSE, sort = !any_group && .op[["sort"]])
      lev <- attr(by, "levels")
      return(setRnDF(.Call(Cpp_varyingl,x,length(lev),by,any_group,FALSE), lev))
    }
    if(is.nmfactor(by)) return(.Call(Cpp_varyingl,x,fnlevels(by),by,any_group,drop))
    by <- qG(by, na.exclude = FALSE, sort = !any_group && .op[["sort"]])
    return(.Call(Cpp_varyingl,x,attr(by,"N.groups"),by,any_group,drop))
  }
  if(!is_GRP(by)) by <- GRP.default(by, sort = !any_group && .op[["sort"]], return.groups = use.g.names && !any_group, call = FALSE)
  if(use.g.names && !any_group && !inherits(x, "data.table") && length(groups <- GRPnames(by)))
    return(setRnDF(.Call(Cpp_varyingl,x,by[[1L]],by[[2L]],any_group,FALSE), groups))
  .Call(Cpp_varyingl,x,by[[1L]],by[[2L]],any_group,drop)
}

varying.list <- function(x, ...) varying.data.frame(x, ...)

varying.pdata.frame <- function(x, effect = 1L, cols = NULL, any_group = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  index <- unclass(findex(x))
  g <- if(length(effect) == 1L) index[[effect]] else finteraction(index[effect], sort = !any_group && .op[["sort"]])
  x <- if(is.null(cols)) fcolsubset(x, attr(x, "names") %!in% names(index[effect])) else colsubset(x, cols)
  if(!any_group && use.g.names) {
    lev <- attr(g, "levels")
    return(setRnDF(.Call(Cpp_varyingl,x,length(lev),g,any_group,FALSE), lev))
  }
  .Call(Cpp_varyingl,x,fnlevels(g),g,any_group,drop)
}

varying.grouped_df <- function(x, any_group = TRUE, use.g.names = FALSE, drop = TRUE, keep.group_vars = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  nam <- attr(x, "names")
  ngn <- nam %!in% g[[5L]]
  if(any_group) {
    if(!all(ngn)) x <- if(drop) .subset(x, ngn) else fcolsubset(x, ngn)
    return(.Call(Cpp_varyingl,x,g[[1L]],g[[2L]],any_group,drop))
  }
  if(is.null(g[[4L]])) keep.group_vars <- FALSE
  ax <- attributes(x)
  ax[["groups"]] <- NULL
  ax[["class"]] <- fsetdiff(ax[["class"]], c("GRP_df", "grouped_df"))
  ax[["row.names"]] <- if(use.g.names) GRPnames(g) else .set_row_names(g[[1L]])
  if(!all(ngn)) {
    if(keep.group_vars) {
      ax[["names"]] <- c(g[[5L]], nam[ngn])
      return(setAttributes(c(g[[4L]],.Call(Cpp_varyingl,.subset(x, ngn),g[[1L]],g[[2L]],FALSE,FALSE)), ax))
    }
    ax[["names"]] <- nam[ngn]
    return(setAttributes(.Call(Cpp_varyingl,.subset(x, ngn),g[[1L]],g[[2L]],FALSE,FALSE), ax))
  } else if(keep.group_vars) {
    ax[["names"]] <- c(g[[5L]], nam)
    return(setAttributes(c(g[[4L]],.Call(Cpp_varyingl,x,g[[1L]],g[[2L]],FALSE,FALSE)), ax))
  } else return(setAttributes(.Call(Cpp_varyingl,x,g[[1L]],g[[2L]],FALSE,FALSE), ax))
}

varying.sf <- function(x, by = NULL, cols = NULL, any_group = TRUE, use.g.names = TRUE, drop = TRUE, ...) {
  clx <- oldClass(x)
  oldClass(x) <- NULL
  x[[attr(x, "sf_column")]] <- NULL
  oldClass(x) <- clx[clx != "sf"]
  if(any(clx == "grouped_df")) return(varying.grouped_df(x, any_group, use.g.names, drop, ...))
  varying.data.frame(x, by, cols, any_group, use.g.names, drop, ...)
}


