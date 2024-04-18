
# Make faster ?
cm <- function(x) if(is.double(x)) x else if(is.character(x) && x == "overall.mean") -Inf else if(isFALSE(x)) Inf else stop("mean must be a number, 'overall.mean' or FALSE")
csd <- function(x) if(is.double(x)) x else if(is.character(x) && x == "within.sd") -Inf else stop("sd must be a number or 'within.sd'")

# TODO: w.type - Implement reliability weights?

fscale <- function(x, ...) UseMethod("fscale") # , x

fscale.default <- function(x, g = NULL, w = NULL, na.rm = .op[["na.rm"]], mean = 0, sd = 1, ...) {
  # if(is.matrix(x) && !inherits(x, "matrix")) return(fscale.matrix(x, g, w, na.rm, mean, sd, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_fscale,x,0L,0L,w,na.rm,cm(mean),csd(sd)))
  g <- G_guo(g)
  .Call(Cpp_fscale,x,g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd))
}

fscale.pseries <- function(x, effect = 1L, w = NULL, na.rm = .op[["na.rm"]], mean = 0, sd = 1, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- group_effect(x, effect)
  res <- if(is.matrix(x))
  .Call(Cpp_fscalem,x,fnlevels(g),g,w,na.rm,cm(mean),csd(sd)) else
  .Call(Cpp_fscale,x,fnlevels(g),g,w,na.rm,cm(mean),csd(sd))
  if(is.double(x)) return(res)
  pseries_to_numeric(res)
}

fscale.matrix <- function(x, g = NULL, w = NULL, na.rm = .op[["na.rm"]], mean = 0, sd = 1, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_fscalem,x,0L,0L,w,na.rm,cm(mean),csd(sd)))
  g <- G_guo(g)
  .Call(Cpp_fscalem,x,g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd))
}

fscale.zoo <- function(x, ...) if(is.matrix(x)) fscale.matrix(x, ...) else fscale.default(x, ...)
fscale.units <- fscale.zoo

fscale.grouped_df <- function(x, w = NULL, na.rm = .op[["na.rm"]], mean = 0, sd = 1, keep.group_vars = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  wsym <- substitute(w)
  nam <- attr(x, "names")
  gn2 <- which(nam %in% g[[5L]])
  gn <- if(keep.group_vars) gn2 else NULL

  if(!is.null(wsym)) {
    w <- eval(wsym, x, parent.frame())
    if(length(wn <- which(nam %in% all.vars(wsym)))) {
      if(any(gn2 %in% wn)) stop("Weights coincide with grouping variables!")
      gn2 <- c(gn2, wn)
      if(keep.w) gn <- c(gn, wn)
    }
  }

  if(length(gn2)) {
    # if(!length(gn)) return(.Call(Cpp_fscalel,x[-gn2],g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd)))
    ax <- attributes(x)
    ax[["names"]] <- c(nam[gn], nam[-gn2]) # first term is removed if !length(gn)
    res <- .Call(Cpp_fscalel, .subset(x, -gn2), g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd))
    if(length(gn)) return(setAttributes(c(.subset(x, gn), res), ax)) else return(setAttributes(res, ax))
  }
  .Call(Cpp_fscalel,x,g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd))
}

fscale.data.frame <- function(x, g = NULL, w = NULL, na.rm = .op[["na.rm"]], mean = 0, sd = 1, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_fscalel,x,0L,0L,w,na.rm,cm(mean),csd(sd)))
  g <- G_guo(g)
  .Call(Cpp_fscalel,x,g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd))
}

fscale.list <- function(x, ...) fscale.data.frame(x, ...)

fscale.pdata.frame <- function(x, effect = 1L, w = NULL, na.rm = .op[["na.rm"]], mean = 0, sd = 1, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- group_effect(x, effect)
  .Call(Cpp_fscalel,x,fnlevels(g),g,w,na.rm,cm(mean),csd(sd))
}


# Standardization Operator

STD <- function(x, ...) UseMethod("STD") # , x

STD.default <- function(x, g = NULL, w = NULL, na.rm = .op[["na.rm"]], mean = 0, sd = 1, ...) {
  # if(is.matrix(x) && !inherits(x, "matrix")) return(STD.matrix(x, g, w, na.rm, mean, sd, ...))
  fscale.default(x, g, w, na.rm, mean, sd, ...)
}

STD.pseries <- function(x, effect = 1L, w = NULL, na.rm = .op[["na.rm"]], mean = 0, sd = 1, ...)
  fscale.pseries(x, effect, w, na.rm, mean, sd, ...)

STD.matrix <- function(x, g = NULL, w = NULL, na.rm = .op[["na.rm"]], mean = 0, sd = 1, stub = .op[["stub"]], ...) {
  res <- fscale.matrix(x, g, w, na.rm, mean, sd, ...)
  if(isTRUE(stub) || is.character(stub)) return(add_stub(res, if(is.character(stub)) stub else "STD."))
  res
}

STD.zoo <- function(x, ...) if(is.matrix(x)) STD.matrix(x, ...) else STD.default(x, ...)
STD.units <- STD.zoo

STD.grouped_df <- function(x, w = NULL, na.rm = .op[["na.rm"]], mean = 0, sd = 1, stub = .op[["stub"]], keep.group_vars = TRUE, keep.w = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- GRP.grouped_df(x, call = FALSE)
  wsym <- substitute(w)
  nam <- attr(x, "names")
  gn2 <- which(nam %in% g[[5L]])
  gn <- if(keep.group_vars) gn2 else NULL

  if(!is.null(wsym)) {
    w <- eval(wsym, x, parent.frame())
    if(length(wn <- which(nam %in% all.vars(wsym)))) {
      if(any(gn2 %in% wn)) stop("Weights coincide with grouping variables!")
      gn2 <- c(gn2, wn)
      if(keep.w) gn <- c(gn, wn)
    }
  }

  if(length(gn2)) {
    ax <- attributes(x)
    ax[["names"]] <- c(nam[gn], do_stub(stub, nam[-gn2], "STD."))
    res <- .Call(Cpp_fscalel, .subset(x, -gn2), g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd))
    if(length(gn)) return(setAttributes(c(.subset(x, gn), res), ax)) else return(setAttributes(res, ax))
  }
  res <- .Call(Cpp_fscalel,x,g[[1L]],g[[2L]],w,na.rm,cm(mean),csd(sd))
  if(isTRUE(stub) || is.character(stub)) return(add_stub(res, if(is.character(stub)) stub else "STD."))
  res
}

# updated (best) version !
STD.pdata.frame <- function(x, effect = 1L, w = NULL, cols = is.numeric,
                            na.rm = .op[["na.rm"]], mean = 0, sd = 1, stub = .op[["stub"]], keep.ids = TRUE,
                            keep.w = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  ax <- attributes(x)
  nam <- ax[["names"]]
  g <- group_effect(x, effect)
  cols_fun <- is.function(cols)

  if(cols_fun && identical(cols, is.numeric)) cols <- which(.Call(C_vtypes, x, 1L))
  else if(length(cols)) cols <- cols2int(cols, x, nam)
  oldClass(x) <- NULL
  if(cols_fun || keep.ids) {
    gn <- which(nam %in% attr(findex(x), "nam")) # Needed for 3+ index variables
    if(length(gn)) {
      if(cols_fun) cols <- fsetdiff(cols, gn)
      else if(is.null(cols)) cols <- seq_along(x)[-gn]
    }
    if(!keep.ids) gn <- NULL
  } else gn <- NULL

  if(is.call(w)) {
    wn <- ckmatch(all.vars(w), nam, "Unknown weight variable:")
    w <- eval(w[[2L]], x, attr(w, ".Environment")) # w <- x[[wn]]
    cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
    if(keep.w) gn <- c(gn, wn)
  }

  if(length(gn) && length(cols)) {
    ax[["names"]] <- c(nam[gn], do_stub(stub, nam[cols], "STD."))
    return(setAttributes(c(x[gn], .Call(Cpp_fscalel,x[cols],fnlevels(g),g,w,na.rm,cm(mean),csd(sd))), ax))
  }
  if(!length(gn)) {
    ax[["names"]] <- do_stub(stub, nam[cols], "STD.")
    return(setAttributes(.Call(Cpp_fscalel,x[cols],fnlevels(g),g,w,na.rm,cm(mean),csd(sd)), ax))
  }
  if(isTRUE(stub) || is.character(stub)) {
    ax[["names"]] <- do_stub(stub, nam, "STD.")
    return(setAttributes(.Call(Cpp_fscalel,x,fnlevels(g),g,w,na.rm,cm(mean),csd(sd)), ax))
  }
  .Call(Cpp_fscalel,`oldClass<-`(x, ax[["class"]]),fnlevels(g),g,w,na.rm,cm(mean),csd(sd))
}

# updated, fast and data.table proof version !
STD.data.frame <- function(x, by = NULL, w = NULL, cols = is.numeric,
                           na.rm = .op[["na.rm"]], mean = 0, sd = 1, stub = .op[["stub"]], keep.by = TRUE,
                           keep.w = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)

  if(is.call(by) || is.call(w)) {
    ax <- attributes(x)
    class(x) <- NULL
    nam <- names(x)

    if(is.call(by)) {
      if(length(by) == 3L) {
        cols <- ckmatch(all.vars(by[[2L]]), nam)
        gn <- ckmatch(all.vars(by[[3L]]), nam)
      } else {
        gn <- ckmatch(all.vars(by), nam)
        cols <- cols2intrmgn(gn, cols, x)
      }
      by <- G_guo(if(length(gn) == 1L) x[[gn]] else x[gn])
      if(!keep.by) gn <- NULL
    } else {
      gn <- NULL
      if(length(cols)) cols <- cols2int(cols, x, nam)
      by <- if(is.null(by)) list(0L, 0L) else G_guo(by)
    }

    if(is.call(w)) {
      wn <- ckmatch(all.vars(w), nam, "Unknown weight variable:")
      w <- eval(w[[2L]], x, attr(w, ".Environment")) # w <- x[[wn]]
      cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
      if(keep.w) gn <- c(gn, wn)
    }

    if(length(gn)) {
      ax[["names"]] <- c(nam[gn], do_stub(stub, nam[cols], "STD."))
      return(setAttributes(c(x[gn], .Call(Cpp_fscalel,x[cols],by[[1L]],by[[2L]],w,na.rm,cm(mean),csd(sd))), ax))
    }
    ax[["names"]] <- do_stub(stub, nam[cols], "STD.")
    return(setAttributes(.Call(Cpp_fscalel,x[cols],by[[1L]],by[[2L]],w,na.rm,cm(mean),csd(sd)), ax))
  } else if(length(cols)) { # Needs to be like this, otherwise subsetting dropps the attributes !!
    ax <- attributes(x)
    class(x) <- NULL
    x <- x[cols2int(cols, x, names(x), FALSE)]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }
  if(isTRUE(stub) || is.character(stub)) attr(x, "names") <- do_stub(stub, attr(x, "names"), "STD.")

  if(is.null(by)) return(.Call(Cpp_fscalel,x,0L,0L,w,na.rm,cm(mean),csd(sd)))
  by <- G_guo(by)
  .Call(Cpp_fscalel,x,by[[1L]],by[[2L]],w,na.rm,cm(mean),csd(sd))
}

STD.list <- function(x, ...) STD.data.frame(x, ...)

