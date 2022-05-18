
ckm <- function(x) if(is.double(x)) x else if(is.character(x) && x == "overall.mean") -Inf else stop("mean must be a number or 'overall.mean'") # better than switch !!

# Note: for principal innovations of this code see fsum.R and fscale.R

fwithin <- function(x, ...) UseMethod("fwithin") # , x

fwithin.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, theta = 1, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fwithin.matrix(x, g, w, na.rm, mean, theta, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_BW,x,0L,0L,NULL,w,na.rm,theta,ckm(mean),FALSE,FALSE))
  g <- G_guo(g)
  .Call(Cpp_BW,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,theta,ckm(mean),FALSE,FALSE)
}

fwithin.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, theta = 1, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- group_effect(x, effect)
  if(is.matrix(x))
  .Call(Cpp_BWm,x,fnlevels(g),g,NULL,w,na.rm,theta,ckm(mean),FALSE,FALSE) else
  .Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,theta,ckm(mean),FALSE,FALSE)
}

fwithin.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, theta = 1, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_BWm,x,0L,0L,NULL,w,na.rm,theta,ckm(mean),FALSE,FALSE))
  g <- G_guo(g)
  .Call(Cpp_BWm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,theta,ckm(mean),FALSE,FALSE)
}

fwithin.data.frame <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, theta = 1, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_BWl,x,0L,0L,NULL,w,na.rm,theta,ckm(mean),FALSE,FALSE))
  g <- G_guo(g)
  .Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,theta,ckm(mean),FALSE,FALSE)
}

fwithin.list <- function(x, ...) fwithin.data.frame(x, ...)

fwithin.pdata.frame <- function(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, theta = 1, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- group_effect(x, effect)
  .Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,theta,ckm(mean),FALSE,FALSE)
}

fwithin.grouped_df <- function(x, w = NULL, na.rm = TRUE, mean = 0, theta = 1,
                               keep.group_vars = TRUE, keep.w = TRUE, ...) {
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
    ax[["names"]] <- c(nam[gn], nam[-gn2]) # first term is removed if !length(gn)
    res <- .Call(Cpp_BWl, .subset(x, -gn2), g[[1L]],g[[2L]],g[[3L]],w,na.rm,theta,ckm(mean),FALSE,FALSE)
    if(length(gn)) return(setAttributes(c(.subset(x, gn), res), ax)) else return(setAttributes(res, ax))
  }
  .Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,theta,ckm(mean),FALSE,FALSE)
}

# Within Operator

W <- function(x, ...) UseMethod("W") # , x

W.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, theta = 1, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(W.matrix(x, g, w, na.rm, mean, theta, ...))
  fwithin.default(x, g, w, na.rm, mean, theta, ...)
}

W.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, theta = 1, ...)
  fwithin.pseries(x, effect, w, na.rm, mean, theta, ...)

W.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, theta = 1, stub = "W.", ...)
  add_stub(fwithin.matrix(x, g, w, na.rm, mean, theta, ...), stub)

W.grouped_df <- function(x, w = NULL, na.rm = TRUE, mean = 0, theta = 1,
                         stub = "W.", keep.group_vars = TRUE, keep.w = TRUE, ...) {
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
    ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[-gn2]) else nam[-gn2])
    res <- .Call(Cpp_BWl, .subset(x, -gn2), g[[1L]],g[[2L]],g[[3L]],w,na.rm,theta,ckm(mean),FALSE,FALSE)
    if(length(gn)) return(setAttributes(c(.subset(x, gn), res), ax)) else return(setAttributes(res, ax))
  }
  add_stub(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,theta,ckm(mean),FALSE,FALSE), stub)
}

W.pdata.frame <- function(x, effect = 1L, w = NULL, cols = is.numeric, na.rm = TRUE, mean = 0, theta = 1,
                          stub = "W.", keep.ids = TRUE, keep.w = TRUE, ...) {

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
    wn <- ckmatch(all.vars(w), nam)
    w <- eval(w[[2L]], x, attr(w, ".Environment")) # w <- x[[wn]]
    cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
    if(keep.w) gn <- c(gn, wn)
  }

  if(length(gn) && length(cols)) {
    ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols])
    return(setAttributes(c(x[gn], .Call(Cpp_BWl,x[cols],fnlevels(g),g,NULL,w,na.rm,theta,ckm(mean),FALSE,FALSE)), ax))
  } else if(!length(gn)) {
    ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
    return(setAttributes(.Call(Cpp_BWl,x[cols],fnlevels(g),g,NULL,w,na.rm,theta,ckm(mean),FALSE,FALSE), ax))
  } else if(is.character(stub)) {
    ax[["names"]] <- paste0(stub, nam)
    return(setAttributes(.Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,theta,ckm(mean),FALSE,FALSE), ax))
  } else return(.Call(Cpp_BWl,`oldClass<-`(x, ax[["class"]]),fnlevels(g),g,NULL,w,na.rm,theta,ckm(mean),FALSE,FALSE))
}

W.data.frame <- function(x, by = NULL, w = NULL, cols = is.numeric, na.rm = TRUE,
                         mean = 0, theta = 1, stub = "W.", keep.by = TRUE, keep.w = TRUE, ...) {

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
      by <- if(is.null(by)) list(0L, 0L, NULL) else G_guo(by)
    }

    if(is.call(w)) {
      wn <- ckmatch(all.vars(w), nam)
      w <- eval(w[[2L]], x, attr(w, ".Environment")) # w <- x[[wn]]
      cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
      if(keep.w) gn <- c(gn, wn)
    }

    if(length(gn)) {
      ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols])
      return(setAttributes(c(x[gn], .Call(Cpp_BWl,x[cols],by[[1L]],by[[2L]],by[[3L]],w,na.rm,theta,ckm(mean),FALSE,FALSE)), ax))
    }
    ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
    return(setAttributes(.Call(Cpp_BWl,x[cols],by[[1L]],by[[2L]],by[[3L]],w,na.rm,theta,ckm(mean),FALSE,FALSE), ax))
  } else if(length(cols)) { # Need to do like this, otherwise list-subsetting drops attributes !
    ax <- attributes(x)
    class(x) <- NULL
    x <- x[cols2int(cols, x, names(x), FALSE)]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }
  if(is.character(stub)) attr(x, "names") <- paste0(stub, attr(x, "names"))

  if(is.null(by)) return(.Call(Cpp_BWl,x,0L,0L,NULL,w,na.rm,theta,ckm(mean),FALSE,FALSE))
  by <- G_guo(by)
  .Call(Cpp_BWl,x,by[[1L]],by[[2L]],by[[3L]],w,na.rm,theta,ckm(mean),FALSE,FALSE)
}

W.list <- function(x, ...) W.data.frame(x, ...)




fbetween <- function(x, ...) UseMethod("fbetween") # , x

fbetween.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fbetween.matrix(x, g, w, na.rm, fill, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_BW,x,0L,0L,NULL,w,na.rm,1,0,TRUE,fill))
  g <- G_guo(g)
  .Call(Cpp_BW,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,1,0,TRUE,fill)
}

fbetween.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- group_effect(x, effect)
  if(is.matrix(x))
  .Call(Cpp_BWm,x,fnlevels(g),g,NULL,w,na.rm,1,0,TRUE,fill) else
  .Call(Cpp_BW,x,fnlevels(g),g,NULL,w,na.rm,1,0,TRUE,fill)
}

fbetween.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_BWm,x,0L,0L,NULL,w,na.rm,1,0,TRUE,fill))
  g <- G_guo(g)
  .Call(Cpp_BWm,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,1,0,TRUE,fill)
}

fbetween.data.frame <- function(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) return(.Call(Cpp_BWl,x,0L,0L,NULL,w,na.rm,1,0,TRUE,fill))
  g <- G_guo(g)
  .Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,1,0,TRUE,fill)
}

fbetween.list <- function(x, ...) fbetween.data.frame(x, ...)

fbetween.pdata.frame <- function(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  g <- group_effect(x, effect)
  .Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,1,0,TRUE,fill)
}

fbetween.grouped_df <- function(x, w = NULL, na.rm = TRUE, fill = FALSE,
                                keep.group_vars = TRUE, keep.w = TRUE, ...) {
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
    ax[["names"]] <- c(nam[gn], nam[-gn2]) # first term is removed if !length(gn)
    res <- .Call(Cpp_BWl, .subset(x, -gn2), g[[1L]],g[[2L]],g[[3L]],w,na.rm,1,0,TRUE,fill)
    if(length(gn)) return(setAttributes(c(.subset(x, gn), res), ax)) else return(setAttributes(res, ax))
  }
  .Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,1,0,TRUE,fill)
}


# Between Operator

B <- function(x, ...) UseMethod("B") # , x

B.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(B.matrix(x, g, w, na.rm, fill, ...))
  fbetween.default(x, g, w, na.rm, fill, ...)
}

B.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, ...)
  fbetween.pseries(x, effect, w, na.rm, fill, ...)

B.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, stub = "B.", ...)
  add_stub(fbetween.matrix(x, g, w, na.rm, fill, ...), stub)

B.grouped_df <- function(x, w = NULL, na.rm = TRUE, fill = FALSE,
                         stub = "B.", keep.group_vars = TRUE, keep.w = TRUE, ...) {
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
    ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[-gn2]) else nam[-gn2])
    res <- .Call(Cpp_BWl, .subset(x, -gn2), g[[1L]],g[[2L]],g[[3L]],w,na.rm,1,0,TRUE,fill)
    if(length(gn)) return(setAttributes(c(.subset(x, gn), res), ax)) else return(setAttributes(res, ax))
  }
  add_stub(.Call(Cpp_BWl,x,g[[1L]],g[[2L]],g[[3L]],w,na.rm,1,0,TRUE,fill), stub)
}

B.pdata.frame <- function(x, effect = 1L, w = NULL, cols = is.numeric, na.rm = TRUE, fill = FALSE,
                          stub = "B.", keep.ids = TRUE, keep.w = TRUE, ...) {
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
    wn <- ckmatch(all.vars(w), nam)
    w <- eval(w[[2L]], x, attr(w, ".Environment")) # w <- x[[wn]]
    cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
    if(keep.w) gn <- c(gn, wn)
  }

  if(length(gn) && length(cols)) {
    ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols])
    return(setAttributes(c(x[gn], .Call(Cpp_BWl,x[cols],fnlevels(g),g,NULL,w,na.rm,1,0,TRUE,fill)), ax))
  } else if(!length(gn)) {
    ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
    return(setAttributes(.Call(Cpp_BWl,x[cols],fnlevels(g),g,NULL,w,na.rm,1,0,TRUE,fill), ax))
  } else if(is.character(stub)) {
      ax[["names"]] <- paste0(stub, nam)
      return(setAttributes(.Call(Cpp_BWl,x,fnlevels(g),g,NULL,w,na.rm,1,0,TRUE,fill), ax))
  } else return(.Call(Cpp_BWl,`oldClass<-`(x, ax[["class"]]),fnlevels(g),g,NULL,w,na.rm,1,0,TRUE,fill))
}

B.data.frame <- function(x, by = NULL, w = NULL, cols = is.numeric, na.rm = TRUE,
                         fill = FALSE, stub = "B.", keep.by = TRUE, keep.w = TRUE, ...) {
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
      by <- if(is.null(by)) list(0L, 0L, NULL) else G_guo(by)
    }

    if(is.call(w)) {
      wn <- ckmatch(all.vars(w), nam)
      w <- eval(w[[2L]], x, attr(w, ".Environment")) # w <- x[[wn]]
      cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
      if(keep.w) gn <- c(gn, wn)
    }

    if(length(gn)) {
      ax[["names"]] <- c(nam[gn], if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols])
      return(setAttributes(c(x[gn], .Call(Cpp_BWl,x[cols],by[[1L]],by[[2L]],by[[3L]],w,na.rm,1,0,TRUE,fill)), ax))
    }
    ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[cols]) else nam[cols]
    return(setAttributes(.Call(Cpp_BWl,x[cols],by[[1L]],by[[2L]],by[[3L]],w,na.rm,1,0,TRUE,fill), ax))
  } else if(length(cols)) { # Necessary, else attributes are dropped by list-subsetting !
    ax <- attributes(x)
    class(x) <- NULL
    x <- x[cols2int(cols, x, names(x), FALSE)]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }
  if(is.character(stub)) attr(x, "names") <- paste0(stub, attr(x, "names"))

  if(is.null(by)) return(.Call(Cpp_BWl,x,0L,0L,NULL,w,na.rm,1,0,TRUE,fill))
  by <- G_guo(by)
  .Call(Cpp_BWl,x,by[[1L]],by[[2L]],by[[3L]],w,na.rm,1,0,TRUE,fill)
}

B.list <- function(x, ...) B.data.frame(x, ...)
