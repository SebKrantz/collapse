
qsu <- function(x, ...) UseMethod("qsu") # , x

qsu.default <- function(x, g = NULL, pid = NULL, w = NULL, higher = FALSE, array = TRUE, stable.algo = TRUE, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(qsu.matrix(x, g, pid, w, higher, array, stable.algo, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) {
    if(is.null(pid)) return(fbstatsCpp(x,higher, w = w, stable.algo = stable.algo))
    pid <- G_guo(pid)
    return(fbstatsCpp(x,higher,0L,0L,pid[[1L]],pid[[2L]],w,stable.algo))
  }
  if(is.atomic(g)) {
    if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
    lev <- attr(g, "levels")
    if(is.null(pid)) return(fbstatsCpp(x,higher,length(lev),g,0L,0L,w,stable.algo,TRUE,TRUE,lev))
    pid <- G_guo(pid)
    return(fbstatsCpp(x,higher,length(lev),g,pid[[1L]],pid[[2L]],w,stable.algo,array,TRUE,lev))
  }
  if(!is_GRP(g)) g <- GRP.default(g, call = FALSE)
  if(is.null(pid)) return(fbstatsCpp(x,higher,g[[1L]],g[[2L]],0L,0L,w,stable.algo,TRUE,TRUE,GRPnames(g)))
  pid <- G_guo(pid)
  fbstatsCpp(x,higher,g[[1L]],g[[2L]],pid[[1L]],pid[[2L]],w,stable.algo,array,TRUE,GRPnames(g))
}

qsu.pseries <- function(x, g = NULL, w = NULL, effect = 1L, higher = FALSE, array = TRUE, stable.algo = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  pid <- group_effect(x, effect)
  if(is.null(g)) return(fbstatsCpp(x,higher,0L,0L,fnlevels(pid),pid,w,stable.algo))
  if(is.atomic(g)) {
    if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
    lev <- attr(g, "levels")
    return(fbstatsCpp(x,higher,length(lev),g,fnlevels(pid),pid,w,stable.algo,array,TRUE,lev))
  }
  if(!is_GRP(g)) g <- GRP.default(g, call = FALSE)
  fbstatsCpp(x,higher,g[[1L]],g[[2L]],fnlevels(pid),pid,w,stable.algo,array,TRUE,GRPnames(g))
}

qsu.matrix <- function(x, g = NULL, pid = NULL, w = NULL, higher = FALSE, array = TRUE, stable.algo = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.null(g)) {
    if(is.null(pid)) return(fbstatsmCpp(x,higher, w = w, stable.algo = stable.algo))
    pid <- G_guo(pid)
    return(fbstatsmCpp(x,higher,0L,0L,pid[[1L]],pid[[2L]],w,stable.algo,array))
  }
  if(is.atomic(g)) {
    if(!is.nmfactor(g)) g <- qF(g, na.exclude = FALSE)
    lev <- attr(g, "levels")
    if(is.null(pid)) return(fbstatsmCpp(x,higher,length(lev),g,0L,0L,w,stable.algo,array,lev))
    pid <- G_guo(pid)
    return(fbstatsmCpp(x,higher,length(lev),g,pid[[1L]],pid[[2L]],w,stable.algo,array,lev))
  }
  if(!is_GRP(g)) g <- GRP.default(g, call = FALSE)
  if(is.null(pid)) return(fbstatsmCpp(x,higher,g[[1L]],g[[2L]],0L,0L,w,stable.algo,array,GRPnames(g)))
  pid <- G_guo(pid)
  fbstatsmCpp(x,higher,g[[1L]],g[[2L]],pid[[1L]],pid[[2L]],w,stable.algo,array,GRPnames(g))
}

qsu.data.frame <- function(x, by = NULL, pid = NULL, w = NULL, cols = NULL, higher = FALSE, array = TRUE, vlabels = FALSE, stable.algo = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  formby <- is.call(by)
  formpid <- is.call(pid)
  formw <- is.call(w)

  # fastest solution!! (see checks below !!)
  if(formby || formpid || formw) {
    v <- NULL
    class(x) <- NULL
    nam <- names(x)
    if(formby) {
      if(length(by) == 3L) {
        v <- ckmatch(all.vars(by[[2L]]), nam)
        byn <- ckmatch(all.vars(by[[3L]]), nam)
      } else byn <- ckmatch(all.vars(by), nam)
      by <- if(length(byn) == 1L) x[[byn]] else GRP.default(x[byn], call = FALSE)
    } else byn <- NULL
    if(formpid) {
      if(length(pid) == 3L) {
        v <- ckmatch(all.vars(pid[[2L]]), nam)
        pidn <- ckmatch(all.vars(pid[[3L]]), nam)
      } else pidn <- ckmatch(all.vars(pid), nam)
      pid <- if(length(pidn) == 1L) x[[pidn]] else GRP.default(x[pidn], return.groups = FALSE, call = FALSE)
    } else pidn <- NULL
    if(formw) {
      widn <- ckmatch(all.vars(w), nam)
      w <- eval(w[[2L]], x, attr(w, ".Environment")) # w <- x[[widn]]
    } else widn <- NULL
    if(is.null(v)) {
      x <- if(is.null(cols)) x[-c(byn, pidn, widn)] else x[cols2int(cols, x, nam, FALSE)]
    } else x <- x[v]
  } else if(length(cols)) x <- .subset(x, cols2int(cols, x, attr(x, "names"), FALSE))

  # Get vlabels
  if(is.function(vlabels) || vlabels)
    attr(x, "names") <- if(is.function(vlabels)) vlabels(x) else
      paste(attr(x, "names"), setv(vlabels(x, use.names = FALSE), NA, ""), sep = ": ")

  # original code:
  if(is.null(by)) {
    if(is.null(pid)) return(fbstatslCpp(x,higher, w = w, stable.algo = stable.algo))
    pid <- G_guo(pid)
    return(drop(fbstatslCpp(x,higher,0L,0L,pid[[1L]],pid[[2L]],w,stable.algo,array)))
  }
  if(is.atomic(by)) {
    if(!is.nmfactor(by)) by <- qF(by, na.exclude = FALSE)
    lev <- attr(by, "levels")
    if(is.null(pid)) return(drop(fbstatslCpp(x,higher,length(lev),by,0L,0L,w,stable.algo,array,lev)))
    pid <- G_guo(pid)
    return(drop(fbstatslCpp(x,higher,length(lev),by,pid[[1L]],pid[[2L]],w,stable.algo,array,lev)))
  }
  if(!is_GRP(by)) by <- GRP.default(by, call = FALSE)
  if(is.null(pid)) return(drop(fbstatslCpp(x,higher,by[[1L]],by[[2L]],0L,0L,w,stable.algo,array,GRPnames(by))))
  pid <- G_guo(pid)
  drop(fbstatslCpp(x,higher,by[[1L]],by[[2L]],pid[[1L]],pid[[2L]],w,stable.algo,array,GRPnames(by)))
}

qsu.list <- function(x, ...) qsu.data.frame(x, ...)

qsu.sf <- function(x, by = NULL, pid = NULL, w = NULL, cols = NULL, higher = FALSE, array = TRUE, vlabels = FALSE, stable.algo = TRUE, ...) {
  oldClass(x) <- NULL
  x[[attr(x, "sf_column")]] <- NULL
  qsu.data.frame(x, by, pid, w, cols, higher, array, vlabels, stable.algo, ...)
}

qsu.pdata.frame <- function(x, by = NULL, w = NULL, cols = NULL, effect = 1L, higher = FALSE, array = TRUE, vlabels = FALSE, stable.algo = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  pid <- group_effect(x, effect)
  x <- unindex(x)

  formby <- is.call(by)
  formw <- is.call(w)

  # fastest solution
  if(formby || formw) {
    v <- NULL
    class(x) <- NULL
    nam <- names(x)
    if(formby) {
      if(length(by) == 3L) {
        v <- ckmatch(all.vars(by[[2L]]), nam)
        byn <- ckmatch(all.vars(by[[3L]]), nam)
      } else byn <- ckmatch(all.vars(by), nam)
      by <- if(length(byn) == 1L) x[[byn]] else GRP.default(x[byn])
    } else byn <- NULL
    if(formw) {
      widn <- ckmatch(all.vars(w), nam)
      w <- eval(w[[2L]], x, attr(w, ".Environment")) # w <- x[[widn]]
    } else widn <- NULL
    if(is.null(v)) {
      x <- if(is.null(cols)) x[-c(byn, widn)] else x[cols2int(cols, x, nam, FALSE)]
    } else x <- x[v]
  } else if(length(cols)) x <- .subset(x, cols2int(cols, x, attr(x, "names"), FALSE))

  if(is.function(vlabels) || vlabels)
    attr(x, "names") <- if(is.function(vlabels)) vlabels(x) else
       paste(attr(x, "names"), setv(vlabels(x, use.names = FALSE), NA, ""), sep = ": ")

  if(is.null(by)) return(drop(fbstatslCpp(x,higher,0L,0L,fnlevels(pid),pid,w,stable.algo,array)))
  if(is.atomic(by)) {
    if(!is.nmfactor(by)) by <- qF(by, na.exclude = FALSE)
    lev <- attr(by, "levels")
    return(drop(fbstatslCpp(x,higher,length(lev),by,fnlevels(pid),pid,w,stable.algo,array,lev)))
  }
  if(!is_GRP(by)) by <- GRP.default(by, call = FALSE)
  drop(fbstatslCpp(x,higher,by[[1L]],by[[2L]],fnlevels(pid),pid,w,stable.algo,array,GRPnames(by)))
}


# Try to speed up ! Printing Takes 100 milliseconds on WDI !
print.qsu <- function(x, digits = 4, nonsci.digits = 9, na.print = "-", return = FALSE, print.gap = 2, ...) {
  vec2mat <- function(x) if(is.array(x)) x else  # outer(1, x) # for variable spacing in vector printing...
    `attributes<-`(x, list(dim = c(1L, length(x)), dimnames = list("", names(x)))) # faster and better !!
  formatfun <- function(x) { # , drop0trailing = FALSE redundat ??
    class(x) <- NULL
    xx <- formatC(vec2mat(round(x, digits)), format = "g", flag = "#",
                  digits = nonsci.digits, big.mark = "'", big.interval = 6, # "\u2009": https://stackoverflow.com/questions/30555232/using-a-half-space-as-a-big-mark-for-knitr-output
                  drop0trailing = TRUE, preserve.width = "individual") # format(unclass(round(x,2)), digits = digits, drop0trailing = TRUE, big.mark = ",", big.interval = 6, scientific = FALSE)
    if(any(ina <- is.na(x))) xx[ina] <- na.print
    xx <- gsub(" ", "", xx, fixed = TRUE) # remove some weird white space (qsu(GGDS10S))
    return(xx)
  }
  xx <- if(is.atomic(x)) formatfun(x) else rapply(x, formatfun, how = "list") # No longer necessary, but keep, maybe you want to print lists using print.qsu.
  if(return) return(xx) else print.default(xx, quote = FALSE, right = TRUE, print.gap = print.gap, ...)
  invisible(x)
}
# View.qsu <- function(x) View(unclass(x))


aperm.qsu <- function(a, perm = NULL, resize = TRUE, keep.class = TRUE, ...) {
  r <- aperm.default(a, perm, resize = resize)
  if(keep.class) oldClass(r) <- oldClass(a)
  r
}


`[.qsu` <- function(x, i, j, ..., drop = TRUE) `oldClass<-`(NextMethod(), oldClass(x))




# testing formula inputs:

# best!!! (subsetting unclassed objects is better), and reassigning x does not take more memry than deleting columns, it is often faster !!, also GRP does not mind unclassed objects !!
# formtest2 <- function(x, by = NULL, xt = NULL, cols = NULL) {
#   formby <- is.call(by)
#   formxt <- is.call(xt)
#
#   # fastest solution: (check: is reassigning x memory efficient ?? should not rather delete columns ??)
#   if(formby || formxt) {
#     v <- NULL
#     class(x) <- NULL # this is faster !!
#     if(formby) {
#       if(length(by) == 3L) {
#         v <- all.vars(by[[2L]])
#         namby <- all.vars(by[[3L]])
#       } else namby <- all.vars(by)
#       by <- if(length(namby) == 1L) x[[namby]] else GRP(x, namby)
#     } else namby <- NULL
#     if(formxt) {
#       if(length(xt) == 3L) {
#         v <- all.vars(xt[[2L]])
#         namxt <- all.vars(xt[[3L]])
#       } else namxt <- all.vars(xt)
#       xt <- if(length(namxt) == 1L) x[[namxt]] else GRP(x, namxt)
#     } else namxt <- NULL
#     if(is.null(v)) { # reassign ?? or set NULL ??? what is more memory efficient ??
#       x <- if(is.null(cols)) x[-match(c(namby,namxt), names(x))] else if(is.function(cols))
#             x[vapply(x, cols, TRUE)] else x[cols]
#     } else x <- x[v]
#   } else if(length(cols)) {
#     # class(x) <- NULL # but unclass is faster !!
#     x <- if(is.function(cols)) unclass(x)[vapply(x, cols, TRUE)] else .subset(x, cols)
#   }
#   return(list(x, by, xt))
# }
#
# formtest1 <- function(x, by = NULL, xt = NULL, cols = NULL) {
#   formby <- is.call(by)
#   formxt <- is.call(xt)
#
#   # fastest solution: (check: is reassigning x memory efficient ?? should not rather delete columns ??)
#   if(formby || formxt) {
#     v <- NULL
#     if(formby) {
#       if(length(by) == 3L) {
#         v <- all.vars(by[[2L]])
#         namby <- all.vars(by[[3L]])
#       } else namby <- all.vars(by)
#       by <- if(length(namby) == 1L) x[[namby]] else GRP(x, namby)
#     } else namby <- NULL
#     if(formxt) {
#       if(length(xt) == 3L) {
#         v <- all.vars(xt[[2L]])
#         namxt <- all.vars(xt[[3L]])
#       } else namxt <- all.vars(xt)
#       xt <- if(length(namxt) == 1L) x[[namxt]] else GRP(x, namxt)
#     } else namxt <- NULL
#     if(is.null(v)) { # reassign ?? or set NULL ??? what is more memory efficient ??
#       x <- if(is.null(cols)) unclass(x)[-match(c(namby,namxt), names(x))] else if(is.function(cols))
#         unclass(x)[vapply(x, cols, TRUE)] else .subset(x, cols)
#     } else x <- unclass(x)[v]
#   } else if(length(cols)) {
#     x <- if(is.function(cols)) unclass(x)[vapply(x, cols, TRUE)] else .subset(x, cols)
#   }
#   return(list(x, by, xt))
# }
#
# formtest3 <- function(x, by = NULL, xt = NULL, cols = NULL) {
#   formby <- is.call(by)
#   formxt <- is.call(xt)
#
#   # fastest solution: (check: is reassigning x memory efficient ?? should not rather delete columns ??)
#   if(formby || formxt) {
#     v <- NULL
#     class(x) <- NULL
#     if(formby) {
#       if(length(by) == 3L) {
#         v <- all.vars(by[[2L]])
#         namby <- all.vars(by[[3L]])
#       } else namby <- all.vars(by)
#       by <- if(length(namby) == 1L) x[[namby]] else GRP(x, namby)
#     } else namby <- NULL
#     if(formxt) {
#       if(length(xt) == 3L) {
#         v <- all.vars(xt[[2L]])
#         namxt <- all.vars(xt[[3L]])
#       } else namxt <- all.vars(xt)
#       xt <- if(length(namxt) == 1L) x[[namxt]] else GRP(x, namxt)
#     } else namxt <- NULL
#     if(is.null(v)) { # reassign ?? or set NULL ??? what is more memory efficient ??
#       if(is.null(cols)) x[match(c(namby,namxt), names(x))] <- NULL else if(is.function(cols))
#         x[!vapply(x, cols, TRUE)] <- NULL else x[-cols] <- NULL
#     } else x[-v] <- NULL
#   } else if(length(cols)) {
#     if(is.function(cols)) x[!vapply(x, cols, TRUE)] <- NULL else x[-cols] <- NULL
#   }
#   return(list(x, by, xt))
# }
#



