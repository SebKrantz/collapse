
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


as.data.frame.qsu <- function(x, ..., gid = "Group", stringsAsFactors = TRUE) {
  d <- dim(x)
  dn <- dimnames(x)
  stnam <- dn[[2L]]
  if(is.null(d)) {
    res <- as.vector(x, "list")
    attr(res, "row.names") <- 1L
    # res <- list(Statistic = names(x), Value = unattrib(x))
    # attr(res, "row.names") <- .set_row_names(length(x))
  } else if(length(d) == 2L) {
    varl <- if(stringsAsFactors) list(`attributes<-`(seq_len(d[1L]), list(levels = dn[[1L]], class = "factor"))) else dn[1L]
    res <- c(varl,  mctl(x))
    names(res) <- c(if(stnam[1L] == "N") "Variable" else "Trans", stnam)
    attr(res, "row.names") <- .set_row_names(d[1L])
  } else if(length(d) == 3L) {
    dimnames(x) <- NULL
    # Special case: qsu(wlddev, PCGDP ~ region, ~ iso3c)
    if(d[3L] == 3L && dn[[3L]][1L] == "Overall") {
      res <- aperm.default(x, c(3L,1L,2L))
      d[c(1L, 3L)] <- d[c(3L, 1L)]
      dn[c(1L, 3L)] <- dn[c(3L, 1L)]
      vn <- gid
    } else {
      vn <- "Variable"
      res <- aperm.default(x, c(1L,3L,2L))
    }
    attributes(res) <- NULL
    dim(res) <- c(d[1L]*d[3L], d[2L])
    varsl <- if(stringsAsFactors) list(`attributes<-`(rep(seq_len(d[3L]), each = d[1L]), list(levels =  dn[[3L]], class = "factor")),
                                `attributes<-`(rep(seq_len(d[1L]), d[3L]), list(levels = dn[[1L]], class = "factor"))) else
                           list(rep(dn[[3L]], each = d[1L]), rep(dn[[1L]], d[3L]))
    res <- c(varsl, mctl(res))
    names(res) <- c(vn, if(stnam[1L] == "N") gid else "Trans", stnam)
    attr(res, "row.names") <- .set_row_names(d[1L]*d[3L])
  } else {
    dimnames(x) <- NULL
    res <- aperm.default(x, c(3L,1L,4L,2L))
    attributes(res) <- NULL
    nr <- d[1L]*3L*d[4L]
    dim(res) <- c(nr, d[2L])
    varsl <- if(stringsAsFactors)
              list(`attributes<-`(rep(seq_len(d[4L]), each = 3L*d[1L]), list(levels = dn[[4L]], class = "factor")),
                   `attributes<-`(rep(seq_len(d[1L]), d[4L], each = 3L), list(levels = dn[[1L]], class = "factor")),
                   `attributes<-`(rep(seq_len(d[3L]), d[1L]*d[4L]), list(levels = dn[[3L]], class = "factor"))) else
              list(rep(dn[[4L]], each = 3L*d[1L]), rep(dn[[1L]], d[4L], each = 3L), rep(dn[[3L]], d[1L]*d[4L]))
    res <- c(varsl,  mctl(res))
    names(res) <- c("Variable", gid, "Trans", stnam)
    attr(res, "row.names") <- .set_row_names(nr)
  }
  class(res) <- "data.frame"
  res
}

