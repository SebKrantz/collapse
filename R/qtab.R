

qtab <- function(..., w = NULL, wFUN = NULL, wFUN.args = NULL,
                   dnn = "auto", sort = .op[["sort"]],
                   na.exclude = TRUE, drop = FALSE, method = "auto") {
  ll <- ...length() == 1L && is.list(..1)
  l <- if(ll) unclass(..1) else list(...)
  n <- length(l)
  dn <- vector("list", n)
  dm <- integer(n)

  names(dn) <- if(is.character(dnn)) {
    if(length(dnn) > 1L) dnn else {
      nam <- names(l)
      nam <- switch(dnn, auto =, namlab =
               if(ll) nam else if(is.null(nam)) .c(...) else
               if(all(has_nam <- nzchar(nam))) nam else
                 `[<-`(nam, !has_nam, value = .c(...)[!has_nam]), dnn)
      if(dnn != "namlab") nam else paste(nam, setv(vlabels(l, use.names = FALSE), NA, ""), sep = ": ")
    }
  } else if(is.function(dnn)) dnn(l) else unlist(dnn, use.names = FALSE)

  # tofact <- function(g) {
  #   if(is.factor(g)) {
  #     if(!na.exclude && !inherits(g, "na.included")) return(addNA2(g))
  #     return(g)
  #   }
  #   groupfact(g, ord = FALSE, fact = TRUE, naincl = !na.exclude, keep = FALSE)
  # }
  g <- qF(l[[1L]], sort = sort, na.exclude = na.exclude, drop = drop, method = method)
  lev <- attr(g, "levels")
  dn[[1L]] <- lev
  dm[1L] <- ngp <- length(lev)
  attributes(g) <- NULL

  if(n > 1L) for (i in 2:n) {
    gi <- qF(l[[i]], sort = sort, na.exclude = na.exclude, drop = drop, method = method)
    lev <- attr(gi, "levels")
    dn[[i]] <- lev
    dm[i] <- length(lev)
    # attributes(gi) <- NULL
    # unattrib(x) + (unattrib(y) - 1L) * fnlevels(x)
    # NA values cause integer overflows...
    # gi %-=% 1L
    # gi %*=% ngp
    # g %+=% gi
    # TODO: what if g is not a deep copy?? -> seems to work so far. I guess qF() or attributes(g) <- NULL creates a deep copy?
    .Call(C_fcrosscolon, g, ngp, gi, na.exclude)
    ngp <- ngp * length(lev)
  }

  if(is.null(w) || is.null(wFUN))
     tab <- .Call(C_fwtabulate, g, w, ngp, na.exclude) # tabulate(g, nbins = ngp)
  else {
    if(is.function(wFUN)) {
       wf <- l1orlst(as.character(substitute(wFUN)))
    } else if (is.character(wFUN)) {
       wf <- wFUN
       wFUN <- match.fun(wFUN)
    } else stop("wFUN needs to be a function or function-string")
    if(na.exclude && anyNA(g)) {
      nna <- whichNA(g, invert = TRUE)
      w <- Csv(w, nna)
      g <- Csv(g, nna)
    }
    attr(g, "N.groups") <- ngp
    oldClass(g) <- c("qG", "na.included")
    if(is.null(wFUN.args)) {
      tab <- if(any(wf == .FAST_STAT_FUN)) wFUN(w, g = g, use.g.names = FALSE) else
             splaplfun(w, g, wFUN)
    } else {
      tab <- if(any(wf == .FAST_STAT_FUN)) do.call(wFUN, c(list(x = w, g = g, use.g.names = FALSE), wFUN.args)) else
             do.call(splaplfun, c(list(x = w, g = g, FUN = wFUN), wFUN.args))
    }
  }

  dim(tab) <- dm
  dimnames(tab) <- dn
  oldClass(tab) <- c("qtab", "table")
  attr(tab, "sorted") <- sort
  attr(tab, "weighted") <- !is.null(w)
  tab
}

qtable <- function(...) qtab(...)
