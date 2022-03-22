qtable <- function(...) {
  r <- GRP(list(...))
  tab <- r$group.sizes
  dn <- lapply(unattrib(r$groups), function(x) as.character(funique.default(x)))
  dim <- lengths(dn)
  pdim <- prod(dim)
  if(pdim == r$N.groups) {
    dim(tab) <- dim
    dimnames(tab) <- dn
    return(tab)
  }
  tabi <- integer(pdim)
  tabi[match(GRPnames(r), interact_names(dn))] <- tab
  dim(tabi) <- dim
  dimnames(tabi) <- dn
  tabi
}

qtable2 <- function(...) {
  r <- lapply(list(...), qF, drop = TRUE) # , keep.attr = FALSE
  dn <- lapply(r, attr, "levels")
  faci <- if(length(dn) == 2L) r[[1L]]:r[[2L]] else Reduce(`:`, r)
  tab <- fnobs.default(faci, faci)
  d <- lengths(dn)
  dim(tab) <- d[length(d):1]
  tab <- aperm(tab)
  dimnames(tab) <- dn
  tab
}

qtable3 <- function(..., na.exclude = TRUE) {
  l <- list(...)
  n <- length(l)
  dn <- vector("list", n)
  dm <- integer(n)
  tofact <- function(g) {
    if(is.factor(g)) {
      if(!na.exclude && !inherits(g, "na.included")) return(addNA2(g))
      return(g)
    }
    groupfact(g, ord = FALSE, fact = TRUE, naincl = !na.exclude, keep = FALSE)
  }
  g <- tofact(l[[1L]])
  lev <- attr(g, "levels")
  dn[[1L]] <- lev
  dm[1L] <- ngp <- length(lev)
  attributes(g) <- NULL
  if(n > 1L) for (i in 2:n) {
    gi <- tofact(l[[i]])
    lev <- attr(gi, "levels")
    dn[[i]] <- lev
    dm[i] <- length(lev)
    attributes(gi) <- NULL
    gi %*=% ngp
    g %+=% gi
    ngp <- ngp * length(lev)
  }
  tab <- tabulate(g, nbins = ngp)
  dim(tab) <- dm
  dimnames(tab) <- dn
  if(n == 1L) drop(tab) else tab
}
