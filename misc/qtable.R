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
