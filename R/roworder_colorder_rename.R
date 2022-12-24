
roworder <- function(X, ..., na.last = TRUE) {
  ovars <- .c(...)
  if(!length(ovars)) stop("... needs to be comma-separated column names, optionally with a '-' prefix for descending order.")
  dec <- startsWith(ovars, "-")
  if(any(dec)) ovars[dec] <- substr(ovars[dec], 2L, 1000000L)
  z <- as.pairlist(.subset(X, ckmatch(ovars, attr(X, "names"))))
  o <- .Call(C_radixsort, na.last, dec, FALSE, FALSE, TRUE, z)
  if(!is.na(na.last) && attr(o, "sorted")) return(condalc(X, inherits(X, "data.table")))
  rn <- attr(X, "row.names")
  res <- .Call(C_subsetDT, X, o, seq_along(unclass(X)), FALSE)
  if(!(is.numeric(rn) || is.null(rn) || rn[1L] == "1")) attr(res, "row.names") <- Csv(rn, o)
  if(inherits(X, "pdata.frame")) {
    index <- findex(X)
    index_o <- .Call(C_subsetDT, index, o, seq_along(unclass(index)), FALSE)
    if(inherits(X, "indexed_frame")) return(reindex(res, index_o))
    attr(res, "index") <- index_o
  }
  res
}

posord <- function(sq, o, pos) switch(pos,
                                      front = c(o, sq[-o]),
                                      end = c(sq[-o], o),
                                      exchange = `[<-`(sq, o[forder.int(o)], value = o),
                                      after = {
                                        if(length(o) == 1L) stop('Need o supply at least 2 columns if pos = "after"')
                                        om1 <- o[-1L]
                                        smo <- sq[-om1]
                                        w1 <- whichv(smo, o[1L])
                                        c(smo[1L:w1], om1, smo[(w1+1L):length(smo)])
                                      },
                                      stop("pos must be 'front', 'end', 'exchange' or 'after'."))

roworderv <- function(X, cols = NULL, neworder = NULL, decreasing = FALSE, na.last = TRUE, pos = "front") {
  if(is.null(neworder)) {
    if(is.null(cols)) {
      if(inherits(X, "sf")) {
        Xo <- X
        oldClass(Xo) <- NULL
        Xo[[attr(Xo, "sf_column")]] <- NULL
        neworder <- radixorderv(Xo, na.last, decreasing)
      } else neworder <- radixorderv(X, na.last, decreasing)
    } else neworder <- radixorderv(colsubset(X, cols), na.last, decreasing)
    if(!is.na(na.last) && attr(neworder, "sorted")) return(condalc(X, inherits(X, "data.table")))
  } else {
    if(!is.integer(neworder)) neworder <- if(is.numeric(neworder)) as.integer(neworder) else if(is.logical(neworder))
                                          which(neworder) else stop("neworder should be integer or logical.")
    if(length(neworder) != fnrow(X)) neworder <- posord(seq_row(X), neworder, pos)
  }
  rn <- attr(X, "row.names")
  res <- .Call(C_subsetDT, X, neworder, seq_along(unclass(X)), FALSE)
  if(!(is.numeric(rn) || is.null(rn) || rn[1L] == "1")) attr(res, "row.names") <- Csv(rn, neworder)
  if(inherits(X, "pdata.frame")) {
    index <- findex(X)
    index_neworder <- .Call(C_subsetDT, index, neworder, seq_along(unclass(index)), FALSE)
    if(inherits(X, "indexed_frame")) return(reindex(res, index_neworder)) # pdata.frame cannot be data.table...
    attr(res, "index") <- index_neworder
  }
  res
}

colorder <- function(.X, ..., pos = "front") { # This also takes names and indices ....
  ax <- attributes(.X)
  oldClass(.X) <- NULL # attributes ?
  nam <- names(.X)
  iX <- seq_along(.X)
  nl <- `names<-`(as.vector(iX, "list"), nam)
  vars <- eval(substitute(c(...)), nl, parent.frame())
  if(!is.integer(vars)) stop(paste0("Unknown columns: ", .c(...)))
  if(length(names(vars))) { # Allow renaming during selection
    nam_vars <- names(vars)
    nonmiss <- nzchar(nam_vars)
    nam[vars[nonmiss]] <- nam_vars[nonmiss]
  }
  if(length(vars) != length(iX)) vars <- posord(iX, vars, pos)
  return(condalcSA(.X[vars], `[[<-`(ax, "names", nam[vars]),
                   any(ax[["class"]] == "data.table")))
}

colorderv <- function(X, neworder = radixorder(names(X)), pos = "front", regex = FALSE, ...) { # This also takes names and indices ....
  ax <- attributes(X)
  oldClass(X) <- NULL # attributes ?
  nam <- names(X)
  if(regex) vars <- rgrep(neworder, nam, ..., sort = FALSE) else {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    vars <- cols2int(neworder, X, nam)
  }
  if(length(vars) != length(X)) vars <- posord(seq_along(X), vars, pos)
  return(condalcSA(X[vars], `[[<-`(ax, "names", nam[vars]),
                   any(ax[["class"]] == "data.table")))
}


frename_core <- function(.x, cols, .nse, ...) {
  args <- if(.nse) substitute(c(...))[-1L] else c(...)
  nam <- attr(.x, "names")
  namarg <- names(args)
  if(is.null(namarg) || !all(nzchar(namarg))) { # The second condition is needed for a function with additional arguments to be passed.
    arg1 <- ..1
    if(length(cols)) ind <- cols2int(cols, .x, nam)
    if(is.function(arg1)) {
      FUN <- if(...length() == 1L) arg1 else # could do special case if ...length() == 2L
        function(x) do.call(arg1, c(list(x), list(...)[-1L]))
      if(is.null(cols)) return(FUN(nam))
      nam[ind] <- FUN(nam[ind])
    } else if(is.character(arg1)) {
      if(is.null(cols)) {
        if(length(arg1) != length(nam)) stop(sprintf("If cols = NULL, the vector or names length = %i must match the object names length = %i.", length(arg1), length(nam)))
        return(arg1)
      }
      if(length(arg1) != length(ind)) stop(sprintf("The vector of names length = %s does not match the number of columns selected = %s.", length(arg1), length(ind)))
      nam[ind] <- arg1
    } else stop("... needs to be expressions colname = newname, a function to apply to the names of columns in cols, or a suitable character vector of names.")
  } else nam[ckmatch(namarg, nam)] <- as.character(args)
  return(nam)
}

frename <- function(.x, ..., cols = NULL, .nse = TRUE) {
  attr(.x, "names") <- frename_core(.x, cols, .nse, ...)
  condalc(.x, inherits(.x, "data.table"))
}

rnm <- frename # rnm clashes with 2 packages.., rme would work but is inconsistent

setrename <- function(.x, ..., cols = NULL, .nse = TRUE) {
  nam <- frename_core(.x, cols, .nse, ...)
  # No longer needed, as also calling setselfref() in C now.
  # if(inherits(.x, "data.table")) {
  #   # Need to allocate here, because the named are captured in ".internal.selfref", so modification be reference still produces an error.
  #   res <- alc(`attr<-`(.x, "names", nam))
  #   assign(as.character(substitute(.x)), res, envir = parent.frame())
  #   return(invisible(res))
  # }
  invisible(.Call(C_setnames, .x, nam))
}

# setrnm <- setrename

relabel <- function(.x, ..., cols = NULL, attrn = "label") { # , sc = TRUE
  args <- list(...)
  nam <- attr(.x, "names")
  namarg <- names(args)
  if(is.null(namarg) || !all(nzchar(namarg))) { # The second condition is needed for a function with additional arguments to be passed.
    arg1 <- args[[1L]]
    if(length(cols)) ind <- cols2int(cols, .x, nam)
    if(is.function(arg1)) {
      lab <- vlabels(.x, attrn, FALSE)
      FUN <- if(length(args) == 1L) arg1 else
        function(x) do.call(arg1, c(list(x), args[-1L]))
      if(is.null(cols)) return(.Call(C_setvlabels, .x, attrn, FUN(lab), NULL))
      args <- FUN(lab[ind])
    } else if(is.character(arg1)) {
      if(is.null(cols)) ind <- if(length(names(arg1))) ckmatch(names(arg1), nam) else NULL
      args <- arg1
    } else stop("... needs to be expressions colname = 'New Label', a function to apply to the names of columns in cols, or a suitable character vector of labels.")
  } else ind <- ckmatch(namarg, nam)
  .Call(C_setvlabels, .x, attrn, args, ind)
}

setrelabel <- function(.x, ..., cols = NULL, attrn = "label")
  invisible(relabel(.x, ..., cols = cols, attrn = attrn))

