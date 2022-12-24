
unlist2d <- function(l, idcols = ".id", row.names = FALSE, recursive = TRUE, id.factor = FALSE, DT = FALSE) {

  if (!is.list(l)) return(l) # stop("l is not a list")
  makeids <- !isFALSE(idcols)
  if(makeids) id.names <- if(isTRUE(idcols)) ".id" else idcols[1L]
  keeprn <- !isFALSE(row.names)
  if(keeprn) row.names <- switch(row.names, `TRUE` = "row.names", row.names)
  idfac <- !isFALSE(id.factor)
  if(idfac) fcclass <- switch(id.factor, `TRUE` = "factor", ordered = c("ordered", "factor"),
                              stop('id.factor needs to be FALSE, TRUE or "ordered"'))
  DATAclass <- if(DT) c("data.table", "data.frame") else "data.frame"

  DFDTl <- function(l) {
    attr(l, "row.names") <- .set_row_names(.Call(C_fnrow, l))
    `oldClass<-`(l, DATAclass)
  }
  # idf <- function(x) if(inherits(x, "data.frame")) 2L else if (!length(x)) 1L else 3L*is.atomic(x) # was if(is.null(x)) 1L -> disregards empty list, bug reported # faster way ? : This is not faster:   2L*inherits(x, "data.frame") + is.null(x) + 3L*is.atomic(x)
  addrn <- function(x) if(any(attr(x, "names") == row.names)) x else c(`names<-`(list(attr(x, "row.names")), row.names), x) # faster way ?
  attol <- function(x) {
    # class(x) <- NULL # tables are also arrays, although only 1D, not because of the class but because they have a dimension attribute.
    if (length(d <- dim(x)) > 1L) { # is.array(x) # length could also be 0... not NULL
      if (length(d) > 2L) { # breaking down HDA
        dn <- dimnames(x)
        dim(x) <- c(d[1L], bprod(d[-1L]))
        if (length(dn)) {
          for (i in 2L:length(d)) if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
          dimnames(x) <- list(dn[[1L]], interact_names(dn[-1L])) # Good ?
        }
      }
      if(keeprn) {
        dn <- dimnames(x)
        x <- `names<-`(c(list(if(is.null(dn[[1L]])) seq_len(d[1L]) else dn[[1L]]), .Call(Cpp_mctl, x, FALSE, 0L)), c(row.names, dn[[2L]]))
      } else x <- .Call(Cpp_mctl, x, TRUE, 0L)
    } else x <- as.vector(x, "list")
    if (is.null(names(x))) names(x) <- paste0("V", seq_along(x))     # it seems this is not yet working for all (i.e. model objects..), also perhaps not start at V1, depending on what other columsn there are.. i.e start at the right position ?
    return(x)
  }
  ul2d <- function(y) {
    if(inherits(y, "data.frame") || is.atomic(y)) return(y)
    if(is.object(y)) oldClass(y) <- NULL # perhaps unclassing y would put more safety ? -> yes !
    ident <- .Call(C_vtypes, y, 6L)  # vapply(`attributes<-`(y, NULL), idf, 1L) # removes names ?
    if(is.list(y) && all(ident > 0L)) {
      if(any(at <- ident == 3L)) y[at] <- lapply(y[at], attol)
      if(keeprn && any(df <- ident == 2L)) y[df] <- lapply(y[df], addrn) # better cbind for data.table ? or x[["row.names"]] =.. and the sort later ?
        if(makeids) {
          if(idfac) {
            y <- y[ident != 1L] # better way ? y[ident!=1L] = NULL ?
            nam <- names(y)
            if(length(nam)) names(y) <- NULL else nam <- as.character(seq_along(y))
            y <- DFDTl(.Call(C_rbindlist, y, TRUE, TRUE, id.names))
            setattributes(.subset2(y, 1L), pairlist(levels = nam, class = fcclass))
            return(y)
          } else return(DFDTl(.Call(C_rbindlist, y[ident != 1L], TRUE, TRUE, id.names)))
        } else return(DFDTl(.Call(C_rbindlist, y[ident != 1L], TRUE, TRUE, NULL)))
    } else lapply(y, ul2d)
  }

  l <- ul2d(l)
  if(recursive) {
    while(!inherits(l, "data.frame")) l <- ul2d(l)
    if(makeids) {
      nams <- attr(l, "names")
      ids <- whichv(nams, id.names)
      nid <- length(ids)
      if(nid > 1L) {
        nids <- seq_len(nid)
        attr(l, "names")[ids] <- if(length(idcols) == nid) idcols else paste(id.names, nids, sep = ".")
        if(keeprn) {
          rn <- whichv(nams, row.names) # with more id's, row.names are automatically generated from the sub-data.frames..
          if(!all(ids == nids) || rn != nid + 1L) .Call(C_setcolorder, l, c(ids, rn, seq_along(nams)[-c(ids, rn)]))
        } else if (!all(ids == nids)) .Call(C_setcolorder, l, c(ids, seq_along(nams)[-ids]))
      } else if(keeprn) { # makes sure row.names comes after ids, even if only one id!
        rn <- whichv(nams, row.names) # length(rn) needed when only vectors... no row names column...
        if(length(rn) && rn != 2L) .Call(C_setcolorder, l,  c(ids, rn, seq_along(nams)[-c(ids, rn)]))
      }
    } else if (keeprn) {
      nams <- attr(l, "names")
      rn <- whichv(nams, row.names)
      if(length(rn) && rn != 1L) .Call(C_setcolorder, l, c(rn, seq_along(nams)[-rn]))
    }
    if(DT) return(alc(l))
  }
  # attr(l, ".internal.selfref") <- NULL
  l
}
