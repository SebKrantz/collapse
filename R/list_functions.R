rapply2d <- function(l, FUN, ..., classes = "data.frame") {
  aply2d <- function(y) if(is.list(y) && !inherits(y, classes)) lapply(y, aply2d) else FUN(y, ...) # is.null(dim(y)) # qsu output shows list of DF can have dim attr.
  aply2d(l) # lapply(x,aply2d) # if this is enabled, rapply2d takes apart data.frame if passed
}

get_elem_FUN <- function(x, FUN, return = "sublist", keep_class = FALSE)
  switch(return, sublist = if(keep_class) fcolsubset(x, vapply(`attributes<-`(x, NULL), FUN, TRUE)) else .subset(x, vapply(`attributes<-`(x, NULL), FUN, TRUE)),
         names = attr(x, "names")[vapply(`attributes<-`(x, NULL), FUN, TRUE)],
         indices = which(vapply(`attributes<-`(x, NULL), FUN, TRUE)),
         named_indices = which(`names<-`(vapply(`attributes<-`(x, NULL), FUN, TRUE), attr(x, "names"))),
         logical = vapply(`attributes<-`(x, NULL), FUN, TRUE),
         named_logical = `names<-`(vapply(`attributes<-`(x, NULL), FUN, TRUE), attr(x, "names")),
         stop("Unknown return option!"))

list_elem <- function(l, return = "sublist", keep.class = FALSE)
    get_elem_FUN(l, is.list, return, keep.class)

atomic_elem <- function(l, return = "sublist", keep.class = FALSE)
    get_elem_FUN(l, is.atomic, return, keep.class)


"list_elem<-" <- function(l, value) {
  al <- attributes(l)
  ilv <- is.list(value)
  len <- if(ilv) length(value) else 1L
  attributes(l) <- NULL # vapply without attributes is faster !
  ind <- which(vapply(l, is.list, TRUE))
  if(len != length(ind)) stop("length(value) must match length(list_elem(l))")
  if(ilv) l[ind] <- value else l[[ind]] <- value
  if(ilv && length(nam <- names(value))) al[["names"]][ind] <- nam
  setAttributes(l, al)
}

"atomic_elem<-" <- function(l, value) {
  al <- attributes(l)
  ilv <- is.list(value)
  len <- if(ilv) length(value) else 1L
  attributes(l) <- NULL
  ind <- which(vapply(l, is.atomic, TRUE))
  if(len != length(ind)) stop("length(value) must match length(list_elem(l))")
  if(ilv) l[ind] <- value else l[[ind]] <- value
  if(ilv && length(nam <- names(value))) al[["names"]][ind] <- nam
  setAttributes(l, al)
}

is_regular <- function(x) is.list(x) || is.atomic(x) # fastest way?

is_unlistable <- function(l, DF.as.list = FALSE) if(DF.as.list) all(unlist(rapply(l, is.atomic, how = "list"), use.names = FALSE)) else
  all(unlist(rapply2d(l, is_regular), use.names = FALSE)) # fastest way?

is.unlistable <- function(l, DF.as.list = FALSE) {
  message("Note that 'is.unlistable' was renamed to 'is_unlistable'. It will not be removed anytime soon, but please use updated function names in new code, see help('collapse-renamed')")
  is_unlistable(l, DF.as.list)
}

# If data.frame, search all, otherwise, make optional counting df or not, but don't search them.
ldepth <- function(l, DF.as.list = FALSE) {
  if (inherits(l, "data.frame")) { # fast defining different functions in if-clause ?
    ld <- function(y,i) if(is.list(y)) lapply(y,ld,i+1L) else i
  } else if(DF.as.list) {
    ld <- function(y,i) {
      df <- inherits(y, "data.frame")
      if(is.list(y) && !df) lapply(y,ld,i+1L) else i+df
    }
  } else {
    ld <- function(y,i) if(is.list(y) && !inherits(y, "data.frame")) lapply(y,ld,i+1L) else i
  }
  base::max(unlist(ld(l, 0L), use.names = FALSE))
}

has_elem <- function(l, elem, recursive = TRUE, DF.as.list = FALSE, regex = FALSE, ...) {
  if(is.function(elem)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(recursive) {
     if(DF.as.list) return(any(unlist(rapply(l, elem, how = "list"), use.names = FALSE)))
     return(any(unlist(rapply2d(l, elem), use.names = FALSE)))
    }
    return(any(vapply(l, elem, TRUE, USE.NAMES = FALSE)))
  } else if(is.character(elem)) {
    if(!regex && !missing(...)) unused_arg_action(match.call(), ...)
    if(recursive) {
      oldClass(l) <- NULL # in case [ behaves weird
      is.subl <- if(DF.as.list) is.list else function(x) is.list(x) && !inherits(x, "data.frame") # could do without, but it seems to remove data.frame attributes, and more speed!
      namply <- function(y) if(any(subl <- vapply(y, is.subl, TRUE)))
        c(names(subl), unlist(lapply(.subset(y, subl), namply), use.names = FALSE)) else names(y) # also overall subl names are important, and .subset for DT subsetting ! # names(which(!subl)) # names(y)[!subl] # which is faster?
      if(regex) return(length(rgrep(elem, namply(l), ...)) > 0L) else return(any(namply(l) %in% elem))
    } else if(regex) return(length(rgrep(elem, names(l), ...)) > 0L) else return(any(names(l) %in% elem))
  } else stop("elem must be a function or character vector of element names or regular expressions")
}

# Experimental:
# elem_names <- function(l, how = c("list", "unlist"), DF.as.list = TRUE) { # need right order for method how = list !!
#   namply <- function(y) if(any(subl <- vapply(y, is.subl, TRUE))) c(names(subl), lapply(.subset(y, subl), namply)) else names(subl)
#   switch(how[1L],
#     unlist = names(rapply(l, function(x) NA)),
#     list =
#   ) rapply(l, function(x) NULL)
#
# }

# General note: What about lists containing data.tables ? '[' subsetting will be wrong !
list_extract_FUN <- function(l, FUN, is.subl, keep.tree = FALSE) {
 regsearch <- function(x) {
  if(any(subl <- vapply(x, is.subl, TRUE, USE.NAMES = FALSE))) { # is.list(x) && a
    wsubl <- which(subl)
    wnsubl <- whichv(subl, FALSE)
    matches <- vapply(x[wnsubl], FUN, TRUE, USE.NAMES = FALSE)
    a <- lapply(x[wsubl], regsearch)
    wa <- vlengths(a, FALSE) > 0L # note that this also gets rid of null elements! could make it length or is.null! # vapply(a, length, 1L, USE.NAMES = FALSE)
    x <- c(x[wnsubl][matches], a[wa]) # The problem here: If all elements in a sublist are atomic, it still retains the sublist itself with NULL inside!
    if(keep.tree || length(x) != 1L)
      return(x[forder.int(c(wnsubl[matches], wsubl[wa]))]) else return(x[[1L]]) # fastest way?
  } else if(length(x)) { # This ensures correct behavior in the final nodes: if (length(x)) because problem encountered in get.elem(V, is.matrix) -> empty xlevels list, the lapply below does not execute
    matches <- which(vapply(x, FUN, TRUE, USE.NAMES = FALSE))
    if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]]) # needs to be !=
  }
 }
 regsearch(l)
}

list_extract_regex <- function(l, exp, is.subl, keep.tree = FALSE, ...) {
  regsearch <- function(x) {
  if(any(subl <- vapply(x, is.subl, TRUE, USE.NAMES = FALSE))) {
      matches <- rgrepl(exp, names(x), ...)
      wres <- which(matches)  #   wres <- rgrep(exp, names(x), ...)
      wnressubl <- if(length(wres)) which(subl & !matches) else which(subl) # fsetdiff(which(subl), wres)
    if(length(wnressubl)) { # faster way?
      a <- lapply(x[wnressubl], regsearch) # is this part still necessary?, or only for keep.tree
      wa <- vlengths(a, FALSE) > 0L # note that this also gets rid of null elements!! could make it length or is.null!, length is better for length 0 lists !! #  vapply(a, length, 1L)
      x <- c(x[wres], a[wa])
      if(keep.tree || length(x) != 1L)
        return(x[forder.int(c(wres, wnressubl[wa]))]) else return(x[[1L]])
    } else if(keep.tree || length(wres) != 1L) return(x[wres]) else return(x[[wres]])
  } else { # This ensures correct behavior in the final nodes:
    matches <- rgrep(exp, names(x), ...)
    if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]]) # needs to be !=
  }
 }
 regsearch(l)
}

list_extract_names <- function(l, nam, is.subl, keep.tree = FALSE) {
 regsearch <- function(x) {
  if(any(subl <- vapply(x, is.subl, TRUE, USE.NAMES = FALSE))) {
      matches <- names(x) %in% nam
      wres <- which(matches) # match(nam, names(x), 0L) # better bcause gives integer(0) -> necessary as cannot do l[[0L]]
      wnressubl <- if(length(wres)) which(subl & !matches) else which(subl) # fsetdiff(which(subl), wres)  # old solution: faster but does not work well if parent list is unnamed ! (i.e. l = list(lm1, lm1))
    if(length(wnressubl)) {
      a <- lapply(x[wnressubl], regsearch)
      wa <- vlengths(a, FALSE) > 0L # vapply(a, length, 1L)
      x <- c(x[wres], a[wa])
      if(keep.tree || length(x) != 1L)
        return(x[forder.int(c(wres, wnressubl[wa]))]) else return(x[[1L]])
    } else if(keep.tree || length(wres) != 1L) return(x[wres]) else return(x[[wres]])
  } else {
    matches <- which(names(x) %in% nam)
    if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]]) # needs to be !=, because interger(0) goes in first..
  }
 }
 regsearch(l)
}

# Idea: Also use indices and logical vectors ? i.e. get first two columns of alist of data.frames ?
# This behaves a bit differently (not find elements everywhere, but also subset inside the list)
list_extract_ind <- function(l, ind, is.subl, keep.tree = FALSE) {
  if(is.logical(ind)) ind <- which(ind)
  if(length(ind) > 1L || keep.tree) {
    regsearch <- function(x) if(is.subl(x)) lapply(x, regsearch) else x[ind]
  } else {
    regsearch <- function(x) if(is.subl(x)) lapply(x, regsearch) else x[[ind]]
  }
  regsearch(l)
}

# Note: all functions currently remove empty list elements !
# keep.tree argument still issues wih xlevels

get_elem <- function(l, elem, recursive = TRUE, DF.as.list = FALSE,
                     keep.tree = FALSE, keep.class = FALSE, regex = FALSE, ...) {
  if(recursive) { # possibly if is.list(x) is redundant, because you check above! -> nah, recursive?
    is.subl <- if(DF.as.list) is.list else function(x) is.list(x) && !inherits(x, "data.frame") # could do without, but it seems to remove data.frame attributes
    if(keep.class) al <- attributes(l) # cll <- class(l) # perhaps generalize to other attributes?
    if(is.function(elem)) {
      if(!missing(...)) unused_arg_action(match.call(), ...)
      l <- list_extract_FUN(l, elem, is.subl, keep.tree)
    } else if(is.character(elem)) {
      if(regex) l <- list_extract_regex(l, elem, is.subl, keep.tree, ...) else {
         if(!missing(...)) unused_arg_action(match.call(), ...)
         l <- list_extract_names(l, elem, is.subl, keep.tree)
      }
    } else {
      if(!missing(...)) unused_arg_action(match.call(), ...)
      l <- list_extract_ind(l, elem, is.subl, keep.tree)
    }
    if(keep.class && is.list(l)) {
      al[["names"]] <- names(l)
      return(setAttributes(l, al)) # class(l) <- cll # when drop.tree is proper, l might not be a list
    } else return(l)
  } else {
    if(is.function(elem)) {
      if(!missing(...)) unused_arg_action(match.call(), ...)
      elem <- which(vapply(l, elem, TRUE, USE.NAMES = FALSE))
    } else if(is.character(elem)) {
      if(regex) elem <- rgrep(elem, names(l), ...) else {
        if(!missing(...)) unused_arg_action(match.call(), ...)
        elem <- which(names(l) %in% elem)
      }
    } else if(is.logical(elem)) {
      if(!missing(...)) unused_arg_action(match.call(), ...)
      elem <- which(elem) # else stop("elem must be a function, character vector or vector of regular expressions!")
    }
    if(keep.tree || length(elem) != 1L) {
      if(keep.class) return(fcolsubset(l, elem)) else return(.subset(l, elem)) # <- # base::Filter(elem, l)
    } else return(.subset2(l, elem))
  }
}


# there is base::getElement

reg_elem <- function(l, recursive = TRUE, keep.tree = FALSE, keep.class = FALSE) {
  if(keep.class) al <- attributes(l)
  # if(inherits(l, "data.frame")) if(keep.class) return(l) else return(unattrib(l))
  if(recursive) {
    is.subl <- function(x) is.list(x) && !inherits(x, "data.frame")
    l <- list_extract_FUN(l, is_regular, is.subl, keep.tree)
    if(keep.class && is.list(l)) {
      al[["names"]] <- names(l)
      return(setAttributes(l, al))
    } else return(l)
  } else {
    matches <- which(vapply(l, is_regular, TRUE, USE.NAMES = FALSE)) # l <- base::Filter(is_regular,l)
    if(keep.tree || length(matches) != 1L) {
      if(keep.class) return(fcolsubset(l, matches)) else return(.subset(l, matches))
    } else return(.subset2(l, matches))
  }
}

irreg_elem <- function(l, recursive = TRUE, keep.tree = FALSE, keep.class = FALSE) {
  is.irregular <- function(x) !(is.list(x) || is.atomic(x)) # is.irregular fastest way?
  if(keep.class) al <- attributes(l)
  # if(inherits(l, "data.frame")) stop("A data.frame is a regular object!")
  if(recursive) {
    is.subl <- function(x) is.list(x) && !inherits(x, "data.frame")
    l <- list_extract_FUN(l, is.irregular, is.subl, keep.tree)
    if(keep.class && is.list(l)) {
      al[["names"]] <- names(l)
      return(setAttributes(l, al))
    } else return(l)
  } else {
    matches <- which(vapply(l, is.irregular, TRUE, USE.NAMES = FALSE)) # l <- base::Filter(is_regular,l)
    if(keep.tree || length(matches) != 1L) {
      if(keep.class) return(fcolsubset(l, matches)) else return(.subset(l, matches))
    } else return(.subset2(l, matches))
  }
}

# TODO: See about big objects!
#microbenchmark(all(rapply(lm,is.atomic)),!is.list(unlist(lm, use.names = FALSE)),all(unlist(rapply2d(lm,is.std), use.names = FALSE)))
#microbenchmark(all(rapply(GGDC,is.atomic)),!is.list(unlist(GGDC, use.names = FALSE)),all(unlist(rapply2d(GGDC,is.std), use.names = FALSE)))

