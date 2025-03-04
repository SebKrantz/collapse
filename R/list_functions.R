rapply2d <- function(l, FUN, ..., classes = "data.frame") {
  aply2d <- function(y) if(is.list(y) && !inherits(y, classes)) lapply(y, aply2d) else FUN(y, ...) # is.null(dim(y)) # qsu output shows list of DF can have dim attr.
  aply2d(l) # lapply(x,aply2d) # if this is enabled, rapply2d takes apart data.frame if passed
}

get_elem_indl <- function(x, indl, return = "sublist", keep_class = FALSE)
  switch(return, sublist = if(keep_class) fcolsubset(x, indl) else .subset(x, indl),
         names = attr(x, "names")[indl],
         indices = which(indl),
         named_indices = which(`names<-`(indl, attr(x, "names"))),
         logical = indl,
         named_logical = `names<-`(indl, attr(x, "names")),
         stop("Unknown return option!"))

list_elem <- function(l, return = "sublist", keep.class = FALSE) {
  if(!is.list(l)) stop("l needs to be a list")
  get_elem_indl(l, .Call(C_vtypes, l, 3L), return, keep.class)
}

atomic_elem <- function(l, return = "sublist", keep.class = FALSE) {
  if(!is.list(l)) stop("l needs to be a list")
  get_elem_indl(l, .Call(C_vtypes, l, 7L), return, keep.class)
}


"list_elem<-" <- function(l, value) {
  if(!is.list(l)) stop("l needs to be a list")
  al <- attributes(l)
  ilv <- is.list(value)
  len <- if(ilv) length(value) else 1L
  attributes(l) <- NULL # vapply without attributes is faster !
  ind <- which(.Call(C_vtypes, l, 3L))
  if(len != length(ind)) stop("length(value) must match length(list_elem(l))")
  if(ilv) l[ind] <- value else l[[ind]] <- value
  if(ilv && length(nam <- names(value))) al[["names"]][ind] <- nam
  setAttributes(l, al)
}

"atomic_elem<-" <- function(l, value) {
  if(!is.list(l)) stop("l needs to be a list")
  al <- attributes(l)
  ilv <- is.list(value)
  len <- if(ilv) length(value) else 1L
  attributes(l) <- NULL
  ind <- which(.Call(C_vtypes, l, 7L))
  if(len != length(ind)) stop("length(value) must match length(list_elem(l))")
  if(ilv) l[ind] <- value else l[[ind]] <- value
  if(ilv && length(nam <- names(value))) al[["names"]][ind] <- nam
  setAttributes(l, al)
}

is_unlistable <- function(l, DF.as.list = FALSE) {
  if(!is.list(l)) return(TRUE)
  if(DF.as.list) return(all(unlist(rapply(l, is.atomic, how = "list"), use.names = FALSE)))
  checkisul <- function(x) if(is.atomic(x) || inherits(x, "data.frame")) TRUE else if(is.list(x)) lapply(x, checkisul) else FALSE
  all(unlist(checkisul(l), use.names = FALSE)) # fastest way?
}

# is.unlistable <- function(l, DF.as.list = FALSE) {
#   .Deprecated(msg = "'is.unlistable' was renamed to 'is_unlistable'. It will be removed end of 2023, see help('collapse-renamed').")
#   is_unlistable(l, DF.as.list)
# }

# If data.frame, search all, otherwise, make optional counting df or not, but don't search them.
ldepth <- function(l, DF.as.list = FALSE) {
  if(!is.list(l)) return(0L)
  if(inherits(l, "data.frame")) { # fast defining different functions in if-clause ?
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
  if(!is.list(l)) stop("l needs to be a list")
  if(is.function(elem)) {
    if(recursive) {
     if(DF.as.list) {
       raply2 <- function(y) if(elem(y, ...)) TRUE else if(is.list(y)) lapply(y, raply2) else FALSE
       return(any(unlist(raply2(l), use.names = FALSE)))
     }
     aply2de <- function(y) if(elem(y, ...)) TRUE else if(is.list(y) && !inherits(y, "data.frame")) lapply(y, aply2de) else FALSE
     return(any(unlist(aply2de(l), use.names = FALSE)))
    }
    return(any(vapply(l, elem, TRUE, ..., USE.NAMES = FALSE)))
  } else if(is.character(elem)) {
    if(!regex && !missing(...)) unused_arg_action(match.call(), ...)
    if(recursive) {
      oldClass(l) <- NULL # in case [ behaves weird
      ret <- 4L - as.logical(DF.as.list) # is.subl <- if(DF.as.list) is.list else function(x) is.list(x) && !inherits(x, "data.frame") # could do without, but it seems to remove data.frame attributes, and more speed!
      namply <- function(y) if(any(subl <- .Call(C_vtypes, y, ret))) # vapply(y, is.subl, TRUE)
        c(names(y), unlist(lapply(.subset(y, subl), namply), use.names = FALSE)) else names(y) # also overall subl names are important, and .subset for DT subsetting ! # names(which(!subl)) # names(y)[!subl] # which is faster?
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


list_extract_FUN <- function(l, FUN, ret, keep.tree = FALSE, nkeep_class = TRUE, invert = FALSE, ...) {
 if(invert) {
   regsearch <- function(x) {
     if(nkeep_class && is.object(x)) oldClass(x) <- NULL
     if(any(subl <- .Call(C_vtypes, x, ret))) {
       matches <- !vapply(x, FUN, TRUE, ..., USE.NAMES = FALSE)
       wsubl <- which(matches & subl)
       if(length(wsubl)) {
         wres <- which(matches & !subl)
         a <- lapply(x[wsubl], regsearch)
         wa <- vlengths(a, FALSE) > 0L
         x <- c(x[wres], a[wa])
         if(keep.tree || length(x) != 1L)
           return(x[forder.int(c(wres, wsubl[wa]))]) else return(x[[1L]])
       } else {
         wres <- which(matches)
         if(keep.tree || length(wres) != 1L) return(x[wres]) else return(x[[wres]])
       }
     } else {
       matches <- whichv(vapply(x, FUN, TRUE, ..., USE.NAMES = FALSE), FALSE)
       if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]])
     }
   }
 } else {
   regsearch <- function(x) {
     if(nkeep_class && is.object(x)) oldClass(x) <- NULL
     if(any(subl <- .Call(C_vtypes, x, ret))) {
       matches <- vapply(x, FUN, TRUE, ..., USE.NAMES = FALSE)
       wres <- which(matches)
       wnressubl <- which(if(length(wres)) subl & !matches else subl)
       if(length(wnressubl)) {
         a <- lapply(x[wnressubl], regsearch)
         wa <- vlengths(a, FALSE) > 0L
         x <- c(x[wres], a[wa])
         if(keep.tree || length(x) != 1L)
           return(x[forder.int(c(wres, wnressubl[wa]))]) else return(x[[1L]])
       } else if(keep.tree || length(wres) != 1L) return(x[wres]) else return(x[[wres]])
     } else {
       matches <- which(vapply(x, FUN, TRUE, ..., USE.NAMES = FALSE))
       if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]])
     }
   }
   ## Previous Version: Does not check the sublists, so cannot find objects through inherits()
   # if(invert) {
   #   # This is rather simple, just negate the vapply calls. could also simple invert the function.. but this is faster...
   #   regsearch <- function(x) {
   #     if(nkeep_class && is.object(x)) oldClass(x) <- NULL
   #     if(any(subl <- .Call(C_vtypes, x, ret))) {
   #       wsubl <- which(subl)
   #       wnsubl <- whichv(subl, FALSE)
   #       matches <- !vapply(x[wnsubl], FUN, TRUE, USE.NAMES = FALSE)
   #       a <- lapply(x[wsubl], regsearch)
   #       wa <- vlengths(a, FALSE) > 0L
   #       x <- c(x[wnsubl][matches], a[wa])
   #       if(keep.tree || length(x) != 1L)
   #         return(x[forder.int(c(wnsubl[matches], wsubl[wa]))]) else return(x[[1L]])
   #     } else if(length(x)) {
   #       matches <- whichv(vapply(x, FUN, TRUE, USE.NAMES = FALSE), FALSE)
   #       if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]])
   #     }
   #   }
   # } else {
   #   regsearch <- function(x) {
   #     if(nkeep_class && is.object(x)) oldClass(x) <- NULL
   #     if(any(subl <- .Call(C_vtypes, x, ret))) { # vapply(x, is.subl, TRUE, USE.NAMES = FALSE) # is.list(x) && a
   #       wsubl <- which(subl)
   #       wnsubl <- whichv(subl, FALSE)
   #       matches <- vapply(x[wnsubl], FUN, TRUE, USE.NAMES = FALSE)
   #       a <- lapply(x[wsubl], regsearch)
   #       wa <- vlengths(a, FALSE) > 0L # note that this also gets rid of null elements! could make it length or is.null! # vapply(a, length, 1L, USE.NAMES = FALSE)
   #       x <- c(x[wnsubl][matches], a[wa]) # The problem here: If all elements in a sublist are atomic, it still retains the sublist itself with NULL inside! -> but c() removes it!!
   #       if(keep.tree || length(x) != 1L)
   #         return(x[forder.int(c(wnsubl[matches], wsubl[wa]))]) else return(x[[1L]]) # fastest way?
   #     } else if(length(x)) { # This ensures correct behavior in the final nodes: if (length(x)) because problem encountered in get.elem(V, is.matrix) -> empty xlevels list, the lapply below does not execute
   #       matches <- which(vapply(x, FUN, TRUE, USE.NAMES = FALSE))
   #       if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]]) # needs to be !=
   #     }
   #   }
   # }
 }
 regsearch(l)
}

list_extract_regex <- function(l, exp, ret, keep.tree = FALSE, nkeep_class = TRUE, invert = FALSE, ...) {
  if(invert) {
    regsearch <- function(x) {
      if(nkeep_class && is.object(x)) oldClass(x) <- NULL
      if(any(subl <- .Call(C_vtypes, x, ret))) {
        matches <- if(is.null(names(x))) rep(TRUE, length(x)) else !rgrepl(exp, names(x), ...) # rgrep with invert??
        wsubl <- which(matches & subl)
        if(length(wsubl)) {
          wres <- which(matches & !subl)
          a <- lapply(x[wsubl], regsearch)
          wa <- vlengths(a, FALSE) > 0L
          x <- c(x[wres], a[wa])
          if(keep.tree || length(x) != 1L)
            return(x[forder.int(c(wres, wsubl[wa]))]) else return(x[[1L]])
        } else {
          wres <- which(matches)
          if(keep.tree || length(wres) != 1L) return(x[wres]) else return(x[[wres]])
        }
      } else {
        matches <- !rgrepl(exp, names(x), ...)
        if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]])
      }
    }
  } else {
    regsearch <- function(x) {
      if(nkeep_class && is.object(x)) oldClass(x) <- NULL
      if(any(subl <- .Call(C_vtypes, x, ret))) {
        matches <- rgrepl(exp, names(x), ...)
        wres <- which(matches)
        # wres <- rgrep(exp, names(x), ...)
        wnressubl <- which(if(length(wres)) subl & !matches else subl)
        # wnressubl <- if(length(wres)) fsetdiff(which(subl), wres) else which(subl)
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
  }
 regsearch(l)
}

list_extract_names <- function(l, nam, ret, keep.tree = FALSE, nkeep_class = TRUE, invert = FALSE) {
 if(invert) {
   regsearch <- function(x) {
     if(nkeep_class && is.object(x)) oldClass(x) <- NULL
     if(any(subl <- .Call(C_vtypes, x, ret))) {
       matches <- if(is.null(names(x))) rep(TRUE, length(x)) else names(x) %!in% nam
       wsubl <- which(matches & subl)
       if(length(wsubl)) {
         wres <- which(matches & !subl)
         a <- lapply(x[wsubl], regsearch)
         wa <- vlengths(a, FALSE) > 0L
         x <- c(x[wres], a[wa])
         if(keep.tree || length(x) != 1L)
           return(x[forder.int(c(wres, wsubl[wa]))]) else return(x[[1L]])
       } else {
         wres <- which(matches)
         if(keep.tree || length(wres) != 1L) return(x[wres]) else return(x[[wres]])
       }
     } else {
       matches <- which(names(x) %!in% nam)
       if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]])
     }
   }
 } else {
   regsearch <- function(x) {
     if(nkeep_class && is.object(x)) oldClass(x) <- NULL
     if(any(subl <- .Call(C_vtypes, x, ret))) {
       matches <- names(x) %in% nam
       wres <- which(matches) # match(nam, names(x), 0L) # better because gives integer(0) -> necessary as cannot do l[[0L]]
       wnressubl <- which(if(length(wres)) subl & !matches else subl) # fsetdiff(which(subl), wres)  # old solution: faster but does not work well if parent list is unnamed ! (i.e. l = list(lm1, lm1))
       if(length(wnressubl)) {
         a <- lapply(x[wnressubl], regsearch)
         wa <- vlengths(a, FALSE) > 0L # vapply(a, length, 1L)
         x <- c(x[wres], a[wa])
         if(keep.tree || length(x) != 1L)
           return(x[forder.int(c(wres, wnressubl[wa]))]) else return(x[[1L]])
       } else if(keep.tree || length(wres) != 1L) return(x[wres]) else return(x[[wres]])
     } else {
       matches <- which(names(x) %in% nam)
       if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]]) # needs to be !=, because integer(0) goes in first..
     }
   }
 }
 regsearch(l)
}

# Idea: Also use indices and logical vectors ? i.e. get first two columns of alist of data.frames ?
# This behaves a bit differently (not find elements everywhere, but also subset inside the list)
list_extract_ind <- function(l, ind, is.subl, keep.tree = FALSE, nkeep_class = TRUE) {
  if(is.logical(ind)) ind <- which(ind)
  if(length(ind) > 1L || keep.tree) {
    regsearch <- function(x) if(is.subl(x)) lapply(x, regsearch) else if(nkeep_class) .subset(x, ind) else x[ind]
  } else {
    # if(ind[1L] < 1L) stop("Cannot subset with single negative indices") # .subset2 throws error...
    regsearch <- function(x) if(is.subl(x)) lapply(x, regsearch) else .subset2(x, ind)
  }
  regsearch(l)
}

# Note: all functions currently remove empty list elements !
# keep.tree argument still issues with xlevels

get_elem <- function(l, elem, recursive = TRUE, DF.as.list = FALSE,
                     keep.tree = FALSE, keep.class = FALSE,
                     regex = FALSE, invert = FALSE, ...) {
  if(!is.list(l)) stop("l needs to be a list")
  if(recursive) {
    ret <- 4L - as.logical(DF.as.list)
    if(keep.class) al <- attributes(l)
    if(is.function(elem)) {
      l <- list_extract_FUN(l, elem, ret, keep.tree, !keep.class, invert, ...)
    } else if(is.character(elem)) {
      if(regex) {
        l <- list_extract_regex(l, elem, ret, keep.tree, !keep.class, invert, ...)
      } else {
        if(!missing(...)) unused_arg_action(match.call(), ...)
        l <- list_extract_names(l, elem, ret, keep.tree, !keep.class, invert)
      }
    } else {
      if(!missing(...)) unused_arg_action(match.call(), ...)
      if(invert) {
        if(is.numeric(elem)) stop("Cannot use option invert = TRUE if elem is indices")
        elem <- !elem
      }
      is.subl <- if(DF.as.list) is.list else function(x) is.list(x) && !inherits(x, "data.frame")
      l <- list_extract_ind(l, elem, is.subl, keep.tree, !keep.class)
    }
    if(keep.class && is.list(l)) {
      al[["names"]] <- names(l)
      return(setAttributes(l, al)) # class(l) <- cll # when drop.tree is proper, l might not be a list
    } else return(l)
  } else {
    if(is.function(elem)) {
      elem <- whichv(vapply(l, elem, TRUE, ..., USE.NAMES = FALSE), TRUE, invert)
    } else if(is.character(elem)) {
      if(regex) elem <- rgrep(elem, names(l), invert = invert, ...) else {
        if(!missing(...)) unused_arg_action(match.call(), ...)
        elem <- which(if(invert) names(l) %!in% elem else names(l) %in% elem)
      }
    } else if(is.logical(elem)) {
      if(!missing(...)) unused_arg_action(match.call(), ...)
      elem <- whichv(elem, TRUE, invert) # else stop("elem must be a function, character vector or vector of regular expressions!")
    }
    if(keep.tree || length(elem) != 1L) {
      if(keep.class) return(fcolsubset(l, elem)) else return(.subset(l, elem))
    } else return(.subset2(l, elem))
  }
}


# there is base::getElement

# 'regular' (is.atomic(x) || is.list(x)) elements, the check now implements in C_vtypes with option 5L.
is_regular_vec <- function(x) .Call(C_vtypes, x, 5L)
is_irregular_vec <- function(x) !.Call(C_vtypes, x, 5L)

# A variant of list_extract_FUN for FUN that can take a list as input and check the elements
list_extract_FUN_vec <- function(l, FUN, ret, keep.tree = FALSE, nkeep_class = TRUE) {
  regsearch <- function(x) {
    if(nkeep_class && is.object(x)) oldClass(x) <- NULL
    if(any(subl <- .Call(C_vtypes, x, ret))) {
      wsubl <- which(subl)
      wnsubl <- whichv(subl, FALSE)
      matches <- FUN(x[wnsubl])
      a <- lapply(x[wsubl], regsearch)
      wa <- vlengths(a, FALSE) > 0L
      x <- c(x[wnsubl][matches], a[wa])
      if(keep.tree || length(x) != 1L)
        return(x[forder.int(c(wnsubl[matches], wsubl[wa]))]) else return(x[[1L]])
    } else if(length(x)) {
      matches <- which(FUN(x))
      if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]])
    }
  }
  regsearch(l)
}

reg_elem <- function(l, recursive = TRUE, keep.tree = FALSE, keep.class = FALSE) {
  if(!is.list(l)) stop("l needs to be a list")
  if(keep.class) al <- attributes(l)
  # if(inherits(l, "data.frame")) if(keep.class) return(l) else return(unattrib(l))
  if(recursive) {
    l <- list_extract_FUN_vec(l, is_regular_vec, 4L, keep.tree, !keep.class)
    if(keep.class && is.list(l)) {
      al[["names"]] <- names(l)
      return(setAttributes(l, al))
    } else return(l)
  } else {
    matches <- which(is_regular_vec(l))
    if(keep.tree || length(matches) != 1L) {
      if(keep.class) return(fcolsubset(l, matches)) else return(.subset(l, matches))
    } else return(.subset2(l, matches))
  }
}

irreg_elem <- function(l, recursive = TRUE, keep.tree = FALSE, keep.class = FALSE) {
  if(!is.list(l)) stop("l needs to be a list")
  if(keep.class) al <- attributes(l)
  if(recursive) {
    l <- list_extract_FUN_vec(l, is_irregular_vec, 4L, keep.tree, !keep.class)
    if(keep.class && is.list(l)) {
      al[["names"]] <- names(l)
      return(setAttributes(l, al))
    } else return(l)
  } else {
    matches <- which(is_irregular_vec(l))
    if(keep.tree || length(matches) != 1L) {
      if(keep.class) return(fcolsubset(l, matches)) else return(.subset(l, matches))
    } else return(.subset2(l, matches))
  }
}

# TODO: See about big objects!
#microbenchmark(all(rapply(lm,is.atomic)),!is.list(unlist(lm, use.names = FALSE)),all(unlist(rapply2d(lm,is.std), use.names = FALSE)))
#microbenchmark(all(rapply(GGDC,is.atomic)),!is.list(unlist(GGDC, use.names = FALSE)),all(unlist(rapply2d(GGDC,is.std), use.names = FALSE)))

