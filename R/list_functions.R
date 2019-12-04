rapply2d <- function(l, FUN, ...) {
  aply2d <- function(y) if (is.list(y) && !inherits(y,  "data.frame")) lapply(y, aply2d) else FUN(y, ...) #is.null(dim(y)) # qsu output shows list of DF can have dim attr.
  #lapply(x,aply2d) # if this is enabled, rapply2d takes apart data.frame is passed
  aply2d(l)
}

list_elem <- function(l, return = c("sublist","names","indices","named_indices"), keep.class = FALSE) {
    switch(return[1L], sublist = if(keep.class) colsubset(l, is.list) else
              unclass(l)[vapply(l, is.list, TRUE, USE.NAMES = FALSE)],
           names = names(which(vapply(l, is.list, TRUE))),
           indices = which(vapply(l, is.list, TRUE, USE.NAMES = FALSE)),
           named_indices = which(vapply(l, is.list, TRUE)),
           stop("Unknown return option!"))
}
"list_elem<-" <- function(l, value) {
  al <- attributes(l)
  ilv <- is.list(value)
  len <- if(ilv) length(value) else 1L
  attributes(l) <- NULL # vapply without attributes is faster !!
  ind <- which(vapply(l, is.list, TRUE, USE.NAMES = FALSE))
  if(len != length(ind)) stop("length(value) must match length(list_elem(l))")
  if(ilv) l[ind] <- value else l[[ind]] <- value
  if(ilv && length(nam <- names(value)) == length(ind)) al[["names"]][ind] <- nam
  return(setAttributes(l, al))
}

atomic_elem <- function(l, return = c("sublist","names","indices","named_indices"), keep.class = FALSE) {
  switch(return[1L], sublist = if(keep.class) colsubset(l, is.atomic) else
    unclass(l)[vapply(l, is.atomic, TRUE, USE.NAMES = FALSE)],
    names = names(which(vapply(l, is.atomic, TRUE))),
    indices = which(vapply(l, is.atomic, TRUE, USE.NAMES = FALSE)),
    named_indices = which(vapply(l, is.atomic, TRUE)),
    stop("Unknown return option!"))
}
"atomic_elem<-" <- function(l, value) {
  al <- attributes(l)
  ilv <- is.list(value)
  len <- if(ilv) length(value) else 1L
  attributes(l) <- NULL
  ind <- which(vapply(l, is.atomic, TRUE, USE.NAMES = FALSE))
  if(len != length(ind)) stop("length(value) must match length(list_elem(l))")
  if(ilv) l[ind] <- value else l[[ind]] <- value
  if(ilv && length(nam <- names(value)) == length(ind)) al[["names"]][ind] <- nam
  return(setAttributes(l, al))
}

is.regular <- function(x) is.list(x) || is.atomic(x) # fastest way??

is.unlistable <- function(l) all(unlist(rapply2d(l, is.regular), use.names = FALSE)) # fastest way??

ldepth <- function(l, DF.as.list = TRUE) { # or list.depth !! (problem here is with list.elem -> I call the argument l, so ldepth is good) / depth / depthl  / l.depth Faster way?? what about data.frames??
  if (inherits(l,  "data.frame")) { # fast defining different functions in if-clause ??
    ld <- function(y,i) if (is.list(y)) lapply(y,ld,i+1L) else i
  } else if (DF.as.list) {
    ld <- function(y,i) {
      df <- inherits(y, "data.frame")
      if (is.list(y) && !df) lapply(y,ld,i+1L) else i+df
    }
  } else {
    ld <- function(y,i) if (is.list(y) && !inherits(y, "data.frame")) lapply(y,ld,i+1L) else i
  }
  return(max(unlist(ld(l, 0L), use.names = FALSE)))
} # If data.frame, search all, otherwise, make optional counting df or not, but don't search them.

has_elem <- function(l, elem, recursive = TRUE, DF.as.list = TRUE, regex = FALSE, ...) { # or data.frame.l ,search.data.frame, unpack.data.frame, as.list.data.frame
  if(is.function(elem)) {
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    if(recursive) {
     if(DF.as.list) return(any(unlist(rapply(l, elem, how = "list"), use.names = FALSE))) else
                    return(any(unlist(rapply2d(l, elem), use.names = FALSE)))
    } else return(any(vapply(l, elem, TRUE, USE.NAMES = FALSE)))
  } else if(is.character(elem)) {
    if(!regex && !missing(...)) stop("Unknown argument ", dotstostr(...))
    if(recursive) {
      class(l) <- NULL # in case [ behaves weird
      is.subl <- if(DF.as.list) is.list else function(x) is.list(x) && !inherits(x, "data.frame") # could do without, but it seems to remove data.frame attributes, and more speed!
      namply <- function(y) if(any(subl <- vapply(y, is.subl, TRUE)))
        c(names(subl), unlist(lapply(unclass(y)[subl], namply), use.names = FALSE)) else names(y) # also overall subl names are important, and unclass for DT subsetting !! !!! # names(which(!subl)) # names(y)[!subl] # which is faster !!
      if(regex) return(length(rgrep(elem, namply(l), ...)) > 0L) else return(any(namply(l) %in% elem))
    } else if(regex) return(length(rgrep(elem, names(l), ...)) > 0L) else return(any(names(l) %in% elem))
  } else stop("elem must be a function or character vector of element names or regular expressions")
}

# Experimental !!
# elem_names <- function(l, how = c("list", "unlist"), DF.as.list = TRUE) { # need right order for method how = list !!
#   namply <- function(y) if(any(subl <- vapply(y, is.subl, TRUE))) c(names(subl), lapply(unclass(y)[subl], namply)) else names(subl)
#   switch(how[1L],
#     unlist = names(rapply(l, function(x) NA)),
#     list =
#   ) rapply(l, function(x) NULL)
#
# }

# General note: What about lists containing data.tables ?? '[' subsetting will be wrong !!
list_extract_FUN <- function(l, FUN, is.subl, keep.tree = FALSE) {
 regsearch <- function(x) {
  if(any(subl <- vapply(x, is.subl, TRUE, USE.NAMES = FALSE))) { # is.list(x) && a
    wsubl <- which(subl)
    wnsubl <- which(!subl)
    matches <- vapply(x[wnsubl], FUN, TRUE, USE.NAMES = FALSE)
    a <- lapply(x[wsubl], regsearch)
    wa <- vapply(a, length, 1L, USE.NAMES = FALSE) > 0L # note that this also gets rid of null elements!! could make it length or is.null!!!
    x <- c(x[wnsubl][matches], a[wa]) # The problem here: If all elements in a sublist are atomic, it still retains the sublist itself with NULL inside!!
    if(keep.tree || length(x) != 1L)
      return(x[order(c(wnsubl[matches], wsubl[wa]))]) else return(x[[1L]]) # fastest way??
  } else if(length(x)) { # This ensures correct behavior in the final nodes: if (length(x)) because problem encountered in get.elem(V, is.matrix) -> empty xlevels list, the lapply below does not execute
    matches <- which(vapply(x, FUN, TRUE, USE.NAMES = FALSE))
    if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]]) # needs to be != !!!!
  }
 }
 return(regsearch(l))
}
list_extract_regex <- function(l, exp, is.subl, keep.tree = FALSE, ...) {
  regsearch <- function(x) {
  if(any(subl <- vapply(x, is.subl, TRUE, USE.NAMES = FALSE))) {
      wres <- rgrep(exp, names(x), ...)
      wnressubl <- setdiff(which(subl), wres)
    if(length(wnressubl)) { # faster way??
      a <- lapply(x[wnressubl], regsearch) # is this part still necessary??, or only for keep.tree
      wa <- vapply(a, length, 1L) > 0L # note that this also gets rid of null elements!! could make it length or is.null!!!
      x <- c(x[wres], a[wa])
      if(keep.tree || length(x) != 1L)
        return(x[order(c(wres, wnressubl[wa]))]) else return(x[[1L]])
    } else if(keep.tree || length(wres) != 1L) return(x[wres]) else return(x[[wres]])
  } else { # This ensures correct behavior in the final nodes:
    matches <- rgrep(exp, names(x), ...) # what if no matches here !!
    if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]]) # needs to be != !!!!
  }
 }
 return(regsearch(l))
}
list_extract_names <- function(l, nam, is.subl, keep.tree = FALSE) {
 regsearch <- function(x) {
  if(any(subl <- vapply(x, is.subl, TRUE, USE.NAMES = FALSE))) {
    wres <- match(nam, names(x))
    if(anyNA(wres)) wres <- integer(0)
    wnressubl <- setdiff(which(subl), wres)
      # matches <- names(x) %in% nam # old solution: faster but does not work well if parent list is unnamed !! (i.e. l = list(lm1, lm1))
      # wres <- which(matches)
      # wnressubl <- which(subl & !matches)
    if(length(wnressubl)) {
      a <- lapply(x[wnressubl], regsearch)
      wa <- vapply(a, length, 1L) > 0L
      x <- c(x[wres], a[wa])
      if(keep.tree || length(x) != 1L)
        return(x[order(c(wres, wnressubl[wa]))]) else return(x[[1L]])
    } else if(keep.tree || length(wres) != 1L) return(x[wres]) else return(x[[wres]])
  } else {
    matches <- match(nam, names(x))  # which(names(x) %in% nam)
    if(anyNA(matches)) matches <- integer(0)
    if(keep.tree || length(matches) != 1L) return(x[matches]) else return(x[[matches]]) # needs to be != !!!!
  }
 }
 return(regsearch(l))
}
# Idea: Also use indices and logical vectors ?? i.e. get first two columns of alist of data.frames ??
# This behaves a bit differently !! (not find elements everywhere, but also subset inside the list !!!)
list_extract_ind <- function(l, ind, is.subl, keep.tree = FALSE) {
  if(is.logical(ind)) ind <- which(ind)
  if(length(ind) > 1L || keep.tree) {
    regsearch <- function(x) if(is.subl(x)) lapply(x, regsearch) else x[ind]
  } else {
    regsearch <- function(x) if(is.subl(x)) lapply(x, regsearch) else x[[ind]]
  }
  return(regsearch(l))
}

# Note: all currently remove empty list elements !!
# keep.tree argument!! !! still issues wih xlevels

get_elem <- function(l, elem, recursive = TRUE, DF.as.list = TRUE,
                     keep.tree = FALSE, keep.class = FALSE, regex = FALSE, ...) {
  if(recursive) { # possibly if is.list(x) is redundant, because you check above!! -> nah, recursive??
    is.subl <- if(DF.as.list) is.list else function(x) is.list(x) && !inherits(x, "data.frame") # could do without, but it seems to remove data.frame attributes
    if(keep.class) al <- attributes(l) # cll <- class(l) # perhaps generalize to other attributes??
    if(is.function(elem)) {
      if(!missing(...)) stop("Unknown argument ", dotstostr(...))
      l <- list_extract_FUN(l, elem, is.subl, keep.tree)
    } else if(is.character(elem)) {
      if(regex) l <- list_extract_regex(l, elem, is.subl, keep.tree, ...) else {
         if(!missing(...)) stop("Unknown argument ", dotstostr(...))
         l <- list_extract_names(l, elem, is.subl, keep.tree)
      }
    } else {
      if(!missing(...)) stop("Unknown argument ", dotstostr(...))
      l <- list_extract_ind(l, elem, is.subl, keep.tree)
    }
    if(keep.class && is.list(l)) {
      al[["names"]] <- names(l)
      return(setAttributes(l, al)) # class(l) <- cll # when drop.tree is proper, l might not be a list
    } else return(l)
  } else {
    if(is.function(elem)) {
      if(!missing(...)) stop("Unknown argument ", dotstostr(...))
      elem <- which(vapply(l, elem, TRUE, USE.NAMES = FALSE))
    } else if(is.character(elem)) {
      if(regex) elem <- rgrep(elem, names(l), ...) else {
        if(!missing(...)) stop("Unknown argument ", dotstostr(...))
        elem <- which(names(l) %in% elem)
      }
    } else if(is.logical(elem)) {
      if(!missing(...)) stop("Unknown argument ", dotstostr(...))
      elem <- which(elem) # else stop("elem must be a function, character vector or vector of regular expressions!")
    }
    if(keep.tree || length(elem) != 1L) {
      if(keep.class) return(colsubset(l, elem)) else return(unclass(l)[elem]) # <- # base::Filter(elem, l)
    } else return(l[[elem]])
  }
}

# Neat example: qDF(get_elem(V,"residuals"))

# there is base::getElement, but still, to be consistent with vars call this get.elem(), and lsubset

reg_elem <- function(l, recursive = TRUE, keep.tree = FALSE, keep.class = FALSE) { # regular.elem, reg.elem, std.elem, data.elem, unl.elem?? data.table methods needed??, add recursive option!!!
  if(keep.class) al <- attributes(l)
  if(inherits(l, "data.frame")) return(l)
  if(recursive) {
    is.subl <- function(x) is.list(x) && !inherits(x, "data.frame")
    l <- list_extract_FUN(l, is.regular, is.subl, keep.tree)
    if(keep.class && is.list(l)) {
      al[["names"]] <- names(l)
      return(setAttributes(l, al))
    } else return(l)
  } else {
    matches <- which(vapply(l, is.regular, TRUE, USE.NAMES = FALSE)) # l <- base::Filter(is.regular,l)
    if(keep.tree || length(matches) != 1L) {
      if(keep.class) return(colsubset(l, matches)) else return(unclass(l)[matches])
    } else return(l[[matches]])
  }
}
irreg_elem <- function(l, recursive = TRUE, keep.tree = FALSE, keep.class = FALSE) { # irreg.elem or better nonstd.elem # add recursive option!!
  is.irregular <- function(x) !(is.list(x) || is.atomic(x)) # is.irregular fastest way??
  if(keep.class) al <- attributes(l)
  if(inherits(l, "data.frame")) stop("A data.frame is a regular object!")
  if(recursive) {
    is.subl <- function(x) is.list(x) && !inherits(x, "data.frame")
    l <- list_extract_FUN(l, is.irregular, is.subl, keep.tree)
    if(keep.class && is.list(l)) {
      al[["names"]] <- names(l)
      return(setAttributes(l, al))
    } else return(l)
  } else {
    matches <- which(vapply(l, is.irregular, TRUE, USE.NAMES = FALSE)) # l <- base::Filter(is.regular,l)
    if(keep.tree || length(matches) != 1L) {
      if(keep.class) return(colsubset(l, matches)) else return(unclass(l)[matches])
    } else return(l[[matches]])
  }
}

# See about big objects!!!
#microbenchmark(all(rapply(lm,is.atomic)),!is.list(unlist(lm, use.names = FALSE)),all(unlist(rapply2d(lm,is.std), use.names = FALSE)))
#microbenchmark(all(rapply(GGDC,is.atomic)),!is.list(unlist(GGDC, use.names = FALSE)),all(unlist(rapply2d(GGDC,is.std), use.names = FALSE)))


# Old versions: --------------------------------
# list.elem <- function(l) base::Filter(is.list,l)
# "list.elem<-" <- function(l, value) {
#   if (inherits(l, "data.table"))
#     l[,which(vapply(l, is.list, TRUE))] = value else
#       l[vapply(l, is.list, TRUE)] = value
#     l # right??
# }
# atomic.elem <- function(l) base::Filter(is.atomic,l) # Problem: identical(md,atomic.elem(md)) is not true!!
# "atomic.elem<-" <- function(l, value) {
#   if (inherits(l, "data.table"))
#     l[,which(vapply(l, is.atomic, TRUE))] = value else
#       l[vapply(l, is.atomic, TRUE)] = value
#     l
# } # https://stackoverflow.com/questions/25130531/how-to-select-only-numeric-columns-from-a-data-table
