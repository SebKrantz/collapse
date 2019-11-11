rapply2d <- function(X, FUN, ...) {
  aply2d <- function(y) if (is.list(y) && !is.data.frame(y)) lapply(y, aply2d) else FUN(y, ...) #is.null(dim(y)) # qsu output shows list of DF can have dim attr.
  #lapply(x,aply2d) # if this is enabled, rapply2d takes apart data.frame is passed
  aply2d(X)
}

# Use _ instead of . (because of classes) or no gap at all ??
# _ is good !! looks like dplyr and good practice !! (maybe there will be a class 'vars' one day )

list_elem <- function(l, return = c("sublist","names","indices","named_indices"), keep.class = FALSE) {
    switch(return[1L], sublist = if(keep.class) colsubset(l, is.list) else
              unclass(l)[vapply(l, is.list, TRUE, USE.NAMES = FALSE)],
           names = names(which(vapply(l, is.list, TRUE))),
           indices = which(vapply(l, is.list, TRUE, USE.NAMES = FALSE)),
           named.indices = which(vapply(l, is.list, TRUE)),
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
    named.indices = which(vapply(l, is.atomic, TRUE)),
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
  if (is.data.frame(l)) { # fast defining different functions in if-clause ??
    ld <- function(y,i) if (is.list(y)) lapply(y,ld,i+1L) else i
  } else if (DF.as.list) {
    ld <- function(y,i) {
      df <- is.data.frame(y)
      if (is.list(y) && !df) lapply(y,ld,i+1L) else i+df
    }
  } else {
    ld <- function(y,i) if (is.list(y) && !is.data.frame(y)) lapply(y,ld,i+1L) else i
  }
  return(max(unlist(ld(l, 0L), use.names = FALSE)))
} # If data.frame, search all, otherwise, make optional counting df or not, but don't search them.
has_elem <- function(l, FoR, recursive = TRUE, DF.as.list = TRUE, regex = FALSE, ...) { # or data.frame.l ,search.data.frame, unpack.data.frame, as.list.data.frame
 # if(!is.list(l)) stop("l is not a list")
  if(is.function(FoR)) {
    if(recursive) {
     if(DF.as.list) return(any(unlist(rapply(l, FoR, how = "list"), use.names = FALSE))) else
                    return(any(unlist(rapply2d(l, FoR), use.names = FALSE)))
    } else return(any(vapply(l, FoR, TRUE, USE.NAMES = FALSE)))
  } else if(is.character(FoR)) {
    if(recursive) {
      class(l) <- NULL # in case [ behaves weird
      is.subl <- if(DF.as.list) is.list else function(x) is.list(x) && !is.data.frame(x) # could do without, but it seems to remove data.frame attributes
      namply <- function(y) if(any(subl <- vapply(y, is.subl, TRUE)))
        c(names(y)[!subl], unlist(lapply(y[subl], namply), use.names = FALSE)) else names(y)
      if(regex) length(rgrep(FoR, namply(l), ...)) > 0L else any(namply(l) %in% FoR)
    } else if(regex) length(rgrep(FoR, names(l), ...)) > 0L else any(names(l) %in% FoR)
  } else stop("FoR must be a function or character vector of element names or regular expressions")
}

get_elem <-function(l, FoR, recursive = TRUE, DF.as.list = TRUE,
                    keep.tree = FALSE, keep.class = FALSE, regex = FALSE, ...) { # FUNorl See if this is implemented in some other package already!!
  # if (!is.list(l)) stop("l is not a list")
  if (DF.as.list) is.subl <- function(x) is.list(x) else is.subl <- function(x) is.list(x) && !is.data.frame(x) # could do without, but it seems to remove data.frame attributes
  if (keep.class) cll <- class(l) # perhaps generalize to other attributes??
  if (is.function(FoR)) {
    if (recursive) { # possibly if is.list(x) is redundant, because you check above!! -> nah, recursive??
      regsearch <- function(x) if (any(subl <- vapply(x, is.subl, TRUE))) { # is.list(x) && a
        wsubl <- which(subl)
        wnsubl <- which(!subl)
        matches <- vapply(x[wnsubl], FoR, TRUE)
        a <- lapply(x[wsubl], regsearch)
        wa <- vapply(a, length, 1L)>0 # note that this also gets rid of null elements!! could make it length or is.null!!!
        x <- c(x[wnsubl][matches], a[wa]) # The problem here: If all elements in a sublist are atomic, it still retains the sublist itself with NULL inside!!
        if (keep.tree || length(x)!=1)
        x[order(c(wnsubl[matches], wsubl[wa]))] else x[[1]] # fastest way??
      } else if (length(x)) { # This ensures correct behavior in the final nodes:
        # if (length(x)) because problem encountered in get.elem(V, is.matrix) -> empty xlevels list, the lapply below does not execute
        matches <- which(vapply(x, FoR, TRUE))
        if (keep.tree || length(matches)!=1) x[matches] else x[[matches]] # needs to be != !!!!
      }
      l <- regsearch(l)
    } else l <- base::Filter(FoR, l)
  } else if (is.character(FoR)) {
    if (recursive) { # possibly if is.list(x) is redundant is.list(x) && any(
      regsearch <- function(x) if (any(subl <- vapply(x, is.subl, TRUE))) {
        if (regex) {
          wres <- rgrep(FoR, names(x), ...)
          wnressubl <- setdiff(which(subl), wres)
        } else {
          matches <- names(x) %in% FoR
          wres <- which(matches)
          wnressubl <- which(subl & !matches) # faster way?? could use setdiff
        }
        if (length(wnressubl)) { # faster way??
          a <- lapply(x[wnressubl], regsearch) # is this part still necessary??, or only for keep.tree
          wa <- vapply(a, length, 1L)>0 # note that this also gets rid of null elements!! could make it length or is.null!!!
          x <- c(x[wres], a[wa])
          if (keep.tree || length(x)!=1)
          x[order(c(wres, wnressubl[wa]))] else x[[1]]
        } else if (keep.tree || length(wres)!=1) x[wres] else x[[wres]]
      } else { # This ensures correct behavior in the final nodes:
        matches <- if(regex) rgrep(FoR, names(x), ...) else which(names(x) %in% FoR)
        if (keep.tree || length(matches)!=1) x[matches] else x[[matches]] # needs to be != !!!!
      }
      l <- regsearch(l)
    } else {
      matches <- if(regex) rgrep(FoR, names(l), ...) else which(names(l) %in% FoR)
      l <- if (keep.tree || length(matches)!=1) l[matches] else l[[matches]]
    }
  } else stop("FoR must be a function or regular expression")
  if (keep.class && is.list(l)) class(l) <- cll # when drop.tree is proper, l might not be a list
  return(l)
} # Does not delete parent tree yet!, make option keep.tree!! also make regex option

# Neat example: plot.ts(do.call(cbind,get.elem(V,"residuals")))
# Better: quickDF(get.elem(V,"residuals"))
# make qDF() and qDT()

# This works well with functions:
  # regseach <- function(x) if (is.list(x) && any(subl <- unlist(lapply(x,is.subl), use.names = FALSE))) {
  #   wsubl <- which(subl)
  #   wnsubl <- which(!subl)
  #   matches <- unlist(lapply(x[wnsubl], FoR), use.names = FALSE)
  #   #x = c(x[wnsubl][matches],lapply(x[wsubl],regseach))
  #   #x[order(c(wnsubl[matches],wsubl))]
  #   a = lapply(x[wsubl],regseach)
  #   wa = unlist(lapply(a, length), use.names = FALSE)>0 # note that this also gets rid of null elements!! could make it length or is.null!!!
  #   x = c(x[wnsubl][matches],a[wa]) # The problem here: If all elements in a sublist are atomic, it still retains the sublist itself with NULL inside!!
  #   x[order(c(wnsubl[matches],wsubl[wa]))]
  # } else x[unlist(lapply(x, FoR), use.names = FALSE)]


# there is base::getElement, but still, to be consistent with vars call this get.elem(), and lsubset
# is only for subsetting a list with another list?
reg_elem <- function(l, recursive = TRUE, keep.tree = FALSE, keep.class = FALSE) { # regular.elem, reg.elem, std.elem, data.elem, unl.elem?? data.table methods needed??, add recursive option!!!
  if(is.list(l)) {
    if(keep.class) cll <- class(l) # perhaps generalize to other attributes??
    if(is.data.frame(l)) return(l)
    if(recursive) {
      is.subl <- function(x) is.list(x) && !is.data.frame(x) # could do without, but it seems to remove data.frame attributes
      is.atordf <- function(x) is.atomic(x) || is.data.frame(x) # generally see waht happens with attributes!!
      regsearch <- function(x) if (any(subl <- vapply(x, is.subl, TRUE))) {
        wsubl <- which(subl)
        wnsubl <- which(!subl)
        atom <- vapply(x[wnsubl], is.atordf, TRUE)
        a <- lapply(x[wsubl], regsearch)
        wa <- vapply(a, length, 1L)>0 # note that this also gets rid of null elements!! could make it length or is.null!!!
        x <- c(x[wnsubl][atom], a[wa]) # The problem here: If all elements in a sublist are atomic, it still retains the sublist itself with NULL inside!!
        if (keep.tree || length(x)!=1)
         x[order(c(wnsubl[atom], wsubl[wa]))] else x[[1]] # fastest way??
    } else if (length(x)) { # This ensures correct behavior in the final nodes:
      # if (length(x)) because problem encountered in get.elem(V, is.matrix) -> empty xlevels list, the lapply below does not execute
      atom <- which(vapply(x, is.atordf, TRUE))
      if (keep.tree || length(atom)!=1) x[atom] else x[[atom]] # needs to be != !!!!
    }
      l <- regsearch(l)
    } else l <- base::Filter(is.regular,l)
    if (keep.class && is.list(l)) class(l) <- cll
    return(l)
  } else if (is.regular(l)) return(l)
} # keep.tree argument!!
# "reg.elem<-" <- function(X, value) { # need data.table options for those?? can a data.table even contain non.atomic stuff??
#   if (is.list(X)) {
#     if (any(class(X) == "data.table"))
#       X[,which(unlist(lapply(X,is.regular), use.names = FALSE, recursive = FALSE))] = value else
#       X[unlist(lapply(X,is.regular), use.names = FALSE, recursive = FALSE)] = value
#   } else if (is.regular(X)) X = value
#   X # right??
# }
irreg_elem <- function(l, recursive = TRUE, keep.tree = FALSE, keep.class = FALSE) { # irreg.elem or better nonstd.elem # add recursive option!!
  is.irregular <- function(x) !(is.list(x) || is.atomic(x)) # is.irregular fastest way??
  if (is.list(l) && !is.data.frame(l)) {
    if (keep.class) cll <- class(l) # perhaps generalize to other attributes??
    if (recursive) {
      is.subl <- function(x) is.list(x) && !is.data.frame(x) # could do without, but it seems to remove data.frame attributes
      nis.atordf <- function(x) !(is.atomic(x) || is.list(x)) # generally see waht happens with attributes!!
      irregsearch <- function(x) if (any(subl <- vapply(x, is.subl, TRUE))) {
        wsubl <- which(subl)
        wnsubl <- which(!subl)
        natom <- vapply(x[wnsubl], nis.atordf, TRUE)
        a <- lapply(x[wsubl], irregsearch)
        wa <- vapply(a, length, 1L)>0
        x <- c(x[wnsubl][natom], a[wa]) # The problem here: If all elements in a sublist are atomic, it still retains the sublist itself with NULL inside!!
        if (keep.tree || length(x)!=1)
          x[order(c(wnsubl[natom], wsubl[wa]))] else x[[1]] # fastest way??
      } else if (length(x)) { # This ensures correct behavior in the final nodes:
        # if (length(x)) because problem encountered in get.elem(V, is.matrix) -> empty xlevels list, the lapply below does not execute
        natom <- which(vapply(x, nis.atordf, TRUE))
        if (keep.tree || length(natom)!=1) x[natom] else x[[natom]] # needs to be != !!!!
      }
      l = irregsearch(l)
    } else l = base::Filter(is.irregular,l)
    if (keep.class && is.list(l)) class(l) <- cll
    return(l)
  } else if (is.irregular(l)) return(l)
} # keep.tree argument!! !! still issues wih xlevels and qr in lm object!!
# "irreg.elem<-" <- function(X, value) {
#   is.irregular <- function(x) !(is.list(x) || is.atomic(x)) # fastest way??
#   if (is.list(X)) {
#     if (any(class(X) == "data.table"))
#       X[,which(unlist(lapply(X,is.irregular), use.names = FALSE, recursive = FALSE))] = value else
#       X[unlist(lapply(X,is.irregular), use.names = FALSE, recursive = FALSE)] = value
#   } else if (is.irregular(X)) X = value
#   X # right??
# }
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
