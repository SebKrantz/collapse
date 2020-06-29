# Cuniqlengths <- data.table:::Cuniqlengths
# Cfrank <- data.table:::Cfrank
# forderv <- data.table:::forderv

GRP <- function(X, ...) UseMethod("GRP") # , X

radixorder <- function(..., na.last = TRUE, decreasing = FALSE, starts = FALSE, group.sizes = FALSE, sort = TRUE) {
  z <- pairlist(...)
  decreasing <- rep_len(as.logical(decreasing), length(z))
  .Call(C_radixsort, na.last, decreasing, starts, group.sizes, sort, z)
}

radixorderv <- function(x, na.last = TRUE, decreasing = FALSE, starts = FALSE, group.sizes = FALSE, sort = TRUE) {
  z <- if(is.atomic(x)) pairlist(x) else as.pairlist(unclass(x))
  decreasing <- rep_len(as.logical(decreasing), length(z))
  .Call(C_radixsort, na.last, decreasing, starts, group.sizes, sort, z)
}

# Added... could also do in GRP.default... but this is better, no match.call etc... match.call takes 4 microseconds. could do both ?? think about possible applications...
GRP.GRP <- function(x) x

GRP.default <- function(X, by = NULL, sort = TRUE, decreasing = FALSE, na.last = TRUE,
                        return.groups = TRUE, return.order = FALSE, call = TRUE, ...) { # , gs = TRUE # o

  if(!missing(...)) {
    args <- list(...)
    if(any(names(args) == "order")) { # all
      decreasing <- args[["order"]] == 1L # ... == 1L
      warning("'order' has been replaced with 'decreasing' and now takes logical arguments. 'order' can still be used but may be removed at some point.")
    } # else unused_arg_action(match.call(), ...)  could also be "group.sizes" ...
  }

  if(is.list(X)) {
    if(inherits(X, "GRP")) return(X) # keep ??
    if(is.null(by)) {
      by <- seq_along(unclass(X))
      namby <- attr(X, "names")
      if(is.null(namby)) attr(X, "names") <- namby <- paste0("Group.", by)
      o <- radixorderv(X, na.last, decreasing, TRUE, TRUE, sort)
    } else {
      if(is.call(by)) {
        namby <- all.vars(by)
        by <- ckmatch(namby, attr(X, "names"))
      } else if(is.character(by)) {
        namby <- by
        by <- ckmatch(by, attr(X, "names"))
      } else if(is.numeric(by)) {
        by <- as.integer(by)
        namby <- attr(X, "names")[by]
        if(is.null(namby)) {
          namby <- paste0("Group.", seq_along(by))
          attr(X, "names") <- paste0("Group.", seq_along(unclass(X))) # best ?
        }
      } else stop("by needs to be either a one-sided formula, character column names or column indices!")
      o <- radixorderv(.subset(X, by), na.last, decreasing, TRUE, TRUE, sort)
    }
  } else {
   if(!is.null(by)) stop("by can only be used to subset list / data.frame columns")
   namby <- l1orlst(as.character(substitute(X))) # paste(all.vars(call), collapse = ".") # good in all circumstances ?
   o <- radixorderv(X, na.last, decreasing, TRUE, TRUE, sort)
  }

  st <- attr(o, "starts")
  gs <- attr(o, "group.sizes")
  sorted <- attr(o, "sorted")

  if(return.groups) {
      ust <- if(sorted) st else o[st]
      groups <- if(is.list(X)) .Call(C_subsetDT, X, ust, by) else
        `names<-`(list(.Call(C_subsetVector, X, ust)), namby) # subsetVector preserves attributes (such as "label")
  } else groups <- NULL

  return(`oldClass<-`(list(N.groups = length(st),
                        group.id = .Call(C_frankds, o, st, gs, TRUE),
                        group.sizes = gs,
                        groups = groups,
                        group.vars = namby,
                        ordered = c(GRP.sort = sort, initially.ordered = sorted),
                        order = if(return.order) `attr<-`(o, "group.sizes", NULL) else NULL,
                        call = if(call) match.call() else NULL), "GRP"))
}

is.GRP <- function(x) inherits(x, "GRP")

GRPnames <- function(x, force.char = TRUE) { # , ...
  groups <- x[[4L]]
  if(is.null(groups)) return(NULL)
  if(length(unclass(groups)) > 1L) return(do.call(paste, c(groups, list(sep = "."))))
  if(force.char) tochar(.subset2(groups, 1L)) else .subset2(groups, 1L) # paste0(groups[[1L]]) prints "NA" but is slow, if assign with rownames<-, cannot have duplicate row names. But, attr<- "row.names" is fine !!
}

group_names.GRP <- GRPnames

print.GRP <- function(x, n = 6, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  ord <- x[[6L]]
  cat(paste("collapse grouping object of length",length(x[[2L]]),"with",
            x[[1L]],ifelse(any(ord),"ordered","unordered"),"groups"), fill = TRUE)
  cat("\nCall: ", paste0(deparse(x[[8L]]),", ",ifelse(ord[2L],"ordered","unordered")), "\n\n", sep = "")
  cat("Distribution of group sizes: ", fill = TRUE)
  print.summaryDefault(summary.default(x[[3L]]))
  if(!is.null(x[[4L]])) {
    ug <- unattrib(x[[4L]])
    cat("\nGroups with sizes: ", fill = TRUE)
    if(length(ug) == 1L) {
      ug <- ug[[1L]]
      if(length(ug) > 2L*n) {
        ind <- seq.int(x[[1L]]-n+1L, x[[1L]])
        print.default(setNames(x[[3L]][1:n], ug[1:n]))
        cat("  ---", fill = TRUE)
        print.default(setNames(x[[3L]][ind], ug[ind]))
      } else print.default(setNames(x[[3L]], ug))
    } else {
      if(length(ug[[1L]]) > 2L*n) {
        ind <- seq.int(x[[1L]]-n+1L, x[[1L]])
        print.default(setNames(x[[3L]][1:n], do.call(paste, c(lapply(ug, function(x)x[1:n]), list(sep = ".")))))
        cat("  ---", fill = TRUE)
        print.default(setNames(x[[3L]][ind], do.call(paste, c(lapply(ug, function(x)x[ind]), list(sep = ".")))))
      } else print.default(setNames(x[[3L]], do.call(paste, c(ug, list(sep = ".")))))
    }
  }
}

plot.GRP <- function(x, breaks = "auto", type = "s", horizontal = FALSE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  oldpar <- par(mfrow = if(horizontal) 1:2 else 2:1, mar = c(3.9,4.1,2.1,1), mgp = c(2.5,1,0))
  on.exit(par(oldpar))
  if(breaks == "auto") {
    ugs <- length(funique(x[[3L]]))
    breaks <- if(ugs > 80) 80 else ugs
  }
  plot(seq_len(x[[1L]]), x[[3L]], type = type, xlab = "Group id", ylab = "Group Size",
       main = paste0("Sizes of ",x[[1L]]," ",ifelse(any(x[[6L]]),"Ordered","Unordered")," Groups"), frame.plot = FALSE)
  if(breaks == 1L) plot(x[[3L]][1L], x[[1L]], type = "h", ylab = "Frequency", xlab = "Group Size",
                        main = "Histogram of Group Sizes", frame.plot = FALSE) else
  hist(x[[3L]], breaks, xlab = "Group Size", main = "Histogram of Group Sizes")
}

as.factor.GRP <- function(x, ordered = FALSE) { # , ...
  if(is.factor(x)) return(x)
  if(!is.GRP(x)) stop("x must be a 'GRP' object")
  f <- x[[2L]]
  gr <- unclass(x[[4L]])
  if(is.null(gr)) {
    attr(f, "levels") <- as.character(seq_len(x[[1L]]))
  } else {
    if(length(gr) == 1L) {
      attr(f, "levels") <- tochar(gr[[1L]]) # or formatC ?
    } else {
      attr(f, "levels") <- do.call(paste, c(gr, list(sep = ".")))
    }
  }
  oldClass(f) <- if(ordered) c("ordered","factor","na.included") else c("factor","na.included") # previously if any(x[[6L]])
  f
}

finteraction <- function(..., ordered = FALSE, sort = TRUE) { # does it drop levels ? -> Yes !
  if(...length() == 1L && is.list(...)) return(as.factor.GRP(GRP.default(..., sort = sort, call = FALSE), ordered))
  as.factor.GRP(GRP.default(list(...), sort = sort, call = FALSE), ordered)
}

GRP.qG <- function(X, ..., group.sizes = TRUE, return.groups = TRUE, call = TRUE) {
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  gvars <- l1orlst(as.character(substitute(X))) # paste(all.vars(call), collapse = ".") # good in all circumstances ?
  ng <- attr(X, "N.groups")
  grl <- return.groups && !is.null(groups <- attr(X, "groups"))
  if(!inherits(X, "na.included")) if(anyNA(unclass(X))) {
    ng <- ng + 1L
    X[is.na(X)] <- ng
    if(grl) groups <- c(groups, NA)
  }
  ordered <- if(is.ordered(X)) c(TRUE,TRUE) else c(FALSE,FALSE)
  attributes(X) <- NULL
  return(`oldClass<-`(list(N.groups = ng,
                        group.id = X,
                        group.sizes = if(group.sizes) tabulate(X, ng) else NULL, # .Internal(tabulate(X, ng))
                        groups = if(grl) `names<-`(list(groups), gvars) else NULL,
                        group.vars = gvars,
                        ordered = ordered,
                        order = NULL,
                        call = if(call) match.call() else NULL), "GRP"))
}

GRP.factor <- function(X, ..., group.sizes = TRUE, return.groups = TRUE, call = TRUE) {
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  nam <- l1orlst(as.character(substitute(X))) # paste(all.vars(call), collapse = ".") # good in all circumstances ?
  if(!inherits(X, "na.included")) X <- addNA2(X)
  lev <- attr(X, "levels")
  nl <- length(lev)
  ordered <- if(is.ordered(X)) c(TRUE, TRUE) else c(FALSE, FALSE)
  attributes(X) <- NULL
  return(`oldClass<-`(list(N.groups = nl,
                        group.id = X,
                        group.sizes = if(group.sizes) tabulate(X, nl) else NULL, # .Internal(tabulate(X, nl))
                        groups = if(return.groups) `names<-`(list(lev), nam) else NULL,
                        group.vars = nam,
                        ordered = ordered,
                        order = NULL,
                        call = if(call) match.call() else NULL), "GRP"))
}

GRP.pseries <- function(X, effect = 1L, ..., group.sizes = TRUE, return.groups = TRUE, call = TRUE) {
  g <- unclass(attr(X, "index")) # index cannot be atomic since plm always adds a time variable !
  if(length(effect) > 1L) return(GRP.default(g[effect], ...))
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  # if(length(g) > 2L) {
  #   mlg <- -length(g)
  #   nam <- paste(names(g)[mlg], collapse = ".")
  #   g <- interaction(g[mlg], drop = TRUE)
  # } else {
    nam <- names(g)[effect]
    g <- g[[effect]] # Fastest way to do this ?
  # }
  lev <- attr(g, "levels")
  nl <- length(lev)
  ordered <- if(is.ordered(g)) c(TRUE,TRUE) else c(FALSE,FALSE)
  attributes(g) <- NULL
  return(`oldClass<-`(list(N.groups = nl,
                        group.id = g,
                        group.sizes = if(group.sizes) tabulate(g, nl) else NULL, # .Internal(tabulate(g, nl))
                        groups = if(return.groups) `names<-`(list(lev), nam) else NULL,
                        group.vars = nam,
                        ordered = ordered,
                        order = NULL,
                        call = if(call) match.call() else NULL), "GRP"))
}
GRP.pdata.frame <- function(X, effect = 1L, ...) GRP.pseries(X, effect, ...)

fgroup_by <- function(X, ..., sort = TRUE, decreasing = FALSE, na.last = TRUE, return.order = FALSE) {      #   e <- substitute(list(...)) # faster but does not preserve attributes of unique groups !!
  attr(X, "groups") <- GRP.default(fselect(X, ...), NULL, sort, decreasing, na.last, TRUE, return.order, FALSE) # `names<-`(eval(e, X, parent.frame()), all.vars(e))
  add_cl <- c("tbl_df", "tbl", "grouped_df")
  oldClass(X) <- c(add_cl, fsetdiff(oldClass(X), add_cl)) # necesssary to avoid printing errors... (i.e. wrong group object etc...)
  X
}

fgroup_vars <- function(X, return = "data") {
  g <- attr(X, "groups")
  vars <- if(is.GRP(g)) g[[5L]] else attr(g, "names")[-length(unclass(g))]
  switch(return[1L],
    data = {
      ax <- attributes(X)
      ax[["class"]] <- ax[["class"]][ax[["class"]] != "grouped_df"]
      ind <- ckmatch(vars, ax[["names"]])
      ax[["names"]] <- vars
      setAttributes(.subset(X, ind), ax[names(ax) != "groups"])
    },
    unique = if(is.GRP(g)) g[[4L]] else fcolsubset(g, -length(unclass(g))), # what about attr(*, ".drop") ??
    names = vars,
    indices = ckmatch(vars, attr(X, "names")),
    named_indices = `names<-`(ckmatch(vars, attr(X, "names")), vars),
    logical = `[<-`(logical(length(unclass(X))), ckmatch(vars, attr(X, "names")), TRUE),
    named_logical = {
      nam <- attr(X, "names")
      `names<-`(`[<-`(logical(length(nam)), ckmatch(vars, nam), TRUE), nam)
    },
    stop("Unknown return option!"))
}

GRP.grouped_df <- function(X, ..., call = TRUE) {
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  # g <- unclass(attr(X, "groups"))
  g <- attr(X, "groups")
  if(is.GRP(g)) return(g)
  oldClass(g) <- NULL
  lg <- length(g)
  gr <- g[[lg]]
  ng <- length(gr)
  gs <- lengths(gr, FALSE)
  return(`oldClass<-`(list(N.groups = ng, # The cpp here speeds up things a lot !!
                        group.id = .Call(Cpp_groups2GRP, gr, fnrow2(X), gs),  # Old: rep(seq_len(ng), gs)[order(unlist(gr, FALSE, FALSE))], # .Internal(radixsort(TRUE, FALSE, FALSE, TRUE, .Internal(unlist(gr, FALSE, FALSE))))
                        group.sizes = gs,
                        groups = g[-lg], # better reclass afterwards ?
                        group.vars = names(g)[-lg],
                        ordered = c(TRUE, TRUE),
                        order = NULL,
                        call = if(call) match.call() else NULL), "GRP"))
}

is.qG <- function(x) inherits(x, "qG")

# TODO: what about NA last option ?
# TODO: More efficient remove missing values ??
radixfact <- function(x, sort, ord, fact, naincl, retgrp = FALSE) {
  o <- .Call(C_radixsort, TRUE, FALSE, fact || naincl || retgrp, naincl, sort, pairlist(x))
  st <- attr(o, "starts")
  f <- if(naincl) .Call(C_frankds, o, st, attr(o, "group.sizes"), TRUE) else # Fastest? -> Seems so..
        .Call(Cpp_groupid, x, o, 1L, TRUE, FALSE)
  if(fact) {
    duplattributes(f, x)
    if(naincl) {
      attr(f, "levels") <- if(attr(o, "sorted")) unattrib(tochar(.Call(C_subsetVector, x, st))) else
            unattrib(tochar(.Call(C_subsetVector, x, o[st]))) # use C_subsetvector ?
    } else {
      attr(f, "levels") <- if(attr(o, "sorted")) unattrib(tochar(na_rm(.Call(C_subsetVector, x, st)))) else
            unattrib(tochar(na_rm(.Call(C_subsetVector, x, o[st]))))
    }
    oldClass(f) <- c(if(ord) "ordered", "factor", if(naincl) "na.included")
  } else {
    if(naincl) attr(f, "N.groups") <- length(st) # the order is important, this before retgrp !!
    if(retgrp) {
      if(naincl) {
         attr(f, "groups") <- if(attr(o, "sorted")) .Call(C_subsetVector, x, st) else .Call(C_subsetVector, x, o[st])
      } else {
         attr(f, "groups") <- if(attr(o, "sorted")) na_rm(.Call(C_subsetVector, x, st)) else na_rm(.Call(C_subsetVector, x, o[st]))
      }
    }
    oldClass(f) <- c(if(ord) "ordered", "qG", if(naincl) "na.included")
  }
  f
}

qF <- function(x, ordered = FALSE, na.exclude = TRUE, sort = TRUE, method = c("auto", "radix", "hash")) {
  if(is.factor(x)) {
    if(na.exclude || inherits(x, "na.included")) {
      if(ordered && !is.ordered(x)) oldClass(x) <- c("ordered", oldClass(x)) # can set unordered ??
      return(x)
    }
    return(`oldClass<-`(addNA2(x), c(if(ordered) "ordered", "factor", "na.included")))
  } else if(is.qG(x)) {
    groups <- attr(x, "groups")
    groups <- if(is.null(groups)) as.character(seq_len(attr(x, "N.groups"))) else tochar(groups)
    nainc <- inherits(x, "na.included")
    if(na.exclude || nainc) {
      clx <- c(if(ordered) "ordered", "factor", if(nainc) "na.included") # can set unordered ??
    } else {
      if(anyNA(unclass(x))) {
        x[is.na(x)] <- attr(x, "N.groups") + 1L
        groups <- c(groups, NA_character_) # faster doing groups[length(groups)+1] <- NA? -> Nope, what you have is fastest !
      }
      clx <- c(if(ordered) "ordered", "factor", "na.included")
    }
    return(`attributes<-`(x, list(levels = groups, class = clx)))
  }
  switch(method[1L], # if((is.character(x) && !na.exclude) || (length(x) < 500 && !(is.character(x) && na.exclude)))
         auto  = if(is.character(x) || is.logical(x) || length(x) < 500L) .Call(Cpp_qF, x, sort, ordered, na.exclude) else
           radixfact(x, sort, ordered, TRUE, !na.exclude),
         radix = radixfact(x, sort, ordered, TRUE, !na.exclude),
         hash = .Call(Cpp_qF, x, sort, ordered, na.exclude),
         stop("Unknown method"))
}

# TODO: Keep if(ordered) "ordered" ?
qG <- function(x, ordered = FALSE, na.exclude = TRUE, sort = TRUE, return.groups = FALSE, method = c("auto", "radix", "hash")) {
  if(inherits(x, c("factor", "qG"))) {
    nainc <- inherits(x, "na.included")
    if(na.exclude || nainc || !anyNA(unclass(x))) {
      newclx <- c(if(ordered) "ordered", "qG", if(nainc || !na.exclude) "na.included")
      if(is.factor(x)) {
        ax <- if(return.groups) list(N.groups = fnlevels(x), groups = attr(x, "levels"), class = newclx) else
          list(N.groups = fnlevels(x), class = newclx)
      } else {
        ax <- if(return.groups) list(N.groups = attr(x, "N.groups"), groups = attr(x, "groups"), class = newclx) else
          list(N.groups = attr(x, "N.groups"), class = newclx)
      }
      return(`attributes<-`(x, ax))
    }
    newclx <- c(if(ordered) "ordered", "qG", "na.included")
    if(is.factor(x)) {
      lev <- attr(x, "levels")
      if(anyNA(lev)) ng <- length(lev) else {
        ng <- length(lev) + 1L
        if(return.groups) lev <- c(lev, NA_character_)
      }
      attributes(x) <- NULL # factor method seems faster, however cannot assign integer, must assign factor level...
    } else {
      if(return.groups && !is.null(lev <- attr(x, "groups"))) lev <- c(lev, NA)
      ng <- attr(x, "N.groups") + 1L
    }
    ax <- if(return.groups) list(N.groups = ng, groups = lev, class = newclx) else
      list(N.groups = ng, class = newclx)
    x[is.na(x)] <- ng
    return(`attributes<-`(x, ax))
  }
  switch(method[1L], # if((is.character(x) && !na.exclude) || (length(x) < 500 && !(is.character(x) && na.exclude)))
         auto  = if(is.character(x) || is.logical(x) || length(x) < 500L) .Call(Cpp_qG, x, sort, ordered, na.exclude, return.groups) else
           radixfact(x, sort, ordered, FALSE, !na.exclude, return.groups),
         radix = radixfact(x, sort, ordered, FALSE, !na.exclude, return.groups),
         hash =  .Call(Cpp_qG, x, sort, ordered, na.exclude, return.groups),
         stop("Unknown method"))
}


radixuniquevec <- function(x, sort) {
  o <- .Call(C_radixsort, TRUE, FALSE, TRUE, FALSE, sort, pairlist(x))
  if(attr(o, "sorted")) .Call(C_subsetVector, x, attr(o, "starts")) else
    .Call(C_subsetVector, x, o[attr(o, "starts")])
}

funique <- function(x, ...) UseMethod("funique")

funique.default <- function(x, sort = TRUE, method = c("auto", "radix", "hash"), ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.array(x)) stop("funique currently only supports atomic vectors and data.frames")
  switch(method[1L],
         auto = if(is.character(x) || is.logical(x) || length(x) < 500L)
           .Call(Cpp_funique, x, sort) else radixuniquevec(x, sort),
         radix = radixuniquevec(x, sort),
         hash = .Call(Cpp_funique, x, sort))
}

# TODO: could make faster still... not using colsubset but something more simple... no attributes needed...
# TODO: Enable by formula use ?? by or cols ?? -> cols is clearer !! also with na_omit, by could imply by-group uniqueness check...
funique.data.frame <- function(x, cols = NULL, sort = TRUE, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  o <- if(is.null(cols)) radixorderv(x, starts = TRUE, sort = sort) else
       radixorderv(colsubset(x, cols), starts = TRUE, sort = sort) # if(is.call(by)) .subset(x, ckmatch(attr(x, "names"), all.vars(by)))
  rn <- attr(x, "row.names")
  if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") {
     if(attr(o, "sorted")) return(.Call(C_subsetDT, x, attr(o, "starts"), seq_along(unclass(x))))
     return(.Call(C_subsetDT, x, o[attr(o, "starts")], seq_along(unclass(x))))
  }
  st <- if(attr(o, "sorted")) attr(o, "starts") else o[attr(o, "starts")]
  res <- .Call(C_subsetDT, x, st, seq_along(unclass(x)))
  attr(res, "row.names") <- rn[st]
  res
}

funique.list <- function(x, cols = NULL, sort = TRUE, ...) funique.data.frame(x, cols, sort, ...)
