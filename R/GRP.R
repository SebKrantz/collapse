# Cuniqlengths <- data.table:::Cuniqlengths
# Cfrank <- data.table:::Cfrank
# forderv <- data.table:::forderv

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

switchGRP <- function(x, na.last = TRUE, decreasing = FALSE, starts = FALSE,
                      group.sizes = FALSE, sort = TRUE, use.group = FALSE) {
  if(use.group) return(.Call(C_group, x, starts, group.sizes))
  z <- if(is.atomic(x)) pairlist(x) else as.pairlist(unclass(x))
  decreasing <- rep_len(as.logical(decreasing), length(z))
  .Call(C_radixsort, na.last, decreasing, starts, group.sizes, sort, z)
}

group <- function(x, starts = FALSE, group.sizes = FALSE) {
  g <- .Call(C_group, x, starts, group.sizes)
  oldClass(g) <- c("qG", "na.included")
  g
}

gsplit <- function(x = NULL, g, use.g.names = FALSE, ...) {
  if(!(is.list(g) && inherits(g, "GRP"))) g <- GRP(g, return.groups = use.g.names, call = FALSE, ...)
  res <- if(is.null(x)) .Call(C_gsplit, 1L, g, TRUE) else if(length(unclass(x)) == length(g[[2L]]))
    .Call(C_gsplit, x, g, FALSE) else if(is.object(x))
      lapply(.Call(C_gsplit, 1L, g, TRUE), function(i) x[i]) else
        stop("length(x) must match length(g)")
  if(use.g.names) names(res) <- GRPnames(g, FALSE)
  res
}

greorder <- function(x, g, ...) {
  if(!(is.list(g) && inherits(g, "GRP"))) g <- GRP(g, return.groups = FALSE, call = FALSE, ...)
  .Call(C_greorder, x, g)
}

G_guo <- function(g) {
  if(is.atomic(g)) {
    if(inherits(g, c("factor", "qG"))) {
      if(inherits(g, "na.included") || !anyNA(unclass(g)))
        return(list(if(is.factor(g)) fnlevels(g) else attr(g, "N.groups"), unattrib(g), NULL))
      if(is.factor(g)) {
        ng <- if(anyNA(lev <- attr(g, "levels"))) length(lev) else length(lev) + 1L
      } else ng <- attr(g, "N.groups") + 1L
      return(list(ng, copyv(unattrib(g), NA_integer_, ng), NULL))
    }
    g <- .Call(C_group, g, FALSE, FALSE)
    return(list(attr(g,"N.groups"), g, NULL))
  }
  if(inherits(g, "GRP")) return(g)
  g <- .Call(C_group, g, FALSE, FALSE)
  return(list(attr(g,"N.groups"), g, NULL))
}

G_t <- function(x) {
  if(is.null(x)) return(x)
  # If integer time variable contains NA, does not break C++ code..
  if(is.atomic(x)) {
    if(is.object(x)) {
      if(inherits(x, c("factor", "qG"))) return(x)
      if(is.numeric(unclass(x))) return(timeid(x, factor = FALSE))
    } else if(is.numeric(x)) {
      # if(is.double(x)) message("Time variable is of type double, but not a date/time object. It is therefore coerced to integer and assumed to represent unitary timesteps. If this is not desired pass timeid(t). To silence this message pass as.integer(t).")
      return(as.integer(x))
    }
    return(qG(x, na.exclude = FALSE, sort = TRUE, method = "hash")) # make sure it is sorted !
  }
  # if(is_GRP(x)) return(x[[2L]]) # Not necessary because GRP.default also returns it..
  return(GRP.default(x, return.groups = FALSE, return.order = FALSE, sort = TRUE, call = FALSE)[[2L]])
}


GRP <- function(X, ...) UseMethod("GRP") # , X

# Added... could also do in GRP.default... but this is better, no match.call etc... match.call takes 4 microseconds. could do both ?? think about possible applications...
GRP.GRP <- function(X, ...) X

GRP.default <- function(X, by = NULL, sort = .op[["sort"]], decreasing = FALSE, na.last = TRUE,
                        return.groups = TRUE, return.order = sort, method = "auto",
                        call = TRUE, ...) {

  use.group <- switch(method, auto = !sort, hash = TRUE, radix = FALSE, stop("method needs to be 'auto', 'hash' or 'radix'."))

  if(is.na(na.last)) stop("here na.last needs to be TRUE or FALSE, otherwise the GRP object does not match the data dimensions.")

  if(is.list(X)) {
    if(inherits(X, "GRP")) return(X)
    if(is.null(by)) {
      by <- seq_along(unclass(X))
      namby <- attr(X, "names")
      if(is.null(namby)) attr(X, "names") <- namby <- paste0("Group.", by)
      o <- switchGRP(X, na.last, decreasing, return.groups || !use.group, TRUE, sort, use.group)
    } else {
      if(is.call(by)) {
        namby <- all.vars(by, unique = FALSE)
        by <- ckmatch(namby, attr(X, "names"))
      } else if(is.character(by)) {
        namby <- by
        by <- ckmatch(by, attr(X, "names"))
      } else {
        by <- if(is.numeric(by)) as.integer(by) else if(is.logical(by)) which(by) else if(is.function(by)) which(vapply(unattrib(X), by, TRUE)) else
              stop("by needs to be either a one-sided formula, character column names, column indices, a logical vector or selector function!")
        namby <- attr(X, "names")[by]
        if(is.null(namby)) {
          namby <- paste0("Group.", seq_along(by))
          attr(X, "names") <- paste0("Group.", seq_along(unclass(X))) # best ?
        }
      }
      o <- switchGRP(.subset(X, by), na.last, decreasing, return.groups || !use.group, TRUE, sort, use.group)
    }
  } else {
   if(length(by)) stop("by can only be used to subset list / data.frame columns")
   namby <- l1orlst(as.character(substitute(X))) # paste(all.vars(call), collapse = ".") # good in all circumstances ?
   o <- switchGRP(X, na.last, decreasing, return.groups || !use.group, TRUE, sort, use.group)
  }

  st <- attr(o, "starts")
  gs <- attr(o, "group.sizes")
  sorted <- if(use.group) NA else attr(o, "sorted")
  if(return.order && !use.group) ao <- attributes(o)[-2L]
  attributes(o) <- NULL

  if(return.groups) {
      # if unit groups, don't subset rows...
      if(length(gs) == length(o) && (use.group || sorted)) {
        ust <- st
        groups <- if(is.list(X)) .Call(C_subsetCols, X, by, FALSE) else `names<-`(list(X), namby)
      } else {
        ust <- if(use.group || sorted) st else if(length(gs) == length(o)) o else .Call(C_subsetVector, o, st, FALSE) # o[st]
        groups <- if(is.list(X)) .Call(C_subsetDT, X, ust, by, FALSE) else
          `names<-`(list(.Call(C_subsetVector, X, ust, FALSE)), namby) # subsetVector preserves attributes (such as "label")
      }
  } else {
    groups <- NULL
    ust <- NULL
  }

  return(`oldClass<-`(list(N.groups = length(gs),
                        group.id = if(use.group) o else .Call(C_frankds, o, st, gs, sorted),
                        group.sizes = gs,
                        groups = groups,
                        group.vars = namby,
                        ordered = c(ordered = sort, sorted = sorted),
                        order = if(return.order && !use.group) `attributes<-`(o, ao) else NULL, # `attributes<-`(o, attributes(o)[-2L]) This does a shallow copy on newer R versions # `attr<-`(o, "group.sizes", NULL): This deep-copies it..
                        group.starts = ust, # Does not need to be computed by group()
                        call = if(call) match.call() else NULL), "GRP"))
}

is_GRP <- function(x) inherits(x, "GRP")
is.GRP <- function(x) {
  .Deprecated(msg = "'is.GRP' was renamed to 'is_GRP'. It will be removed end of 2023, see help('collapse-renamed').")
  inherits(x, "GRP")
}

length.GRP <- function(x) length(x[[2L]])

GRPnames <- function(x, force.char = TRUE, sep = ".") { # , ...
  groups <- x[[4L]]
  if(is.null(groups)) return(NULL)
  if(length(unclass(groups)) > 1L) return(do.call(paste, c(groups, list(sep = sep))))
  if(force.char) tochar(.subset2(groups, 1L)) else .subset2(groups, 1L) # paste0(groups[[1L]]) prints "NA" but is slow, if assign with rownames<-, cannot have duplicate row names. But, attr<- "row.names" is fine !!
}

GRPid <- function(x, sort = FALSE, ...) {
  if(!missing(...) && any(names(dots <- list(...)) == "g")) {
    g <- dots$g
    if(!inherits(g, "GRP")) stop("g must be a 'GRP' object")
    res <- g$group.id
    if(!missing(x) && is.list(x)) return(lapply(x, function(y) res))
    return(res)
  }
  return(GRP(x, sort = sort, return.groups = FALSE, return.order = FALSE, call = FALSE, ...)$group.id)
}

GRPN <- function(x, expand = TRUE, ...) {
  if(!missing(...) && any(names(dots <- list(...)) == "g")) {
    g <- dots$g
    if(!inherits(g, "GRP")) stop("g must be a 'GRP' object")
    res <- if(any(names(dots) == "TRA")) .Call(C_subsetVector, g$group.sizes, g$group.id, FALSE) else g$group.sizes
    if(!missing(x) && is.list(x)) return(lapply(x, function(y) res))
    return(res)
  }
  g <- GRP(x, sort = FALSE, return.groups = FALSE, return.order = FALSE, call = FALSE, ...)
  if(expand) .Call(C_subsetVector, g$group.sizes, g$group.id, FALSE) else g$group.sizes
}

# dplyr-style n(): only for masking if collapse_mask = "all"
n_internal <- function(x, g, TRA, ...) {
  if(missing(g)) {
    if(missing(x)) stop("if data is not grouped need to call n() on a column")
    return(if(is.list(x)) fnrow(x) else length(x))
  }
  if(!inherits(g, "GRP")) stop("g must be a 'GRP' object")
  if(missing(TRA)) return(g$group.sizes)
  .Call(C_subsetVector, g$group.sizes, g$group.id, FALSE)
}
# group_names.GRP <- function(x, force.char = TRUE) {
#   .Deprecated("GRPnames")
#   GRPnames(x, force.char)
# }

print.GRP <- function(x, n = 6, ...) {
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  ord <- x[[6L]]
  cat(paste0("collapse grouping object of length ", length(x[[2L]]), " with ",
            x[[1L]], if(isTRUE(any(ord))) " ordered" else if(anyNA(ord)) "" else " unordered", " groups"), fill = TRUE)
  cat("\nCall: ", paste0(deparse(x[["call"]]), if(is.na(ord[2L])) "" else if(ord[2L]) ", X is sorted" else ", X is unsorted"), "\n\n", sep = "")
  cat("Distribution of group sizes: ", fill = TRUE)
  print.summaryDefault(summary.default(x[[3L]]), ...)
  if(!is.null(x[[4L]])) {
    ug <- unattrib(x[[4L]])
    cat("\nGroups with sizes: ", fill = TRUE)
    if(length(ug) == 1L) {
      ug <- ug[[1L]]
      if(length(ug) > 2L*n) {
        ind <- seq.int(x[[1L]]-n+1L, x[[1L]])
        print.default(setNames(x[[3L]][1:n], ug[1:n]), ...)
        cat("  ---", fill = TRUE)
        print.default(setNames(x[[3L]][ind], ug[ind]), ...)
      } else print.default(setNames(x[[3L]], ug), ...)
    } else {
      if(length(ug[[1L]]) > 2L*n) {
        ind <- seq.int(x[[1L]]-n+1L, x[[1L]])
        print.default(setNames(x[[3L]][1:n], do.call(paste, c(lapply(ug, function(x) x[1:n]), list(sep = ".")))), ...)
        cat("  ---", fill = TRUE)
        print.default(setNames(x[[3L]][ind], do.call(paste, c(lapply(ug, function(x) x[ind]), list(sep = ".")))), ...)
      } else print.default(setNames(x[[3L]], do.call(paste, c(ug, list(sep = ".")))), ...)
    }
  }
}

plot.GRP <- function(x, breaks = "auto", type = "l", horizontal = FALSE, ...) {
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  if(x[[1L]] <= 1e4) {
    oldpar <- par(mfrow = if(horizontal) 1:2 else 2:1,
                  mar = c(3.9,4.1,2.1,1), mgp = c(2.5,1,0))
    on.exit(par(oldpar))
  }
  if(breaks == "auto") {
    ugs <- fndistinct.default(x[[3L]])
    breaks <- if(ugs > 80) 80 else ugs
  }
  if(x[[1L]] <= 1e4) plot(seq_len(x[[1L]]), x[[3L]], type = type, xlab = "Group id", ylab = "Group Size",
                          xlim = c(1L, x[[1L]]), ylim = c(0L, bmax(x[[3L]])),
                          main = paste0("Sizes of ", x[[1L]], if(isTRUE(any(x[[6L]]))) " Ordered" else if(anyNA(x[[6L]])) "" else " Unordered", " Groups"), frame.plot = FALSE, ...)
  # grid()
  if(breaks == 1L) plot(x[[3L]][1L], x[[1L]], type = "h", ylab = "Frequency", xlab = "Group Size",
                        main = "Histogram of Group Sizes", frame.plot = FALSE, ...) else
  hist(x[[3L]], breaks, xlab = "Group Size",
       main = paste0("Histogram of Group Sizes", if(x[[1L]] > 1e4) paste0(" (N = ", x[[1L]], ")") else ""), ...)
}

as_factor_GRP <- function(x, ordered = FALSE, sep = ".") { # , ...
  # if(is.factor(x)) return(x)
  # if(!is_GRP(x)) stop("x must be a 'GRP' object")
  f <- x[[2L]]
  gr <- unclass(x[[4L]])
  if(is.null(gr)) {
    attr(f, "levels") <- as.character(seq_len(x[[1L]]))
  } else {
    if(length(gr) == 1L) {
      attr(f, "levels") <- tochar(gr[[1L]]) # or formatC ?
    } else {
      attr(f, "levels") <- do.call(paste, c(gr, list(sep = sep)))
    }
  }
  oldClass(f) <- if(ordered) c("ordered","factor","na.included") else c("factor","na.included") # previously if any(x[[6L]])
  f
}

as.factor_GRP <- function(x, ordered = FALSE) {
  .Deprecated(msg = "'as.factor_GRP' was renamed to 'as_factor_GRP'. It will be removed end of 2023, see help('collapse-renamed').")
  as_factor_GRP(x, ordered)
}

finteraction <- function(..., factor = TRUE, ordered = FALSE, sort = factor && .op[["sort"]], method = "auto", sep = ".") { # does it drop levels ? -> Yes !
  X <- if(...length() == 1L && is.list(..1)) ..1 else list(...)
  if(factor) return(as_factor_GRP(GRP.default(X, sort = sort, return.order = FALSE, method = method, call = FALSE), ordered, sep))
  if(sort || method == "radix") {
    g <- GRP.default(X, sort = sort, return.groups = FALSE, return.order = FALSE, method = method, call = FALSE)
    res <- g[[2L]]
    attr(res, "N.groups") <- g[[1L]]
  } else res <- .Call(C_group, X, FALSE, FALSE)
  oldClass(res) <- c(if(ordered) "ordered", "qG", "na.included")
  res
}

itn <- function(...) finteraction(...)

GRP.qG <- function(X, ..., group.sizes = TRUE, return.groups = TRUE, call = TRUE) {
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  gvars <- l1orlst(as.character(substitute(X))) # paste(all.vars(call), collapse = ".") # good in all circumstances ?
  ng <- attr(X, "N.groups")
  grl <- return.groups && length(groups <- attr(X, "groups"))
  if(!inherits(X, "na.included")) if(anyNA(unclass(X))) {
    ng <- ng + 1L
    X <- .Call(C_setcopyv, X, NA, ng, FALSE, FALSE, FALSE) # X[is.na(X)] <- ng
    if(grl) groups <- c(groups, NA)
  }
  st <- attr(X, "starts")
  ordered <- is.ordered(X)
  attributes(X) <- NULL

  return(`oldClass<-`(list(N.groups = ng,
                        group.id = X,
                        group.sizes = if(group.sizes) .Call(C_fwtabulate, X, NULL, ng, FALSE) else NULL, # tabulate(X, ng)  # .Internal(tabulate(X, ng))
                        groups = if(grl) `names<-`(list(groups), gvars) else NULL,
                        group.vars = gvars,
                        ordered = c(ordered = if(ordered) TRUE else NA, sorted = issorted(X)),
                        order = NULL, # starts = NULL, maxgrpn = NULL,
                        group.starts = st,
                        call = if(call) match.call() else NULL), "GRP"))
}

GRP.factor <- function(X, ..., group.sizes = TRUE, drop = FALSE, return.groups = TRUE, call = TRUE) {
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  nam <- l1orlst(as.character(substitute(X))) # paste(all.vars(call), collapse = ".") # good in all circumstances ?
  if(!inherits(X, "na.included")) X <- addNA2(X)
  if(drop) X <- .Call(Cpp_fdroplevels, X, FALSE)
  lev <- attr(X, "levels")
  nl <- length(lev)
  ordered <- is.ordered(X)
  attributes(X) <- NULL
  return(`oldClass<-`(list(N.groups = nl,
                        group.id = X,
                        group.sizes = if(group.sizes) .Call(C_fwtabulate, X, NULL, nl, FALSE) else NULL, # tabulate(X, nl) # .Internal(tabulate(X, nl))
                        groups = if(return.groups) `names<-`(list(lev), nam) else NULL,
                        group.vars = nam,
                        ordered = c(ordered = if(ordered) TRUE else NA, sorted = issorted(X)),
                        order = NULL, # starts = NULL, maxgrpn = NULL,
                        group.starts = NULL,
                        call = if(call) match.call() else NULL), "GRP"))
}

GRP.pseries <- function(X, effect = 1L, ..., group.sizes = TRUE, return.groups = TRUE, call = TRUE) {
  g <- unclass(findex(X)) # index cannot be atomic since plm always adds a time variable !
  if(length(effect) > 1L) return(GRP.default(g[effect], ...))
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  # if(length(g) > 2L) {
  #   mlg <- -length(g)
  #   nam <- paste(names(g)[mlg], collapse = ".")
  #   g <- interaction(g[mlg], drop = TRUE)
  # } else {
    nam <- if(is.character(effect)) effect else names(g)[effect]
    g <- g[[effect]] # Fastest way to do this ?
  # }
  lev <- attr(g, "levels")
  nl <- length(lev)
  ordered <- is.ordered(g)
  attributes(g) <- NULL

  return(`oldClass<-`(list(N.groups = nl,
                        group.id = g,
                        group.sizes = if(group.sizes) .Call(C_fwtabulate, g, NULL, nl, FALSE) else NULL, # tabulate(g, nl) # .Internal(tabulate(g, nl))
                        groups = if(return.groups) `names<-`(list(lev), nam) else NULL,
                        group.vars = nam,
                        ordered = c(ordered = if(ordered) TRUE else NA, sorted = issorted(g)),
                        order = NULL, # starts = NULL, maxgrpn = NULL,
                        group.starts = NULL,
                        call = if(call) match.call() else NULL), "GRP"))
}
GRP.pdata.frame <- function(X, effect = 1L, ..., group.sizes = TRUE, return.groups = TRUE, call = TRUE)
  GRP.pseries(X, effect, ..., group.sizes = group.sizes, return.groups = return.groups, call = call)

fgroup_by <- function(.X, ..., sort = .op[["sort"]], decreasing = FALSE, na.last = TRUE, return.groups = TRUE, return.order = sort, method = "auto") {          #   e <- substitute(list(...)) # faster but does not preserve attributes of unique groups !
  clx <- oldClass(.X)
  oldClass(.X) <- NULL
  m <- match(c("GRP_df", "grouped_df", "data.frame"), clx, nomatch = 0L)
  dots <- substitute(list(...))
  # vars <- all.vars(dots, unique = FALSE)
  # In case sequences of columns are passed... Think: can enable fgroup_by(mtcars, 1:cyl)
  if(any(all_funs(dots) == ":")) { # length(vars)+1L != length(dots) && any(all.names(dots) == ":")
  # Note that fgroup_by(mtcars, bla = round(mpg / cyl), vs:am) only groups by vs, and am. fselect(mtcars, bla = round(mpg / cyl), vs:am) also does the wrong thing.
    nl <- `names<-`(as.vector(seq_along(.X), "list"), names(.X))
    vars <- eval(substitute(c(...)), nl, parent.frame())
    e <- .X[vars]
    # This allows renaming...
    if(length(nam_vars <- names(vars))) {
      nonmiss <- nzchar(nam_vars)
      names(e)[nonmiss] <- nam_vars[nonmiss]
    }
    # e <- fselect(if(m[2L]) fungroup(.X) else .X, ...)
  } else {
    e <- eval(dots, .X, parent.frame())
    name <- names(e)
    vars <- all.vars(dots, unique = FALSE)
    # If something else than NSE cols is supplied, see https://github.com/SebKrantz/collapse/issues/320
    # Note: doesn't support fgroup_by(mtcars, cyl / vs), but ok, this should be named...
    # fgroup_by(mtcars, c("cyl", "vs")) gives vars == character(0)
    if(length(e) == 1L && is.null(name) && (length(vars) != 1L || !anyv(names(.X), vars))) { # !is.symbol(dots[[2L]]) || length(e[[1L]]) != length(.X[[1L]]) || is.function(e[[1L]] # Fixes #320
      e <- .X[cols2int(e[[1L]], .X, names(.X), FALSE)]
    } else {
      if(length(name)) {  # fgroup_by(mtcars, bla = round(mpg / cyl), vs, am)
        nonmiss <- nzchar(name) # -> using as.character(dots[-1L]) instead of vars
        if(!all(nonmiss)) names(e) <- `[<-`(as.character(dots[-1L]), nonmiss, value = name[nonmiss])
      } else names(e) <- vars
    }
  }
  attr(.X, "groups") <- GRP.default(e, NULL, sort, decreasing, na.last, return.groups, return.order, method, FALSE)
  # if(any(clx == "sf")) oldClass(.X) <- clx[clx != "sf"]
  # attr(.X, "groups") <- GRP.default(fselect(if(m[2L]) fungroup(.X) else .X, ...), NULL, sort, decreasing, na.last, TRUE, return.order, method, FALSE)
    # Needed: wlddev %>% fgroup_by(country) gives error if dplyr is loaded. Also sf objects etc..
    # .rows needs to be list(), NULL won't work !! Note: attaching a data.frame class calls data frame methods, even if "list" in front! -> Need GRP.grouped_df to restore object !
    # attr(.X, "groups") <- `oldClass<-`(c(g, list(.rows = list())), c("GRP", "data.frame")) # `names<-`(eval(e, .X, parent.frame()), all.vars(e))
  oldClass(.X) <- c("GRP_df",  if(length(mp <- m[m != 0L])) clx[-mp] else clx, "grouped_df", if(m[3L]) "data.frame") # clx[-m] doesn't work if clx is only "data.table" for example
    # simplest, but .X is coerced to data.frame. Through the above solution it can be a list and only receive the 'grouped_df' class
    # add_cl <- c("grouped_df", "data.frame")
    # oldClass(.X) <- c(fsetdiff(oldClass(.X), add_cl), add_cl)
  if(any(clx == "data.table")) return(alc(.X))
  .X
}

gby <- fgroup_by

group_by_vars <- function(X, by = NULL, ...) {
  clx <- oldClass(X)
  oldClass(X) <- NULL
  m <- match(c("GRP_df", "grouped_df", "data.frame"), clx, nomatch = 0L)
  if(length(by)) by <- cols2int(by, X, names(X), FALSE)
  attr(X, "groups") <- GRP.default(X[by], NULL, ..., call = FALSE) # Need to unclass because of sf and regrouping! (and some functions expect unclassed)
  oldClass(X) <- c("GRP_df",  if(length(mp <- m[m != 0L])) clx[-mp] else clx, "grouped_df", if(m[3L]) "data.frame")
  if(any(clx == "data.table")) return(alc(X))
  X
}

print.GRP_df <- function(x, ...) {
  print(fungroup(x), ...) # better !! (the method could still print groups attribute etc. ) And can also get rid of .rows() in fgroup_by and other fuzz..
  # but better keep for now, other functions in dplyr might check this and only preserve attributes if they exist. -> Nah. select(UGA_sf, addr_cname) doesn't work anyway..
  # NextMethod()
  g <- attr(x, "groups")
  if(is_GRP(g)) { # Issue Patrice flagged !
    # oldClass(g) <- NULL # could get rid of this if get rid of "data.frame" class.
    if(length(g[[3L]])) {
      su <- unclass(qsu.default(g[[3L]], stable.algo = FALSE))
      stats <- if(su[4L] == su[5L]) paste0(" [", g[[1L]], " | ", round(su[2L]), " (", round(su[3L], 1L), ")]") else
        paste0(" [", g[[1L]], " | ", round(su[2L]), " (", round(su[3L], 1L), ") ", su[4L], "-", su[5L], "]")
    } else
      stats <- paste0(" [", g[[1L]], " | ", round(length(g[[2L]]) / g[[1L]]), "]")
    # Groups: # if(any(g[[6L]])) "ordered groups" else "unordered groups", -> ordered 99% of times...
    cat("\nGrouped by: ", paste(g[[5L]], collapse = ", "), stats, "\n")
    if(inherits(x, "pdata.frame"))
      message("\nNote: 'pdata.frame' methods for flag, fdiff, fgrowth, fcumsum, fbetween, fwithin, fscale, qsu and varying\n      take precedence over the 'grouped_df' methods for these functions.")
  }
}

print.invisible <- function(x, ...) cat("")

# Still solve this properly for data.table...
`[.GRP_df` <- function(x, ...) {
  clx <- oldClass(x)
  if(any(clx == "data.table")) {
    res <- NextMethod()
    if(any(clx == "invisible")) { # for chaining...
      clx <- clx[clx != "invisible"]
      oldClass(res) <- clx # in case of early return (reduced rows)...
    }
    if(any(grepl(":=", .c(...)))) {
      eval.parent(substitute(x <- res))
      oldClass(res) <- c("invisible", clx) # return(invisible(res)) -> doesn't work here for some reason
    } else {
      if(!(is.list(res) && fnrow(res) == fnrow(x))) return(fungroup(res))
      if(is.null(attr(res, "groups"))) attr(res, "groups") <- attr(x, "groups")
      oldClass(res) <- clx
    }
  } else {
    res <- `[`(fungroup(x), ...) # does not respect data.table properties, but better for sf data frame and others which check validity of "groups" attribute
    if(!(is.list(res) && fnrow(res) == fnrow(x))) return(res)
    attr(res, "groups") <- attr(x, "groups")
    oldClass(res) <- clx
  }
  res
}

# missing doesn't work, its invisible return...
# `[.GRP_df` <- function(x, ...) {
#   tstop <- function(x) if(missing(x)) NULL else x
#   res <- tstop(NextMethod()) # better than above (problems with data.table method, but do further checks...)
#   if(is.null(res)) return(NULL)
#   if(!(is.list(res) && fnrow(res) == fnrow(x))) return(fungroup(res))
#   if(is.null(g <- attr(res, "groups"))) attr(res, "groups") <- g
#   oldClass(res) <- oldClass(x)
#   return(res)
# }

# also needed to sort out errors with dplyr ...
`[[.GRP_df` <-  function(x, ...) UseMethod("[[", fungroup(x)) # function(x, ..., exact = TRUE) .subset2(x, ..., exact = exact)
`[<-.GRP_df` <- function(x, ..., value) UseMethod("[<-", fungroup(x))
`[[<-.GRP_df` <- function(x, ..., value) UseMethod("[[<-", fungroup(x))

# Produce errors...
# print_GRP_df_core <- function(x) {
#   g <- attr(x, "groups")
#   cat("\nGrouped by: ", paste(g[[5L]], collapse = ", "),
#       # if(any(g[[6L]])) "ordered groups" else "unordered groups", -> ordered 99% of times...
#       paste0(" [", g[[1L]], " | ", round(length(g[[2L]]) / g[[1L]]), " (", round(fsd.default(g[[3L]]), 1), ")]"))
#   if(inherits(x, "pdata.frame"))
#     message("\nNote: 'pdata.frame' methods for flag, fdiff, fgrowth, fbetween, fwithin and varying\n      take precedence over the 'grouped_df' methods for these functions.")
# }
#
# head.GRP_df <- function(x, ...) {
#   NextMethod()
#   print_GRP_df_core(x)
# }
#
# tail.GRP_df <- function(x, ...) {
#   NextMethod()
#   print_GRP_df_core(x)
# }


fungroup <- function(X, ...) {
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  clx <- oldClass(X)
  attr(X, "groups") <- NULL
  oldClass(X) <- fsetdiff(clx, c("GRP_df", "grouped_df"))  # clx[clx != "grouped_df"]
  if(any(clx == "data.table")) return(alc(X))
  X
}


condCopyAttrib <- function(x, d) {
  if(is.object(x)) return(x)
  cld <- oldClass(d)
  condalcSA(x, list(names = attr(x, "names"),
                    row.names = .set_row_names(.Call(C_fnrow, x)),
                    class = cld[cld %!in% c("GRP_df", "grouped_df", "sf", "pdata.frame", "indexed_frame")]),
            any(cld == "data.table"))
  # attr(d, "groups") <- NULL
  # attr(d, "row.names") <- NULL
  # x <- copyMostAttributes(x, d)
  # attr(x, "row.names") <- rn
  # oldClass(x) <- fsetdiff(cld, c("GRP_df", "grouped_df", "sf"))
  # if(any(cld == "data.table")) return(alc(x))
  # x
}

fgroup_vars <- function(X, return = "data") {
  g <- attr(X, "groups")
  if(!is.list(g)) stop("attr(X, 'groups') is not a grouping object")
  vars <- if(is_GRP(g)) g[[5L]] else attr(g, "names")[-length(unclass(g))]
  switch(return,
    data = .Call(C_subsetCols, fungroup(X), ckmatch(vars, attr(X, "names")), TRUE),
    unique = if(is_GRP(g)) condCopyAttrib(g[[4L]], X) else .Call(C_subsetCols, g, -length(unclass(g)), FALSE), # what about attr(*, ".drop") ??
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

GRP.grouped_df <- function(X, ..., return.groups = TRUE, call = TRUE) {
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  # g <- unclass(attr(X, "groups"))
  g <- attr(X, "groups")
  if(is_GRP(g)) return(g) # return(`oldClass<-`(.subset(g, 1:8), "GRP")) # To avoid data.frame methods being called
  if(!is.list(g)) stop("attr(X, 'groups') is not a grouping object")
  oldClass(g) <- NULL
  lg <- length(g)
  gr <- g[[lg]]
  ng <- length(gr)
  gs <- vlengths(gr, FALSE)
  id <- .Call(C_groups2GRP, gr, fnrow(X), gs)
  return(`oldClass<-`(list(N.groups = ng, # The C code here speeds up things a lot !!
                        group.id = id,  # Old: rep(seq_len(ng), gs)[order(unlist(gr, FALSE, FALSE))], # .Internal(radixsort(TRUE, FALSE, FALSE, TRUE, .Internal(unlist(gr, FALSE, FALSE))))
                        group.sizes = gs,
                        groups = if(return.groups) g[-lg] else NULL, # better reclass afterwards ? -> Nope, this is only used in internal codes...
                        group.vars = names(g)[-lg],
                        ordered = c(ordered = TRUE, sorted = issorted(id)), # Important to have NA here, otherwise wrong result in gsplit (wrong optimization)
                        order = NULL, # starts = NULL, maxgrpn = NULL,
                        group.starts = NULL,
                        call = if(call) match.call() else NULL), "GRP"))
}

is_qG <- function(x) is.integer(x) && inherits(x, "qG")
is.qG <- function(x) {
  .Deprecated(msg = "'is.qG' was renamed to 'is_qG'. It will be removed end of 2023, see help('collapse-renamed').")
  inherits(x, "qG")
}

na_rm2 <- function(x, sort) {
  if(sort) return(if(is.na(x[length(x)])) x[-length(x)] else x)
  na_rm(x) # if(anyNA(x)) x[!is.na(x)] else x # use na_rm here when speed fixed..
}

Csv <- function(x, i) .Call(C_subsetVector, x, i, FALSE)

# What about NA last option to radixsort ? -> Nah, vector o becomes too short...

radixfact <- function(x, sort, ord, fact, naincl, keep, retgrp = FALSE) {
  o <- .Call(C_radixsort, TRUE, FALSE, fact || naincl || retgrp, naincl, sort, pairlist(x))
  st <- attr(o, "starts")
  sorted <- attr(o, "sorted")
  f <- if(naincl) .Call(C_frankds, o, st, attr(o, "group.sizes"), sorted) else # Fastest? -> Seems so..
        .Call(Cpp_groupid, x, if(sorted) NULL else o, 1L, TRUE, FALSE)
  if(fact) {
    if(keep) duplAttributes(f, x) else attributes(f) <- NULL
    rawlev <- Csv(x, if(sorted) st else Csv(o, st))
    attr(f, "levels") <- unattrib(tochar(if(naincl) rawlev else na_rm2(rawlev, sort)))
    oldClass(f) <- c(if(ord) "ordered", "factor", if(naincl) "na.included")
  } else {
    if(naincl) attr(f, "N.groups") <- length(st) # the order is important, this before retgrp !!
    if(retgrp) {
      rawlev <- Csv(x, if(sorted) st else Csv(o, st))
      attr(f, "groups") <-  if(naincl) rawlev else na_rm2(rawlev, sort)
    }
    oldClass(f) <- c(if(ord) "ordered", "qG", if(naincl) "na.included")
  }
  f
}

# TODO: Why is numeric to character conversion so slow?...
groupfact <- function(x, ord, fact, naincl, keep, retgrp = FALSE) {
  g <- .Call(C_groupat, x, fact || retgrp, naincl)
  if(fact) {
    st <- attr(g, "starts")
    if(keep) duplAttributes(g, x) else attributes(g) <- NULL
    attr(g, "levels") <- unattrib(tochar(if(length(st) == length(g)) x else Csv(x, st)))
    oldClass(g) <- c(if(ord) "ordered", "factor", if(naincl) "na.included")
  } else {
    if(retgrp) {
      st <- attr(g, "starts")
      attributes(g) <- NULL
      attr(g, "N.groups") <- length(st)
      attr(g, "groups") <- if(length(st) == length(g)) x else Csv(x, st)
    }
    oldClass(g) <- c(if(ord) "ordered", "qG", if(naincl) "na.included")
  }
  g
}

# TODO: Why is numeric to character conversion so slow?... this really does away with the added speed...
groupfact_sorted <- function(x, ord, fact, naincl, keep, retgrp = FALSE) {
  g <- .Call(C_groupat, x, TRUE, naincl)
  st <- attr(g, "starts")
  ng <- length(st)
  lev <- if(ng == length(x)) x else Csv(x, st)
  o <- forder.int(lev)
  # TODO: keep always add class na.included?? -> Could add anyNA attribute as output from groupat... also for groupfact...
  if(!attr(o, "sorted")) {
    if(fact || retgrp) lev <- Csv(lev, o)
    o <- forder.int(o) # This is necessary. Can optimize??
    g <- if(naincl) Csv(unattrib(o), g) else o[g] # [ propagates NA's
  }
  if(fact) {
    if(keep) duplAttributes(g, x) else attributes(g) <- NULL
    attr(g, "levels") <- unattrib(tochar(lev))
    oldClass(g) <- c(if(ord) "ordered", "factor", if(naincl) "na.included")
  } else {
    attributes(g) <- NULL
    attr(g, "N.groups") <- ng
    if(retgrp) attr(g, "groups") <- lev
    oldClass(g) <- c(if(ord) "ordered", "qG", if(naincl) "na.included")
  }
  g
}

hashfact <- function(x, sort, ord, fact, naincl, keep, retgrp = FALSE) {
  if(sort) return(groupfact_sorted(x, ord, fact, naincl, keep, retgrp)) # return(.Call(Cpp_qF, x, ord, !naincl, keep, if(fact) 1L else 2L+retgrp))
  groupfact(x, ord, fact, naincl, keep, retgrp)
}

as_factor_qG <- function(x, ordered = FALSE, na.exclude = TRUE) {
  groups <- if(is.null(attr(x, "groups"))) as.character(seq_len(attr(x, "N.groups"))) else tochar(attr(x, "groups"))
  nainc <- inherits(x, "na.included")
  if(na.exclude || nainc) {
    clx <- c(if(ordered) "ordered", "factor", if(nainc) "na.included") # can set unordered ??
  } else {
    if(anyNA(unclass(x))) {
      x <- .Call(C_setcopyv, x, NA, attr(x, "N.groups") + 1L, FALSE, FALSE, FALSE) # x[is.na(x)] <- attr(x, "N.groups") + 1L
      groups <- c(groups, NA_character_) # faster doing groups[length(groups)+1] <- NA? -> Nope, what you have is fastest !
    }
    clx <- c(if(ordered) "ordered", "factor", "na.included")
  }
  return(`attributes<-`(x, list(levels = groups, class = clx)))
}

as.factor_qG <- function(x, ordered = FALSE, na.exclude = TRUE) {
  .Deprecated(msg = "'as.factor_qG' was renamed to 'as_factor_qG'. It will be removed end of 2023, see help('collapse-renamed').")
  as_factor_qG(x, ordered, na.exclude)
}

qF <- function(x, ordered = FALSE, na.exclude = TRUE, sort = .op[["sort"]], drop = FALSE,
               keep.attr = TRUE, method = "auto") {
  if(is.factor(x) && sort) {
    if(!keep.attr && !all(names(ax <- attributes(x)) %in% c("levels", "class")))
      attributes(x) <- ax[c("levels", "class")]
    if(na.exclude || inherits(x, "na.included")) {
      clx <- oldClass(x)
      if(ordered && !any(clx == "ordered")) oldClass(x) <- c("ordered", clx) else
      if(!ordered && any(clx == "ordered")) oldClass(x) <- clx[clx != "ordered"]
      if(drop) return(.Call(Cpp_fdroplevels, x, !inherits(x, "na.included"))) else return(x)
    }
    x <- addNA2(x)
    oldClass(x) <- c(if(ordered) "ordered", "factor", "na.included")
    if(drop) return(.Call(Cpp_fdroplevels, x, FALSE)) else return(x)
  }
  if(is_qG(x)) return(as_factor_qG(x, ordered, na.exclude)) #  && sort??
  switch(method, # if((is.character(x) && !na.exclude) || (length(x) < 500 && !(is.character(x) && na.exclude)))
         auto  = if(is.double(x) && sort) # is.character(x) || is.logical(x) || !sort || length(x) < 500L
                   radixfact(x, sort, ordered, TRUE, !na.exclude, keep.attr) else if(sort && length(x) < 100000L && !is.object(x))
                     .Call(Cpp_qF, x, ordered, na.exclude, keep.attr, 1L) else
                 hashfact(x, sort, ordered, TRUE, !na.exclude, keep.attr),
         radix = radixfact(x, sort, ordered, TRUE, !na.exclude, keep.attr),
         hash  = hashfact(x, sort, ordered, TRUE, !na.exclude, keep.attr), # .Call(Cpp_qF, x, sort, ordered, na.exclude, keep.attr, 1L),
         rcpp_hash = .Call(Cpp_qF, x, ordered, na.exclude, keep.attr, 1L),
         stop("Unknown method:", method))
}

qG <- function(x, ordered = FALSE, na.exclude = TRUE, sort = .op[["sort"]],
               return.groups = FALSE, method = "auto") {
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
      if(return.groups && length(lev <- attr(x, "groups"))) lev <- c(lev, NA)
      ng <- attr(x, "N.groups") + 1L
    }
    ax <- if(return.groups) list(N.groups = ng, groups = lev, class = newclx) else
      list(N.groups = ng, class = newclx)
     # x[is.na(x)] <- ng
    return(`attributes<-`(.Call(C_setcopyv, x, NA, ng, FALSE, FALSE, FALSE), ax))
  }
  switch(method, # if((is.character(x) && !na.exclude) || (length(x) < 500 && !(is.character(x) && na.exclude)))
         auto  = if(is.double(x) && sort) # is.character(x) || is.logical(x) || !sort || length(x) < 500L
           radixfact(x, sort, ordered, FALSE, !na.exclude, FALSE, return.groups) else if(sort && length(x) < 100000L)
             .Call(Cpp_qF, x, ordered, na.exclude, FALSE, 2L+return.groups) else
           hashfact(x, sort, ordered, FALSE, !na.exclude, FALSE, return.groups),
         radix = radixfact(x, sort, ordered, FALSE, !na.exclude, FALSE, return.groups),
         hash  = hashfact(x, sort, ordered, FALSE, !na.exclude, FALSE, return.groups), # .Call(Cpp_qF, x, sort, ordered, na.exclude, FALSE, 2L+return.groups),
         rcpp_hash = .Call(Cpp_qF, x, ordered, na.exclude, FALSE, 2L+return.groups),
         stop("Unknown method:", method))
}


radixuniquevec <- function(x, sort, na.last = TRUE, decreasing = FALSE) {
  o <- .Call(C_radixsort, na.last, decreasing, TRUE, FALSE, sort, pairlist(x))
  if(attr(o, "maxgrpn") <= 1L && (!sort || attr(o, "sorted"))) return(x)
  Csv(x, if(attr(o, "sorted")) attr(o, "starts") else Csv(o, attr(o, "starts")))
}

funique <- function(x, ...) UseMethod("funique")

funique.default <- function(x, sort = FALSE, method = "auto", ...) {
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.array(x)) stop("funique currently only supports atomic vectors and data.frames")
  switch(method,
         auto = if(sort && is.numeric(x) && length(x) > 500L) radixuniquevec(x, sort, ...) else
                if(sort) .Call(Cpp_sortunique, x) else .Call(C_funique, x),
         radix = radixuniquevec(x, sort, ...),
         hash = if(sort) .Call(Cpp_sortunique, x) else .Call(C_funique, x),
         stop("method needs to be 'auto', 'hash' or 'radix'.")) # , ... adding dots gives error message too strict, package default is warning..
}

# could make faster still... not using colsubset but something more simple... no attributes needed...
# Enable by formula use ?? by or cols ?? -> cols is clearer !! also with na_omit, by could imply by-group uniqueness check...
funique.data.frame <- function(x, cols = NULL, sort = FALSE, method = "auto", ...) {
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  use.group <- switch(method, auto = !sort, hash = TRUE, radix = FALSE, stop("method needs to be 'auto', 'hash' or 'radix'."))
  o <- switchGRP(if(is.null(cols)) x else colsubset(x, cols), starts = TRUE, sort = sort, use.group = use.group, ...)
  if((use.group && length(o) == attr(o, "N.groups")) || (!use.group && attr(o, "maxgrpn") <= 1L && (!sort || attr(o, "sorted")))) # return(x)
     return(if(inherits(x, "data.table")) alc(x) else x)
  st <- if(use.group || attr(o, "sorted")) attr(o, "starts") else Csv(o, attr(o, "starts"))
  rn <- attr(x, "row.names")
  res <- .Call(C_subsetDT, x, st, seq_along(unclass(x)), FALSE)
  if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(res)
  attr(res, "row.names") <- Csv(rn, st)
  res
}

## Problem: could be confused to mean unique values within groups. Also can use ffirst() to achieve something similar
# funique.grouped_df <- function(x, ...) {
#   g <- GRP.grouped_df(x, call = FALSE)
#   if(g[[1L]] == length(g[[2L]])) return(fungroup(x))
#   st <- if(length(g$group.starts)) g$group.starts else .Call(C_ffirst, seq_along(g[[2L]]), g[[1L]], g[[2L]], NULL, FALSE)
#   rn <- attr(x, "row.names")
#   attr(x, "groups") <- NULL
#   oldClass(x) <- fsetdiff(oldClass(x), c("GRP_df", "grouped_df"))
#   res <- .Call(C_subsetDT, x, st, seq_along(unclass(x)), FALSE)
#   if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(res)
#   attr(res, "row.names") <- Csv(rn, st)
#   res
# }


funique.list <- function(x, cols = NULL, sort = FALSE, method = "auto", ...) funique.data.frame(x, cols, sort, method, ...)

funique.sf <- function(x, cols = NULL, sort = FALSE, method = "auto", ...) {
  use.group <- switch(method, auto = !sort, hash = TRUE, radix = FALSE, stop("method needs to be 'auto', 'hash' or 'radix'."))
  cols <- if(is.null(cols)) whichv(attr(x, "names"), attr(x, "sf_column"), TRUE) else
                            cols2int(cols, x, attr(x, "names"), FALSE)
  o <- switchGRP(.subset(x, cols), starts = TRUE, sort = sort, use.group = use.group, ...)
  if((use.group && length(o) == attr(o, "N.groups")) || (!use.group && attr(o, "maxgrpn") <= 1L && (!sort || attr(o, "sorted")))) return(x)
  st <- if(use.group || attr(o, "sorted")) attr(o, "starts") else Csv(o, attr(o, "starts"))
  rn <- attr(x, "row.names")
  res <- .Call(C_subsetDT, x, st, seq_along(unclass(x)), FALSE)
  if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(res)
  attr(res, "row.names") <- Csv(rn, st)
  res
}

funique.pseries <- function(x, sort = FALSE, method = "auto", drop.index.levels = "id", ...) {
  if(is.array(x)) stop("funique currently only supports atomic vectors and data.frames")
  use.group <- switch(method, auto = !sort, hash = TRUE, radix = FALSE, stop("method needs to be 'auto', 'hash' or 'radix'."))
  o <- switchGRP(x, starts = TRUE, sort = sort, use.group = use.group, ...)
  if((use.group && length(o) == attr(o, "N.groups")) || (!use.group && attr(o, "maxgrpn") <= 1L && (!sort || attr(o, "sorted")))) return(x)
  st <- if(use.group || attr(o, "sorted")) attr(o, "starts") else Csv(o, attr(o, "starts"))
  res <- Csv(x, st)
  if(length(names(x))) names(res) <- Csv(names(x), st)
  index <- findex(x)
  index_ss <- droplevels_index(.Call(C_subsetDT, index, st, seq_along(unclass(index)), FALSE), drop.index.levels)
  attr(res, if(inherits(x, "indexed_series")) "index_df" else "index") <- index_ss
  res
}

funique.pdata.frame <- function(x, cols = NULL, sort = FALSE, method = "auto", drop.index.levels = "id", ...) {
  use.group <- switch(method, auto = !sort, hash = TRUE, radix = FALSE, stop("method needs to be 'auto', 'hash' or 'radix'."))
  o <- switchGRP(if(is.null(cols)) x else colsubset(x, cols), starts = TRUE, sort = sort, use.group = use.group, ...)
  if((use.group && length(o) == attr(o, "N.groups")) || (!use.group && attr(o, "maxgrpn") <= 1L && (!sort || attr(o, "sorted")))) # return(x)
    return(if(inherits(x, "data.table")) alc(x) else x)
  st <- if(use.group || attr(o, "sorted")) attr(o, "starts") else Csv(o, attr(o, "starts"))
  rn <- attr(x, "row.names")
  res <- .Call(C_subsetDT, x, st, seq_along(unclass(x)), FALSE)
  if(!(is.numeric(rn) || is.null(rn) || rn[1L] == "1")) attr(res, "row.names") <- Csv(rn, st)
  index <- findex(x)
  index_ss <- droplevels_index(.Call(C_subsetDT, index, st, seq_along(unclass(index)), FALSE), drop.index.levels)
  if(inherits(x, "indexed_frame")) return(reindex(res, index_ss))
  attr(res, "index") <- index_ss
  res
}

fnunique <- function(x) {
  if(is.list(x) && length(unclass(x)) == 1L) x <- .subset2(x, 1L)
  if(is.atomic(x) && !is.complex(x)) .Call(C_fndistinct, x, NULL, FALSE, 1L) else
    attr(.Call(C_group, x, FALSE, FALSE), "N.groups")
}

any_duplicated <- function(x) fnunique(x) < (if(is.atomic(x)) length(x) else .Call(C_fnrow, x))

fduplicated <- function(x, all = FALSE) {
  if(all) {
    g <- .Call(C_group, x, FALSE, FALSE)
    ng <- attr(g, "N.groups")
    if(ng == length(g)) return(.Call(C_alloc, FALSE, length(g), TRUE))
    gs <- .Call(C_fwtabulate, g, NULL, ng, FALSE)
    return(.Call(C_subsetVector, gs != 1L, g, FALSE))
  }
  g <- .Call(C_group, x, TRUE, FALSE)
  starts <- attr(g, "starts")
  if(length(starts) == length(g)) return(.Call(C_alloc, FALSE, length(g), TRUE))
  .Call(C_setcopyv, .Call(C_alloc, TRUE, length(g), TRUE), starts, FALSE, FALSE, TRUE, TRUE)
}

fdroplevels <- function(x, ...) UseMethod("fdroplevels")

fdroplevels.default <- function(x, ...) {
  message("Trying to drop levels from an unsupported object: returning object")
  x
}

fdroplevels.factor <- function(x, ...) {
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  clx <- class(x)
  if(!any(clx == "factor")) stop("x needs to be a factor")
  .Call(Cpp_fdroplevels, x, !any(clx == "na.included"))
}

fdroplevels.data.frame <- function(x, ...) {
  # if(!missing(...)) unused_arg_action(match.call(), ...)
  res <- duplAttributes(lapply(unattrib(x), function(y)
    if(is.factor(y)) .Call(Cpp_fdroplevels, y, !inherits(y, "na.included")) else y), x)
  if(inherits(x, "data.table")) return(alc(res))
  res
}

fdroplevels.list <- function(x, ...) {
  duplAttributes(lapply(unattrib(x), function(y)
    if(is.factor(y)) .Call(Cpp_fdroplevels, y, !inherits(y, "na.included")) else y), x)
}


