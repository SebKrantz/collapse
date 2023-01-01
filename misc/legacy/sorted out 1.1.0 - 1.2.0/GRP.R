# Cuniqlengths <- data.table:::Cuniqlengths
# Cfrank <- data.table:::Cfrank
# forderv <- data.table:::forderv

GRP <- function(X, ...) UseMethod("GRP") # , X

forderv <- function(x, by = seq_along(x), retGrp = FALSE, sort = TRUE, order = 1L, na.last = FALSE) {
  if(is.atomic(x)) {
    if(!missing(by) && !is.null(by)) stop("x is a single vector, non-NULL 'by' doesn't make sense")
    by <- NULL
  } else {
    if(!length(unclass(x))) return(integer(0L))
    if(length(order) == 1L) order <- rep(order, length(by))
  }
  .Call(C_forder, x, by, retGrp, sort, as.integer(order), na.last)
}

GRP.default <- function(X, by = NULL, sort = TRUE, order = 1L, na.last = TRUE,
                        return.groups = TRUE, return.order = FALSE, ...) { # , gs = TRUE # o

  if(!missing(...)) unused_arg_action(match.call(), ...)
  call <- match.call()

  if(is.list(X)) {
    if(is.null(by)) {
      by <- seq_along(unclass(X))
      namby <- attr(X, "names")
      if(is.null(namby)) attr(X, "names") <- namby <- paste0("Group.", by)
    } else if(is.call(by)) {
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
 } else namby <- paste(all.vars(call), collapse = ".")  # deparse(substitute(X)), all.vars is faster !

  o <- forderv(X, by, TRUE, sort, order, na.last)
  f <- attr(o, "starts")

  if(length(o)) { # if ordered, returns 0
    len <- .Call(C_uniqlengths, f, length(o))
    grpuo <- .Call(C_frank, o, f, len, "dense")
    ordered <- c(GRP.sort = sort, initially.ordered = FALSE)
    if(return.groups) { # subsetVector preserves variable labels !
        groups <- if(is.list(X)) .Call(C_subsetDT, X, o[f], by) else `names<-`(list(.Call(C_subsetVector, X, o[f])), namby)
    } else groups <- NULL
  } else {
    lx <- if(is.list(X)) fnrow2(X) else length(X)
    len <- .Call(C_uniqlengths, f, lx) # data.table:::Cuniqlengths # or cumsubtract?
    grpuo <- rep.int(seq_along(len), len) # rep.int fastest ?? -> about same speed as rep
    ordered <- c(GRP.sort = sort, initially.ordered = TRUE)
    if(return.groups) {
      groups <- if(is.list(X)) .Call(C_subsetDT, X, f, by) else `names<-`(list(.Call(C_subsetVector, X, f)), namby)
    } else groups <- NULL
  }
  return(`class<-`(list(N.groups = length(f),
                        group.id = grpuo,
                        group.sizes = len,
                        groups = groups,
                        group.vars = namby,
                        ordered = ordered,
                        order = if(return.order) o else NULL,
                        call = call), "GRP"))
}

is.GRP <- function(x) inherits(x, "GRP")

group_names.GRP <- function(x, force.char = TRUE) { # , ...
  if(is.null(x[[4L]])) return(NULL)
  groups <- x[[4L]]
  if(length(groups) == 1L) {
   if(force.char && !is.character(groups[[1L]])) paste0(groups[[1L]]) else groups[[1L]]
  } else do.call(paste, c(groups, list(sep = ".")))
}

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

as.factor.GRP <- function(x) { # , ...
  if(is.factor(x)) return(x)
  if(!is.GRP(x)) stop("x must be a 'GRP' object")
  f <- x[[2L]]
  gr <- unclass(x[[4L]])
  if(is.null(gr)) {
    attr(f, "levels") <- as.character(seq_len(x[[1L]]))
  } else {
    if(length(gr) == 1L) {
      attr(f, "levels") <- if(is.character(gr[[1L]])) gr[[1L]] else as.character(gr[[1L]]) # or formatC ?
    } else {
      attr(f, "levels") <- do.call(paste, c(gr, list(sep = ".")))
    }
  }
  class(f) <- if(any(x[[6L]])) c("ordered","factor","na.included") else c("factor","na.included") # NA included ?
  return(f)
}

GRP.qG <- function(X, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  ng <- attr(X, "N.groups")
  if(!inherits(X, "na.included") || anyNA(unclass(X))) X[is.na(X)] <- ng
  ordered <- if(is.ordered(X)) c(TRUE,TRUE) else c(FALSE,FALSE)
  attributes(X) <- NULL
  call <- match.call()
  return(`class<-`(list(N.groups = ng,
                        group.id = X,
                        group.sizes = tabulate(X, ng), # .Internal(tabulate(X, ng))
                        groups = NULL,
                        group.vars = paste(all.vars(call), collapse = "."),
                        ordered = ordered,
                        order = NULL,
                        call = call), "GRP"))
}

GRP.factor <- function(X, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(!inherits(X, "na.included")) X <- addNA2(X)
  lev <- attr(X, "levels")
  nl <- length(lev)
  ordered <- if(is.ordered(X)) c(TRUE,TRUE) else c(FALSE,FALSE)
  attributes(X) <- NULL
  call <- match.call()
  nam <- paste(all.vars(call), collapse = ".")
  return(`class<-`(list(N.groups = nl,
                        group.id = X,
                        group.sizes = tabulate(X, nl), # .Internal(tabulate(X, nl))
                        groups = `names<-`(list(lev), nam),
                        group.vars = nam,
                        ordered = ordered,
                        order = NULL,
                        call = call), "GRP"))
}

GRP.pseries <- function(X, effect = 1L, ...) {
  g <- unclass(attr(X, "index")) # index cannot be atomic since plm always adds a time variable !
  if(length(effect) > 1L) return(GRP.default(g[effect], ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
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
  return(`class<-`(list(N.groups = nl,
                        group.id = g,
                        group.sizes = tabulate(g, nl), # .Internal(tabulate(g, nl))
                        groups = `names<-`(list(lev), nam),
                        group.vars = nam,
                        ordered = ordered,
                        order = NULL,
                        call = match.call()), "GRP"))
}
GRP.pdata.frame <- function(X, effect = 1L, ...) GRP.pseries(X, effect, ...)

fgroup_by <- function(X, ..., sort = TRUE, order = 1L, na.last = TRUE, return.order = FALSE) {      #   e <- substitute(list(...)) # faster but does not preserve attributes of unique groups !!
  attr(X, "groups") <- GRP.default(fselect(X, ...), NULL, sort, order, na.last, TRUE, return.order) # `names<-`(eval(e, X, parent.frame()), all.vars(e))
  add_cl <- c("tbl_df", "tbl", "grouped_df")
  class(X) <- c(add_cl, fsetdiff(class(X), add_cl)) # necesssary to avoid printing errors... (i.e. wrong group object etc...)
  X
}

fgroup_vars <- function(X, return = c("data", "unique")) {
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

GRP.grouped_df <- function(X, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  # g <- unclass(attr(X, "groups"))
  g <- attr(X, "groups")
  if(is.GRP(g)) return(g)
  class(g) <- NULL
  lg <- length(g)
  gr <- g[[lg]]
  ng <- length(gr)
  gs <- lengths(gr, FALSE)
  return(`class<-`(list(N.groups = ng, # The cpp here speeds up things a lot !!
                        group.id = .Call(Cpp_groups2GRP, gr, fnrow2(X), gs),  # Old: rep(seq_len(ng), gs)[order(unlist(gr, FALSE, FALSE))], # .Internal(radixsort(TRUE, FALSE, FALSE, TRUE, .Internal(unlist(gr, FALSE, FALSE))))
                        group.sizes = gs,
                        groups = g[-lg], # better reclass afterwards ?
                        group.vars = names(g)[-lg],
                        ordered = c(TRUE, TRUE),
                        order = NULL,
                        call = match.call()), "GRP"))
}

