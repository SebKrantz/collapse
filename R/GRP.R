# Cuniqlengths <- data.table:::Cuniqlengths
# Cfrank <- data.table:::Cfrank
# forderv <- data.table:::forderv

# or getGRP: more speed using tabulate ?? or code your own??
GRP <- function(X, ...) { # by = NULL, sort = TRUE, order = 1L, na.last = FALSE, return.groups = TRUE, return.order = FALSE
  UseMethod("GRP", X)
}

forderv <- function(x, by = seq_along(x), retGrp = FALSE, sort = TRUE, order = 1L, na.last = FALSE) {
  if(is.atomic(x)) {
    if(!missing(by) && !is.null(by)) stop("x is a single vector, non-NULL 'by' doesn't make sense")
    by = NULL
  } else {
    if(!length(x)) return(integer(0L))
    if(length(order) == 1L) order = rep(order, length(by))
  }
  order = as.integer(order)
  .Call(forder, x, by, retGrp, sort, order, na.last)
}

GRP.default <- function(X, by = NULL, sort = TRUE, order = 1L, na.last = FALSE,
                        return.groups = TRUE, return.order = FALSE, ...) { # , gs = TRUE # o

  call <- match.call()

  if(is.list(X) && is.null(by)) {
    by <- seq_along(X)
    namby <- names(X)
  } else if(is.call(by)) {
    namby <- all.vars(by)
    by <- anyNAerror(match(namby, names(X)), "Unknown column names!")
  } else if(is.null(by)) {
    namby <- paste(all.vars(call), collapse = ".")  # deparse(substitute(X)) # all.vars is faster !!!
  } else if(is.character(by)) {
    namby <- by
    by <- anyNAerror(match(by, names(X)), "Unknown column names!")
  } else {
    by <- as.integer(by)
    namby <- names(X)[by]
  }

  o <- forderv(X, by, TRUE, sort, order, na.last)
  f <- attr(o, "starts")

  if(length(o)) { # if ordered, returns 0
    len <- .Call(uniqlengths, f, length(o))
    grpuo <- .Call(frank, o, f, len, "dense")
    ordered <- c(GRP.sort = sort, initially.ordered = FALSE)
    if(return.groups) { # subsetDT is faster than base, and subsetVector preserves variable labels !!
        groups <- if(is.list(X)) .Call(subsetDT, X, o[f], by) else `names<-`(list(.Call(subsetVector, X, o[f])), namby)
    } else groups <- NULL
  } else {
    lx <- if(inherits(X, "data.frame")) nrow(X) else if(is.atomic(X)) length(X) else length(X[[1L]])
    len <- .Call(uniqlengths, f, lx) # data.table:::Cuniqlengths # or cumsubtract !!!
    grpuo <- rep.int(seq_along(len), len) # rep.int fastest ?? -> about same speed as rep !!
    ordered <- c(GRP.sort = sort, initially.ordered = TRUE)
    if(return.groups) {
      groups <- if(is.list(X)) .Call(subsetDT, X, f, by) else `names<-`(list(.Call(subsetVector, X, f)), namby)
    } else groups <- NULL
  }
  if(!return.order) o <- NULL
  return(`class<-`(list(N.groups = length(f),
                        group.id = grpuo,
                        group.sizes = len,
                        groups = groups,
                        group.vars = namby,
                        ordered = ordered,
                        order = o,
                        call = call), "GRP"))
}

is.GRP <- function(x) inherits(x, "GRP")

group_names.GRP <- function(g, force.char = TRUE) {
  if(is.null(g[[4L]])) return(NULL)
  groups <- g[[4L]]
  if(length(groups) == 1L) {
   if(force.char && !is.character(groups[[1L]])) paste0(groups[[1L]]) else groups[[1L]]
  } else do.call(paste, c(groups, list(sep = ".")))
}

print.GRP <- function(g, n = 6, ...) {
  ord <- g[[6L]]
  cat(paste("collapse grouping object of length",length(g[[2L]]),"with",
            g[[1L]],ifelse(any(ord),"ordered","unordered"),"groups"), fill = TRUE)
  cat("\nCall: ", paste0(deparse(g[[8L]]),", ",ifelse(ord[2L],"ordered","unordered")), "\n\n", sep = "")
  cat("Distribution of group sizes: ", fill = TRUE)
  print.summaryDefault(summary.default(g[[3L]]))
  if(!is.null(g[[4L]])) {
    ug <- g[[4L]]
    cat("\nGroups with sizes: ", fill = TRUE)
    if(length(ug) == 1L) {
      ug <- ug[[1L]]
      if(length(ug) > 2L*n) {
        ind <- seq.int(g[[1L]]-n+1L, g[[1L]])
        print.default(setNames(g[[3L]][1:n], ug[1:n]))
        cat("  ---", fill = TRUE)
        print.default(setNames(g[[3L]][ind], ug[ind]))
      } else print.default(setNames(g[[3L]], ug))
    } else {
      if(length(ug[[1L]]) > 2L*n) {
        ind <- seq.int(g[[1L]]-n+1L, g[[1L]])
        print.default(setNames(g[[3L]][1:n], do.call(paste, c(lapply(ug, function(x)x[1:n]), list(sep = ".")))))
        cat("  ---", fill = TRUE)
        print.default(setNames(g[[3L]][ind], do.call(paste, c(lapply(ug, function(x)x[ind]), list(sep = ".")))))
      } else print.default(setNames(g[[3L]], do.call(paste, c(ug, list(sep = ".")))))
    }
  }
}

plot.GRP <- function(g, breaks = "auto", type = "s", horizontal = FALSE, ...) {
  settings <- par(c("mfrow","mar","mgp"))
  par(mfrow = if(horizontal) 1:2 else 2:1, mar = c(3.9,4.1,2.1,1), mgp = c(2.5,1,0))
  if(breaks == "auto") {
    ugs <- length(unique.default(g[[3L]]))
    breaks <- if(ugs > 80) 80 else ugs
  }
  plot(seq_len(g[[1L]]), g[[3L]], type = type, xlab = "Group id", ylab = "Group Size",
       main = paste0("Sizes of ",g[[1L]]," ",ifelse(any(g[[6L]]),"Ordered","Unordered")," Groups"), frame.plot = FALSE)
  if(breaks == 1L) plot(g[[3L]][1L], g[[1L]], type = "h", ylab = "Frequency", xlab = "Group Size",
                        main = "Histogram of Group Sizes", frame.plot = FALSE) else
  hist(g[[3L]], breaks, xlab = "Group Size", main = "Histogram of Group Sizes")
  par(settings)
}

as.factor.GRP <- function(g) {
  if(is.factor(g)) return(g)
  if(!is.GRP(g)) stop("g must be a 'GRP' object")
  f <- g[[2L]]
  gr <- g[[4L]]
  if(is.null(gr)) {
    attr(f, "levels") <- as.character(seq_len(g[[1L]]))
  } else {
    if(length(gr) == 1L) {
      attr(f, "levels") <- if(is.character(gr[[1L]])) gr[[1L]] else as.character(gr[[1L]]) # or formatC ??
    } else {
      attr(f, "levels") <- do.call(paste, c(gr, list(sep = ".")))
    }
  }
  class(f) <- if(any(g[[6L]])) c("ordered","factor") else "factor"
  return(f)
}

GRP.qG <- function(X, ...) {
  # nam <- deparse(substitute(X)) # takes 9 microseconds !!, all vars on call is faster !!
  ng <- attr(X, "N.groups")
  ordered <- if(is.ordered(X)) c(TRUE,TRUE) else c(FALSE,FALSE)
  attributes(X) <- NULL
  call <- match.call()
  return(`class<-`(list(N.groups = ng,
                        group.id = X,
                        group.sizes = tabulate(X, ng), # .Internal(tabulate(X, ng)),
                        groups = NULL,
                        group.vars = paste(all.vars(call), collapse = "."),
                        ordered = ordered,
                        order = NULL,
                        call = call), "GRP"))
}

GRP.factor <- function(X, ...) {
  # nam <- deparse(substitute(X)) # takes 9 microseconds !!
  lev <- attr(X, "levels")
  nl <- length(lev)
  ordered <- if(is.ordered(X)) c(TRUE,TRUE) else c(FALSE,FALSE)
  attributes(X) <- NULL
  call <- match.call()
  nam <- paste(all.vars(call), collapse = ".")
  return(`class<-`(list(N.groups = nl,
                        group.id = X,
                        group.sizes = tabulate(X, nl), # .Internal(tabulate(X, nl)), # faster tabulating factors ??? -> Nope, same speed of nl is supplied !!
                        groups = `names<-`(list(lev), nam),
                        group.vars = nam,
                        ordered = ordered,
                        order = NULL,
                        call = call), "GRP"))
}

GRP.pseries <- function(X, ...) {
  g <- attr(X, "index") # index cannot be atomic since plm always adds a time variable !!
  if(length(g) > 2L) {
    mlg <- -length(g)
    nam <- paste(names(g)[mlg], collapse = ".")
    g <- interaction(g[mlg], drop = TRUE)
  } else {
    nam <- names(g)[1L]
    g <- g[[1L]]
  }
  lev <- attr(g, "levels")
  nl <- length(lev)
  ordered <- if(is.ordered(g)) c(TRUE,TRUE) else c(FALSE,FALSE)
  attributes(g) <- NULL
  return(`class<-`(list(N.groups = nl,
                        group.id = g,
                        group.sizes = tabulate(g, nl), # .Internal(tabulate(g, nl)), # faster tabulating factors ??? -> Nope, same speed of nl is supplied !!
                        groups = `names<-`(list(lev), nam),
                        group.vars = nam,
                        ordered = ordered,
                        order = NULL,
                        call = match.call()), "GRP"))
}
GRP.pdata.frame <- function(X, ...) GRP.pseries(X, ...)

GRP.grouped_df <- function(X, ...) {
  g <- unclass(attr(X, "groups"))
  lg <- length(g)
  gr <- g[[lg]]
  ng <- length(gr)
  gs <- lengths(gr) # faster sorting still ?? qsort ?? or data.table forder -> nope, slower than order !!
  return(`class<-`(list(N.groups = ng,
                        group.id = rep(seq_len(ng), gs)[sort.list(unlist(gr, FALSE, FALSE))], # .Internal(radixsort(TRUE, FALSE, FALSE, TRUE, .Internal(unlist(gr, FALSE, FALSE))))  #
                        group.sizes = gs,
                        groups = g[-lg],
                        group.vars = names(g)[-lg],
                        ordered = c(TRUE, TRUE),
                        order = NULL,
                        call = match.call()), "GRP"))
}

