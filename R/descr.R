
# Super fast tabulation of a single atomic vector, with various sorting options
fsorttable <- function(x, srt, w = NULL) {
  if(is.factor(x)) {
    lev <- attr(x, "levels")
    t <- .Call(C_fwtabulate, x, w, length(lev), !inherits(x, "na.included")) # tabulate(x, nbins = length(lev)) # skips missing values !!
    names(t) <- lev
    sorted <- TRUE
  } else {
    sorted <- FALSE
    g <- .Call(C_groupat, x, TRUE, FALSE) # FALSE = keeps NA
    t <- .Call(C_fwtabulate, g, w, attr(g, "N.groups"), TRUE) # TRUE = check for NA's and skip them
    names(t) <- Csv(x, attr(g, "starts"))
    # This seems is slightly faster with not too many distinct values, but less straightforward
    # g <- .Call(C_group, x, TRUE, is.null(w))
    # t <- if(is.null(w)) attr(g, "group.sizes") else
    #   .Call(C_fwtabulate, g, w, attr(g, "N.groups"), FALSE)
    # nam <- Csv(x, attr(g, "starts"))
    # names(t) <- nam
    # if(anyNA(nam)) t <- t[-whichNA(nam)]
  }
  switch(srt,
    value = if(sorted || attr(o <- forder.int(names(t)), "sorted")) t else t[o],
    # "quick" sort seems best, based on multiple datasets, but "radix" (second best) keeps ties in order...
    # sort.int(t, method = "radix", decreasing = TRUE, na.last = TRUE)
    freq = if(attr(o <- forder.int(t, decreasing = TRUE), "sorted")) t else t[o],
    none = t,
    stop("sort.table must be one of 'value', 'freq' or 'none'"))
}

# Same for grouped data, building on qtab()
sorttable2D <- function(x, f, srt, w = NULL) {
  if(is.factor(x)) sorted <- TRUE
  else {
    sorted <- switch(srt, value = TRUE, FALSE)
    x <- qF(x, sort = sorted)
  }
  t <- qtab(x, f, w = w, dnn = NULL)
  switch(srt,
         value = if(sorted || attr(o <- forder.int(dimnames(t)[[1L]]), "sorted")) t else t[o, , drop = FALSE],
         freq = if(attr(o <- forder.int(frowSums(t), decreasing = TRUE), "sorted")) t else t[o, , drop = FALSE],
         none = t,
         stop("sort.table must be one of 'value', 'freq' or 'none'"))
}
# Extended version including totals and transpose option: better do that in print!
# sorttable2D <- function(x, f, srt, w = NULL, transpose = FALSE) {
#   if(is.factor(x)) sorted <- TRUE
#   else {
#     sorted <- switch(srt, value = TRUE, FALSE)
#     x <- qF(x, sort = sorted)
#   }
#   if(transpose) {
#     t <- qtab(f, x, w = w, dnn = NULL)
#     tot <- unattrib(fsum.matrix(t))
#     t <- rbind(t, Total = tot)
#   } else {
#     t <- qtab(x, f, w = w, dnn = NULL)
#     tot <- if(is.double(w)) frowSums(t) else as.integer(frowSums(t))
#     t <- cbind(t, Total = tot)
#   }
#   switch(srt,
#          value = if(sorted || attr(o <- forder.int(dimnames(t)[[1L+transpose]]), "sorted")) t else if(transpose) t[, o, drop = FALSE] else t[o, , drop = FALSE],
#          freq = if(attr(o <- forder.int(tot, decreasing = TRUE), "sorted")) t else if(transpose) t[, o, drop = FALSE] else t[o, , drop = FALSE],
#          none = t,
#          stop("sort.table must be one of 'value', 'freq' or 'none'"))
# }


# X = wlddev; by = ~ income; w = ~ POP;
# cols = NULL; Ndistinct = TRUE; higher = TRUE; table = TRUE; sort.table = "freq"
# Qprobs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99); Qtype = 7L
# label.attr = 'label'; stepwise = FALSE; nam = "wlddev"; dotsok = TRUE
# fndistinctC = collapse:::fndistinctC; fsumC = collapse:::fsumC;
# fsorttable = collapse:::fsorttable; frowSums = collapse:::frowSums

# Expects X to be a plain list and nam the name of the dataset
descr_core <- function(X, nam, by = NULL, w = NULL, Ndistinct = TRUE, higher = TRUE, table = TRUE, sort.table = "freq",
                       Qprobs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), Qtype = 7L,
                       label.attr = 'label', stepwise = FALSE, ...) {

  dotsok <- if(missing(...)) TRUE else names(substitute(c(...))[-1L]) %!in% c('pid','g')

  # Checking for numeric data
  num <- .Call(C_vtypes, X, 1L) # vapply(unattrib(X), is.numeric, TRUE)
  Nnum <- bsum(num)

  # Define functions to process numeric data
  if(Nnum > 0L) {

    if(Ndistinct && dotsok) {
      armat <- if(is.null(by)) function(x, y) c(x[1L], Ndist = y, x[-1L]) else
                               function(x, y) cbind(x[, 1L, drop = FALSE], Ndist = y, x[, -1L])
      numstats <- function(x, ...) armat(qsu.default(x, by, w = w, higher = higher, ...), fndistinctC(x, by))
    } else numstats <- function(x, ...) qsu.default(x, by, w = w, higher = higher, ...)

    quantiles <- if(is.null(by)) function(x, ...) .quantile(x, Qprobs, w, type = Qtype, names = TRUE, ...) else
                                 function(x, ...) BY.default(x, by, .quantile, probs = Qprobs, w = w, type = Qtype, names = TRUE, expand.wide = TRUE, ...)

    # This function will be applied to different columns.
    descrnum <- if(is.numeric(Qprobs)) function(x, ...) list(Class = class(x), Label = attr(x, label.attr), Stats = numstats(x, ...),
                                                             Quant = quantiles(x)) else
                                       function(x, ...) list(Class = class(x), Label = attr(x, label.attr), Stats = numstats(x, ...))
  }

  # Non-numeric data, assumed to have at least some categorical variables (could also be date)
  if(Nnum != length(num)) {

    if(table && !is.null(by)) {
      f <- as_factor_GRP(by)
      descrcat <- function(x) {
        tab <- sorttable2D(x, f, sort.table, w)
        list(Class = class(x), Label = attr(x, label.attr),
             Stats = if(Ndistinct) cbind(N = fsum.matrix(tab), Ndist = fsum.matrix(tab > 0L)) else cbind(N = fsum.matrix(tab)), # Good ? Faster ???
             Table = tab)
      }
    } else if(table) {
      descrcat <- function(x) {
        tab <- fsorttable(x, sort.table, w)
        list(Class = class(x), Label = attr(x, label.attr),
             Stats = if(Ndistinct) c(N = fsumC(tab), Ndist = length(tab)) else `names<-`(fsumC(tab), "N"),
             Table = tab)
      }
    } else {
      descrcat <- function(x) list(Class = class(x), Label = attr(x, label.attr),
                                   Stats = if(Ndistinct) c(N = fnobsC(x), Ndist = fndistinctC(x)) else `names<-`(fnobsC(x), "N"))
    }

  }

  descrdate <- if(is.null(by)) function(x) list(Class = class(x), Label = attr(x, label.attr),
                                                Stats = c(if(Ndistinct) c(N = fnobsC(x), Ndist = fndistinctC(x)) else `names<-`(fnobsC(x), "N"), `names<-`(frange(x), c("Min", "Max")))) else
                               function(x) list(Class = class(x), Label = attr(x, label.attr),
                                                Stats = cbind(N = fnobs.default(x, by), Ndist = if(Ndistinct) fndistinctC(x, by) else NULL,
                                                              Min = fmin.default(x, by, use.g.names = FALSE), Max = fmax.default(x, by, use.g.names = FALSE)))

  # Result vector and attributes
  res <- vector('list', length(X))
  ares <- list(names = names(X), name = nam, N = fnrow(X),
               arstat = !dotsok, table = table, groups = by, weights = w, class = "descr")

  # Computation
  if(stepwise) { # This means we compute one by one, mainly for printing...
    attributes(res) <- ares
    cat('Dataset: ', nam,', ', length(res), ' Variables, N = ', ares$N, "\n", sep = "")
    cat(paste(rep("-", .Options$width), collapse = ""), "\n", sep = "")
    for(i in seq_along(X)) {
      invisible(readline(prompt = sprintf("Press [enter] for variable %s/%s or [esc] to exit", i, length(res))))
      xi <- X[[i]]
      res[[i]] <- if(is.numeric(xi)) descrnum(xi, ...) else if(is_date(xi)) descrdate(xi) else descrcat(xi)
      print(res[i], noheader = TRUE)
    }
  } else {
    res[num] <- lapply(X[num], descrnum, ...)
    if(Nnum != length(num)) {
      date <- vapply(unattrib(X), is_date, TRUE)
      if(any(date)) {
        res[date] <- lapply(X[date], descrdate)
        cat <- !(num | date)
      } else cat <- !num
      res[cat] <- lapply(X[cat], descrcat)
    }
    attributes(res) <- ares
  }
  return(if(stepwise) invisible(res) else res)
}

# Since v1.9.0, descr() is generic, with a grouped_df method
descr <- function(X, ...) UseMethod("descr")

descr.default <- function(X, by = NULL, w = NULL, cols = NULL, Ndistinct = TRUE, higher = TRUE, table = TRUE, sort.table = "freq",
                          Qprobs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), Qtype = 7L,
                          label.attr = 'label', stepwise = FALSE, ...) {

  # Getting input information
  nam <- l1orlst(as.character(substitute(X)))

  # Unclassing and (if necessary) transforming X
  if(is.list(X)) {
    is_sf <- inherits(X, "sf")
    # if(inherits(X, "POSIXlt")) X <- list(X = as.POSIXct(X))
    if(inherits(X, "pdata.frame")) X <- unindex(X)
    class(X) <- NULL
    if(is_sf) X[[attr(X, "sf_column")]] <- NULL
  } else {
    if(inherits(X, "pseries")) X <- unindex(X)
    is_1D <- is.null(dim(X))
    X <- unclass(qDF(X))
    if(is_1D) names(X) <- nam
  }

  # Processing by and w arguments: inspired by qsu()
  if(is.call(by) || is.call(w)) {
    v <- NULL
    if(is.call(by)) {
      if(length(by) == 3L) {
        v <- ckmatch(all.vars(by[[2L]]), names(X))
        byn <- ckmatch(all.vars(by[[3L]]), names(X))
      } else byn <- ckmatch(all.vars(by), names(X))
      by <- GRP.default(X, byn, call = FALSE) # , ...
    } else {
      if(!is.null(by)) by <- GRP.default(by, call = FALSE) # , ...
      byn <- NULL
    }
    if(is.call(w)) {
      widn <- ckmatch(all.vars(w), names(X))
      w <- eval(w[[2L]], X, attr(w, ".Environment"))
    } else widn <- NULL
    X <- X[if(length(v)) v else if(is.null(cols)) -c(byn, widn) else cols2int(cols, X, names(X), FALSE)]
  } else {
    if(!is.null(by)) by <- GRP.default(by, call = FALSE) # , ...
    if(length(cols)) X <- X[cols2int(cols, X, names(X), FALSE)]
  }

  descr_core(X, nam, by, w, Ndistinct, higher, table, sort.table, Qprobs, Qtype, label.attr, stepwise, ...)
}

# Benefit of grouped_df method: better control on how data is grouped with fgroup_by(), selection with fselect() etc.
descr.grouped_df <- function(X, w = NULL, Ndistinct = TRUE, higher = TRUE, table = TRUE, sort.table = "freq",
                             Qprobs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), Qtype = 7L,
                             label.attr = 'label', stepwise = FALSE, ...) {

  # Getting input information
  nam <- l1orlst(as.character(substitute(X)))
  wsym <- substitute(w)
  by <- GRP.grouped_df(X, call = FALSE)

  # Unclassing and (if necessary) transforming X
  is_sf <- inherits(X, "sf")
  if(inherits(X, "pdata.frame")) X <- unindex(X)
  class(X) <- NULL
  if(is_sf) X[[attr(X, "sf_column")]] <- NULL

  # Getting group indices
  byn <- which(names(X) %in% by[[5L]])
  # Processing weights and combining indices with group indices
  if(!is.null(wsym)) {
    w <- eval(wsym, X, parent.frame()) # This allows w to be a function of multiple variables
    if(length(wn <- which(names(X) %in% all.vars(wsym)))) {
      if(any(byn %in% wn)) stop("Weights coincide with grouping variables!")
      byn <- c(byn, wn)
    }
  }

  if(length(byn)) X <- X[-byn] # Subsetting X

  descr_core(X, nam, by, w, Ndistinct, higher, table, sort.table, Qprobs, Qtype, label.attr, stepwise, ...)
}

# Methods ----------------------------------------------------------

`[.descr` <- function(x, ...) copyMostAttributes(.subset(x, ...), x)

# Idea: compact = TRUE: combines stats and quantiles, and omits summaries of table frequencies. enable by default for grouped summaries ?
print.descr <- function(x, n = 14, perc = TRUE, digits = 2, t.table = TRUE, summary = TRUE, reverse = FALSE, stepwise = FALSE, ...) {
  noheader <- !missing(...) && isTRUE(list(...)$noheader)
  oldClass(x) <- NULL
  w <- paste(rep("-", .Options$width), collapse = "")
  arstat <- attr(x, "arstat")
  DSname <- attr(x, "name")
  DSN <- attr(x, "N")
  cb <- function(...) if(t.table) cbind(...) else formatC(rbind(...), drop0trailing = TRUE)
  ct <- function(z) if(t.table) cbind(Freq = z) else z
  if(reverse) x <- rev.default(x) else if(!noheader) {
    cat('Dataset: ', DSname,', ',length(x), ' Variables, N = ', DSN, "\n", sep = "")
    cat(w, "\n", sep = "")
  }
  nam <- names(x) # Needs to be here
  for(i in seq_along(x)) {
    if(stepwise) invisible(readline(prompt = sprintf("Press [enter] for variable %s/%s or [esc] to exit", i, length(x))))
    xi <- x[[i]]
    cat(nam[i]," (",strclp(xi[[1L]]),"): ",xi[[2L]], "\n", sep = "")
    stat <- xi[[3L]]
    if(stat[[1L]] != DSN) cat("Statistics (", round((1-stat[[1L]]/DSN)*100, 2), "% NAs)\n", sep = "")
    else cat("Statistics\n")
    if(any(xi[[1L]] %in% c("Date", "POSIXct")))
      print.default(c(stat[1:2], as.character(`oldClass<-`(stat[3:4], xi[[1L]]))),
                     quote = FALSE, right = TRUE, print.gap = 2)
    else print.qsu(stat, digits)
    if(length(xi) > 3L) {
      if(arstat) cat("\n")
      if(names(xi)[4L] == "Table") {
        cat("Table\n")
        t <- unclass(xi[[4L]])
        if(length(t) <= n) {
          if(perc) print.default(cb(Freq = t, Perc = round(t/bsum(t)*100, digits)), right = TRUE, print.gap = 2, quote = FALSE) else
            print.table(ct(t))
        } else {
          t1 <- t[seq_len(n)]
          st <- bsum(t)
          rem <- `names<-`(st-bsum(t1), sprintf("... %s Others", length(t)-n))
          if(perc) {
            pct <- unattrib(t1)/st*100
            print.default(cb(Freq = c(t1, rem), Perc = round(c(pct, 100-bsum(pct)), digits)), right = TRUE, print.gap = 2, quote = FALSE)
            # cat("...\n")
          } else {
            print.table(ct(c(t1, rem)))
            # cat("...\n")
          }
          if(summary) {
            cat("\nSummary of Table Frequencies\n")
            print.summaryDefault(summary.default(t), digits)
          }
        }
      } else {
        cat("Quantiles\n")
        print.qsu(xi[[4L]], digits)
      }
    }
    cat(w, "\n", sep = "") # More compressed -> better !
    # cat("\n", w, "\n", sep = "")
  }
  if(reverse && !noheader) cat('Dataset: ', DSname,', ',length(x), ' Variables, N = ', DSN, "\n", sep = "")
  invisible(x)
}

# Note: This does not work for array stats (using g or pid.. )
as.data.frame.descr <- function(x, ..., gid = "Group", wsum = FALSE, stringsAsFactors = !is.null(x$group)) {
   if(attr(x, "arstat")) stop("Cannot handle arrays of statistics!")
  has_tables <- attr(x, "table")
  nam <- attr(x, "names")
  attributes(x) <- NULL # faster lapply
   if(has_tables) {
     r <- lapply(x, function(z) c(list(Class = strclp(z[[1L]]), Label = null2NA(z[[2L]])),
          unlist(`names<-`(lapply(z[names(z) != "Table"][-(1:2)], as.vector, "list"), NULL), recursive = FALSE)))
   } else {
     r <- lapply(x, function(z) c(list(Class = strclp(z[[1L]]), Label = null2NA(z[[2L]])),
          unlist(`names<-`(lapply(z[-(1:2)], as.vector, "list"), NULL), recursive = FALSE)))
   }
  names(r) <- nam
   r <- .Call(C_rbindlist, r, TRUE, TRUE, "Variable")
   if(allNA(r[["Label"]])) r[["Label"]] <- NULL
   attr(r, "row.names") <- .set_row_names(length(r[[1L]]))
   class(r) <- "data.frame"
   r
}
