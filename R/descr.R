
fsorttable <- function(x, srt) {
  if(is.factor(x)) {
    lev <- attr(x, "levels")
    t <- .Call(C_fwtabulate, x, NULL, length(lev), TRUE) # tabulate(x, nbins = length(lev)) # skips missing values !!
    names(t) <- lev
    sorted <- TRUE
  } else {
    sorted <- FALSE
    g <- .Call(C_group, x, TRUE, TRUE)
    t <- attr(g, "group.sizes")
    nam <- Csv(x, attr(g, "starts"))
    names(t) <- nam
    if(anyNA(nam)) t <- t[-whichNA(nam)] # TODO: speef up using fwtabulate and groupat?
  }
  switch(srt,
    value = if(sorted || attr(o <- forder.int(names(t)), "sorted")) t else t[o],
    # "quick" sort seems best, based on multiple datasets, but "radix" (second best) keeps ties in order...
    # sort.int(t, method = "radix", decreasing = TRUE, na.last = TRUE)
    freq = if(attr(o <- forder.int(t, decreasing = TRUE), "sorted")) t else t[o],
    none = t,
    stop("sort.table must be one of 'value', 'freq' or 'none'"))
}


descr <- function(X, Ndistinct = TRUE, higher = TRUE, table = TRUE, sort.table = "freq",
                  Qprobs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), cols = NULL,
                  label.attr = 'label', ...) {
  nam <- l1orlst(as.character(substitute(X)))

  armat <- function(x, y) c(x[1L], Ndist = y, x[-1L])
  # natrm <- function(x) if(is.na(names(x)[length(x)])) x[-length(x)] else x # Remove NA from table !

  dotsok <- if(missing(...)) TRUE else names(substitute(c(...))[-1L]) %!in% c('pid','g')

  numstats <- if(Ndistinct && dotsok) function(x, ...) armat(qsu.default(x, higher = higher, ...), fndistinctC(x)) else function(x, ...) qsu.default(x, higher = higher, ...)

  descrnum <- if(is.numeric(Qprobs)) function(x, ...) list(Class = class(x), Label = attr(x, label.attr), Stats = numstats(x, ...),
                                                           Quant = quantile(na_rm(x), probs = Qprobs)) else
                                     function(x, ...) list(Class = class(x), Label = attr(x, label.attr), Stats = numstats(x, ...))

  descrcat <- if(table) function(x) {
                                    tab <- fsorttable(x, sort.table)
                                    list(Class = class(x), Label = attr(x, label.attr),
                                         Stats = if(Ndistinct) c(N = fsumC(tab), Ndist = length(tab)) else `names<-`(fsumC(tab), 'Nobs'),
                                         Table = tab) # natrm(fnobs.default(x, x)) # table(x). fnobs is a lot Faster, but includes NA as level !
                        } else
                        function(x) list(Class = class(x), Label = attr(x, label.attr),
                                         Stats = if(Ndistinct) c(N = fnobsC(x), Ndist = fndistinctC(x)) else `names<-`(fnobsC(x), 'Nobs'))

  descrdate <- function(x) list(Class = class(x), Label = attr(x, label.attr),
                                Stats = c(if(Ndistinct) c(N = fnobsC(x), Ndist = fndistinctC(x)) else `names<-`(fnobsC(x), 'Nobs'),
                                          `names<-`(frange(x), c("Min", "Max"))))

  if(is.list(X)) {
    is_sf <- inherits(X, "sf")
    # if(inherits(X, "POSIXlt")) X <- list(X = as.POSIXct(X))
    if(inherits(X, "indexed_frame")) X <- unindex(X)
    class(X) <- NULL
    if(is_sf) X[[attr(X, "sf_column")]] <- NULL
  } else X <- unclass(qDF(X))
  if(length(cols)) X <- X[cols2int(cols, X, names(X), FALSE)]
  res <- vector('list', length(X))
  num <- .Call(C_vtypes, X, 1L) # vapply(unattrib(X), is.numeric, TRUE)
  res[num] <- lapply(X[num], descrnum, ...)
  if(!all(num)) {
    date <- vapply(unattrib(X), is_date, TRUE)
    if(any(date)) {
      res[date] <- lapply(X[date], descrdate)
      cat <- !(num | date)
    } else cat <- !num
    res[cat] <- lapply(X[cat], descrcat)
  }
  attributes(res) <- list(names = names(X), name = nam, N = length(X[[1L]]),
                          arstat = !dotsok, table = table, class = "descr")
  res
}

`[.descr` <- function(x, ...) copyMostAttributes(.subset(x, ...), x)

print.descr <- function(x, n = 14, perc = TRUE, digits = 2, t.table = TRUE, summary = TRUE, reverse = FALSE, stepwise = FALSE, ...) {
  oldClass(x) <- NULL
  w <- paste(rep("-", .Options$width), collapse = "")
  arstat <- attr(x, "arstat")
  DSname <- attr(x, "name")
  DSN <- attr(x, "N")
  cb <- function(...) if(t.table) cbind(...) else formatC(rbind(...), drop0trailing = TRUE)
  ct <- function(z) if(t.table) cbind(Freq = z) else z
  if(reverse) x <- rev.default(x) else {
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
  if(reverse) cat('Dataset: ', DSname,', ',length(x), ' Variables, N = ', DSN, "\n", sep = "")
  invisible(x)
}

# Note: This does not work for array stats (using g or pid.. )
as.data.frame.descr <- function(x, ...) {
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
