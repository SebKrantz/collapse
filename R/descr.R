

# Could make label attribute optional !
descr <- function(X, Ndistinct = TRUE, higher = TRUE, table = TRUE,
                  Qprobs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), cols = NULL,
                  label.attr = 'label', ...) {
  nam <- l1orlst(as.character(substitute(X)))

  armat <- function(x, y) c(x[1L], Ndist = y, x[-1L])
  natrm <- function(x) if(is.na(names(x)[length(x)])) x[-length(x)] else x # Remove NA from table !

  dotsok <- if(!missing(...)) names(substitute(c(...))[-1L]) %!in% c('pid','g') else TRUE

  numstats <- if(Ndistinct && dotsok) function(x, ...) armat(qsu.default(x, higher = higher, ...), fndistinctCpp(x)) else function(x, ...) qsu.default(x, higher = higher, ...)

  descrnum <- if(is.numeric(Qprobs)) function(x, ...) list(Class = class(x), Label = attr(x, label.attr), Stats = numstats(x, ...),
                                                           Quant = quantile(na_rm(x), probs = Qprobs)) else
                                         function(x, ...) list(Class = class(x), Label = attr(x, label.attr), Stats = numstats(x, ...))
  # Could make this more efficient ?
  descrcat <- function(x, tab = table) if(tab) list(Class = class(x), Label = attr(x, label.attr),
                                                    Stats = if(Ndistinct) c(N = fnobsC(x), Ndist = fndistinctCpp(x)) else `names<-`(fnobsC(x), 'Nobs'),
                                                    Table = natrm(fnobs.default(x, x))) else # table(x). fnobs is a lot Faster, but includes NA as level !
                                                      list(Class = class(x), Label = attr(x, label.attr),
                                                           Stats = if(Ndistinct) c(N = fnobsC(x), Ndist = fndistinctCpp(x)) else `names<-`(fnobsC(x), 'Nobs'))

  descrdate <- function(x) list(Class = class(x), Label = attr(x, label.attr),
                                Stats = c(if(Ndistinct) c(N = fnobsC(x), Ndist = fndistinctCpp(x)) else `names<-`(fnobsC(x), 'Nobs'),
                                          `names<-`(range(x, na.rm = TRUE), c("Min", "Max"))))


  if(is.list(X)) {
    is_sf <- inherits(X, "sf")
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

print.descr <- function(x, n = 7, perc = TRUE, digits = 2, t.table = TRUE, summary = TRUE, ...) {
  w <- paste(rep("-", .Options$width), collapse = "")
  nam <- names(x)
  arstat <- attr(x, "arstat")
  cb <- function(...) if(t.table) cbind(...) else formatC(rbind(...), drop0trailing = TRUE)
  ct <- function(z) if(t.table) cbind(Freq = z) else z
  cat('Dataset: ', attr(x,"name"),', ',length(x), ' Variables, N = ', attr(x, "N"), "\n", sep = "")
  cat(w, "\n", sep = "")
  for(i in seq_along(x)) {
    xi <- x[[i]]
    namxi <- names(xi)
    cat(nam[i]," (",strclp(xi[[1L]]),"): ",xi[[2L]], "\n", sep = "")
    cat(namxi[3L], ": \n", sep = "")
    print.qsu(xi[[3L]], digits)
    if(length(xi) > 3L) {
      if(arstat) cat("\n")
      cat(namxi[4L], ": \n", sep = "")
      if(namxi[4L] == "Table") {
        t <- unclass(xi[[4L]])
        if(length(t) <= 2*n) {
          if(perc) print.default(cb(Freq = t, Perc = round(t/bsum(t)*100, digits)), right = TRUE, print.gap = 2, quote = FALSE) else
            print.table(ct(t))
        } else {
          lt <- length(t)
          t1 <- t[seq_len(n)]
          t2 <- t[seq(lt-n, lt)]
          if(perc) {
            st <- bsum(t)
            print.default(cb(Freq = t1, Perc = round(t1/st*100, digits)), right = TRUE, print.gap = 2, quote = FALSE)
            cat("  ---\n")
            print.default(cb(Freq = t2, Perc = round(t2/st*100, digits)), right = TRUE, print.gap = 2, quote = FALSE)
          } else {
            print.table(ct(t1))
            cat("  ---\n")
            print.table(ct(t2))
          }
          if(summary) {
            cat("\nSummary of Table: \n")
            print.summaryDefault(summary.default(t), digits)
          }
        }
      } else print.qsu(xi[[4L]], digits)
    }
    cat(w, "\n", sep = "") # More compressed -> better !
    # cat("\n", w, "\n", sep = "")
  }
  invisible(x)
}
# Not pasteclass !!

# Note: This does not work for array stats (using g or pid.. )
as.data.frame.descr <- function(x, ...) {
   if(attr(x, "arstat")) stop("Cannot handle arrays of statistics!")
   if(attr(x, "table")) {
     r <- lapply(x, function(z) c(list(Class = strclp(z[[1L]]), Label = null2NA(z[[2L]])),
          unlist(`names<-`(lapply(z[names(z) != "Table"][-(1:2)], as.vector, "list"), NULL), recursive = FALSE)))
   } else {
     r <- lapply(x, function(z) c(list(Class = strclp(z[[1L]]), Label = null2NA(z[[2L]])),
          unlist(`names<-`(lapply(z[-(1:2)], as.vector, "list"), NULL), recursive = FALSE)))
   }
   r <- .Call(C_rbindlist, r, TRUE, TRUE, "Variable")
   if(allNA(r[["Label"]])) r[["Label"]] <- NULL
   attr(r, "row.names") <- .set_row_names(length(r[[1L]]))
   class(r) <- "data.frame"
   r
}
