

# Could make label attribute optional !!
descr <- function(X, Ndistinct = TRUE, higher = TRUE, table = TRUE,
                  Qprobs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99),
                  label.attr = 'label', ...) {
  nam <- deparse(substitute(X))

  armat <- function(x, y) c(x[1L], Ndist = y, x[-1L])

  dotsok <- if(!missing(...)) names(substitute(c(...))[-1L]) %!in% c('pid','g') else TRUE

  numstats <- if(Ndistinct && dotsok) function(x, ...) armat(qsu.default(x, higher = higher, ...), fNdistinctCpp(x)) else function(x, ...) qsu.default(x, higher = higher, ...)

  descrnum <- if(is.numeric(Qprobs)) function(x, ...) list(Class = class(x), Label = attr(x, label.attr), Stats = numstats(x, ...),
                                                          Quant = quantile(x, probs = Qprobs, na.rm = TRUE)) else
                                         function(x, ...) list(Class = class(x), Label = attr(x, label.attr), Stats = numstats(x, ...))
  # Could make this more efficient ???
  descrcat <- function(x, tab = table) if(tab) list(Class = class(x), Label = attr(x, label.attr),
                                                    Stats = if(Ndistinct) c(N = fNobsCpp(x), Ndist = fNdistinctCpp(x)) else `names<-`(fNobsCpp(x), 'Nobs'),
                                                    Table = table(x)) else
                                                      list(Class = class(x), Label = attr(x, label.attr),
                                                           Stats = if(Ndistinct) c(N = fNobsCpp(x), Ndist = fNdistinctCpp(x)) else `names<-`(fNobsCpp(x), 'Nobs'))
  class(X) <- NULL
  res <- vector('list', length(X))
  num <- vapply(X, is.numeric, TRUE, USE.NAMES = FALSE)
  res[num] <- lapply(X[num], descrnum, ...)
  if(!all(num)) {
    date <- vapply(X, is.Date, TRUE, USE.NAMES = FALSE)
    if(any(date)) {
      res[date] <- lapply(X[date], descrcat, FALSE)
      cat <- !(num | date)
    } else cat <- !num
    res[cat] <- lapply(X[cat], descrcat)
  }
  attributes(res) <- list(names = names(X), name = nam, N = length(X[[1L]]),
                          arstat = !dotsok, class = "descr")
  return(res)
}

print.descr <- function(x, n = 6, perc = TRUE, summary = TRUE, ...) {
  tabbperc <- if(perc) function(t) `names<-`(paste0(t," (",round(t/sum(t)*100,1),"%)"), names(t)) else
    function(t) `names<-`(unclass(t), names(t))
  w <- paste(rep("-", .Options$width), collapse = "")
  nam <- names(x)
  arstat <- attr(x, "arstat")
  cat('Dataset: ', attr(x,"name"),', ',length(x), ' Variables, N = ', attr(x, "N"), "\n", sep = "")
  cat(w, "\n", sep = "")
  for(i in seq_along(x)) {
    xi <- x[[i]]
    namxi <- names(xi)
    cat(nam[i]," (",xi[[1L]],"): ",xi[[2L]], "\n", sep = "")
    cat(namxi[3L], ": \n", sep = "")
    print.qsu(xi[[3L]])
    if(length(xi) > 3L) {
      if(arstat) cat("\n")
      cat(namxi[4L], ": \n", sep = "")
      if(namxi[4L] == "Table") {
        if(length(xi[[4L]]) <= 2*n) print.table(tabbperc(xi[[4L]])) else {
          t <- unclass(xi[[4L]])
          lt <- length(t)
          print.table(tabbperc(t[seq_len(n)]))
          cat("  ---\n")
          print.table(tabbperc(t[seq(lt-n, lt)]))
          if(summary) {
            cat("\nSummary of Table: \n")
            print.summaryDefault(summary.default(t))
          }
        }
      } else print.qsu(xi[[4L]])
    }
    cat(w, "\n", sep = "") # More compressed -> better !
    # cat("\n", w, "\n", sep = "")
  }
  invisible(x)
}


# Note: This does not work for array stats (using g or pid.. )
as.data.frame.descr <- function(x, ...) {
   if(attr(x, "arstat")) stop("Cannot handle arrays of statistics!")
   r <- lapply(x, function(z) unlist(`names<-`(lapply(z[names(z) != "Table"], as.vector, "list"), NULL), recursive = FALSE))
   r <- .Call(C_rbindlist, r, TRUE, TRUE, "Variable")
   if(names(r)[3L] == "") names(r)[2:3] <- c("Class", "Label") else names(r)[2L] <- "Class"
   attr(r, "row.names") <- .set_row_names(length(r[[1L]]))
   class(r) <- "data.frame"
   return(r)
}
