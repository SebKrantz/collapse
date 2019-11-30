# Note: do.call(cbind, lapply(X, FUN)) is always faster than for loop !!

comp <- function(x, val) do.call(cbind, lapply(x, `==`, val))
comp_grepl <- function(x, val) do.call(cbind, lapply(x, function(y) grepl(val, y)))

# possibly even faster by converting df to matrix and back for multiple comp ??
Recode <- function(X, ..., copy = FALSE, reserve.na.nan = TRUE, regex = FALSE) {
  if(missing(...)) stop("Recode requires arguments of the form: value = replacement")
  if(is.list(X) && !inherits(X, "data.frame")) stop("Recode only works with atomic objects or data.frames")
  args <- list(...)
  nam <- names(args)

  if(reserve.na.nan && any(wm <- !is.na(ma <- match(c("NaN","NA"), nam)))) { # more efficient way to code this ??
    if(wm[1L]) {  # note: does not give multiple-matching error !!
      if(is.list(X)) {
        X[do.call(cbind, lapply(X, is.nan))] <- args[[ma[1L]]]
      } else X[is.nan(X)] <- args[[ma[1L]]]
    }
    if(wm[2L]) X[is.na(X)] <- args[[ma[2L]]] # is.na already accounts for NaN, so this needs to the the order!!
    if(sum(wm) == length(args)) return(X)
    args <- args[-ma[wm]]
    nam <- names(args)
  }

  arglen <- length(args)
  onearg <- arglen == 1L

  if(regex) {
    if(is.list(X)) {
      if(onearg) {
        X[comp_grepl(X, nam)] <- args[[1L]]
      } else if(copy) {
        for (i in seq_along(X)) {
          Xi <- X[[i]]
          for (j in seq_len(arglen)) X[[i]][grepl(nam[j], Xi)] <- args[[j]]
        }
      } else for (j in seq_len(arglen)) X[comp_grepl(X, nam[j])] <- args[[j]]
    } else {
      if(onearg) X[grepl(nam, X)] <- args[[1L]] else if(copy) {
        Y <- X
        for (j in seq_len(arglen)) X[grepl(nam[j], Y)] <- args[[j]]
      } else for (j in seq_len(arglen)) X[grepl(nam[j], X)] <- args[[j]]
    }
  } else {
    if(is.list(X)) {
      options(warn = -1) # faster than suppressWarnings !!
      numnam <- as.numeric(nam) # suppressWarnings(numnam <- as.numeric(nam)) -> this line takes most time !!
      if(onearg && !is.na(numnam)) nam <- numnam else {
        nam <- as.vector(nam, "list")
        if(any(numnaml <- !is.na(numnam))) nam[numnaml] <- numnam[numnaml]
      }
      options(warn = 1)
      if(onearg) X[comp(X, nam)] <- args[[1L]] else if(copy) {
        Y <- X
        for (j in seq_len(arglen)) X[comp(Y, nam[[j]])] <- args[[j]]
      } else for (j in seq_len(arglen)) X[comp(X, nam[[j]])] <- args[[j]]
    } else {
      if(!is.character(X)) nam <- as.numeric(nam)
      if(onearg) X[X == nam] <- args[[1L]] else if(copy) {
        Y <- X
        for (j in seq_len(arglen)) X[Y == nam[j]] <- args[[j]]
      } else for (j in seq_len(arglen)) X[X == nam[j]] <- args[[j]]
    }
  }
  return(X)
}

# profvis({
#   v <- c("a","b","c")
#   Recode(v, a = "b")
# })

# Remove Inf (Infinity) and NaN (Not a number) from vectors or data frames:
replace_non_finite <- function(X, value = NA, replace.nan = TRUE) {
  if(is.list(X)) {
    if(!inherits(X, "data.frame")) stop("replace_non_finite only works with atomic objects or data.frames")
    if(replace.nan) {
      infornan <- function(x) is.infinite(x) | is.nan(x)
      X[do.call(cbind, lapply(X, infornan))] <- value       # for (i in seq_along(X)) X[[i]][!is.finite(X[[i]])] <- value #   is.infinite(X[[i]]) | is.nan(X[[i]])
    } else X[do.call(cbind, lapply(X, is.infinite))] <- value       # for (i in seq_along(X)) X[[i]][is.infinite(X[[i]])] <- value
  } else {
    if(replace.nan) {
      X[is.infinite(X) | is.nan(X)] <- value #  !is.finite(X) also replaces NA
    } else X[is.infinite(X)] <- value
  }
  return(X)
}


replace_outliers <- function(X, limits, value = NA, single.limit = c("SDs","min","max")) {
  if(is.list(X) && !inherits(X, "data.frame")) stop("replace_outliers only works with atomic objects or data.frames")
  if(length(limits) > 1L) return(replace(X, if(is.atomic(X)) X < limits[1L] | X > limits[2L] else do.call(cbind, lapply(X, function(x) x < limits[1L] | x > limits[2L])), value)) # could use data.table::between -> but not faster !!
  switch(single.limit[1L],
    SDs = replace(X, if(is.atomic(X)) abs(fscale(X)) > limits else abs(do.call(cbind, fscale(X))) > limits , value), # this is slightly faster, but doesn't allow grouped scaling: do.call(cbind, lapply(X, function(x) abs(fscale.default(x)) > limits))
    min = replace(X, if(is.atomic(X)) X < limits else do.call(cbind, lapply(X, `<`, limits)), value),
    max = replace(X, if(is.atomic(X)) X > limits else do.call(cbind, lapply(X, `>`, limits)), value),
    stop("Unknown single.limit option"))
}


# replace_outliers <- function(X, limits, value = NA, single.limit = c("SDs","Min","Max")) {
#   if(length(limits) > 1L) {
#     betw <- function(x, lim) x < lim[1L] | x > lim[2L]
#     if(is.list(X))  {
#       for (i in seq_along(X)) X[[i]][betw(X[[i]], limits)] <- value
#     } else X[betw(X, limits)] <- value
#     return(X)
#   } else {
#     switch(single.limit[1L],
#            SDs = {
#              if(is.list(X))  {
#                gsd <- function(x, lim) abs(x) > lim*fsdCpp(x)
#                for (i in seq_along(X)) X[[i]][gsd(X[[i]], limits)] <- value
#              } else X[abs(X) > lim*fsd(X)] <- value
#              return(X)
#            },
#            Min = {
#              if(is.list(X))  {
#                for (i in seq_along(X)) X[[i]][X[[i]] < limits] <- value
#              } else X[X < limits] <- value
#              return(X)
#            },
#            Max = {
#              if(is.list(X))  {
#                for (i in seq_along(X)) X[[i]][X[[i]] > limits] <- value
#              } else X[X > limits] <- value
#              return(X)
#            },
#            stop("Unknown single.limit option"))
#   }
# }
