# Note: don't change the order of these arguments !!!
scv <- function(x, v, r, set = FALSE, inv = FALSE) .Call(C_setcopyv, x, v, r, inv, set, FALSE)

# inspired by ?dplyr::recode
# Think about adopting this code for as_numeric_factor and as_character_factor
recode_num <- function(X, ..., default = NULL, missing = NULL, set = FALSE) {
  if(missing(...)) stop("recode_num requires arguments of the form: value = replacement")
  args <- list(...)
  nam <- as.numeric(names(args))
  # nzchar(names(args)) ... check non-empty names ? -> nah, this package is not for dummies
  if(anyNA(nam)) stop(paste("Non-numeric arguments:", paste(names(args)[is.na(nam)], collapse = ", ")))
  arglen <- length(args)
  missingl <- !is.null(missing)
  if(missingl && any(nam == missing))  warning(paste0("To improve performance missing values are replaced prior to recode, so this replaces all missing values with ",
                                               missing, " and those are then again replaced with ", args[[which(nam == missing)]], ". If this is not desired, call replace_NA after recode with missing = NULL."))
  if(arglen == 1L) {
    args <- args[[1L]]
    if(is.null(default)) {
      if(missingl) {
        repfun <- function(y) if(is.numeric(y)) {
          z <- scv(y, NA, missing, set) # y[is.na(y)] <- missing
          scv(z, nam, args, TRUE) # `[<-`(y, y == nam, value = args)
          } else y
      } else {
        repfun <- function(y) if(is.numeric(y)) scv(y, nam, args, set) else y # `[<-`(y, y == nam, value = args)
      }
    } else {
      nr <- if(is.atomic(X)) NROW(X) else fnrow2(X)
      if(missingl) {
        repfun <- function(y) if(is.numeric(y)) {
          z <- scv(y, NA, default, set, TRUE) # duplAttributes(alloc(default, nr), y)
          scv(z, NA, missing, TRUE) # z[is.na(y)] <- missing # could put behind -> better but inconsistent
          `[<-`(z, whichv(y, nam), value = args) # y == nam
        } else y
      } else {
        repfun <- function(y) if(is.numeric(y)) scv(scv(y, nam, default, set, TRUE), nam, args, TRUE) else y # `[<-`(duplAttributes(alloc(default, nr), y), y == nam, value = args)
      }
    }
  } else {
    seqarg <- seq_len(arglen)
    if(is.null(default)) {
      if(missingl) {
        repfun <- function(y) if(is.numeric(y)) {
          y <- scv(y, NA, missing, set) # y[is.na(y)] <- missing
          z <- y
          for(i in seqarg) z[whichv(y, nam[i])] <- args[[i]]
          z
        } else y
      } else {
        repfun <- function(y) if(is.numeric(y)) {
          z <- y
          for(i in seqarg) z[whichv(y, nam[i])] <- args[[i]]
          z
        } else y
      }
    } else {
      nr <- if(is.atomic(X)) NROW(X) else fnrow2(X)
      if(missingl) {
        repfun <- function(y) if(is.numeric(y)) {
          z <- scv(y, NA, default, set, TRUE) # duplAttributes(alloc(default, nr), y)
          scv(z, NA, missing, TRUE)     # z[is.na(y)] <- missing # could put behind -> better but inconsistent
          for(i in seqarg) z[whichv(y, nam[i])] <- args[[i]]
          z
        } else y
      } else {
        repfun <- function(y) if(is.numeric(y)) {
          z <- scv(y, nam[1L], default, set, TRUE) # duplAttributes(alloc(default, nr), y)
          scv(z, nam[1L], args[[1L]], TRUE)
          for(i in seqarg[-1L]) z[whichv(y, nam[i])] <- args[[i]]
          z
        } else y
      }
    }
  }
  if(is.list(X)) {
    res <- duplAttributes(lapply(unattrib(X), repfun), X)
    return(if(inherits(X, "data.table")) alc(res) else res)
  }
  if(!is.numeric(X)) stop("X needs to be numeric or a list")
  repfun(X)
}

recode_char <- function(X, ..., default = NULL, missing = NULL, regex = FALSE,
                        ignore.case = FALSE, fixed = FALSE, set = FALSE) {
  if(missing(...)) stop("recode_char requires arguments of the form: value = replacement")
  args <- list(...)
  nam <- names(args)
  arglen <- length(args)
  missingl <- !is.null(missing)
  if(missingl && any(nam == missing))  warning(paste0("To improve performance missing values are replaced prior to recode, so this replaces all missing values with ",
                                                      missing, " and those are then again replaced with ", args[[which(nam == missing)]], ". If this is not desired, call replace_NA after recode with missing = NULL."))
  if(regex) {
    if(arglen == 1L) {
      args <- args[[1L]]
      if(is.null(default)) {
        if(missingl) {
          repfun <- function(y) if(is.character(y)) {
             y <- scv(y, NA, missing, set)  # y[is.na(y)] <- missing
            `[<-`(y, grepl(nam, y, ignore.case, FALSE, fixed), value = args)
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) `[<-`(y, grepl(nam, y, ignore.case, FALSE, fixed), value = args) else y
        }
      } else {
        nr <- if(is.atomic(X)) NROW(X) else fnrow2(X)
        if(missingl) {
          repfun <- function(y) if(is.character(y)) {
            z <- scv(y, NA, default, set, TRUE)  # duplAttributes(alloc(default, nr), y)
            scv(z, NA, missing, TRUE) # z[is.na(y)] <- missing # could put behind -> better but inconsistent
            `[<-`(z, grepl(nam, y, ignore.case, FALSE, fixed), value = args)
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) `[<-`(duplAttributes(alloc(default, nr), y), grepl(nam, y, ignore.case, FALSE, fixed), value = args) else y
        }
      }
    } else {
      seqarg <- seq_len(arglen)
      if(is.null(default)) {
        if(missingl) {
          repfun <- function(y) if(is.character(y)) {
            y <- scv(y, NA, missing, set) # y[is.na(y)] <- missing
            z <- y
            for(i in seqarg) z[grepl(nam[i], y, ignore.case, FALSE, fixed)] <- args[[i]]
            z
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) {
            z <- y
            for(i in seqarg) z[grepl(nam[i], y, ignore.case, FALSE, fixed)] <- args[[i]]
            z
          } else y
        }
      } else {
        nr <- if(is.atomic(X)) NROW(X) else fnrow2(X)
        if(missingl) {
          repfun <- function(y) if(is.character(y)) {
            z <- scv(y, NA, default, set, TRUE)  # duplAttributes(alloc(default, nr), y)
            scv(z, NA, missing, TRUE)  # z[is.na(y)] <- missing # could put behind -> better but inconsistent
            for(i in seqarg) z[grepl(nam[i], y, ignore.case, FALSE, fixed)] <- args[[i]]
            z
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) {
            z <- duplAttributes(alloc(default, nr), y)
            for(i in seqarg) z[grepl(nam[i], y, ignore.case, FALSE, fixed)] <- args[[i]]
            z
          } else y
        }
      }
    }
  } else {
    if(arglen == 1L) {
      args <- args[[1L]]
      if(is.null(default)) {
        if(missingl) {
          repfun <- function(y) if(is.character(y)) {
            y <- scv(y, NA, missing, set) # y[is.na(y)] <- missing
            scv(y, nam, args, TRUE)  # `[<-`(y, y == nam, value = args)
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) scv(y, nam, args, set) else y # `[<-`(y, y == nam, value = args)
        }
      } else {
        nr <- if(is.atomic(X)) NROW(X) else fnrow2(X)
        if(missingl) {
          repfun <- function(y) if(is.character(y)) {
            z <- scv(y, NA, default, set, TRUE)  # duplAttributes(alloc(default, nr), y)
            scv(z, NA, missing, TRUE) #  z[is.na(y)] <- missing # could put behind -> better but inconsistent
            `[<-`(z, whichv(y, nam), value = args)
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) scv(scv(y, nam, default, set, TRUE), nam, args, TRUE) else y # `[<-`(duplAttributes(alloc(default, nr), y), y == nam, value = args)
        }
      }
    } else {
      seqarg <- seq_len(arglen)
      if(is.null(default)) {
        if(missingl) {
          repfun <- function(y) if(is.character(y)) {
            y <- scv(y, NA, missing, set) # y[is.na(y)] <- missing
            z <- y
            for(i in seqarg) z[whichv(y, nam[i])] <- args[[i]]
            z
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) {
            z <- y
            for(i in seqarg) z[whichv(y, nam[i])] <- args[[i]]
            z
          } else y
        }
      } else {
        nr <- if(is.atomic(X)) NROW(X) else fnrow2(X)
        if(missingl) {
          repfun <- function(y) if(is.character(y)) {
            z <- scv(y, NA, default, set, TRUE)  # duplAttributes(alloc(default, nr), y)
            scv(z, NA, missing, TRUE) # z[is.na(y)] <- missing # could put behind -> better but inconsistent
            for(i in seqarg) z[whichv(y, nam[i])] <- args[[i]]
            z
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) {
            z <- scv(y, nam[1L], default, set, TRUE)  # duplAttributes(alloc(default, nr), y)
            scv(z, nam[1L], args[[1L]], TRUE)
            for(i in seqarg[-1L]) z[whichv(y, nam[i])] <- args[[i]]
            z
          } else y
        }
      }
    }
  }
  if(is.list(X)) {
    res <- duplAttributes(lapply(unattrib(X), repfun), X)
    return(if(inherits(X, "data.table")) alc(res) else res)
  }
  if(!is.character(X)) stop("X needs to be character or a list")
  repfun(X)
}

replace_NA <- function(X, value = 0L, cols = NULL, set = FALSE) {
  if(set) {
    if(is.list(X)) {
      if(is.null(cols)) {
        lapply(unattrib(X), scv, NA, value, TRUE)
      } else if(is.function(cols)) {
        lapply(unattrib(X), function(y) if(cols(y)) scv(y, NA, value, TRUE) else y)
      } else {
        cols <- cols2int(cols, X, attr(X, "names"), FALSE)
        lapply(unattrib(X)[cols], scv, NA, value, TRUE)
      }
    } else scv(X, NA, value, TRUE) # `[<-`(X, is.na(X), value = value)
    return(invisible(X))
  }
  if(is.list(X)) {
    if(is.null(cols)) return(condalc(duplAttributes(lapply(unattrib(X), scv, NA, value), X), inherits(X, "data.table"))) # function(y) `[<-`(y, is.na(y), value = value)
    if(is.function(cols)) return(condalc(duplAttributes(lapply(unattrib(X),
      function(y) if(cols(y)) scv(y, NA, value) else y), X), inherits(X, "data.table")))
    clx <- oldClass(X)
    oldClass(X) <- NULL
    cols <- cols2int(cols, X, names(X), FALSE)
    X[cols] <- lapply(unattrib(X[cols]), scv, NA, value) #  function(y) `[<-`(y, is.na(y), value = value)
    return(condalc(`oldClass<-`(X, clx), any(clx == "data.table")))
  }
  scv(X, NA, value) # `[<-`(X, is.na(X), value = value)
}

# Remove Inf (Infinity) and NaN (Not a number) from vectors or data frames:
replace_Inf <- function(X, value = NA, replace.nan = FALSE) {
  if(is.list(X)) {
    # if(!inherits(X, "data.frame")) stop("replace_non_finite only works with atomic objects or data.frames")
    res <- duplAttributes(lapply(unattrib(X),
             if(replace.nan) (function(y) if(is.numeric(y)) `[<-`(y, is.infinite(y) | is.nan(y), value = value) else y) else
                             (function(y) if(is.numeric(y)) `[<-`(y, is.infinite(y), value = value) else y)), X)
    return(if(inherits(X, "data.table")) alc(res) else res)
  }
  if(!is.numeric(X)) stop("Infinite values can only be replaced in numeric objects!")
  if(replace.nan) return(`[<-`(X, is.infinite(X) | is.nan(X), value = value)) #  !is.finite(X) also replaces NA
  `[<-`(X, is.infinite(X), value = value)
}

# replace_non_finite <- function(X, value = NA, replace.nan = TRUE) {
#   .Deprecated("replace_Inf")
#   replace_Inf(X, value, replace.nan)
# }

replace_outliers <- function(X, limits, value = NA, single.limit = c("SDs", "min", "max", "overall_SDs")) {
  ll <- length(limits)
  if(lg1 <- ll > 1L) {
    if(ll > 2L) stop("length(limits) must be 1 or 2")
    l1 <- limits[1L]
    l2 <- limits[2L]
  }
  if(is.list(X)) {
    # if(!inherits(X, "data.frame")) stop("replace_outliers only works with atomic objects or data.frames")
    if(lg1) {
      res <- duplAttributes(lapply(unattrib(X), function(y) if(is.numeric(y)) `[<-`(y, y < l1 | y > l2, value = value) else y), X) # could use data.table::between -> but it seems not faster !
    } else {
      res <- switch(single.limit[1L], # Allows grouped scaling if X is a grouped_df, but requires extra memory equal to X ... extra argument gSDs ?
           SDs = {
             if(inherits(X, c("grouped_df", "pdata.frame"))) {
              num <- .Call(C_vtypes, X, 1L) # vapply(unattrib(X), is.numeric, TRUE)
              num <- if(inherits(X, "grouped_df")) num & !fgroup_vars(X, "logical") else
                      num & attr(findex(X), "names") %!in% attr(X, "names")
              clx <- oldClass(X)
              STDXnum <- fscale(fcolsubset(X, num))
              oldClass(X) <- NULL
              X[num] <- mapply(function(z, y) `[<-`(z, abs(y) > limits, value = value), unattrib(X[num]), unattrib(STDXnum), SIMPLIFY = FALSE)
              `oldClass<-`(X, clx)
             } else duplAttributes(lapply(unattrib(X), function(y) if(is.numeric(y)) `[<-`(y, abs(fscaleCpp(y)) > limits, value = value) else y), X)
           },
           min = duplAttributes(lapply(unattrib(X), function(y) if(is.numeric(y)) `[<-`(y, y < limits, value = value) else y), X),
           max = duplAttributes(lapply(unattrib(X), function(y) if(is.numeric(y)) `[<-`(y, y > limits, value = value) else y), X),
           overall_SDs = duplAttributes(lapply(unattrib(X), function(y) if(is.numeric(y)) `[<-`(y, abs(fscaleCpp(y)) > limits, value = value) else y), X),
           stop("Unknown single.limit option"))
    }
    return(if(inherits(res, "data.table")) alc(res) else res)
  }
  if(!is.numeric(X)) stop("Outliers can only be replaced in numeric objects!")
  if(lg1) return(`[<-`(X, X < l1 | X > l2, value = value))
  switch(single.limit[1L],
    SDs =, overall_SDs = `[<-`(X, abs(fscale(X)) > limits, value = value),
    min = `[<-`(X, X < limits, value = value),
    max = `[<-`(X, X > limits, value = value),
    stop("Unknown single.limit option"))
}



# pad or fpad? x is vector, matrix or data.frame
pad_atomic <- function(x, i, n, value) {
  ax <- attributes(x)
  tx <- typeof(x)
  if(typeof(value) != tx) value <- as.vector(value, tx)
  if(is.matrix(x)) {
    k <- dim(x)[2L]
    m <- .Call(C_alloc, value, n * k)  # matrix(value, n, k)
    dim(m) <- c(n, k)
    m[i, ] <- x
    if(length(ax) == 1L) return(m)
    ax[["dim"]] <- c(n, k)
    # Could also pad row-names? perhaps with names of i ??
    if(length(ax[["dimnames"]][[1L]])) ax[["dimnames"]] <- list(NULL, ax[["dimnames"]][[2L]])
    if(is.object(x)) ax[["class"]] <- NULL
    return(`attributes<-`(m, ax)) # fastest ??
  }
  r <- .Call(C_alloc, value, n) # matrix(value, n) # matrix is faster than rep_len !!!!
  r[i] <- x
  if(is.null(ax)) return(r)
  if(length(names(x))) {
    if(length(ax) == 1L) return(r)
    ax[["names"]] <- NULL
  }
  return(`attributes<-`(r, ax))
}

# microbenchmark::microbenchmark(x[-i] <- ri, x[i2] <- ri)
# Unit: milliseconds
# expr       min       lq     mean   median       uq       max neval cld
# x[-i] <- ri 255.16654 420.7083 491.7369 446.0340 476.3324 1290.7396   100   b
# x[i2] <- ri  80.18755 136.8012 157.0027 146.8156 166.7158  311.5526   100  a
# microbenchmark::microbenchmark(seq_along(x)[-i])
# Unit: milliseconds
# expr      min       lq     mean   median       uq      max neval
# seq_along(x)[-i] 506.0745 541.7975 605.0245 567.8115 585.8384 1341.035   100

pad <- function(X, i, value = NA, method = c("auto", "xpos", "vpos")) { # 1 - i is same length as X, fill missing, 2 - i is positive: insert missing values in positions
  ilog <- is.logical(i)
  ineg <- i[1L] < 0L
  n <- if(is.list(X)) length(.subset2(X, 1L)) else if(is.matrix(X)) dim(X)[1L] else length(X)
  xpos <- switch(method[1L], auto = if(ilog) bsum(i) == n else if(ineg) FALSE else length(i) == n,
                 xpos = TRUE, vpos = FALSE, stop("Unknown method: ", method[1L]))
  n <- if(ilog) length(i) else if(xpos && !ineg) bmax(i) else n + length(i)
  if(is.atomic(X)) return(pad_atomic(X, if(xpos || ineg) i else if(ilog) !i else -i, n, value))
  if(!is.list(X)) stop("X must be atomic or a list")
  if(ilog) {
    i <- if(xpos) which(i) else whichv(i, FALSE)
  } else if(!xpos) {
    i <- seq_len(n)[if(ineg) i else -i]
  }
  ax <- attributes(X)
  attributes(X) <- NULL
  res <- lapply(X, pad_atomic, i, n, value)
  if(length(ax[["row.names"]])) ax[["row.names"]] <- .set_row_names(n)
  return(condalcSA(res, ax, any(ax[["class"]] == "data.table")))
}

# Something like this already exists?? -> should work with lists as well...


# # Previous version of Recode (Until collapse 1.1.0), Now depreciated in favor or recode_num and recode_char
# comp <- function(x, val) do.call(cbind, lapply(x, `==`, val))
# comp_grepl <- function(x, val) do.call(cbind, lapply(x, function(y) grepl(val, y)))
#
# Recode <- function(X, ..., copy = FALSE, reserve.na.nan = TRUE, regex = FALSE) {
#   .Deprecated("recode_num, recode_char")
#   if(missing(...)) stop("Recode requires arguments of the form: value = replacement")
#   if(is.list(X) && !inherits(X, "data.frame")) stop("Recode only works with atomic objects or data.frames")
#   args <- list(...)
#   nam <- names(args)
#
#   if(reserve.na.nan && any(wm <- !is.na(ma <- match(c("NaN","NA"), nam)))) { # more efficient way to code this ??
#     if(wm[1L]) {  # note: does not give multiple-matching error !!
#       if(is.list(X)) {
#         X[do.call(cbind, lapply(X, is.nan))] <- args[[ma[1L]]]
#       } else X[is.nan(X)] <- args[[ma[1L]]]
#     }
#     if(wm[2L]) X[is.na(X)] <- args[[ma[2L]]] # is.na already accounts for NaN, so this needs to the the order!!
#     if(bsum(wm) == length(args)) return(X)
#     args <- args[-ma[wm]]
#     nam <- names(args)
#   }
#
#   arglen <- length(args)
#   onearg <- arglen == 1L
#
#   if(regex) {
#     if(is.list(X)) {
#       if(onearg) {
#         X[comp_grepl(X, nam)] <- args[[1L]]
#       } else if(copy) {
#         for (i in seq_along(X)) {
#           Xi <- X[[i]]
#           for (j in seq_len(arglen)) X[[i]][grepl(nam[j], Xi)] <- args[[j]]
#         }
#       } else for (j in seq_len(arglen)) X[comp_grepl(X, nam[j])] <- args[[j]]
#     } else {
#       if(onearg) X[grepl(nam, X)] <- args[[1L]] else if(copy) {
#         Y <- X
#         for (j in seq_len(arglen)) X[grepl(nam[j], Y)] <- args[[j]]
#       } else for (j in seq_len(arglen)) X[grepl(nam[j], X)] <- args[[j]]
#     }
#   } else {
#     if(is.list(X)) {
#       oldopts <- options(warn = -1) # faster than suppressWarnings !!
#       on.exit(options(oldopts))
#       numnam <- as.numeric(nam) # suppressWarnings(numnam <- as.numeric(nam)) -> this line takes most time !!
#       if(onearg && !is.na(numnam)) nam <- numnam else {
#         nam <- as.vector(nam, "list")
#         if(any(numnaml <- !is.na(numnam))) nam[numnaml] <- numnam[numnaml]
#       }
#       if(onearg) X[comp(X, nam)] <- args[[1L]] else if(copy) {
#         Y <- X
#         for (j in seq_len(arglen)) X[comp(Y, nam[[j]])] <- args[[j]]
#       } else for (j in seq_len(arglen)) X[comp(X, nam[[j]])] <- args[[j]]
#     } else {
#       if(!is.character(X)) nam <- as.numeric(nam)
#       if(onearg) X[X == nam] <- args[[1L]] else if(copy) {
#         Y <- X
#         for (j in seq_len(arglen)) X[Y == nam[j]] <- args[[j]]
#       } else for (j in seq_len(arglen)) X[X == nam[j]] <- args[[j]]
#     }
#   }
#   return(X)
# }


# Experimental:
# recode_num2 is slightly slower than recode_num above .... (because 2 times apply) but almost identical..
# recode_num2 <- function(X, ...) { # , regex = FALSE # , default = NULL
#   if(missing(...)) stop("Recode requires arguments of the form: value = replacement")
#   args <- list(...)
#   nam <- as.numeric(names(args))
#   if(anyNA(nam)) stop(paste("Non-numeric arguments:", paste(names(args)[is.na(nam)], collapse = ", ")))
#   arglen <- length(args)
#
#   if(arglen == 1L) {
#     args <- args[[1L]]
#     repfun <- function(y) `[<-`(y, y == nam, value = args)
#   } else {
#     seqarg <- seq_len(arglen)
#     repfun <- function(y) {
#       z <- y
#       for(i in seqarg) z[y == nam[i]] <- args[[i]]
#       z
#     }
#   }
#   if(is.list(X)) {
#     num <- vapply(unattrib(X), is.numeric, TRUE)
#     clx <- class(X)
#     class(X) <- NULL
#     X[num] <- lapply(X[num], repfun)
#     return(`oldClass<-`(X, clx))
#   }
#   if(!is.numeric(X)) stop("X needs to be numeric or a list")
#   return(repfun(X))
# }
#
#
#
# # possibly even faster by converting df to matrix and back for multiple comp ??
# Recode2 <- function(X, ...) { # , regex = FALSE # , default = NULL
#   if(missing(...)) stop("Recode requires arguments of the form: value = replacement")
#   args <- list(...)
#   nam <- names(args)
#   arglen <- length(args)
#   onearg <- arglen == 1L
#
#   if(onearg) {
#     args <- args[[1L]]
#     repfun <- function(y) {
#       if(is.numeric(y)) {
#         v <- as.numeric(nam)
#         if(is.na(v)) {
#           warning("Trying to replace a non-numeric expressiom in a numeric column: Skipping this column")
#           return(y)
#         }
#         `[<-`(y, y == v, value = args)
#       } else `[<-`(y, y == nam, value = args)
#     }
#   } else {
#     seqarg <- seq_len(arglen)
#     repfun <- function(y) {
#       z <- y
#       if(is.numeric(y)) {
#         if(is.null(v)) {
#           v <<- as.numeric(nam)
#           if(anyNA(v)) {
#             warning("Trying to replace a non-numeric expressiom in a numeric column: Skipping those expressions")
#             v <- na_rm(v)
#             if(!length(v)) return(y)
#           }
#         }
#         for(i in seqarg) z[y == v[i]] <- args[[i]]
#       } else {
#         for(i in seqarg) z[y == nam[i]] <- args[[i]]
#       }
#       z
#     }
#   }
#
#   if(is.list(X)) {
#     return(duplAttributes(lapply(unattrib(X), repfun), X))
#   }
#   if(!is.character(X)) storage.mode(nam) <- storage.mode(X)
#   return(repfun(X))
# }
#
#
# rec1 <- function(x, ...) {
#   args <- list(...)
#   nam <- as.numeric(names(args))
#   ax <- attributes(x)
#   attributes(x) <- NULL
#   for(i in seq_along(args)) {
#     ni <- nam[i] # Faster!
#     argi <- args[[i]]
#     x <- lapply(x, function(y) `[<-`(y, y == ni, value = argi))
#   }
#   setAttributes(x, ax)
# }
# # More memory efficient !!
# rec2 <- function(x, ...) {
#   args <- list(...)
#   nam <- as.numeric(names(args))
#   if(length(args) > 1L) repfun <- function(y) {
#     for(i in seq_along(args)) y[y == nam[i]] <- args[[i]]
#     y
#   } else repfun <- function(y) `[<-`(y, y == nam, value = args[[1L]])
#   dapply(x, repfun)
# }
# # Even more memory efficient and faster !!
# rec3 <- function(x, ...) {
#   args <- list(...)
#   nam <- as.numeric(names(args))
#   if(length(args) > 1L) repfun <- function(y) {
#     z <- y
#     for(i in seq_along(args)) z[y == nam[i]] <- args[[i]]
#     z
#   } else repfun <- function(y) `[<-`(y, y == nam, value = args[[1L]])
#   dapply(x, repfun)
# }
#
#
# tf <- function(x, ...) list(...)
#
# # profvis({
# #   v <- c("a","b","c")
# #   Recode(v, a = "b")
# # })




# Previous Versions (until collapse 1.1.0)

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
