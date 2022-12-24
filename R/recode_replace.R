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
      nr <- if(is.atomic(X)) NROW(X) else fnrow(X)
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
      nr <- if(is.atomic(X)) NROW(X) else fnrow(X)
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
        nr <- if(is.atomic(X)) NROW(X) else fnrow(X)
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
        nr <- if(is.atomic(X)) NROW(X) else fnrow(X)
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
        nr <- if(is.atomic(X)) NROW(X) else fnrow(X)
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
        nr <- if(is.atomic(X)) NROW(X) else fnrow(X)
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
              X[num] <- .mapply(function(z, y) `[<-`(z, abs(y) > limits, value = value), list(unattrib(X[num]), unattrib(STDXnum)), NULL)
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
  n <- if(is.list(X) || is.matrix(X)) fnrow(X) else length(X)
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
