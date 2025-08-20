# Note: don't change the order of these arguments !!!
scv <- function(x, v, r, set = FALSE, inv = FALSE, vind1 = FALSE) .Call(C_setcopyv, x, v, r, inv, set, vind1)

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
                                               missing, " and those are then again replaced with ", args[[which(nam == missing)]], ". If this is not desired, call replace_na after recode with missing = NULL."))
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
          nas <- is.na(y)
          z <- scv(y, nas, missing, set, vind1 = TRUE)
          ind <- whichv(z, nam)
          scv(z, nas, default, TRUE, TRUE, vind1 = TRUE) # duplAttributes(alloc(default, nr), y)
          scv(z, ind, args, TRUE, vind1 = TRUE) # y == nam
        } else y
      } else {
        repfun <- function(y) if(is.numeric(y)) scv(scv(y, nam, default, set, TRUE), nam, args, TRUE) else y # `[<-`(duplAttributes(alloc(default, nr), y), y == nam, value = args)
      }
    }
  } else {
    seqarg <- seq_len(arglen)
    if(is.null(default)) {
      repfun <- function(y) if(is.numeric(y)) {
        if(missingl) y <- scv(y, NA, missing, set) # y[is.na(y)] <- missing
        else if(!set) y <- scv(y, 1L, y[1L], vind1 = TRUE) # copy
        z <- scv(y, 1L, y[1L], vind1 = TRUE) # copy
        for(i in seqarg) scv(y, whichv(z, nam[i]), args[[i]], TRUE, vind1 = TRUE)
        y
      } else y
      # repfun <- function(y) if(is.numeric(y)) {
      #   if(missingl) y <- scv(y, NA, missing, set) # y[is.na(y)] <- missing
      #   if(set) { # Note: not strictly the way this should work...
      #     for(i in seqarg) scv(y, nam[i], args[[i]], TRUE)
      #     return(y)
      #   }
      #   z <- scv(y, 1L, y[1L], vind1 = TRUE) # copy
      #   for(i in seqarg) scv(z, whichv(y, nam[i]), args[[i]], TRUE, vind1 = TRUE)
      #   z
      # } else y
    } else {
      nr <- if(is.atomic(X)) NROW(X) else fnrow(X)
      if(missingl) {
        repfun <- function(y) if(is.numeric(y)) {
          nas <- is.na(y)
          y <- scv(y, nas, missing, set, vind1 = TRUE)
          z <- scv(y, 1L, y[1L], vind1 = TRUE) # Copy
          scv(y, nas, default, TRUE, TRUE, vind1 = TRUE)
          for(i in seqarg) scv(y, whichv(z, nam[i]), args[[i]], TRUE, vind1 = TRUE)
          y
        } else y
      } else {
        repfun <- function(y) if(is.numeric(y)) {
          z <- scv(y, 1L, y[1L], vind1 = TRUE) # Copy
          y <- scv(y, nam[1L], default, set, TRUE) # duplAttributes(alloc(default, nr), y)
          scv(y, nam[1L], args[[1L]], TRUE)
          for(i in seqarg[-1L]) scv(y, whichv(z, nam[i]), args[[i]], TRUE, vind1 = TRUE)
          y
        } else y
      }
    }
  }

  if(is.list(X)) {
    if(set) {
      lapply(unattrib(X), repfun)
      return(invisible(X))
    }
    res <- duplAttributes(lapply(unattrib(X), repfun), X)
    return(if(inherits(X, "data.table")) alc(res) else res)
  }
  if(!is.numeric(X)) stop("X needs to be numeric or a list")
  res <- repfun(X)
  return(if(set) invisible(res) else res)
}

recode_char <- function(X, ..., default = NULL, missing = NULL, regex = FALSE,
                        ignore.case = FALSE, fixed = FALSE, set = FALSE) {
  if(missing(...)) stop("recode_char requires arguments of the form: value = replacement")
  args <- list(...)
  nam <- names(args)
  arglen <- length(args)
  missingl <- !is.null(missing)
  if(missingl && any(nam == missing))  warning(paste0("To improve performance missing values are replaced prior to recode, so this replaces all missing values with ",
                                                      missing, " and those are then again replaced with ", args[[which(nam == missing)]], ". If this is not desired, call replace_na after recode with missing = NULL."))
  if(regex) {
    if(arglen == 1L) {
      args <- args[[1L]]
      if(is.null(default)) {
        if(missingl) {
          repfun <- function(y) if(is.character(y)) {
            y <- scv(y, NA, missing, set)  # y[is.na(y)] <- missing
            scv(y, grepl(nam, y, ignore.case, FALSE, fixed), args, TRUE, vind1 = TRUE)
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) scv(y, grepl(nam, y, ignore.case, FALSE, fixed), args, set, vind1 = TRUE) else y
        }
      } else {
        nr <- if(is.atomic(X)) NROW(X) else fnrow(X)
        if(missingl) {
          repfun <- function(y) if(is.character(y)) {
            nas <- is.na(y)
            z <- scv(y, nas, missing, set, vind1 = TRUE)
            ind <- grepl(nam, z, ignore.case, FALSE, fixed)
            scv(z, nas, default, TRUE, TRUE, vind1 = TRUE) # duplAttributes(alloc(default, nr), y)
            scv(z, ind, args, TRUE, vind1 = TRUE)
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) {
            ind <- grepl(nam, y, ignore.case, FALSE, fixed)
            scv(scv(y, ind, default, set, TRUE, vind1 = TRUE), ind, args, TRUE, vind1 = TRUE)
          } else y
        }
      }
    } else {
      seqarg <- seq_len(arglen)
      if(is.null(default)) {
        repfun <- function(y) if(is.character(y)) {
          if(missingl) y <- scv(y, NA, missing, set) # y[is.na(y)] <- missing
          else if(!set) y <- scv(y, 1L, y[1L], vind1 = TRUE) # copy
          z <- scv(y, 1L, y[1L], vind1 = TRUE) # copy
          for(i in seqarg) scv(y, grepl(nam[i], z, ignore.case, FALSE, fixed), args[[i]], TRUE, vind1 = TRUE)
          y
        } else y
      } else {
        nr <- if(is.atomic(X)) NROW(X) else fnrow(X)
        if(missingl) {
          repfun <- function(y) if(is.character(y)) {
            nas <- is.na(y)
            y <- scv(y, nas, missing, set, vind1 = TRUE)
            z <- scv(y, 1L, y[1L], vind1 = TRUE) # Copy
            scv(y, nas, default, TRUE, TRUE, vind1 = TRUE)
            for(i in seqarg) scv(y, grepl(nam[i], z, ignore.case, FALSE, fixed), args[[i]], TRUE, vind1 = TRUE)
            y
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) {
            z <- scv(y, 1L, y[1L], vind1 = TRUE) # Copy
            y <- scv(y, seq_along(y), default, set, vind1 = TRUE)  # Initialize all to default
            for(i in seqarg) scv(y, grepl(nam[i], z, ignore.case, FALSE, fixed), args[[i]], TRUE, vind1 = TRUE)
            y
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
            z <- scv(y, NA, missing, set) # y[is.na(y)] <- missing
            scv(z, nam, args, TRUE) # `[<-`(y, y == nam, value = args)
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) scv(y, nam, args, set) else y # `[<-`(y, y == nam, value = args)
        }
      } else {
        nr <- if(is.atomic(X)) NROW(X) else fnrow(X)
        if(missingl) {
          repfun <- function(y) if(is.character(y)) {
            nas <- is.na(y)
            z <- scv(y, nas, missing, set, vind1 = TRUE)
            ind <- whichv(z, nam)
            scv(z, nas, default, TRUE, TRUE, vind1 = TRUE) # duplAttributes(alloc(default, nr), y)
            scv(z, ind, args, TRUE, vind1 = TRUE) # y == nam
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) scv(scv(y, nam, default, set, TRUE), nam, args, TRUE) else y # `[<-`(duplAttributes(alloc(default, nr), y), y == nam, value = args)
        }
      }
    } else {
      seqarg <- seq_len(arglen)
      if(is.null(default)) {
        repfun <- function(y) if(is.character(y)) {
          if(missingl) y <- scv(y, NA, missing, set) # y[is.na(y)] <- missing
          else if(!set) y <- scv(y, 1L, y[1L], vind1 = TRUE) # copy
          z <- scv(y, 1L, y[1L], vind1 = TRUE) # copy
          for(i in seqarg) scv(y, whichv(z, nam[i]), args[[i]], TRUE, vind1 = TRUE)
          y
        } else y
      } else {
        nr <- if(is.atomic(X)) NROW(X) else fnrow(X)
        if(missingl) {
          repfun <- function(y) if(is.character(y)) {
            nas <- is.na(y)
            y <- scv(y, nas, missing, set, vind1 = TRUE)
            z <- scv(y, 1L, y[1L], vind1 = TRUE) # Copy
            scv(y, nas, default, TRUE, TRUE, vind1 = TRUE)
            for(i in seqarg) scv(y, whichv(z, nam[i]), args[[i]], TRUE, vind1 = TRUE)
            y
          } else y
        } else {
          repfun <- function(y) if(is.character(y)) {
            z <- scv(y, 1L, y[1L], vind1 = TRUE) # Copy
            y <- scv(y, nam[1L], default, set, TRUE) # duplAttributes(alloc(default, nr), y)
            scv(y, nam[1L], args[[1L]], TRUE)
            for(i in seqarg[-1L]) scv(y, whichv(z, nam[i]), args[[i]], TRUE, vind1 = TRUE)
            y
          } else y
        }
      }
    }
  }

  if(is.list(X)) {
    if(set) {
      lapply(unattrib(X), repfun)
      return(invisible(X))
    }
    res <- duplAttributes(lapply(unattrib(X), repfun), X)
    return(if(inherits(X, "data.table")) alc(res) else res)
  }
  if(!is.character(X)) stop("X needs to be character or a list")
  res <- repfun(X)
  return(if(set) invisible(res) else res)
}


na_locf <- function(x, set = FALSE) .Call(C_na_locf, x, set)
na_focb <- function(x, set = FALSE) .Call(C_na_focb, x, set)

na_locf_ph <- function(x, ph1, ph2, set = FALSE) .Call(C_na_locf, x, set)
na_focb_ph <- function(x, ph1, ph2, set = FALSE) .Call(C_na_focb, x, set)

replace_na <- function(X, value = 0L, cols = NULL, set = FALSE, type = "const") {
  FUN <- switch(type, const =, value = scv, locf = na_locf_ph, focb = na_focb_ph,
                stop("Unknown type:", type))
  if(set) {
    if(is.list(X)) {
      if(is.null(cols)) {
        lapply(unattrib(X), FUN, NA, value, TRUE)
      } else if(is.function(cols)) {
        lapply(unattrib(X), function(y) if(cols(y)) FUN(y, NA, value, TRUE) else y)
      } else {
        cols <- cols2int(cols, X, attr(X, "names"), FALSE)
        lapply(unattrib(X)[cols], FUN, NA, value, TRUE)
      }
    } else FUN(X, NA, value, TRUE) # `[<-`(X, is.na(X), value = value)
    return(invisible(X))
  }
  if(is.list(X)) {
    if(is.null(cols)) return(condalc(duplAttributes(lapply(unattrib(X), FUN, NA, value), X), inherits(X, "data.table"))) # function(y) `[<-`(y, is.na(y), value = value)
    if(is.function(cols)) return(condalc(duplAttributes(lapply(unattrib(X),
      function(y) if(cols(y)) FUN(y, NA, value) else y), X), inherits(X, "data.table")))
    clx <- oldClass(X)
    oldClass(X) <- NULL
    cols <- cols2int(cols, X, names(X), FALSE)
    X[cols] <- lapply(unattrib(X[cols]), FUN, NA, value) #  function(y) `[<-`(y, is.na(y), value = value)
    return(condalc(`oldClass<-`(X, clx), any(clx == "data.table")))
  }
  FUN(X, NA, value) # `[<-`(X, is.na(X), value = value)
}

replace_NA <- replace_na

# Remove Inf (Infinity) and NaN (Not a number) from vectors or data frames:
replace_inf <- function(X, value = NA, replace.nan = FALSE, set = FALSE) {
  if(set) {
    if(is.list(X)) {
      lapply(unattrib(X), if(replace.nan) (function(y) if(is.numeric(y)) scv(y, is.infinite(y) | is.nan(y), value, TRUE, vind1 = TRUE) else y) else
                                          (function(y) if(is.numeric(y)) scv(y, is.infinite(y), value, TRUE, vind1 = TRUE) else y))
    }
    if(!is.numeric(X)) stop("Infinite values can only be replaced in numeric objects!")
    if(replace.nan) scv(X, is.infinite(X) | is.nan(X), value, TRUE, vind1 = TRUE) else scv(X, is.infinite(X), value, TRUE, vind1 = TRUE)
    return(invisible(X))
  }
  if(is.list(X)) {
    # if(!inherits(X, "data.frame")) stop("replace_non_finite only works with atomic objects or data.frames")
    res <- duplAttributes(lapply(unattrib(X),
             if(replace.nan) (function(y) if(is.numeric(y)) scv(y, is.infinite(y) | is.nan(y), value, vind1 = TRUE) else y) else
                             (function(y) if(is.numeric(y)) scv(y, is.infinite(y), value, vind1 = TRUE) else y)), X)
    return(if(inherits(X, "data.table")) alc(res) else res)
  }
  if(!is.numeric(X)) stop("Infinite values can only be replaced in numeric objects!")
  if(replace.nan) return(scv(X, is.infinite(X) | is.nan(X), value, vind1 = TRUE)) #  !is.finite(X) also replaces NA
  scv(X, is.infinite(X), value, vind1 = TRUE)
}

replace_Inf <- replace_inf

# replace_non_finite <- function(X, value = NA, replace.nan = TRUE) {
#   .Deprecated("replace_Inf")
#   replace_Inf(X, value, replace.nan)
# }

Crepoutl <- function(x, limits, value, single_limit, set = FALSE) .Call(C_replace_outliers, x, limits, value, single_limit, set)

sd_limits <- function(x, limits) {
  st <- fbstatsCpp(x, stable.algo = FALSE, setn = FALSE)
  st[2L] + st[3L] * c(-limits, limits)
}

mad_limits <- function(x, limits) {
  med <- fmedian.default(x)
  mad <- fmedian.default(abs(x - med))
  med + mad * c(-limits, limits)
}

# scaling data using MAD
mad_trans <- function(x) {
  if(inherits(x, c("pseries", "pdata.frame"))) {
    g <- GRP(x)
    tmp <- fmedian(x, g, TRA = "-")
    tmp %/=% fmedian(if(is.list(tmp)) lapply(tmp, abs) else abs(tmp), g, TRA = "fill", set = TRUE)
    return(tmp)
  }
  tmp <- fmedian(x, TRA = "-")
  tmp %/=% fmedian(if(is.list(tmp)) dapply(tmp, abs) else abs(tmp), TRA = "fill", set = TRUE)
  return(tmp)
}

replace_outliers <- function(X, limits, value = NA,
                             single.limit = c("sd", "mad", "min", "max"),
                             ignore.groups = FALSE,
                             set = FALSE) {

  if(length(limits) == 1L) {
   # "overall_" arguments are legacy, now accommodated via the ignore.groups argument
   sl <- switch(single.limit[1L], SDs = 4L, min = 2L, max = 3L,
                overall_SDs = 5L, sd = 4L, mad = 6L,
                MADs = 6L, overall_MADs = 7L, # Just in case
                stop("Unknown single.limit option: ", single.limit[1L]))
   if(sl == 5L || sl == 7L) ignore.groups <- TRUE
  } else sl <- 0L

  if(sl > 3L) { # Outliers according to standard deviation or MAD threshold
    if(is.list(X)) {
      if(!ignore.groups && inherits(X, c("grouped_df", "pdata.frame"))) {
        if(is.character(value)) stop("clipping is not yet supported with grouped/panel data and SDs/MADs thresholds.")
        num <- .Call(C_vtypes, X, 1L) # vapply(unattrib(X), is.numeric, TRUE)
        num <- if(inherits(X, "grouped_df")) num & !fgroup_vars(X, "logical") else
          num & attr(findex(X), "names") %!in% attr(X, "names")
        clx <- oldClass(X)
        STDXnum <- if(sl > 5L) mad_trans(fcolsubset(X, num)) else fscale(fcolsubset(X, num))
        oldClass(X) <- NULL
        res <- .mapply(function(z, y) scv(z, abs(y) > limits, value, set, vind1 = TRUE),
                       list(unattrib(X[num]), unattrib(STDXnum)), NULL)
        if(set) return(invisible(X))
        X[num] <- res
        res <- `oldClass<-`(X, clx)
      } else {
        limit_fun <- if(sl > 5L) mad_limits else sd_limits
        res <- lapply(unattrib(X), function(y) if(is.numeric(y)) Crepoutl(y, limit_fun(y, limits), value, sl, set) else y)
        if(set) return(invisible(X))
        res <- duplAttributes(res, X)
      }
      return(if(inherits(res, "data.table")) alc(res) else res)
    }
    if(is.matrix(X)) {
      if(is.character(value)) stop("clipping is not yet supported with matrices and SDs/MADs thresholds.")
      res <- scv(X, abs(if(sl > 5L) mad_trans(X) else fscale(X)) > limits, value, set, vind1 = TRUE)
    } else {
      res <- Crepoutl(X, if(sl > 5L) mad_limits(X, limits) else sd_limits(X, limits), value, sl, set)
    }
    return(if(set) invisible(res) else res)
  }

  # Standard cases
  if(set) {
    if(is.list(X)) lapply(unattrib(X), function(y) if(is.numeric(y)) Crepoutl(y, limits, value, sl, set) else y) else
      Crepoutl(X, limits, value, sl, set)
    return(invisible(X))
  }

  if(is.list(X)) {
    res <- duplAttributes(lapply(unattrib(X), function(y) if(is.numeric(y)) Crepoutl(y, limits, value, sl, set) else y), X)
    return(if(inherits(res, "data.table")) alc(res) else res)
  }

  Crepoutl(X, limits, value, sl, set)
}



# pad or fpad? x is vector, matrix or data.frame
pad_atomic <- function(x, i, n, value) {
  ax <- attributes(x)
  tx <- typeof(x)
  if(typeof(value) != tx) value <- as.vector(value, tx)
  if(is.matrix(x)) {
    k <- dim(x)[2L]
    m <- .Call(C_alloc, value, n * k, TRUE)  # matrix(value, n, k)
    dim(m) <- c(n, k)
    m[i, ] <- x
    if(length(ax) == 1L) return(m)
    ax[["dim"]] <- c(n, k)
    # Could also pad row-names? perhaps with names of i ??
    if(length(ax[["dimnames"]][[1L]])) ax[["dimnames"]] <- list(NULL, ax[["dimnames"]][[2L]])
    if(is.object(x)) ax[["class"]] <- NULL
    return(`attributes<-`(m, ax)) # fastest ??
  }
  r <- .Call(C_alloc, value, n, TRUE) # matrix(value, n) # matrix is faster than rep_len !!!!
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
