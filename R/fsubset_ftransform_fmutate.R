

fsubset <- function(.x, ...) UseMethod("fsubset")
sbt <- fsubset

# Also not really faster than default for numeric (but a bit faster for factors ...)
fsubset.default <- function(.x, subset, ...) {
  # if(is.matrix(.x) && !inherits(.x, "matrix")) return(fsubset.matrix(.x, subset, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.logical(subset)) return(.Call(C_subsetVector, .x, which(subset), FALSE))
  .Call(C_subsetVector, .x, subset, TRUE)
}

fsubset.matrix <- function(.x, subset, ..., drop = FALSE) {
  if(missing(...)) return(.x[subset, , drop = drop])  # better row subsetting ? (like df, method? use mctl ?)
  nl <- `names<-`(as.vector(1L:ncol(.x), "list"), dimnames(.x)[[2L]])
  vars <- eval(substitute(c(...)), nl, parent.frame())
  if(missing(subset)) return(.x[, vars, drop = drop])
  .x[subset, vars, drop = drop]
}

fsubset.zoo <- function(.x, ...) if(is.matrix(.x)) fsubset.matrix(.x, ...) else fsubset.default(.x, ...)
fsubset.units <- fsubset.zoo

# No lazy eval
ss <- function(x, i, j, check = TRUE) {
  if(is.atomic(x)) if(is.matrix(x)) return(if(missing(j)) x[i, , drop = FALSE] else if(missing(i)) x[, j, drop = FALSE] else x[i, j, drop = FALSE]) else return(x[i])
  mj <- missing(j)
  if(mj) j <- seq_along(unclass(x)) else if(is.integer(j)) { # if(missing(i)) stop("Need to supply either i or j or both")
    if(missing(i)) return(.Call(C_subsetCols, x, j, TRUE))
    if(check && any(j < 0L)) j <- seq_along(unclass(x))[j]
  } else {
    if(is.character(j)) {
      j <- ckmatch(j, attr(x, "names"))
    } else if(is.logical(j)) {
      if(check && length(j) != length(unclass(x))) stop("If j is logical, it needs to be of length ncol(x)")
         j <- which(j)
    } else if(is.numeric(j)) {
     j <- if(check && any(j < 0)) seq_along(unclass(x))[j] else as.integer(j)
    } else stop("j needs to be supplied integer indices, character column names, or a suitable logical vector")
    if(missing(i)) return(.Call(C_subsetCols, x, j, TRUE))
  }
  if(!is.integer(i)) {
    if(is.numeric(i)) i <- as.integer(i) else if(is.logical(i)) {
      nr <- fnrow(x)
      if(check && length(i) != nr) stop("i needs to be integer or logical(nrow(x))") # which(r & !is.na(r)) not needed !
      i <- which(i)
      if(length(i) == nr) return(if(mj) x else .Call(C_subsetCols, x, j, TRUE))
      check <- FALSE
    } else stop("i needs to be integer or logical(nrow(x))")
  }
  rn <- attr(x, "row.names")
  if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(.Call(C_subsetDT, x, i, j, check))
  res <- .Call(C_subsetDT, x, i, j, check)
  attr(res, "row.names") <- .Call(C_subsetVector, rn, i, check)
  res
}

fsubset.data.frame <- function(.x, subset, ...) {
  r <- eval(substitute(subset), .x, parent.frame()) # Needs to be placed above any column renaming
  if(missing(...)) vars <- seq_along(unclass(.x)) else {
    ix <- seq_along(unclass(.x))
    nl <- `names<-`(as.vector(ix, "list"), attr(.x, "names"))
    vars <- eval(substitute(c(...)), nl, parent.frame())
    nam_vars <- names(vars)
    if(is.integer(vars)) {
      if(any(vars < 0L)) vars <- ix[vars]
    } else {
      if(is.character(vars)) vars <- ckmatch(vars, names(nl)) else if(is.numeric(vars)) {
        vars <- if(any(vars < 0)) ix[vars] else as.integer(vars)
      } else stop("... needs to be comma separated column names, or column indices")
    }
    if(length(nam_vars)) {
      nonmiss <- nzchar(nam_vars)
      attr(.x, "names")[vars[nonmiss]] <- nam_vars[nonmiss]
    }
  }
  checkrows <- TRUE
  if(is.logical(r)) {
    nr <- fnrow(.x)
    if(length(r) != nr) stop("subset needs to be an expression evaluating to logical(nrow(.x)) or integer") # which(r & !is.na(r)) not needed !
    r <- which(r)
    if(length(r) == nr) if(missing(...)) return(.x) else return(.Call(C_subsetCols, .x, vars, TRUE))
    checkrows <- FALSE
  } else if(is.numeric(r)) r <- as.integer(r) else
    stop("subset needs to be an expression evaluating to logical(nrow(.x)) or integer")
  rn <- attr(.x, "row.names")
  res <- .Call(C_subsetDT, .x, r, vars, checkrows)
  if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(res)
  attr(res, "row.names") <- .Call(C_subsetVector, rn, r, checkrows)
  res
}


fsubset.pseries <- function(.x, subset, ..., drop.index.levels = "id") {
  if(is.array(.x)) stop("fsubset does not support pseries matrices")
  if(!missing(...)) unused_arg_action(match.call(), ...)
  checkrows <- TRUE
  if(!is.integer(subset)) {
    if(is.numeric(subset)) subset <- as.integer(subset) else if(is.logical(subset)) {
      subset <- which(subset)
      if(length(subset) == length(.x)) return(.x)
      checkrows <- FALSE
    } else stop("subset needs to be integer or logical")
  }
  res <- .Call(C_subsetVector, .x, subset, checkrows)
  if(length(names(.x))) names(res) <- .Call(C_subsetVector, names(.x), subset, checkrows)
  index <- findex(.x)
  index_ss <- droplevels_index(.Call(C_subsetDT, index, subset, seq_along(unclass(index)), checkrows), drop.index.levels)
  attr(res, if(inherits(.x, "indexed_series")) "index_df" else "index") <- index_ss
  res
}

# Exact same code as .data.frame, just adding a block to deal with the index
fsubset.pdata.frame <- function(.x, subset, ..., drop.index.levels = "id") {
  r <- eval(substitute(subset), .x, parent.frame()) # Needs to be placed above any column renaming
  if(missing(...)) vars <- seq_along(unclass(.x)) else {
    ix <- seq_along(unclass(.x))
    nl <- `names<-`(as.vector(ix, "list"), attr(.x, "names"))
    vars <- eval(substitute(c(...)), nl, parent.frame())
    nam_vars <- names(vars)
    if(is.integer(vars)) {
      if(any(vars < 0L)) vars <- ix[vars]
    } else {
      if(is.character(vars)) vars <- ckmatch(vars, names(nl)) else if(is.numeric(vars)) {
        vars <- if(any(vars < 0)) ix[vars] else as.integer(vars)
      } else stop("... needs to be comma separated column names, or column indices")
    }
    if(length(nam_vars)) {
      nonmiss <- nzchar(nam_vars)
      attr(.x, "names")[vars[nonmiss]] <- nam_vars[nonmiss]
    }
  }
  checkrows <- TRUE
  if(is.logical(r)) {
    nr <- fnrow(.x)
    if(length(r) != nr) stop("subset needs to be an expression evaluating to logical(nrow(.x)) or integer") # which(r & !is.na(r)) not needed !
    r <- which(r)
    if(length(r) == nr) if(missing(...)) return(.x) else return(.Call(C_subsetCols, .x, vars, TRUE))
    checkrows <- FALSE
  } else if(is.numeric(r)) r <- as.integer(r) else
    stop("subset needs to be an expression evaluating to logical(nrow(.x)) or integer")
  rn <- attr(.x, "row.names")
  res <- .Call(C_subsetDT, .x, r, vars, checkrows)
  if(!(is.numeric(rn) || is.null(rn) || rn[1L] == "1")) attr(res, "row.names") <- .Call(C_subsetVector, rn, r, checkrows)
  index <- findex(.x)
  index_ss <- droplevels_index(.Call(C_subsetDT, index, r, seq_along(unclass(index)), checkrows), drop.index.levels)
  if(inherits(.x, "indexed_frame")) return(reindex(res, index_ss))
  attr(res, "index") <- index_ss
  res
}

fsubset.grouped_df <- function(.x, subset, ...) stop("fsubset() does not support grouped data: please subset your data before grouping it")

# Example:
# fsubset(GGDC10S, Variable == "VA" & Year > 1990, Country, Year, AGR:SUM)

ftransform_core <- function(X, value) { # value is unclassed, X has all attributes
  ax <- attributes(X) # keep like this ?
  oldClass(X) <- NULL
  nam <- names(value)
  if(!length(nam) || fanyDuplicated(nam)) stop("All replacement expressions have to be uniquely named")
  namX <- names(X) # !length also detects character(0)
  if(!length(namX) || fanyDuplicated(namX)) stop("All columns of .data have to be uniquely named")
  le <- vlengths(value, FALSE)
  nr <- .Call(C_fnrow, X)
  rl <- le == nr # checking if computed values have the right length
  inx <- match(nam, namX) # calling names on a plain list is really fast -> no need to save objects..
  matched <- !is.na(inx)
  if(all(rl)) { # All computed vectors have the right length
    if(any(matched)) X[inx[matched]] <- value[matched]
  } else { # Some do not
    if(any(1L < le & !rl)) stop("Lengths of replacements must be equal to nrow(.data) or 1, or NULL to delete columns")
    if(any(le1 <- le == 1L)) value[le1] <- lapply(value[le1], alloc, nr) # Length 1 arguments. can use TRA ?, or rep_len, but what about date variables ?
    if(any(le0 <- le == 0L)) { # best order -> yes, ftransform(mtcars, bla = NULL) just returns mtcars, but could also put this error message:
      if(any(le0 & !matched)) stop(paste("Can only delete existing columns, unknown columns:", paste(nam[le0 & !matched], collapse = ", ")))
      if(all(le0)) {
        X[inx[le0]] <- NULL
        return(`oldClass<-`(X, ax[["class"]]))
      }
      matched <- matched[!le0]
      value <- value[!le0] # value[le0] <- NULL
      if(any(matched)) X[inx[!le0][matched]] <- value[matched] # index is wrong after first deleting, thus we delete after !
      X[inx[le0]] <- NULL
    } else if(any(matched)) X[inx[matched]] <- value[matched] # NULL assignment ... -> Nope !
  }
  if(all(matched)) return(`oldClass<-`(X, ax[["class"]]))
  ax[["names"]] <- c(names(X), names(value)[!matched])
  setAttributes(c(X, value[!matched]), ax)
}

ftransform <- function(.data, ...) { # `_data` ?
  if(!is.list(.data)) stop(".data needs to be a list of equal length columns or a data.frame")
  e <- eval(substitute(list(...)), .data, parent.frame())
  if(is.null(names(e)) && length(e) == 1L && is.list(e[[1L]])) e <- unclass(e[[1L]]) # support list input -> added in v1.3.0
  return(condalc(ftransform_core(.data, e), inherits(.data, "data.table")))
}

tfm <- ftransform

`ftransform<-` <- function(.data, value) {
  if(!is.list(.data)) stop(".data needs to be a list of equal length columns or a data.frame")
  if(!is.list(value)) stop("value needs to be a named list")
  return(condalc(ftransform_core(.data, unclass(value)), inherits(.data, "data.table")))
}
`tfm<-` <- `ftransform<-`

# Example:
# ftransform(mtcars, cyl = cyl + 10, vs2 = 1, mpg = NULL)

eval_exp <- function(nam, exp, pe) {
  nl <- `names<-`(as.vector(seq_along(nam), "list"), nam)
  eval(exp, nl, pe)
}

ftransformv <- function(.data, vars, FUN, ..., apply = TRUE) {
  if(!is.list(.data)) stop(".data needs to be a list of equal length columns or a data.frame")
  if(!is.function(FUN)) stop("FUN needs to be a function")
  clx <- oldClass(.data)
  vs <- tryCatch(vars, error = function(e) NULL)
  if(apply) {
    oldClass(.data) <- NULL
    if(is.null(vs)) vs <- eval_exp(names(.data), substitute(vars), parent.frame())
    vars <- cols2int(vs, .data, names(.data), FALSE)
    value <- `names<-`(.data[vars], NULL)
    value <- if(missing(...)) lapply(value, FUN) else
      eval(substitute(lapply(value, FUN, ...)), .data, parent.frame())
  } else {
    nam <- attr(.data, "names")
    if(is.null(vs)) vs <- eval_exp(nam, substitute(vars), parent.frame())
    vars <- cols2int(vs, .data, nam, FALSE)
    value <- .Call(C_subsetCols, .data, vars, FALSE)
    value <- if(missing(...)) unclass(FUN(value)) else # unclass needed here ? -> yes for lengths...
      unclass(eval(substitute(FUN(value, ...)), .data, parent.frame()))
    if(!identical(names(value), nam[vars]))
      return(condalc(ftransform_core(.data, value), any(clx == "data.table")))
    oldClass(.data) <- NULL
  }
  le <- vlengths(value, FALSE)
  nr <- .Call(C_fnrow, .data)
  if(allv(le, nr)) .data[vars] <- value else if(allv(le, 1L))
    .data[vars] <- lapply(value, alloc, nr) else {
      if(apply) names(value) <- names(.data)[vars]
      .data <- ftransform_core(.data, value)
  }
  return(condalc(`oldClass<-`(.data, clx), any(clx == "data.table")))
}

tfmv <- ftransformv


settransform <- function(.data, ...) {
  assign(as.character(substitute(.data)), ftransform(.data, ...), envir = parent.frame())
  invisible(.data)
}
# eval.parent(substitute(.data <- get0("ftransform", envir = getNamespace("collapse"))(.data, ...))) # can use `<-`(.data, ftransform(.data,...)) but not faster ..

settfm <- settransform

settransformv <- function(.data, ...) {
  assign(as.character(substitute(.data)), ftransformv(.data, ...), envir = parent.frame())
  invisible(.data)
}
# eval.parent(substitute(.data <- get0("ftransformv", envir = getNamespace("collapse"))(.data, vars, FUN, ..., apply = apply)))

settfmv <- settransformv


fcompute_core <- function(.data, e, keep = NULL) {
  ax <- attributes(.data)
  nam <- ax[["names"]]
  if(!length(nam) || fanyDuplicated(nam)) stop("All columns of .data have to be uniquely named")
  if(length(keep)) {
    keep <- cols2int(keep, .data, nam, FALSE)
    if(any(m <- match(names(e), nam[keep], nomatch = 0L))) {
      temp <- .subset(.data, keep)
      pos <- m > 0L
      temp[m[pos]] <- e[pos]
      e <- c(temp, e[!pos])
    } else e <- c(.subset(.data, keep), e)
  }
  if(inherits(.data, "sf") && !any(names(e) == attr(.data, "sf_column")))
        e <- c(e, .subset(.data, attr(.data, "sf_column")))
  ax[["names"]] <- names(e)
  le <- vlengths(e, FALSE)
  nr <- fnrow(.data)
  rl <- le == nr
  if(all(rl)) return(condalcSA(e, ax, inherits(.data, "data.table"))) # All computed vectors have the right length
  if(any(1L < le & !rl)) stop("Lengths of replacements must be equal to nrow(.data) or 1")
  e[!rl] <- lapply(e[!rl], alloc, nr)
  return(condalcSA(e, ax, inherits(.data, "data.table")))
}


fcompute <- function(.data, ..., keep = NULL) { # within ?
  if(!is.list(.data)) stop(".data needs to be a list of equal length columns or a data.frame")
  e <- eval(substitute(list(...)), .data, parent.frame())
  if(is.null(names(e)) && length(e) == 1L && is.list(e[[1L]])) e <- unclass(e[[1L]]) # support list input -> added in v1.3.0
  return(fcompute_core(.data, e, keep))
}


fcomputev <- function(.data, vars, FUN, ..., apply = TRUE, keep = NULL) {
  if(!is.list(.data)) stop(".data needs to be a list of equal length columns or a data.frame")
  if(!is.function(FUN)) stop("FUN needs to be a function")
  vs <- tryCatch(vars, error = function(e) NULL)
  nam <- attr(.data, "names")
  if(is.null(vs)) vs <- eval_exp(nam, substitute(vars), parent.frame())
  vars <- cols2int(vs, .data, nam, FALSE)
  if(apply) {
    value <- `names<-`(.subset(.data, vars), NULL)
    value <- if(missing(...)) lapply(value, FUN) else
      eval(substitute(lapply(value, FUN, ...)), .data, parent.frame())
    names(value) <- nam[vars]
  } else {
    value <- .Call(C_subsetCols, .data, vars, FALSE)
    value <- if(missing(...)) unclass(FUN(value)) else # unclass needed here ? -> yes for lengths...
      unclass(eval(substitute(FUN(value, ...)), .data, parent.frame()))
  }
  return(fcompute_core(.data, value, keep)) # Note: Need to do this, value could be scalars or vectors
}


# fmutate
fFUN_mutate_add_groups <- function(z) {
  if(!is.call(z)) return(z)
  cz <- as.character(z[[1L]])
  if(length(cz) > 1L) cz <- if(any(cz == "collapse")) cz[length(cz)] else "" # needed if collapse::fmean etc..
  if(any(cz == .FAST_FUN_MOPS)) {
    z$g <- quote(.g_)
    if(any(cz == .FAST_STAT_FUN_POLD) && is.null(z$TRA)) z$TRA <- 1L
    # if(is.null(z$TRA)) z$TRA <- 1L
     # z$use.g.names <- FALSE # Not necessary

  } # This works for nested calls (nothing more required, but need to put at the end..)
  if(length(z) > 2L || is.call(z[[2L]])) return(as.call(lapply(z, fFUN_mutate_add_groups))) # Need because:  mpg - fmean(mpg)
  z
}


gsplit_single_apply <- function(x, g, ex, v, encl, unl = TRUE) {
  funexpr <- quote(function(.x_yz_) .x_yz_)
  funexpr[[3]] <- eval(call("substitute", ex, `names<-`(list(quote(.x_yz_)), v)), NULL, NULL)
  funexpr[[4]] <- NULL
  fun <- eval(funexpr, encl, baseenv())
  res <- lapply(gsplit(x, g), fun)
  if(unl) copyMostAttributes(unlist(res, FALSE, FALSE), x) else res
}

# Old version: more expensive...
# gsplit_single_apply <- function(x, g, ex, v, encl)
#   copyMostAttributes(unlist(lapply(gsplit(x, g), function(i) eval(ex, `names<-`(list(i), v), encl)), FALSE, FALSE), x)

gsplit_multi_apply <- function(x, g, ex, encl, SD = FALSE) {
  sx <- seq_along(x)
  gs <- gsplit(NULL, g)
  if(!SD) return(lapply(gs, function(i) eval(ex, .Call(C_subsetDT, x, i, sx, FALSE), encl)))
  funexpr <- substitute(function(.data) expr, list(expr = ex))
  funexpr[[4]] <- NULL
  fun <- eval(funexpr, encl, baseenv())
  lapply(gs, function(i) fun(.Call(C_subsetDT, x, i, sx, FALSE)))
}

othFUN_compute <- function(x) {
  if(length(x) == 2L) # No additional function arguments
    return(substitute(lapply(.gsplit_(a, .g_), b),
                      list(a = x[[2L]], b = x[[1L]])))
  # With more arguments, things become more complex..
  as.call(c(list(quote(lapply), substitute(.gsplit_(a, .g_), list(a = x[[2L]]))), as.list(x[-2L])))
}

keep_v <- function(d, v) copyMostAttributes(null_rm(.subset(d, unique.default(v))), d)

acr_get_cols <- function(.cols, d, nam, ce) {
  # Note: .cols is passed through substitute() before it enters here. Thus only an explicit NULL is NULL up front
  if(is.null(.cols)) return(if(is.null(d[[".g_"]])) seq_along(nam) else seq_along(nam)[nam %!in% c(".g_", ".gsplit_", d[[".g_"]]$group.vars)])
  nl <- `names<-`(as.vector(seq_along(nam), "list"), nam)
  cols <- eval(.cols, nl, ce)
  # Needed for programming usage, because you can pass a variable that is null
  if(is.null(cols)) return(if(is.null(d[[".g_"]])) seq_along(nam) else seq_along(nam)[nam %!in% c(".g_", ".gsplit_", d[[".g_"]]$group.vars)])
  if(is.logical(cols)) return(which(cols)) # if .g_ etc. is added to data, length check for logical vectors will fail
  if(is.null(d[[".g_"]]) || is.character(cols) || (is.numeric(cols) && cols[1L] > 0)) return(cols2int(cols, d, nam))
  cols2intrmgn(match(c(".g_", ".gsplit_", d[[".g_"]]$group.vars), nam), cols, d)
}

# Also used in collap()
acr_get_funs <- function(.fnsexp, .fns, ...) {

  if(is.function(.fns)) {
    namfun <- l1orlst(as.character(.fnsexp))
    .fns <- `names<-`(list(.fns), namfun)
  } else if(is.list(.fns)) {
    namfun <- names(.fns)
    # In programming usage, could simply pass a list of functions l, in which case this is not a call..
    if(is.call(.fnsexp) && (.fnsexp[[1L]] == quote(list) || .fnsexp[[1L]] == quote(c))) { # or we could have funlist[[i]] which is also sorted out here..
      nf <- all.vars(.fnsexp, unique = FALSE)
      if(length(nf) == length(.fns)) {
        names(.fns) <- nf
        if(is.null(namfun)) namfun <- nf
      } else {
        nf <- vapply(.fnsexp[-1L], function(x) l1orlst(all.vars(x)), character(1L), USE.NAMES = FALSE)
        names(.fns) <- nf
        if(is.null(namfun)) namfun <- as.character(seq_along(.fns))
      }
    } else if(is.null(namfun)) names(.fns) <- namfun <- as.character(seq_along(.fns))
  } else if(is.character(.fns)) {
    namfun <- names(.fns)
    names(.fns) <- .fns
    .fns <- lapply(.fns, ...) # lapply(.fns, match.fun())
    if(is.null(namfun)) namfun <- names(.fns)
  } else stop(".fns must be a function, list of functions or character vector of function names")

  return(list(namfun = namfun, funs = .fns))
}


fungroup2 <- function(X, ocl) {
  attr(X, "groups") <- NULL
  oldClass(X) <- fsetdiff(ocl, c("GRP_df", "grouped_df"))
  X
}


setup_across <- function(.cols, .fnsexp, .fns, .names, .apply, .transpose, .FFUN) {
  pe <- parent.frame(n = 4L)
  d <- unclass(pe$.data) # Safer to unclass here also...
  ce <- parent.frame(n = 5L) # Caller environment
  # return(list(.cols, .fns, .names, d))
  nam <- names(d)
  cols <- acr_get_cols(.cols, d, nam, ce)
  funs <- acr_get_funs(.fnsexp, .fns, get, mode = "function", envir = ce)
  namfun <- funs$namfun
  fun <- funs$funs

  if(length(.names) && !is.logical(.names)) {
    if(is.function(.names)) {
      names <- if(isFALSE(.transpose)) # .names(nam[cols], namfun)
        as.vector(outer(nam[cols], namfun, .names)) else
        as.vector(t(outer(nam[cols], namfun, .names)))
    } else {
      if(length(.names) == 1L && .names == "flip") {
        names <- if(isFALSE(.transpose))
          as.vector(outer(nam[cols], namfun, function(z, f) paste(f, z, sep = "_"))) else
          as.vector(t(outer(nam[cols], namfun, function(z, f) paste(f, z, sep = "_"))))
      } else {
        if(length(.names) != length(namfun) * length(cols)) stop("length(.names) must match length(.fns) * length(.cols)")
        names <- .names
      }
    }
  } else {
    # Third version: .names = FALSE does nothing. Allows fmutate(mtcars, across(cyl:vs, list(L, D, G), n = 1:3))
    # This makes sense, because if .transpose = "auto" and the lengths of generated columns are unequal, you cannot use generated names anyway because they would mismatch..
    names <- if((is.null(.names) && length(namfun) == 1L) || (isFALSE(.names) && length(namfun) > 1L)) NULL else if(isFALSE(.names)) # this allows you to force names false for a single function...
             nam[cols] else if(isFALSE(.transpose))
             as.vector(outer(nam[cols], namfun, paste, sep = "_")) else
             as.vector(t(outer(nam[cols], namfun, paste, sep = "_")))
    # Second version: .names = TRUE auto generates names, .names = FALSE yields default names (no change to names by the function),
    # and .names = NULL (default) yields function names or auto names if multiple functions...
    # names <- if(is.null(.names) && length(namfun) == 1L) NULL else if(!isFALSE(.names))
    #          as.vector(t(outer(nam[cols], namfun, paste, sep = "_"))) else if(length(namfun) == 1L)
    #          nam[cols] else stop("Computed columns need to be uniquely named. If .names = FALSE, can only use one function, or need to supply custom names!")

    # First version: requires .names = FALSE for renaming functions like L, W etc...
    # names <- if(isFALSE(.names)) NULL else
    #   if(length(namfun) == 1L && !isTRUE(.names)) nam[cols] else
    #    as.vector(t(outer(nam[cols], namfun, paste, sep = "_")))
  }

  if(is.logical(.apply)) {
    aplvec <- if(.apply) rep_len(TRUE, length(fun)) else rep_len(FALSE, length(fun))
  } else {
    .apply <- switch(.apply, auto = NA, stop(".apply must be 'auto', TRUE or FALSE"))
    aplvec <- names(fun) %!in% .FFUN
  }
  .data_ <- if(all(aplvec)) d[cols] else .Call(C_subsetCols, if(is.null(d[[".g_"]])) `oldClass<-`(d, pe$cld) else fungroup2(d, pe$cld), cols, FALSE)

  # Note: Keep the order and the names !!!
  list(data = d,
       .data_ = .data_, # cols = cols,
       funs = fun,
       aplvec = aplvec,
       ce = ce,
       names = names)
}

across <- function(.cols = NULL, .fns, ..., .names = NULL, .apply = "auto", .transpose = "auto") {
  stop("across() can only work inside fmutate() and fsummarise()")
}

do_across <- function(.cols = NULL, .fns, ..., .names = NULL, .apply = "auto", .transpose = "auto", .eval_funi, .summ = TRUE) {
  # nodots <- missing(...)
  # return(setup_across(substitute(.cols), substitute(.fns), .fns, .names, .apply, .FAST_FUN_MOPS))
  setup <- setup_across(substitute(.cols), substitute(.fns), .fns, .names, .apply, .transpose, .FAST_FUN_MOPS)
  seqf <- seq_along(setup$funs)
  names <- setup$names
  # return(eval_funi(seqf, ...))
  # return(lapply(seqf, eval_funi, ...))
  if(length(seqf) == 1L) {
    res <- .eval_funi(seqf, setup[[1L]], setup[[2L]], setup[[3L]], setup[[4L]], setup[[5L]], ...)  # eval_funi(seqf, aplvec, funs, nodots, .data_, data, ce, ...)
    # return(res)
  } else {
    # motivated by: fmutate(mtcars, across(cyl:vs, list(L, D, G), n = 1:3))
    r <- lapply(seqf, .eval_funi, setup[[1L]], setup[[2L]], setup[[3L]], setup[[4L]], setup[[5L]], ...) # do.call(lapply, c(list(seqf, eval_funi), setup[1:5], list(...))) # lapply(seqf, eval_funi, aplvec, funs, nodots, .data_, data, ce, ...)
    # return(r)
    if(isFALSE(.transpose) || (is.character(.transpose) && !all_eq(vlengths(r, FALSE)))) {
      # stop("reached here")
      res <- unlist(r, FALSE, use.names = TRUE) # need use.names= TRUE here
      # return(list(res = res, r = r))
    } else {
      res <- unlist(t_list2(r), FALSE, FALSE)
      if(is.null(names(res)) && is.null(names))
        names(res) <- unlist(t_list2(lapply(r, names)), FALSE, FALSE)
    }
  }
  if(.summ) return(if(is.null(names)) res else `names<-`(res, names))
  return(`[<-`(setup$data, if(is.null(names)) names(res) else names, value = res))
}

mutate_funi_simple <- function(i, data, .data_, funs, aplvec, ce, ...) { # g is unused here...
  .FUN_ <- funs[[i]]
  nami <- names(funs)[i]
  if(aplvec[i]) {
    value <- if(missing(...)) lapply(unattrib(.data_), .FUN_) else
      do.call(lapply, c(list(unattrib(.data_), .FUN_), eval(substitute(list(...)), data, ce)), envir = ce) # eval(substitute(lapply(unattrib(.data_), .FUN_, ...)), c(list(.data_ = .data_), data), ce)
    names(value) <- names(.data_)
  } else if(any(nami == .FAST_STAT_FUN_POLD)) {
    if(missing(...)) return(unclass(.FUN_(.data_, TRA = 1L))) # Old way: Not necessary to construct call.. return(unclass(eval(as.call(list(as.name(nami), quote(.data_), TRA = 1L))))) # faster than substitute(.FUN_(.data_, TRA = 1L), list(.FUN_ = as.name(nami)))
    # if(any(...names() == "TRA")) # This down not work because it substitutes setup[[]] from mutate_across !!!
    #   return(unclass(eval(substitute(.FUN_(.data_, ...)), c(list(.data_ = .data_), data), ce)))
    # return(unclass(eval(substitute(.FUN_(.data_, ..., TRA = 1L)), c(list(.data_ = .data_), data), ce)))
    fcal <- as.call(c(list(as.name(nami), quote(.data_)), as.list(substitute(list(...))[-1L])))
    if(is.null(fcal$TRA)) fcal$TRA <- 1L
    return(unclass(eval(fcal, c(list(.data_ = .data_), data), ce)))
  } else {
    value <- if(missing(...)) .FUN_(.data_) else
      do.call(.FUN_, c(list(.data_), eval(substitute(list(...)), data, ce)), envir = ce) # Object setup not found: eval(substitute(.FUN_(.data_, ...)), c(list(.data_ = .data_), data), ce)
    oldClass(value) <- NULL
    if(any(nami == .FAST_FUN_MOPS)) return(value) # small improvement for fast funs...
  }
  # return(unclass(r))
  # fcal <- if(missing(...)) as.call(list(funs[[nami]], quote(.data_))) else
  #         as.call(c(list(funs[[nami]], quote(.data_)), as.list(substitute(list(...))[-1L]))) # , parent.frame()
  #   # substitute(list(...), parent.frame())
  #   # substitute(FUN(.data_, ...), list(FUN = funs[[nami]], ...))
  #   # as.call(substitute(list(funs[[nami]], quote(.data_), ...)))
  #   # substitute(FUN(.data_, ...), list(FUN = funs[[nami]]))  #
  # if(any(nami == .FAST_STAT_FUN_POLD) && is.null(fcal$TRA)) fcal$TRA <- 1L
  # fast functions have a data.frame method, thus can be applied simultaneously to all columns
  # return(fcal)
  # return(eval(fcal, c(list(.data_ = .data_), data), setup$ce))
  lv <- vlengths(value, FALSE)
  nr <- .Call(C_fnrow, data)
  if(allv(lv, nr)) return(value)
  if(allv(lv, 1L)) return(lapply(value, alloc, nr))
  stop("Without groups, NROW(value) must either be 1 or nrow(.data)")
}

dots_apply_grouped <- function(d, g, f, dots) {
  attributes(d) <- NULL
  n <- length(d[[1L]])
  # Arguments same length as data
  if(length(ln <- whichv(vlengths(dots, FALSE), n))) {
    asl <- lapply(dots[ln], gsplit, g)
    if(length(dots) > length(ln)) {
      mord <- dots[-ln]
      if(is.null(names(mord)) && is.null(names(asl))) warning("If some arguments have the same length as the data (vectors) while others have length 1 (scalars), please ensure that at least one of the two have keywords e.g. argname = value. This is because the latter are passed to the 'MoreArgs' argument of .mapply, and thus the order in which arguments are passed to the function might be different from your top-level call. In particular, .mapply will first pass the vector valued arguments followed by the scalar valued ones.")
      FUN <- function(x) .mapply(f, c(list(gsplit(x, g)), asl), mord) # do.call(mapply, c(list(f, gsplit(x, g), SIMPLIFY = FALSE, USE.NAMES = FALSE, MoreArgs = mord), asl))
    } else FUN <- function(x) .mapply(f, c(list(gsplit(x, g)), asl), NULL) # do.call(mapply, c(list(f, gsplit(x, g), SIMPLIFY = FALSE, USE.NAMES = FALSE), asl))
    return(lapply(d, function(y) copyMostAttributes(unlist(FUN(y), FALSE, FALSE), y)))
  }
  # No arguments to be split
  do.call(lapply, c(list(d, copysplaplfun, g, f), dots))
}

dots_apply_grouped_bulk <- function(d, g, f, dots) {
  n <- fnrow(d)
  dsp <- rsplit.data.frame(d, g, simplify = FALSE, flatten = TRUE, use.names = FALSE)
  if(is.null(dots)) return(lapply(dsp, f))
  # Arguments withs ame length as data
  if(length(ln <- whichv(vlengths(dots, FALSE), n))) {
    asl <- lapply(dots[ln], gsplit, g)
    if(length(dots) > length(ln)) {
      mord <- dots[-ln]
      if(is.null(names(mord)) && is.null(names(asl))) warning("If some arguments have the same length as the data (vectors) while others have length 1 (scalars), please ensure that at least one of the two have keywords e.g. argname = value. This is because the latter are passed to the 'MoreArgs' argument of .mapply, and thus the order in which arguments are passed to the function might be different from your top-level call. In particular, .mapply will first pass the vector valued arguments followed by the scalar valued ones.")
    } else mord <- NULL
    return(.mapply(f, c(list(dsp), asl), mord))
  }
  # No arguments to be split
  do.call(lapply, c(list(dsp, f), dots))
}

mutate_grouped_expand <- function(value, g) {
  lv <- vlengths(value, FALSE)
  nr <- length(g[[2L]])
  if(allv(lv, nr)) {
    if(!isTRUE(g$ordered[2L])) {
      if(length(value) < 4L) { # optimal?
        value <- lapply(value, function(x, g) .Call(C_greorder, x, g), g)
      } else {
        ind <- .Call(C_greorder, seq_len(nr), g)
        value <- .Call(C_subsetDT, value, ind, seq_along(value), FALSE)
      }
    }
    return(value)
  }
  if(!allv(lv, g[[1L]])) stop("With groups, NROW(value) must either be ng or nrow(.data)")
  return(.Call(C_subsetDT, value, g[[2L]], seq_along(value), FALSE))
}

mutate_funi_grouped <- function(i, data, .data_, funs, aplvec, ce, ...) {
  g <- data[[".g_"]]
  .FUN_ <- funs[[i]]
  nami <- names(funs)[i]
  apli <- aplvec[i]
  if(apli) {
    value <- if(missing(...)) lapply(unattrib(.data_), copysplaplfun, g, .FUN_) else
             dots_apply_grouped(.data_, g, .FUN_, eval(substitute(list(...)), data, ce)) # Before: do.call(lapply, c(list(unattrib(.data_), copysplaplfun, g, .FUN_), eval(substitute(list(...)), data, ce)), envir = ce)
  } else if(any(nami == .FAST_STAT_FUN_POLD)) {
    if(missing(...)) return(unclass(.FUN_(.data_, g = g, TRA = 1L)))
    fcal <- as.call(c(list(as.name(nami), quote(.data_), g = quote(.g_)), as.list(substitute(list(...))[-1L])))
    if(is.null(fcal$TRA)) fcal$TRA <- 1L
    return(unclass(eval(fcal, c(list(.data_ = .data_), data), ce)))
  } else if(any(nami == .FAST_FUN_MOPS)) {
    if(any(nami == .OPERATOR_FUN)) {
      value <- if(missing(...)) .FUN_(.data_, by = g) else
        do.call(.FUN_, c(list(.data_, by = g), eval(substitute(list(...)), data, ce)), envir = ce)
    } else {
      value <- if(missing(...)) .FUN_(.data_, g = g) else
        do.call(.FUN_, c(list(.data_, g = g), eval(substitute(list(...)), data, ce)), envir = ce)
    }
    oldClass(value) <- NULL
    return(value)
  } else { # stop("In grouped computations, .apply = FALSE only works with .FAST_FUN and .OPERATOR_FUN")
    value <- dots_apply_grouped_bulk(.data_, g, .FUN_, if(missing(...)) NULL else eval(substitute(list(...)), data, ce))
    value <- .Call(C_rbindlist, unclass(value), FALSE, FALSE, NULL)
    oldClass(value) <- NULL
  }
  if(apli) names(value) <- names(.data_)
  return(mutate_grouped_expand(value, g))
}


do_grouped_expr <- function(ei, nfun, .data, g, pe) {
  v <- all.vars(ei) # unique = FALSE -> not needed anymore... can turn expressions into functions...
  if(length(v) > 1L) {
    # Could include global environmental variables e.g. fmutate(data, new = mean(var) + q)
    namd <- names(.data)
    if(length(wv <- na_rm(match(v, namd))) > 1L) return(unlist(gsplit_multi_apply(.data[wv], g, ei, pe), FALSE, FALSE))
    return(gsplit_single_apply(.data[[wv]], g, ei, namd[wv], pe))
  }
  if(nfun == 1L) {
    res <- eval(othFUN_compute(ei), .data, pe)
    return(copyMostAttributes(unlist(res, FALSE, FALSE), .data[[v]]))
  }
  gsplit_single_apply(.data[[v]], g, ei, v, pe)
}

# Same as above, without unlisting...
do_grouped_expr_list <- function(ei, .data, g, pe, .cols, ax, mutate = FALSE) {
  v <- all.vars(ei)
  if(any(v == ".data")) {
    .data[names(.data) %in% c(".g_", ".gsplit_", if(is.null(.cols)) g$group.vars)] <- NULL
    if(is.character(ax)) { # for fmutate
      cld <- ax
      ax <- attributes(.data)
      ax[["groups"]] <- NULL
      # ax[["names"]] <- fsetdiff(ax[["names"]], c(".g_", ".gsplit_")) # Redundant, removed above...
      ax[["class"]] <- fsetdiff(cld, c("GRP_df", "grouped_df"))
    }
    if(length(.cols)) .data <- colsubset(.data, .cols)
    ax[["names"]] <- names(.data)
    setattributes(.data, ax)
    res <- gsplit_multi_apply(.data, g, ei, pe, TRUE)
  } else if(length(v) > 1L) {
    namd <- names(.data)
    res <- if(length(wv <- na_rm(match(v, namd))) > 1L)
      gsplit_multi_apply(.data[wv], g, ei, pe) else
      gsplit_single_apply(.data[[wv]], g, ei, namd[wv], pe, FALSE)
  } else {
    res <- if(length(all_funs(ei)) == 1L) eval(othFUN_compute(ei), .data, pe) else
      gsplit_single_apply(.data[[v]], g, ei, v, pe, FALSE)
  }
  res <- .Call(C_rbindlist, res, FALSE, FALSE, NULL)
  if(mutate) return(mutate_grouped_expand(res, g))
  res
}


fmutate <- function(.data, ..., .keep = "all", .cols = NULL) {
  if(!is.list(.data)) stop(".data needs to be a list of equal length columns or a data.frame")
  e <- substitute(list(...))
  nam <- names(e)
  nullnam <- is.null(nam)
  # if(!length(nam)) stop("All replacement expressions have to be named")
  pe <- parent.frame()
  cld <- oldClass(.data) # This needs to be called cld, because across fetches it from here !!
  oldClass(.data) <- NULL
  nr <- .Call(C_fnrow, .data)
  namdata <- names(.data)
  if(is.null(namdata) || fanyDuplicated(namdata)) stop("All columns of .data have to be uniquely named")
  if(!is.character(.keep)) .keep <- cols2char(.keep, .data, namdata) # allowing .keep to be NULL
  gdfl <- any(cld == "grouped_df")
  if(gdfl) {
    g <- GRP.grouped_df(.data, return.groups = FALSE, call = FALSE)
    .data[c(".g_", ".gsplit_")] <- list(g, gsplit)
    for(i in 2:length(e)) {
      ei <- e[[i]]
      if(nullnam || nam[i] == "") { # Across
        if(ei[[1L]] == quote(across) || ei[[1L]] == quote(acr)) {
          ei[[1L]] <- quote(do_across)
          ei$.eval_funi <- quote(mutate_funi_grouped)
          ei$.summ <- FALSE
          # return(eval(ei, enclos = pe))
          .data <- eval(ei, list(do_across = do_across, mutate_funi_grouped = mutate_funi_grouped), pe) # ftransform_core(.data, eval(ei, pe))
        } else {
          r <- do_grouped_expr_list(ei, .data, g, pe, .cols, cld, TRUE)
          .data[names(r)] <- r
        }
      } else { # Tagged vector expressions
        if(is.null(ei)) {
          .data[[nam[i]]] <- NULL
          next
        }
        eif <- all_funs(ei)
        if(any(eif %in% .FAST_FUN_MOPS)) {
          .data[[nam[i]]] <- eval(fFUN_mutate_add_groups(ei), .data, pe)
        } else if(length(eif)) {
          r <- do_grouped_expr(ei, length(eif), .data, g, pe)
          .data[[nam[i]]] <- if(length(r) == g[[1L]])
               .Call(C_subsetVector, r, g[[2L]], FALSE) else # .Call(C_TRA, .data[[v]], r, g[[2L]], 1L) # Faster than simple subset r[g[[2L]] ??]
               .Call(C_greorder, r, g) # r[forder.int(forder.int(g[[2L]]))] # Seems twice is necessary...
        } else { # something like bla = 1 or mpg = vs
          r <- eval(ei, .data, pe)
          if(length(r) == 1L) r <- alloc(r, nr)
          else if(length(r) != nr) stop("length mismatch")
          .data[[nam[i]]] <- r
        }
      }
    }
    .data[c(".g_", ".gsplit_")] <- NULL
  } else { # Without groups...
    for(i in 2:length(e)) { # This is good and very fast
      ei <- e[[i]]
      if(nullnam || nam[i] == "") { # Across
        if(ei[[1L]] == quote(across) || ei[[1L]] == quote(acr)) { # stop("expressions need to be named or start with across(), or its shorthand acr().")
          ei[[1L]] <- quote(do_across)
          ei$.eval_funi <- quote(mutate_funi_simple)
          ei$.summ <- FALSE
          # return(eval(ei, enclos = pe))
          .data <- eval(ei, list(do_across = do_across, mutate_funi_simple = mutate_funi_simple), pe) # ftransform_core(.data, eval(ei, enclos = pe))
        } else {
          r <- eval(ei, .data, pe)
          .data[names(r)] <- r
        }
      } else { # Tagged vector expressions
        r <- eval(ei, .data, pe)
        if(!is.null(r)) { # don't use length(), because only NULL removes list elements...
          if(length(r) == 1L) r <- alloc(r, nr)
          else if(length(r) != nr) stop("length mismatch")
        }
        .data[[nam[i]]] <- r
      }
    }
  }
  # Implementing .keep argument
  # TODO: Implement .keep with across...
  .data <- if(length(.keep) > 1L) keep_v(.data, c(.keep, nam[-1L])) else
    switch(.keep,
           all = .data,
           used = keep_v(.data, c(namdata[namdata %in% c(if(gdfl) g$group.vars, unlist(lapply(e[-1L], all.vars), FALSE, FALSE), nam[-1L])], nam[-1L])),
           unused = keep_v(.data, c(namdata[namdata %in% c(if(gdfl) g$group.vars, fsetdiff(namdata, unlist(lapply(e[-1L], all.vars), FALSE, FALSE)), nam[-1L])], nam[-1L])),
           none = keep_v(.data, c(if(gdfl) g$group.vars, nam[-1L])), # g$group.vars[g$group.vars %!in% nam[-1L]] -> inconsistent and inefficient...
           keep_v(.data, c(.keep, nam[-1L])))

  oldClass(.data) <- cld
  return(condalc(.data, any(cld == "data.table")))
}

# or mut / mte? () If you need o choose a vowel, u is more distinctive, lut for consistency let's stock with consonants
mtt <- fmutate # Note: see if function(.data, ...) fmutate(.data, ...) is possible (what about objects in global environment?)




