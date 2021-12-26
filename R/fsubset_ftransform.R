

fsubset <- function(x, ...) UseMethod("fsubset")
sbt <- fsubset

# Also not really faster than default for numeric (but a bit faster for factors ...)
fsubset.default <- function(x, subset, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fsubset.matrix(x, subset, ...))
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.logical(subset)) return(.Call(C_subsetVector, x, which(subset), FALSE))
  .Call(C_subsetVector, x, subset, TRUE)
}

fsubset.matrix <- function(x, subset, ..., drop = FALSE) {
  if(missing(...)) return(x[subset, , drop = drop])  # better row subsetting ? (like df, method? use mctl ?)
  nl <- `names<-`(as.vector(1L:ncol(x), "list"), dimnames(x)[[2L]])
  vars <- eval(substitute(c(...)), nl, parent.frame())
  if(missing(subset)) return(x[, vars, drop = drop])
  x[subset, vars, drop = drop]
}

# No lazy eval
ss <- function(x, i, j) {
  if(is.atomic(x)) if(is.array(x)) return(if(missing(j)) x[i, , drop = FALSE] else x[i, j, drop = FALSE]) else return(x[i])
  mj <- missing(j)
  if(mj) j <- seq_along(unclass(x)) else if(is.integer(j)) {
    if(any(j < 0L)) j <- seq_along(unclass(x))[j]
  } else {
    if(is.character(j)) {
      j <- ckmatch(j, attr(x, "names"))
    } else if(is.logical(j)) {
      if(length(j) != length(unclass(x))) stop("If j is logical, it needs to be of length ncol(x)")
         j <- which(j)
    } else if(is.numeric(j)) {
     j <- if(any(j < 0)) seq_along(unclass(x))[j] else as.integer(j)
    } else stop("j needs to be supplied integer indices, character column names, or a suitable logical vector")
  }
  checkrows <- TRUE
  if(!is.integer(i)) {
    if(is.numeric(i)) i <- as.integer(i) else if(is.logical(i)) {
      nr <- fnrow2(x)
      if(length(i) != nr) stop("i needs to be integer or logical(nrow(x))") # which(r & !is.na(r)) not needed !
      i <- which(i)
      if(length(i) == nr) if(mj) return(x) else return(.Call(C_subsetCols, x, j, TRUE))
      checkrows <- FALSE
    } else stop("i needs to be integer or logical(nrow(x))")
  }
  rn <- attr(x, "row.names")
  if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(.Call(C_subsetDT, x, i, j, checkrows))
  return(`attr<-`(.Call(C_subsetDT, x, i, j, checkrows), "row.names", rn[i]))
}

fsubset.data.frame <- function(x, subset, ...) {
  r <- eval(substitute(subset), x, parent.frame()) # Needs to be placed above any column renaming
  if(missing(...)) vars <- seq_along(unclass(x)) else {
    ix <- seq_along(unclass(x))
    nl <- `names<-`(as.vector(ix, "list"), attr(x, "names"))
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
      attr(x, "names")[vars[nonmiss]] <- nam_vars[nonmiss]
    }
  }
  checkrows <- TRUE
  if(is.logical(r)) {
    nr <- fnrow2(x)
    if(length(r) != nr) stop("subset needs to be an expression evaluating to logical(nrow(x)) or integer") # which(r & !is.na(r)) not needed !
    r <- which(r)
    if(length(r) == nr) if(missing(...)) return(x) else return(.Call(C_subsetCols, x, vars, TRUE))
    checkrows <- FALSE
  } else if(is.numeric(r)) r <- as.integer(r) else
    stop("subset needs to be an expression evaluating to logical(nrow(x)) or integer")
  rn <- attr(x, "row.names")
  if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(.Call(C_subsetDT, x, r, vars, checkrows))
  return(`attr<-`(.Call(C_subsetDT, x, r, vars, checkrows), "row.names", rn[r]))
}

# Example:
# fsubset(GGDC10S, Variable == "VA" & Year > 1990, Country, Year, AGR:SUM)

ftransform_core <- function(X, value) { # value is unclassed, X has all attributes
  ax <- attributes(X) # keep like this ?
  oldClass(X) <- NULL
  nam <- names(value)
  if(!length(nam) || fanyDuplicated(nam)) stop("All replacement expressions have to be uniquely named")
  namX <- names(X) # !length also detects character(0)
  if(!length(namX) || fanyDuplicated(namX)) stop("All columns of X have to be uniquely named")
  le <- lengths(value, FALSE)
  nr <- length(X[[1L]])
  rl <- le == nr # checking if computed values have the right length
  inx <- match(nam, namX) # calling names on a plain list is really fast -> no need to save objects..
  matched <- !is.na(inx)
  if(all(rl)) { # All computed vectors have the right length
    if(any(matched)) X[inx[matched]] <- value[matched]
  } else { # Some do not
    if(any(1L < le & !rl)) stop("Lengths of replacements must be equal to nrow(X) or 1, or NULL to delete columns")
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

ftransformv <- function(.data, vars, FUN, ..., apply = TRUE) {
  if(!is.list(.data)) stop(".data needs to be a list of equal length columns or a data.frame")
  if(!is.function(FUN)) stop("FUN needs to be a function")
  clx <- oldClass(.data)
  if(apply) {
    oldClass(.data) <- NULL
    vars <- cols2int(vars, .data, names(.data), FALSE)
    value <- `names<-`(.data[vars], NULL)
    value <- if(missing(...)) lapply(value, FUN) else
      eval(substitute(lapply(value, FUN, ...)), .data, parent.frame())
  } else {
    nam <- attr(.data, "names")
    vars <- cols2int(vars, .data, nam, FALSE)
    value <- .Call(C_subsetCols, .data, vars, FALSE)
    value <- if(missing(...)) unclass(FUN(value)) else # unclass needed here ? -> yes for lengths...
      unclass(eval(substitute(FUN(value, ...)), .data, parent.frame()))
    if(!identical(names(value), nam[vars]))
      return(condalc(ftransform_core(.data, value), any(clx == "data.table")))
    oldClass(.data) <- NULL
  }
  le <- lengths(value, FALSE)
  nr <- length(.data[[1L]])
  if(all(le == nr)) .data[vars] <- value else if(all(le == 1L))
    .data[vars] <- lapply(value, alloc, nr) else {
      if(apply) names(value) <- names(.data)[vars]
      .data <- ftransform_core(.data, value)
  }
  return(condalc(`oldClass<-`(.data, clx), any(clx == "data.table")))
}

tfmv <- ftransformv


settransform <- function(.data, ...)
  assign(as.character(substitute(.data)), ftransform(.data, ...), envir = parent.frame())
# eval.parent(substitute(.data <- get0("ftransform", envir = getNamespace("collapse"))(.data, ...))) # can use `<-`(.data, ftransform(.data,...)) but not faster ..

settfm <- settransform

settransformv <- function(.data, vars, FUN, ..., apply = TRUE)
  assign(as.character(substitute(.data)), ftransformv(.data, vars, FUN, ..., apply = apply), envir = parent.frame())
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
  le <- lengths(e, FALSE)
  nr <- fnrow2(.data)
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
  if(apply) {
    nam <- attr(.data, "names")
    vars <- cols2int(vars, .data, nam, FALSE)
    value <- `names<-`(.subset(.data, vars), NULL)
    value <- if(missing(...)) lapply(value, FUN) else
      eval(substitute(lapply(value, FUN, ...)), .data, parent.frame())
    names(value) <- nam[vars]
  } else {
    vars <- cols2int(vars, .data, attr(.data, "names"), FALSE)
    value <- .Call(C_subsetCols, .data, vars, FALSE)
    value <- if(missing(...)) unclass(FUN(value)) else # unclass needed here ? -> yes for lengths...
      unclass(eval(substitute(FUN(value, ...)), .data, parent.frame()))
  }
  return(fcompute_core(.data, value, keep)) # Note: Need to do this, value could be scalars or vectors
}


# Fmutate
fFUN_mutate_add_groups <- function(z) {
  if(!is.call(z)) return(z)
  cz <- as.character(z[[1L]])
  if(any(cz == .FAST_FUN_MOPS)) {
    z$g <- quote(.g_)
    if(any(cz == .FAST_STAT_FUN_POLD)) {
      if(is.null(z$TRA)) z$TRA <- 1L
      z$use.g.names <- FALSE
    }
  } # This works for nested calls (nothing more required, but need to put at the end..)
  if(is.call(z[[2L]])) return(as.call(lapply(z, fFUN_mutate_add_groups)))
  z
}

# TODO: Improve rsplit...
# gsplit_multi <- function(x, g)
#   lapply(gsplit(1L, g, toint = TRUE), .Call, .NAME = C_subsetDT, seq_along(x), FALSE)

gsplit_single_apply <- function(x, g, ex, encl)
  copyMostAttributes(unlist(lapply(gsplit(x, g), function(i) eval(ex, `names<-`(list(i), v), encl)), FALSE, FALSE), x)

# TODO: enable for fsummarise as well...
gsplit_multi_apply <- function(x, g, ex, encl) {
  sx <- seq_along(x)
  unlist(lapply(gsplit(1L, g, toint = TRUE),
         function(i) eval(ex, .Call(C_subsetDT, x, i, sx, FALSE), encl)), FALSE, FALSE)
}
# Also todo: keep argument...
fmutate <- function(.data, ...) {# , TRA = "replace_fill" # TODO: Implement TRA !!!
  if(!is.list(.data)) stop(".data needs to be a list of equal length columns or a data.frame")
  namdata <- attr(.data, "names")
  if(!length(namdata) || fanyDuplicated(namdata)) stop("All columns of .data have to be uniquely named")
  e <- substitute(list(...))
  nam <- names(e)
  if(!length(nam)) stop("All replacement expressions have to be named")
  nr <- length(.subset2(.data, 1L))
  pe <- parent.frame()
  cld <- oldClass(.data)
  oldClass(.data) <- NULL
  if(any(cld == "grouped_df")) { # What about NULL assignment with grouped data? and what about across??
    g <- GRP.grouped_df(.data, call = FALSE)
    .data[[".g_"]] <- g
    for(i in 2:length(e)) {
      ei <- e[[i]]
      eiv <- all.vars(ei, functions = TRUE)
      if(any(eiv %in% .FAST_FUN_MOPS)) {
        ei <- fFUN_mutate_add_groups(ei)
        .data[[nam[i]]] <- eval(ei, .data, pe)
      } else {
        v <- all.vars(ei) # Note: Still issue with copyMostAttrib for othFUN_compute when collapse is not attached...
        r <- if(length(v) > 1L) gsplit_multi_apply(.data[v], g, ei, pe) else if(length(eiv) == 2L)
             eval(othFUN_compute(ei), .data, pe) else gsplit_single_apply(.data[[v]], g, ei, pe)
        r <- if(length(r) == g[[1L]]) copyMostAttributes(r[g[[2L]]], r) else # .Call(Cpp_TRA, .data[[v]], r, g[[2L]], 1L) # Faster than simple subset r[g[[2L]] ??]
             greorder(r, g) # r[forder.int(forder.int(g[[2L]]))] # Seems twice is necessary...
        .data[[nam[i]]] <- r
      }
    }
    .data[[".g_"]] <- NULL
  } else { # Without groups...
    for(i in 2:length(e)) { # This is good and very fast
      r <- eval(e[[i]], .data, pe)
      if(!is.null(r)) {
        if(length(r) == 1L) r <- alloc(r, nr)
        else if(length(r) != nr) stop("length mismatch")
      }
      .data[[nam[i]]] <- r
    }
  }
  oldClass(.data) <- cld
  return(condalc(.data, any(cld == "data.table")))
}

mte <- fmutate



# OLD versions and experimental stuff:

# fssm <- function(x, subset) { # not faster than native [ !!
#   ax <- attributes(x)
#   d <- dim(x)
#   ax[["dimnames"]][[1L]] <- ax[["dimnames"]][[1L]][subset]
#   ax[["dim"]] <- c(length(subset), d[2L])
#   ic <- seq_len(d[2L]) * d[1L] - d[1L]
#   setAttributes(.Call(C_subsetVector, x, outer(subset, ic, FUN = "+"), TRUE), ax)
# }

# Older version: But classes for [ can also be very useful for certain objects !!
# fsubset.matrix <- function(x, subset, select, drop = FALSE, ...) {
#   if(!missing(...)) stop("Unknown argument ", dotstostr(...))
#   if(missing(select)) {
#     if(is.object(x)) return(`oldClass<-`(unclass(x)[subset, , drop = drop], oldClass(x))) else
#       return(x[subset, , drop = drop])
#   } else {
#     nl <- as.vector(1L:ncol(x), "list")
#     names(nl) <- dimnames(x)[[2L]]
#     vars <- eval(substitute(select), nl, parent.frame())
#     if(is.object(x)) {
#       if(missing(subset)) return(`oldClass<-`(unclass(x)[, vars, drop = drop], class(x))) else
#         return(`oldClass<-`(unclass(x)[subset, vars, drop = drop], oldClass(x)))
#     } else {
#       if(missing(subset)) return(x[, vars, drop = drop]) else
#         return(x[subset, vars, drop = drop])
#     }
#   }
# }

# older version -> more like base::subset
# fsubset.data.frame <- function(x, subset, select, ...) {
#   if(!missing(...)) stop("Unknown argument ", dotstostr(...))
#   if(missing(select)) vars <- seq_along(unclass(x)) else {
#     nl <- `names<-`(as.vector(seq_along(unclass(x)), "list"), attr(x, "names"))
#     vars <- eval(substitute(select), nl, parent.frame())
#     if(!is.integer(vars)) vars <- if(is.character(vars)) ckmatch(vars, names(nl)) else which(vars)
#   }                   # Best solution ??
#   if(missing(subset)) return(colsubset(x, vars)) else { # if(is.atomic(subset))  # rep_len(TRUE, length(x[[1L]])) else {
#     r <- eval(substitute(subset), x, parent.frame()) #     # e <- substitute(subset) # if(e[[1L]] == ":") ... but what about objects? -> just keep this !!
#     if(is.logical(r)) r <- which(r) # which(r & !is.na(r)) is.na not needed !!
#   } # improve qDF !!!
#   rn <- attr(x, "row.names") # || is.integer(rn) # maybe many have character converted integers ??
#   if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(.Call(C_subsetDT, x, r, vars, TRUE))
#   return(`attr<-`(.Call(C_subsetDT, x, r, vars, TRUE), "row.names", rn[r])) # fast ?? scalable ??
# }


# transform(mtcars, newc = cyl > 5, bla = cyl > 3)

# See also with and within. What about keeping attributes ??
