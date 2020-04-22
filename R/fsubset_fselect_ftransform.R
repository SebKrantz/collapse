

# or rows ? get_rows ??
fsubset <- function(x, ...) UseMethod("fsubset")
# ss <- fsubset # or fss ??

fsubset.default <- function(x, subset, ...) { # also not really faster than default for numeric !! (but a bit faster for factors ...)
  if(!missing(...)) unused_arg_warning(match.call(), ...)
  if(is.logical(subset)) .Call(C_subsetVector, x, which(subset)) else
    .Call(C_subsetVector, x, subset)
}

fsubset.matrix <- function(x, subset, ..., drop = FALSE) {
  if(missing(...)) return(x[subset, , drop = drop]) else { # better row subsetting ??? (like df, method? use mctl ??)
    nl <- `names<-`(as.vector(1L:ncol(x), "list"), dimnames(x)[[2L]])
    vars <- eval(substitute(c(...)), nl, parent.frame()) # better than list(...) ?? -> Yes !!! great !!
    if(missing(subset)) return(x[, vars, drop = drop]) else
      return(x[subset, vars, drop = drop])
  }
}

# fsubset(GGDC10S, Variable == "VA" & Year > 1990, Country, Year, AGR:SUM)

# Also make ss (or fss) -> just indices, no lazy eval ???
ss <- function(data, subset, select) {
  if(missing(select)) select <- seq_along(unclass(data)) else if(is.integer(select)) {
    if(any(select < 0L)) select <- seq_along(unclass(data))[select]
  } else {
    select <- if(is.character(select)) ckmatch(select, attr(data, "names")) else if(is.logical(select))
              which(select) else if(any(select < 0)) seq_along(unclass(data))[select] else as.integer(select)
  }
  if(!is.integer(subset)) subset <- if(is.logical(subset)) which(subset) else as.integer(subset)
  rn <- attr(data, "row.names")
  if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(.Call(C_subsetDT, data, subset, select))
  return(`attr<-`(.Call(C_subsetDT, data, subset, select), "row.names", rn[r]))
}

# make fore flexible selection ???? need vector and matrix methods ?? -> yes !! make more like fselect !!!
# fsubset(GGDC10S, 1) -> error: double, not integer
fsubset.data.frame <- function(x, subset, ...) {
  if(missing(...)) vars <- seq_along(unclass(x)) else {
    ix <- seq_along(unclass(x))
    nl <- `names<-`(as.vector(ix, "list"), attr(x, "names"))
    vars <- eval(substitute(c(...)), nl, parent.frame()) # better than list(...) ?? -> Yes !!! great !!
    if(is.integer(vars)) {
      if(any(vars < 0L)) vars <- ix[vars]
    } else {
      vars <- if(is.character(vars)) ckmatch(vars, names(nl)) else if(is.logical(vars))
               which(vars) else if(any(vars < 0)) ix[vars] else as.integer(vars)
    }
  }
  r <- eval(substitute(subset), x, parent.frame()) #     # e <- substitute(subset) # if(e[[1L]] == ":") ... but what about objects? -> just keep this !!
  if(!is.integer(r)) r <- if(is.logical(r)) which(r) else as.integer(r) # which(r & !is.na(r)) is.na not needed !!
  rn <- attr(x, "row.names") # || is.integer(rn) # maybe many have character converted integers ??????????????????
  if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(.Call(C_subsetDT, x, r, vars))
  return(`attr<-`(.Call(C_subsetDT, x, r, vars), "row.names", rn[r])) # fast ?? scalable ??
}
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
#   rn <- attr(x, "row.names") # || is.integer(rn) # maybe many have character converted integers ??????????????????
#   if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(.Call(C_subsetDT, x, r, vars))
#   return(`attr<-`(.Call(C_subsetDT, x, r, vars), "row.names", rn[r])) # fast ?? scalable ??
# }



# Think: Perhaps you have underestimated the ease and value of non-standard evaluation !!! -> eval checks matches automatically...
# Rehink the use of formulas vs. non-standard eval ... what about operators !!, what about fast fun ??
# Make better grouped_df methods !! No deparse(substitute(...)) !!
# also make fgroup_by and fast pipe !!!





#  microbenchmark(ftransform(mtcars, mpg = 2)) -> beat 53 milliseconds...

# or trans? get_trans?
# short tr or ftr
# ftransform(mtcars, cyl = cyl + 10, vs2 = 1, mpg = NULL)

ftransform <- function(X, ...) { # `_data` ??
  if(!is.list(X)) stop("X needs to be a list of equal length columns or a data.frame")
  ax <- attributes(X) # keep like this ??
  class(X) <- NULL
  e <- eval(substitute(list(...)), X, parent.frame()) # a list of computed values. What about attributes ?????
  le <- lengths(e, FALSE)
  nr <- length(X[[1L]])
  rl <- le == nr # checking if computed values have the right length
  inx <- match(names(e), names(X)) # calling names on a plain list is really fast -> no need to save objects..
  matched <- !is.na(inx)
  if(all(rl)) { # All computed vectors have the right length
    if(any(matched)) X[inx[matched]] <- e[matched]
  } else { # Some do not
    if(any(1L < le & !rl)) stop("Lengths of replacements must be equal to nrow(X) or 1, or NULL to delete columns")
    if(any(le1 <- le == 1L)) e[le1] <- lapply(e[le1], rep, nr) # Length 1 arguments. can use TRA ???, or rep_len, but what about date variables ??
    if(any(le0 <- le == 0L)) { # best order -> yes, ftransform(mtcars, bla = NULL) just returns mtcars, but could also put this error message:
      if(any(le0 & !matched)) stop(paste("Can only delete existing columns, columns",paste(names(e)[le0 & !matched], collapse = ", "),"not found in X"))
      if(all(le0)) {
        X[inx[le0]] <- NULL
        return(`oldClass<-`(X, ax[["class"]]))
      }
      matched <- matched[!le0]
      e <- e[!le0] # e[le0] <- NULL
      if(any(matched)) X[inx[!le0][matched]] <- e[matched] # index is wrong after first deleting, thus we delete after !!
      X[inx[le0]] <- NULL
    } else if(any(matched)) X[inx[matched]] <- e[matched] # NULL assignment ... -> Nope !!
  }
  if(all(matched)) return(`oldClass<-`(X, ax[["class"]]))
  ax[["names"]] <- c(names(X), names(e)[!matched])
  return(setAttributes(c(X, e[!matched]), ax))
}

# also make settransform ?? or add_vars <- fcompute(...)
#settrans
settransform <- function(X, ...) eval.parent(substitute(X <- ftransform(X, ...)))

# short cp or fcp
# compute_vars
fcompute <- function(X, ...) { # within ??
  ax <- attributes(X)
  e <- eval(substitute(list(...)), X, parent.frame())
  ax[["names"]] <- names(e)
  le <- lengths(e, FALSE)
  nr <- length(X[[1L]])
  rl <- le == nr
  if(all(rl)) return(setAttributes(e, ax)) # All computed vectors have the right length
  if(any(1L < le & !rl)) stop("Lengths of replacements must be equal to nrow(X) or 1")
  e[!rl] <- lapply(e[!rl], rep, nr)
  return(setAttributes(e, ax))
}




# fssm <- function(x, subset) { # not faster than native [ !!
#   ax <- attributes(x)
#   d <- dim(x)
#   ax[["dimnames"]][[1L]] <- ax[["dimnames"]][[1L]][subset]
#   ax[["dim"]] <- c(length(subset), d[2L])
#   ic <- seq_len(d[2L]) * d[1L] - d[1L]
#   setAttributes(.Call(C_subsetVector, x, outer(subset, ic, FUN = "+")), ax)
# }

# Older version: But classes for [ can also be very useful for certain objects !!
# fsubset.matrix <- function(x, subset, select, drop = FALSE, ...) {
#   if(!missing(...)) stop("Unknown argument ", dotstostr(...))
#   if(missing(select)) {
#     if(is.object(x)) return(`oldClass<-`(unclass(x)[subset, , drop = drop], class(x))) else
#       return(x[subset, , drop = drop])
#   } else {
#     nl <- as.vector(1L:ncol(x), "list")
#     names(nl) <- dimnames(x)[[2L]]
#     vars <- eval(substitute(select), nl, parent.frame())
#     if(is.object(x)) {
#       if(missing(subset)) return(`oldClass<-`(unclass(x)[, vars, drop = drop], class(x))) else
#         return(`oldClass<-`(unclass(x)[subset, vars, drop = drop], class(x)))
#     } else {
#       if(missing(subset)) return(x[, vars, drop = drop]) else
#         return(x[subset, vars, drop = drop])
#     }
#   }
# }


# transform(mtcars, newc = cyl > 5, bla = cyl > 3)

# See also with and within. What about keeping attributes ??
