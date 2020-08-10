

fsubset <- function(x, ...) UseMethod("fsubset")
sbt <- fsubset

# Also not really faster than default for numeric (but a bit faster for factors ...)
fsubset.default <- function(x, subset, ...) {
  if(!missing(...)) unused_arg_action(match.call(), ...)
  if(is.logical(subset)) return(.Call(C_subsetVector, x, which(subset)))
  .Call(C_subsetVector, x, subset)
}

fsubset.matrix <- function(x, subset, ..., drop = FALSE) {
  if(missing(...)) return(x[subset, , drop = drop])  # better row subsetting ? (like df, method? use mctl ?)
  nl <- `names<-`(as.vector(1L:ncol(x), "list"), dimnames(x)[[2L]])
  vars <- eval(substitute(c(...)), nl, parent.frame()) # better than list(...) ? -> Yes, great !
  if(missing(subset)) return(x[, vars, drop = drop])
  x[subset, vars, drop = drop]
}

# Just indices, no lazy eval
ss <- function(data, i, j) {
  if(missing(j)) j <- seq_along(unclass(data)) else if(is.integer(j)) {
    if(any(j < 0L)) j <- seq_along(unclass(data))[j]
  } else {
    j <- if(is.character(j)) ckmatch(j, attr(data, "names")) else if(is.logical(j))
         which(j) else if(any(j < 0)) seq_along(unclass(data))[j] else as.integer(j)
  }
  if(!is.integer(i)) i <- if(is.logical(i)) which(i) else as.integer(i)
  rn <- attr(data, "row.names")
  if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(.Call(C_subsetDT, data, i, j))
  return(`attr<-`(.Call(C_subsetDT, data, i, j), "row.names", rn[i]))
}

fsubset.data.frame <- function(x, subset, ...) {
  if(missing(...)) vars <- seq_along(unclass(x)) else {
    ix <- seq_along(unclass(x))
    nl <- `names<-`(as.vector(ix, "list"), attr(x, "names"))
    vars <- eval(substitute(c(...)), nl, parent.frame())
    if(is.integer(vars)) {
      if(any(vars < 0L)) vars <- ix[vars]
    } else {
      vars <- if(is.character(vars)) ckmatch(vars, names(nl)) else if(is.logical(vars))
               which(vars) else if(any(vars < 0)) ix[vars] else as.integer(vars)
    }
  }
  r <- eval(substitute(subset), x, parent.frame()) # e <- substitute(subset) # if(e[[1L]] == ":") ... but what about objects? -> just keep this
  if(!is.integer(r)) r <- if(is.logical(r)) which(r) else as.integer(r) # which(r & !is.na(r)) not needed !
  rn <- attr(x, "row.names") # || is.integer(rn) # maybe many have character converted integers ?
  if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(.Call(C_subsetDT, x, r, vars))
  return(`attr<-`(.Call(C_subsetDT, x, r, vars), "row.names", rn[r])) # fast ?? scalable ??
}

# Example:
# fsubset(GGDC10S, Variable == "VA" & Year > 1990, Country, Year, AGR:SUM)

ftransform <- function(X, ...) { # `_data` ?
  if(!is.list(X)) stop("X needs to be a list of equal length columns or a data.frame")
  ax <- attributes(X) # keep like this ?
  oldClass(X) <- NULL
  e <- eval(substitute(list(...)), X, parent.frame()) # a list of computed values. What about attributes ?
  if(length(e) == 1L && is.list(e[[1L]])) e <- unclass(e[[1L]]) # support list input -> added in v1.3.0
  nam <- names(e)
  if(is.null(nam) || any(nam == "")) stop("all expressions have to be named")
  le <- lengths(e, FALSE)
  nr <- length(X[[1L]])
  rl <- le == nr # checking if computed values have the right length
  inx <- match(nam, names(X)) # calling names on a plain list is really fast -> no need to save objects..
  matched <- !is.na(inx)
  if(all(rl)) { # All computed vectors have the right length
    if(any(matched)) X[inx[matched]] <- e[matched]
  } else { # Some do not
    if(any(1L < le & !rl)) stop("Lengths of replacements must be equal to nrow(X) or 1, or NULL to delete columns")
    if(any(le1 <- le == 1L)) e[le1] <- lapply(e[le1], rep, nr) # Length 1 arguments. can use TRA ?, or rep_len, but what about date variables ?
    if(any(le0 <- le == 0L)) { # best order -> yes, ftransform(mtcars, bla = NULL) just returns mtcars, but could also put this error message:
      if(any(le0 & !matched)) stop(paste("Can only delete existing columns, columns",paste(names(e)[le0 & !matched], collapse = ", "),"not found in X"))
      if(all(le0)) {
        X[inx[le0]] <- NULL
        return(`oldClass<-`(X, ax[["class"]]))
      }
      matched <- matched[!le0]
      e <- e[!le0] # e[le0] <- NULL
      if(any(matched)) X[inx[!le0][matched]] <- e[matched] # index is wrong after first deleting, thus we delete after !
      X[inx[le0]] <- NULL
    } else if(any(matched)) X[inx[matched]] <- e[matched] # NULL assignment ... -> Nope !
  }
  if(all(matched)) return(`oldClass<-`(X, ax[["class"]]))
  ax[["names"]] <- c(names(X), names(e)[!matched])
  setAttributes(c(X, e[!matched]), ax)
}
tfm <- ftransform # of trfm ? trf is easiest to type... Lets go with consistency and take first, middle and end consonant

# Note: Only edit the code for ftransform, not this one...
# About 4 microseconds faster than above... and more secure because standard eval..
`ftransform<-` <- function(X, value) {
  if(!is.list(X)) stop("X needs to be a list of equal length columns or a data.frame")
  if(!is.list(value)) stop("value needs to be a named list")
  ax <- attributes(X)
  oldClass(X) <- NULL
  oldClass(value) <- NULL
  nam <- names(value)
  if(is.null(nam) || any(nam == "")) stop("all list-elements have to be named")
  le <- lengths(value, FALSE)
  nr <- length(X[[1L]])
  rl <- le == nr
  inx <- match(nam, names(X))
  matched <- !is.na(inx)
  if(all(rl)) {
    if(any(matched)) X[inx[matched]] <- value[matched]
  } else {
    if(any(1L < le & !rl)) stop("Lengths of replacements must be equal to nrow(X) or 1, or NULL to delete columns")
    if(any(le1 <- le == 1L)) value[le1] <- lapply(value[le1], rep, nr)
    if(any(le0 <- le == 0L)) {
      if(any(le0 & !matched)) stop(paste("Can only delete existing columns, columns",paste(names(value)[le0 & !matched], collapse = ", "),"not found in X"))
      if(all(le0)) {
        X[inx[le0]] <- NULL
        return(`oldClass<-`(X, ax[["class"]]))
      }
      matched <- matched[!le0]
      value <- value[!le0]
      if(any(matched)) X[inx[!le0][matched]] <- value[matched]
      X[inx[le0]] <- NULL
    } else if(any(matched)) X[inx[matched]] <- value[matched]
  }
  if(all(matched)) return(`oldClass<-`(X, ax[["class"]]))
  ax[["names"]] <- c(names(X), names(value)[!matched])
  setAttributes(c(X, value[!matched]), ax)
}
`tfm<-` <- `ftransform<-`

# Example:
# ftransform(mtcars, cyl = cyl + 10, vs2 = 1, mpg = NULL)

settransform <- function(X, ...) eval.parent(substitute(X <- ftransform(X, ...))) # can use `<-`(X, ftransform(X,...)) but not faster ..
# settrans settrfm -> settrf is easiest to type
settfm <- settransform

# compute_vars
fcompute <- function(X, ...) { # within ?
  ax <- attributes(X)
  e <- eval(substitute(list(...)), X, parent.frame())
  if(length(e) == 1L && is.list(e[[1L]])) e <- unclass(e[[1L]]) # support list input -> added in v1.3.0 # sensible ??? what application ??
  ax[["names"]] <- names(e)
  le <- lengths(e, FALSE)
  nr <- fnrow2(X)
  rl <- le == nr
  if(all(rl)) return(setAttributes(e, ax)) # All computed vectors have the right length
  if(any(1L < le & !rl)) stop("Lengths of replacements must be equal to nrow(X) or 1")
  e[!rl] <- lapply(e[!rl], rep, nr)
  setAttributes(e, ax)
}



# OLD versions and experimental stuff:

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
#   if(is.numeric(rn) || is.null(rn) || rn[1L] == "1") return(.Call(C_subsetDT, x, r, vars))
#   return(`attr<-`(.Call(C_subsetDT, x, r, vars), "row.names", rn[r])) # fast ?? scalable ??
# }


# transform(mtcars, newc = cyl > 5, bla = cyl > 3)

# See also with and within. What about keeping attributes ??
