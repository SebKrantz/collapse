# library(Rcpp)
# sourceCpp("R/C++/psmat.cpp", rebuild = TRUE) # Todo: What about factors and dates ?? Matrix return fastest ??
# sourceCpp("R/C++/qFqG.cpp", rebuild = TRUE)
# qF <- function(x, ordered = TRUE) {
#   if(is.factor(x)) return(x)
#   qFCpp(x, ordered)
# }

# General Note : if add ... can put extra arguments that are not used. if not, gives unuse argument error !!!

psmat <- function(x, ...) { # g, t = NULL, cols = NULL, transpose = FALSE, simplify = FALSE
  UseMethod("psmat", x)
}
psmat.default <- function(x, g, t = NULL, transpose = FALSE, ...) {
  if(is.matrix(x)) stop("x is already a matrix")
  if(is.atomic(g) && length(g) == 1) {
    if(transpose) matrix(x, ncol = round(g), dimnames =
    list(seq_len(length(x)/round(g)), paste0("GRP.",seq_len(g)))) else
    matrix(x, nrow = round(g), byrow = TRUE,
    dimnames = list(paste0("GRP.",seq_len(g)), seq_len(length(x)/round(g))))
  } else {
  if(!is.factor(g)) if(is.atomic(g)) g <- qF(g) else if(all(class(g) == "GRP"))
                    g <- as.factor.GRP(g) else g <- as.factor.GRP(GRP(g)) # interaction(lapply(g, qF))
  if(is.null(t)) {
    message("No timevar provided: Assuming Balanced Panel")
    psmatCpp(x, g, NULL, transpose)
  } else {
    if(!is.factor(t)) if(is.atomic(t)) t <- qF(t) else if(all(class(t) == "GRP"))
                      t <- as.factor.GRP(t) else t <- as.factor.GRP(GRP(t)) # interaction(lapply(t, qF))
    psmatCpp(x, g, t, transpose)
    }
  }
}
psmat.data.frame <- function(x, by, t = NULL, cols = NULL, transpose = FALSE, array = TRUE, ...) {
  if(is.atomic(by) && length(by) == 1) {
    n <- round(by)
    if(transpose) {
      dn <- list(seq_len(length(x)/n), paste0("GRP.",seq_len(by)))
      res <- lapply(x, matrix, ncol = n, dimnames = dn)
    } else {
      dn <- list(paste0("GRP.",seq_len(by)), seq_len(length(x)/round(by)))
      res <- lapply(x, matrix, nrow = n, byrow = TRUE, dimnames = dn)
    }
  } else {
    if(is.call(by)) {
      nam <- names(x)
      if(length(by) == 3) {
        v <- match(all.vars(by[[2L]]), nam)
        by <- match(all.vars(by[[3L]]), nam)
      } else {
        by <- match(all.vars(by), nam)
        v <- if(is.null(cols)) seq_along(x)[-by] else if(is.function(cols))
             setdiff(which(vapply(x, cols, TRUE)), by) else if(is.character(cols))
             match(cols, nam) else cols
      }
      class(x) <- NULL
      by <- if(length(by) == 1) x[[by]] else GRP(x, by, return.groups = FALSE)
      if(is.call(t)) { # If time-variable supplied !!
        t <- match(all.vars(t), nam)
        v <- setdiff(v, t)
        t <- if(length(t) == 1) x[[t]] else GRP(x, t, return.groups = FALSE)
      }
      x <- x[v]
    } else if(!is.null(cols)) {
      class(x) <- NULL
      x <- if(is.function(cols)) x[vapply(x, cols, TRUE)] else x[cols]
    }
    if(!is.factor(by)) if(is.atomic(by)) by <- qF(by) else if(all(class(by) == "GRP"))
      by <- as.factor.GRP(by) else by <- as.factor.GRP(GRP(by)) # interaction(lapply(by, qF))
      if(is.null(t)) {
        message("No timevar provided: Assuming Balanced Panel")
        res <- lapply(x, psmatCpp, by, NULL, transpose)
      } else {
        if(!is.factor(t)) if(is.atomic(t)) t <- qF(t) else if(all(class(t) == "GRP"))
          t <- as.factor.GRP(t) else t <- as.factor.GRP(GRP(t)) # interaction(lapply(t, qF))
        res <- lapply(x, psmatCpp, by, t, transpose)
      }
  }
  if(array) {
    if(length(res) == 1L) return(res[[1L]]) else return(simplify2array(res))
  } else return(res)
}
psmat.pseries <- function(x, transpose = FALSE, ...) {
  index <- attr(x, "index")
  if(is.matrix(x)) stop("x is already a matrix")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  psmatCpp(x, index[[1L]], index[[2L]], transpose)
}
psmat.pdata.frame <- function(x, transpose = FALSE, array = TRUE, ...) {
  index <- attr(x, "index")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  res <- lapply(x, psmatCpp, index[[1L]], index[[2L]], transpose)
  if(array) {
    if(length(res) == 1L) return(res[[1L]]) else return(simplify2array(res))
  } else return(res)
}


# is.balanced.panel()
# is.unsorted.panel <- function(x, g, t) {
# }
#
# index = attr(dGGDC,"index")
# x = index[[1L]]
# y = index[[2L]]
# if (length(x) != length(y))
#   stop("The length of the two vectors differs\n")
# x <- x[drop = TRUE]
# y <- y[drop = TRUE]
# z <- table(x, y)
# if (any(as.vector(z) == 0)) {
#   balanced <- FALSE
# }
# else {
#   balanced <- TRUE
# }
# if (any(as.vector(z) > 1))
#   warning("duplicate couples (id-time)\n")
# return(balanced)
