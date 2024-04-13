
# TODO: More tests for attribute handling + Optimize linear fitting...

demean <- function(x, fl, weights, ..., means = FALSE) {
  if(length(fl) == 1L && is.null(attr(fl, "slope.flag"))) {
    clx <- oldClass(x) # Need to do this because could call fbetween.grouped_df of fbetween.pseries / pdata.frame
    if(means) return(`oldClass<-`(fbetween(unclass(x), fl[[1L]], weights, na.rm = FALSE), clx)) else
      return(`oldClass<-`(fwithin(unclass(x), fl[[1L]], weights, na.rm = FALSE), clx))
  }
  msg <- "For higher-dimensional centering and projecting out interactions need to install.packages('%s'), then unload [detach('package:collapse', unload = TRUE)] and reload [library(collapse)]."
  res <- getenvFUN("fixest_demean", msg)(x, fl, attr(fl, "slope.vars"), attr(fl, "slope.flag"),
                                         weights = weights, ..., notes = FALSE, im_confident = TRUE)
  if(!means) return(duplAttributes(res, x))
    # if(!is.matrix(x)) dim(res) <- NULL # also need for flmres... e.g. with weights... intercept is no longer always added, so res needs to be a matrix...
    # Need matrix dimensions... for subset in variable.wise... do.call(cbind, fl[!fc]) needs to be preserved... # return(if(means) x - drop(res) else drop(res))
  if(is.atomic(res)) return(duplAttributes(x - res, x))
  duplAttributes(.mapply(`-`, list(unattrib(x), unattrib(res)), NULL), x)
}

myModFrame <- function(f, data) {
  t <- terms.formula(f)
  v <- attr(t, "variables")
  res <- eval(v, data, parent.frame()) # faster than res <- eval(substitute(with(data, e), list(e = v)))
  attributes(res) <- list(names = as.character(v[-1L]),
                          row.names = .set_row_names(fnrow(data)),
                          class = "data.frame",
                          terms = t)
  res
}
# Example:
# mf <- myModFrame( ~ factor(cyl)*poly(carb, 2) + factor(cyl):factor(vs) + factor(cyl):factor(vs):wt + factor(cyl):mpg + factor(am) + factor(hp > 146):qsec + vs + carb:am, data = mtcars)
# mf <- myModFrame( ~ factor(cyl)*poly(carb, 2) + factor(cyl):factor(vs) + factor(cyl):mpg + factor(am) + factor(hp > 146):qsec + vs + carb:am, data = mtcars)

finteract <- function(x, facts, mf) { # x and facts are logical
  f <- which(x & facts)
  if(length(f) == 1L) mf[[f]] else if(length(f) == 2L) do.call(`:`, mf[f]) else
    as_factor_GRP(GRP.default(mf[f], call = FALSE))
}

slinteract <- function(sl, facts, mf) { # sl and facts are  logical
  sl <- which(sl & !facts)
  res <- if(length(sl) == 1L) mf[[sl]] else do.call(`*`, mf[sl])
  if(is.matrix(res)) mctl(res) else list(res)
}

# This is probably the craziest piece of code in the whole package:
# It takes a model.frame as input and computes from it the inputs for both fixest::demean()
# and linear model fitting


getfl <- function(mf) {

  facts <- .Call(C_vtypes, mf, 2L) # vapply(unattrib(mf), is.factor, TRUE)

  # Any factors
  if(any(facts)) {

    terms <- attributes(attr(mf, "terms"))
    clmf <- oldClass(mf)
    oldClass(mf) <- NULL # good ??
    tl <- terms[["term.labels"]]
    factors <- terms[[2L]]
    fctterms <- fcolSums(factors[facts, , drop = FALSE]) > 0
    fctinteract <- fctterms & fcolSums(factors) > 1 # best ??

    # Any interactions involving factors
    if(any(fctinteract)) {
      modelterms <- tl[!fctterms & tl %in% names(which(rowSums(factors) <= 1))]
      single <- rowSums(factors[facts, , drop = FALSE] > 0L) == 1 # These are either single factors or factors only appearing inside an interaction...
      factors <- factors[, fctinteract, drop = FALSE]
      nointeract <- frowSums(factors[facts, , drop = FALSE]) == 0  # These are factors not appearing in interactions
      singlefct <- names(which(single & nointeract)) # better way ??  # tl[fctterms & !fctinteract]
      intterms <-  mctl(factors > 0L, TRUE) # Need names here
      fctfct <- colSums(factors[!facts, , drop = FALSE]) == 0 # These are factor-factor interactions... need names

      fctdat <- NULL # best way to do this ?? or as before with pre-allocation ??
      lsf <- length(singlefct)
      lff <- bsum(fctfct)
      if(lsf) fctdat <- mf[singlefct] # unattrib() -> wrap around at the end... Nah, better with names...
      if(lff) fctdat <- c(fctdat, lapply(intterms[fctfct], finteract, TRUE, mf))

      # Any heterogenous slopes
      if(lff != length(intterms)) {
        intslope <- intterms[!fctfct]
        slflag <- integer(lsf)

        factors <- factors[facts, !fctfct, drop = FALSE]
        dimnames(factors) <- NULL

        # Could have imp:exp and imp:exp:year, so we need to partial match imp:exp in all slope terms...
        imc <- im <- pmatch(names(which(fctfct)), names(intslope), nomatch = 0L) # need names to match here !!
        if(any(im)) { # first the fact:fact in order (only add slopes), then the other ones
          if(!all(im)) im <- im[im > 0L]

          # Check for duplicate factors in interactions (largely independent of the other stuff)
          dupchk <- factors[, -im, drop = FALSE] > 0L # same as intslopes...
          if(any(dupfct <- frowSums(dupchk) > 1)) { # Check for factors with multiple slopes...
            if(bsum(dupfct) > 1L) stop("Cannot currently support multiple factors with multiple slopes...")
            dupfct <- which(dupchk[dupfct, ]) # This accounts for im
            fctdat <- c(fctdat, lapply(c(intslope[-im][dupfct[1L]], intslope[-im][-dupfct]), finteract, facts, mf))
          } else
            fctdat <- c(fctdat, lapply(intslope[-im], finteract, facts, mf)) # only get factors not already in fctfct...

          slopes <- lapply(c(intslope[im], intslope[-im]), slinteract, facts, mf)
          lsl <- lengths(slopes, FALSE) # No names here
          lim <- seq_along(im)
          imc[imc > 0L] <- lsl[lim] # This is ok, these are also included elsewhere
          slflag <- c(slflag, imc)

          if(length(lsl) != length(lim)) { # The other cases... if exist
            othmc <- lsl[-lim]
            if(any(alone <- single & !nointeract)) {
              alone <- fcolSums(factors[alone, -im, drop = FALSE]) > 0 # This finds the terms corresponding to a factor appearing in an interaction but nowhere else..
              othmc[alone] <- -othmc[alone]
            }
            if(any(dupfct)) { # reordering if dupfct... putting it in front..
              slopes[-lim] <- c(slopes[-lim][dupfct], slopes[-lim][-dupfct])
              othmc <- c(bsum(othmc[dupfct]), othmc[-dupfct])
            }
            slflag <- c(slflag, othmc)
          }
          # this shows single factors not interacted... set slflag to negative...
          # what about double interactions only with slope ??? i.e. only imp:exp:year -> also negative flag...

        } else { # No double factor interactions with slopes.. Only simple slopes interactions.. (what about dupfact of two different double interactions with slope, but no factfact?)
          dupchk <- factors > 0L # same as intslopes...
          if(any(dupfct <- frowSums(dupchk) > 1)) { # Check for factors with multiple slopes...
            if(bsum(dupfct) > 1L) stop("Cannot currently support multiple factors with multiple slopes...")
            dupfct <- which(dupchk[dupfct, ])
            fctdat <- c(fctdat, lapply(c(intslope[dupfct[1L]], intslope[-dupfct]), finteract, facts, mf))
          } else fctdat <- c(fctdat, lapply(intslope, finteract, facts, mf))
          slopes <- lapply(intslope, slinteract, facts, mf) # getting slopes, independent of dupfct...
          lsl <- lengths(slopes, FALSE)
          if(any(alone <- single & !nointeract)) { # Any factor occurring only inside an interaction... This is independent of dupfact and thre associated reordering...
            alone <- fcolSums(factors[alone, , drop = FALSE]) > 0
            lsl[alone] <- -lsl[alone]
          }
          if(any(dupfct)) { # reordering if dupfct... putting it in front..
            slopes <- c(slopes[dupfct], slopes[-dupfct])
            lsl <- c(bsum(lsl[dupfct]), lsl[-dupfct])
          }
          slflag <- c(slflag, integer(lff), lsl)
        }
        attr(fctdat, "slope.vars") <- unlist(slopes, recursive = FALSE) # , FALSE, FALSE)
        attr(fctdat, "slope.flag") <- slflag  # c(integer(length(fctdat)-length(intslope)), lengths(slopes)) # what about other slopes (not poly??)
      }
      # drop unused factor levels ??
    } else {
      modelterms <- tl[!fctterms]
      fctdat <- mf[facts]
    }
    slflag <- attr(fctdat, "slope.flag")
    if(length(modelterms)) { # Intercept only needed if facts with only negative slope flag...
      form <- paste0(if(is.null(slflag) || any(slflag > 0L)) "~ -1 + " else "~ ", paste(modelterms, collapse = " + "))
      moddat <- model.matrix.default(as.formula(form), data = `oldClass<-`(mf, clmf))
    } else {
      moddat <- if(is.null(slflag) || any(slflag > 0L)) NULL else
                alloc(1, length(mf[[1L]]))
    }
  } else {
    fctdat <- NULL
    moddat <- model.matrix.default(attr(mf, "terms"), data = mf) # .External2(stats:::C_modelmatrix, attr(mf, "terms"), mf)
  }
  list(fl = fctdat, xmat = moddat)
}

# Keeps attributes ? -> Yes !
# fastest way ? or better use vectors ? -> this is faster than lapply(fl, `[`, cc) !
subsetfl <- function(fl, cc) {
  slopes <- attr(fl, "slope.vars") # fl could be a data.frame, slope vars not (getfl() unclasses)
  if(is.null(names(fl))) names(fl) <- seq_along(unclass(fl))
  if(is.null(slopes)) return(.Call(C_subsetDT, fl, cc, seq_along(unclass(fl)), FALSE))
  attr(fl, "slope.vars") <- NULL
  if(is.null(names(slopes))) names(slopes) <- seq_along(slopes)
  res <- .Call(C_subsetDT, fl, cc, seq_along(fl), FALSE)
  attr(res, "slope.vars") <- .Call(C_subsetDT, slopes, cc, seq_along(slopes), FALSE) # fdroplevels ??
  res
}

# Old version:
# subsetfl <- function(fl, cc) {
#   lapply(fl, function(f) { # use CsubsetDT or CsubsetVector ?? also check NA in regressors ??
#     x <- attr(f, "x")
#     if(is.null(x)) return(.Call(C_subsetVector, f, cc, FALSE)) else
#       return(`attr<-`(.Call(C_subsetVector, f, cc, FALSE), "x",
#                       if(is.matrix(x)) x[cc, , drop = FALSE] else
#                         .Call(C_subsetVector, x, cc, FALSE)))
#   })
# }

# Examples:
# str(getfl(myModFrame( ~ cyl + carb, data = mtcars)))
# str(getfl(myModFrame( ~ factor(cyl)*carb, data = mtcars)))
# str(getfl(myModFrame( ~ factor(cyl) + factor(am), data = mtcars)))
# str(getfl(myModFrame( ~ factor(cyl):factor(am), data = mtcars)))
# str(getfl(myModFrame( ~ mpg + factor(cyl)*carb, data = mtcars)))
# str(getfl(myModFrame( ~ mpg + factor(cyl) + factor(am), data = mtcars)))
# str(getfl(myModFrame( ~ mpg + factor(cyl):factor(am), data = mtcars)))
# str(getfl(myModFrame( ~ mpg + factor(cyl):factor(am):vs, data = mtcars))) # wow !!
# str(getfl(myModFrame( ~ mpg + factor(cyl):factor(am)*vs, data = mtcars))) # wow !!
# str(getfl(myModFrame( ~ mpg + factor(cyl):factor(am) + factor(cyl):factor(am):vs, data = mtcars))) # wow !!
# str(getfl(myModFrame( ~ mpg + factor(cyl):mpg + factor(am):mpg + factor(cyl):factor(am), data = mtcars)))
# str(getfl(model.frame( ~ factor(cyl)*carb + vs + wt:gear + wt:gear:carb, data = mtcars)))

# (Weighted) linear model fitting for vectors and lists...

# Neded to sort out some insufficiencies of base R default functions when dealing with dimensions
`%**%` <- function(x, y) if(length(y) > 1L) x %*% y else x * y
tcrossprod2 <- function(x, y) if(length(x) > 1L) tcrossprod(x, y) else `dim<-`(x * y, c(1L, length(y)))

# y = x; X = xmat; w = w; meth = lm.method
flmres <- function(y, X, w = NULL, meth = "qr", resi = TRUE, ...) {
  # n <- dim(X)[1L]
  # if(n != NROW(y)) stop("NROW(y) must match nrow(X)")
  dimnames(X) <- NULL # faster ??
  if(length(w)) {
    # if(length(w) != n) stop("w must be numeric and length(w) == nrow(X)")
    wts <- sqrt(w)
    if(is.atomic(y)) {
      dimnames(y) <- NULL
      return(drop(switch(meth,
                         qr = {
                           fit <- X %**% qr.coef(qr(X * wts, ...), y * wts) # same as lm...
                           if(resi) y - fit else fit
                         },
                         chol = {
                           fit <- X * wts
                           fit <- X %*% chol2inv(chol(crossprod(fit), ...)) %*% crossprod(fit, y * wts)
                           if(resi) y - fit else fit
                         },
                         stop("Only methods 'qr' and 'chol' are supported"))))
    }
    attributes(y) <- NULL
      return(switch(meth,
                     qr = {
                       calc <- qr(X * wts, ...)
                       if(resi) lapply(y, function(z) drop(z - X %**% qr.coef(calc, z * wts))) else
                                lapply(y, function(z) drop(X %**% qr.coef(calc, z * wts)))
                     },
                     chol = {
                       calc <- X * wts
                       calc <- X %*% tcrossprod2(chol2inv(chol(crossprod(calc), ...)), calc)
                       if(resi) lapply(y, function(z) drop(z - calc %*% (z * wts))) else
                                lapply(y, function(z) drop(calc %*% (z * wts)))
                     },
                     stop("Only methods 'qr' and 'chol' are supported")))
  }
  if(is.atomic(y)) {
    dimnames(y) <- NULL
    return(drop(switch(meth,
           qr = if(resi) qr.resid(qr(X, ...), y) else qr.fitted(qr(X, ...), y),
           chol = {
            fit <- X %*% chol2inv(chol(crossprod(X), ...)) %*% crossprod(X, y)
            if(resi) y - fit else fit
           },
           stop("Only methods 'qr' and 'chol' are supported"))))
  }
  attributes(y) <- NULL
  return(switch(meth,
                qr = {
                  calc <- qr(X, ...)
                  if(resi) lapply(y, function(z) drop(qr.resid(calc, z))) else
                           lapply(y, function(z) drop(qr.fitted(calc, z)))
                },
                chol = {
                  calc <- X %*% tcrossprod2(chol2inv(chol(crossprod(X), ...)), X)
                  if(resi) lapply(y, function(z) drop(z - calc %*% z)) else
                           lapply(y, function(z) drop(calc %*% z))
                },
                stop("Only methods 'qr' and 'chol' are supported")))
}




fhdwithin <- function(x, ...) UseMethod("fhdwithin") # , x

fhdwithin.default <- function(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, lm.method = "qr", ...) {
  # if(is.matrix(x) && !inherits(x, "matrix")) return(fhdwithin.matrix(x, fl, w, na.rm, fill, ...))
  ax <- attributes(x)
  if(na.rm) {
    cc <- complete.cases(x, fl, w) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- whichv(cc, FALSE)
      cc <- which(cc)
      w <- w[cc] # Note this here !!
      if(!fill) {
        if(length(names(x))) ax[["names"]] <- Csv(names(x), cc) # best ??
        x <- Csv(x, cc)
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) {
    fc <- .Call(C_vtypes, fl, 2L) # vapply(unattrib(fl), is.factor, TRUE)
    fcl <- any(fc)
    # if(!fcl && !missing(...)) unused_arg_action(match.call(), ...)
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demean(do.call(cbind, fl[!fc]), fl[fc], w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(alloc(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    # if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(Csv(fl, cc)) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    # if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl = Csv(fl, cc)) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(nallfc || !fcl) {
    if(na.rm && fill) {
      x[-cc] <- NA
      x[cc] <- flmres(if(nallfc) demean(Csv(x, cc), fl, w, ...) else Csv(x, cc), xmat, w, lm.method, ...)
      return(setAttributes(x, ax))
    } else return(setAttributes(flmres(if(nallfc) demean(x, fl, w, ...) else x, xmat, w, lm.method, ...), ax))
  } else if(na.rm && fill) {
    x[-cc] <- NA
    x[cc] <- demean(Csv(x, cc), fl, w, ...)
    return(setAttributes(x, ax))
  } else return(setAttributes(demean(x, fl, w, ...), ax))
}
fhdwithin.pseries <- function(x, effect = "all", w = NULL, na.rm = .op[["na.rm"]], fill = TRUE, ...) {
  if(is.matrix(x)) stop("higher-dimensional centering of matrix pseries is currently not supported. You can use fhdwithin.matrix(x, ix(x), fill = TRUE)")
  ix <- findex(x)
  namix <- attr(ix, "names")
  if(is.character(effect) && length(effect) == 1L && effect == "all") {
    effect <- seq_along(namix)
  } else effect <- cols2int(effect, ix, namix)
  g <- .subset(ix, effect)

  if(na.rm && length(cc <- whichv(x, NA, TRUE)) != length(x)) {
     g <- .Call(C_subsetDT, g, cc, seq_along(g), FALSE) # lapply(g, `[`, cc) -> slower !
    if(fill) {
      x[cc] <- demean(Csv(unattrib(x), cc), g, w[cc], ...) # keeps attributes ?? -> Yes !!
      return(x)
    }
    ax <- attributes(x)
    attributes(x) <- NULL
    xcc <- Csv(x, cc)
    nix <- length(unclass(ix))
    if(nix != length(g)) {
      toss <- seq_len(nix)[-effect]
      reix <- copyMostAttributes(c(.Call(C_subsetDT, ix, cc, toss, FALSE), g)[namix], ix)
    } else reix <- copyMostAttributes(g, ix)
    attr(reix, "row.names") <- .set_row_names(length(cc))
    ax[[if(any(ax$class == "indexed_series")) "index_df" else "index"]] <- reix
    ax$na.rm <- seq_along(x)[-cc]
    if(length(ax$names)) ax$names <- Csv(ax$names, cc)
    res <- setAttributes(demean(xcc, g, w[cc], ...), ax)
  } else res <- demean(x, g, w, ...) # keeps attributes ?? -> Yes !!
  if(is.double(x)) return(res)
  pseries_to_numeric(res)
}

# x = mNA; fl = m; lm.method = "qr"
fhdwithin.matrix <- function(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, lm.method = "qr", ...) {
  ax <- attributes(x)
  if(na.rm) {
    cc <- complete.cases(x, fl, w) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- whichv(cc, FALSE)
      cc <- which(cc)
      w <- w[cc]
      if(!fill) {
        if(length(dimnames(x)[[1L]])) ax[["dimnames"]][[1L]] <- Csv(dimnames(x)[[1L]], cc) # best ??
        ax[["dim"]][1L] <- length(cc)
        x <- x[cc, , drop = FALSE]
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) {
    fc <- .Call(C_vtypes, fl, 2L) # vapply(unattrib(fl), is.factor, TRUE)
    fcl <- any(fc)
    # if(!fcl && !missing(...)) unused_arg_action(match.call(), ...)
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demean(do.call(cbind, fl[!fc]), fl[fc], w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(alloc(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    # if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(Csv(fl, cc)) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    # if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl = Csv(fl, cc)) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(nallfc || !fcl) {
    if(na.rm && fill) {
      x[-cc, ] <- NA # What about weights cc ?????
      x[cc, ] <- flmres(if(nallfc) demean(x[cc, ], fl, w, ...) else x[cc, ], xmat, w, lm.method, ...)
      return(setAttributes(x, ax))
    } else return(setAttributes(flmres(if(nallfc) demean(x, fl, w, ...) else x, xmat, w, lm.method, ...), ax))
  } else if(na.rm && fill) {
    x[-cc, ] <- NA
    x[cc, ] <- demean(x[cc, ], fl, w, ...)
    return(setAttributes(x, ax))
  } else return(setAttributes(demean(x, fl, w, ...), ax))
}

fhdwithin.zoo <- function(x, ...) if(is.matrix(x)) fhdwithin.matrix(x, ...) else fhdwithin.default(x, ...)
fhdwithin.units <- fhdwithin.zoo

# x = collapse:::colsubset(pwlddev, is.numeric)
fhdwithin.pdata.frame <- function(x, effect = "all", w = NULL, na.rm = .op[["na.rm"]], fill = TRUE, variable.wise = TRUE, ...) {
  ix <- findex(x)
  namix <- attr(ix, "names")
  if(is.character(effect) && length(effect) == 1L && effect == "all") {
    effect <- seq_along(namix)
  } else effect <- cols2int(effect, ix, namix)
  g <- .subset(ix, effect)

  if(na.rm && fill && variable.wise) {
    ax <- attributes(x)
    attributes(x) <- NULL
    varwisecomp <- function(x, fl, w, ...) lapply(x, function(y) {
      ycc <- whichv(y, NA, TRUE)
      y[ycc] <- demean(Csv(unattrib(y), ycc), subsetfl(fl, ycc), w[ycc], ...)
      return(y)
    })
    return(setAttributes(varwisecomp(x, g, w, ...), ax))
  } else if(na.rm && any(miss <- missDF(x))) {
    cc <- whichv(miss, FALSE)
    gcc <- .Call(C_subsetDT, g, cc, seq_along(g), FALSE)
    Y <- demean(.Call(C_subsetDT, x, cc, seq_along(unclass(x)), FALSE), gcc, w[cc], ...)
    if(fill) {
      ax <- attributes(x)
      ax[["na.rm"]] <- which(miss)
      return(setAttributes(.Call(C_lassign, Y, fnrow(x), cc, NA_real_), ax))
    }
    attr(Y, "row.names") <- attr(x, "row.names")[cc] # row.names of pdata.frame are special.
    nix <- length(unclass(ix))
    if(nix != length(g)) {
      toss <- seq_len(nix)[-effect]
      reix <- copyMostAttributes(c(.Call(C_subsetDT, ix, cc, toss, FALSE), gcc)[namix], ix)
    } else reix <- copyMostAttributes(gcc, ix)
    attr(reix, "row.names") <- .set_row_names(length(cc))
    attr(Y, if(inherits(x, "indexed_frame")) "index_df" else "index") <- reix
    attr(Y, "na.rm") <- which(miss)
    return(Y)
  } else return(demean(x, g, w, ...)) # setAttributes(, ax) -> Not needed anymore (included in demean())
}

# x = data[5:6]; fl = data[-(5:6)]; variable.wise = TRUE
fhdwithin.data.frame <- function(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, variable.wise = FALSE, lm.method = "qr", ...) {
  ax <- attributes(x)

  if(na.rm) {
    cc <- if(variable.wise) complete.cases(fl, w) else complete.cases(x, fl, w) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- whichv(cc, FALSE)
      cc <- which(cc)
      w <- w[cc]
      if(!variable.wise) {
        if(fill) nrx <- fnrow(x) else if(is.character(ax[["row.names"]]))
          ax[["row.names"]] <- ax[["row.names"]][cc] else ax[["row.names"]] <- .set_row_names(length(cc)) # best ??
        x <- .Call(C_subsetDT, x, cc, seq_along(unclass(x)), FALSE)
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) { # fl is a list !!
    fc <- .Call(C_vtypes, fl, 2L) # vapply(unattrib(fl), is.factor, TRUE)
    fcl <- any(fc)
    # if(!fcl && !missing(...)) unused_arg_action(match.call(), ...)
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demean(do.call(cbind, fl[!fc]), fl[fc], w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(alloc(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
      # if(!missing(...)) unused_arg_action(match.call(), ...)
      xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
      nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(Csv(fl, cc)) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    # if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl = Csv(fl, cc)) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(variable.wise) {
    if(na.rm) { # this means there were mising values in fl, which were already removed!
      return(setAttributes(lapply(unattrib(x), function(y) {
        y[-cc] <- NA # which is not faster !!
        ycc <- whichv(y, NA, TRUE)
        y_cc <- Csv(y, cc)
        YC <- whichv(y_cc, NA, TRUE)
        y_cc <- Csv(y_cc, YC)
        wc <- w[YC]
        y[ycc] <- if(nallfc) flmres(demean(y_cc, subsetfl(fl, YC), wc, ...), xmat[YC, , drop = FALSE], wc, lm.method, ...) else if(fcl)
          demean(y_cc, subsetfl(fl, YC), wc, ...) else flmres(y_cc, xmat[YC, , drop = FALSE], wc, lm.method, ...)
        return(y)
      }), ax))
    }
    return(setAttributes(lapply(unattrib(x), function(y) {
      ycc <- whichv(y, NA, TRUE)
      y_cc <- Csv(y, ycc)
      wc <- w[ycc]
      y[ycc] <- if(nallfc) flmres(demean(y_cc, subsetfl(fl, ycc), wc, ...), xmat[ycc, , drop = FALSE], wc, lm.method, ...) else if(fcl)
        demean(y_cc, subsetfl(fl, ycc), wc, ...) else flmres(y_cc, xmat[ycc, , drop = FALSE], wc, lm.method, ...)
      return(y)
    }), ax)) # Rfast fastlm??
  } else { # at this point missing values are already removed from x and fl !!
    Y <- if(nallfc || !fcl) flmres(if(nallfc) demean(x, fl, w, ...) else x, xmat, w, lm.method, ...) else demean(x, fl, w, ...)
    if(na.rm && fill)  # x[cc, ] <- Y; x[-cc, ] <- NA
      return(setAttributes(.Call(C_lassign, Y, nrx, cc, NA_real_), ax))
    return(setAttributes(Y, ax))
  }
}

fhdwithin.list <- function(x, ...) fhdwithin.data.frame(x, ...)


# Note: could also do Mudlack and add means to second regression -> better than two-times centering ??
HDW <- function(x, ...) UseMethod("HDW") # , x

HDW.default <- function(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, lm.method = "qr", ...) {
  # if(is.matrix(x) && !inherits(x, "matrix")) return(HDW.matrix(x, fl, w, na.rm, fill, lm.method, ...))
  fhdwithin.default(x, fl, w, na.rm, fill, lm.method, ...)
}

HDW.pseries <- function(x, effect = "all", w = NULL, na.rm = .op[["na.rm"]], fill = TRUE, ...)
  fhdwithin.pseries(x, effect, w, na.rm, fill, ...)

HDW.matrix <- function(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, stub = .op[["stub"]], lm.method = "qr", ...) {
  res <- fhdwithin.matrix(x, fl, w, na.rm, fill, lm.method, ...)
  if(isTRUE(stub) || is.character(stub)) return(add_stub(res, if(is.character(stub)) stub else "HDW."))
  res
}

HDW.zoo <- function(x, ...) if(is.matrix(x)) HDW.matrix(x, ...) else HDW.default(x, ...)
HDW.units <- HDW.zoo

# x = mtcars; fl = ~ qF(cyl):carb; w = wdat; stub = FALSE
HDW.data.frame <- function(x, fl, w = NULL, cols = is.numeric, na.rm = .op[["na.rm"]], fill = FALSE,
                           variable.wise = FALSE, stub = .op[["stub"]], lm.method = "qr", ...) {
  if(is.call(fl)) {
    ax <- attributes(x)
    nam <- ax[["names"]]
    if(length(fl) == 3L) {
      fvars <- ckmatch(all.vars(fl[[3L]]), nam)
      Xvars <- ckmatch(all.vars(fl[[2L]]), nam)
      fl[[2L]] <- NULL
    } else {
      fvars <- ckmatch(all.vars(fl), nam)
      Xvars <- cols2intrmgn(fvars, cols, x) # if(length(cols)) fsetdiff(cols2int(cols, x, nam), fvars) else seq_along(unclass(x))[-fvars]
    }
    ax[["names"]] <- do_stub(stub, nam[Xvars], "HDW.")

    if(na.rm) {
      miss <- missDF(x, if(variable.wise) fvars else c(Xvars, fvars))
      if(missw <- length(w) && anyNA(w)) miss <- miss | is.na(w)
      if(missw || any(miss)) {
        ax[["na.rm"]] <- which(miss)
        cc <- whichv(miss, FALSE)
        w <- w[cc]
        if(!variable.wise) if(fill) nrx <- fnrow(x) else if(is.character(ax[["row.names"]]))
          ax[["row.names"]] <- ax[["row.names"]][cc] else ax[["row.names"]] <- .set_row_names(length(cc)) # best ??
      } else na.rm <- FALSE
    }

    xmat <- NULL
    list2env(getfl(myModFrame(fl, if(na.rm) .Call(C_subsetDT, x, cc, fvars, FALSE) else .subset(x, fvars))), envir = environment())
    fcl <- !is.null(fl)
    nallfc <- fcl && !is.null(xmat)
    if(nallfc) xmat <- demean(xmat, fl, w, ...)

    if(variable.wise) {
      if(na.rm) {
        return(setAttributes(lapply(.subset(x, Xvars), function(y) {
          y[-cc] <- NA
          ycc <- whichv(y, NA, TRUE)
          y_cc <- Csv(y, cc)
          YC <- whichv(y_cc, NA, TRUE)
          y_cc <- Csv(y_cc, YC)
          wc <- w[YC]
          y[ycc] <- if(nallfc) flmres(demean(y_cc, subsetfl(fl, YC), wc, ...), xmat[YC, , drop = FALSE], wc, lm.method, ...) else if(fcl)
            demean(y_cc, subsetfl(fl, YC), wc, ...) else  flmres(y_cc, xmat[YC, , drop = FALSE], wc, lm.method, ...)
          return(y)
        }), ax))
      }
      return(setAttributes(lapply(.subset(x, Xvars), function(y) {
        ycc <- whichv(y, NA, TRUE)
        y_cc <- Csv(y, ycc)
        wc <- w[ycc]
        y[ycc] <- if(nallfc) flmres(demean(y_cc, subsetfl(fl, ycc), wc, ...), xmat[ycc, , drop = FALSE], wc, lm.method, ...) else if(fcl)
          demean(y_cc, subsetfl(fl, ycc), wc, ...) else flmres(y_cc, xmat[ycc, , drop = FALSE], wc, lm.method, ...)
        return(y)
      }), ax))
    } else { # at this point missing values are already removed from  fl !!
      Y <- if(na.rm) .Call(C_subsetDT, x, cc, Xvars, FALSE) else .subset(x, Xvars)
      Y <- if(nallfc || !fcl) flmres(if(nallfc) demean(Y, fl, w, ...) else Y, xmat, w, lm.method, ...) else demean(Y, fl, w, ...)
      if(na.rm && fill)  # x[cc, ] <- Y; x[-cc, ] <- NA
        return(setAttributes(.Call(C_lassign, Y, nrx, cc, NA_real_), ax))
      return(setAttributes(Y, ax))
    }
  }  # fl is not a formula !!
 res <- fhdwithin.data.frame(if(is.null(cols)) x else colsubset(x, cols), fl, w, na.rm, fill, variable.wise, lm.method, ...)
 if(isTRUE(stub) || is.character(stub)) return(add_stub(res, if(is.character(stub)) stub else "HDW."))
 res
}

HDW.pdata.frame <- function(x, effect = "all", w = NULL, cols = is.numeric, na.rm = .op[["na.rm"]], fill = TRUE,
                            variable.wise = TRUE, stub = .op[["stub"]], ...) {
  res <- fhdwithin.pdata.frame(fcolsubset(x, cols2intrmgn(which(attr(x, "names") %in% attr(findex(x), "nam")), cols, x)), effect, w, na.rm, fill, variable.wise, ...)
  if(isTRUE(stub) || is.character(stub)) return(add_stub(res, if(is.character(stub)) stub else "HDW."))
  res
}

HDW.list <- function(x, ...) HDW.data.frame(x, ...)

# Theory: y = ?1 x1 + ?2 x2 + e
# FWT: M2 y = ?1 M2 x1 + e so residuals: e = M2 y - ?1 M2 x1 and fitted:
# Now M = I - x(x'x)-1x' = I - P.
# So (I-P2) y = ?1 (I-P2) x1 + e or y - P2 y = ?1 x1 - ?1 P2 x1 + e

# I want y - e = y^ = ?1 x1 + ?2 x2
# so
# P2 y = ?1 P2 x1 + ?2 x2
# Haven't quite figgured it out, but my solution is to just subtract the demeaned data !!

# Note: Only changes to fhdwithin is in the computation part: Perhaps you can combine the code in some better way to reduce code duplication ??


fhdbetween <- function(x, ...) UseMethod("fhdbetween") # , x


fhdbetween.default <- function(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, lm.method = "qr", ...) {
  # if(is.matrix(x) && !inherits(x, "matrix")) return(fhdwithin.matrix(x, fl, w, na.rm, fill, lm.method, ...))
  ax <- attributes(x)
  if(na.rm) {
    cc <- complete.cases(x, fl, w) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- whichv(cc, FALSE)
      cc <- which(cc)
      w <- w[cc] # Note this here !!
      if(!fill) {
        if(length(names(x))) ax[["names"]] <- Csv(names(x), cc) # best ??
        x <- Csv(x, cc)
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) {
    fc <- .Call(C_vtypes, fl, 2L) # vapply(unattrib(fl), is.factor, TRUE)
    fcl <- any(fc)
    # if(!fcl && !missing(...)) unused_arg_action(match.call(), ...)
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demean(do.call(cbind, fl[!fc]), fl[fc], w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(alloc(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    # if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(Csv(fl, cc)) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    # if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl = Csv(fl, cc)) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }
  # Only this part of the code is different from fhdwithin...
  if(nallfc || !fcl) {
    if(na.rm && fill) {
      x[-cc] <- NA
      xcc <- Csv(x, cc)
      x[cc] <- if(nallfc) xcc - flmres(demean(xcc, fl, w, ...), xmat, w, lm.method, ...) else
        flmres(xcc, xmat, w, lm.method, FALSE, ...)
      return(setAttributes(x, ax))
    } else return(setAttributes(if(nallfc)
      x - flmres(demean(x, fl, w, ...), xmat, w, lm.method, ...) else
        flmres(x, xmat, w, lm.method, FALSE, ...), ax))
  } else if(na.rm && fill) {
    x[-cc] <- NA
    x[cc] <- demean(Csv(x, cc), fl, w, ..., means = TRUE)
    return(setAttributes(x, ax))
  } else return(setAttributes(demean(x, fl, w, ..., means = TRUE), ax))
}


fhdbetween.pseries <- function(x, effect = "all", w = NULL, na.rm = .op[["na.rm"]], fill = TRUE, ...)
  fhdwithin.pseries(x, effect, w, na.rm, fill, ..., means = TRUE)

fhdbetween.matrix <- function(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, lm.method = "qr", ...) {
  ax <- attributes(x)
  if(na.rm) {
    cc <- complete.cases(x, fl, w) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- whichv(cc, FALSE)
      cc <- which(cc)
      w <- w[cc]
      if(!fill) {
        if(length(dimnames(x)[[1L]])) ax[["dimnames"]][[1L]] <- dimnames(x)[[1L]][cc] # best ??
        ax[["dim"]][1L] <- length(cc)
        x <- x[cc, , drop = FALSE]
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) {
    fc <- .Call(C_vtypes, fl, 2L) # vapply(unattrib(fl), is.factor, TRUE)
    fcl <- any(fc)
    # if(!fcl && !missing(...)) unused_arg_action(match.call(), ...)
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demean(do.call(cbind, fl[!fc]), fl[fc], w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(alloc(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    # if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(Csv(fl, cc)) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    # if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl = Csv(fl, cc)) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }
  # Only this part of the code is different from fhdwithin...
  if(nallfc || !fcl) {
    if(na.rm && fill) {
      x[-cc, ] <- NA
      xcc <- x[cc, ] # What about weights cc ? -> done above...
      x[cc, ] <- if(nallfc) xcc - flmres(demean(xcc, fl, w, ...), xmat, w, lm.method, ...) else
        flmres(xcc, xmat, w, lm.method, FALSE, ...)
      return(setAttributes(x, ax))
    } else return(setAttributes(if(nallfc)
      x - flmres(demean(x, fl, w, ...), xmat, w, lm.method, ...) else
        flmres(x, xmat, w, lm.method, FALSE, ...), ax))
  } else if(na.rm && fill) {
    x[-cc, ] <- NA
    x[cc, ] <- demean(x[cc, ], fl, w, ..., means = TRUE)
    return(setAttributes(x, ax))
  } else return(setAttributes(demean(x, fl, w, ..., means = TRUE), ax))
}

fhdbetween.zoo <- function(x, ...) if(is.matrix(x)) fhdbetween.matrix(x, ...) else fhdbetween.default(x, ...)
fhdbetween.units <- fhdbetween.zoo

fhdbetween.pdata.frame <- function(x, effect = "all", w = NULL, na.rm = .op[["na.rm"]], fill = TRUE, variable.wise = TRUE, ...)
  fhdwithin.pdata.frame(x, effect, w, na.rm, fill, variable.wise, ..., means = TRUE)


fhdbetween.data.frame <- function(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, variable.wise = FALSE, lm.method = "qr", ...) {
  ax <- attributes(x)

  if(na.rm) {
    cc <- if(variable.wise) complete.cases(fl, w) else complete.cases(x, fl, w) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- whichv(cc, FALSE)
      cc <- which(cc)
      w <- w[cc]
      if(!variable.wise) {
        if(fill) nrx <- fnrow(x) else if(is.character(ax[["row.names"]]))
          ax[["row.names"]] <- ax[["row.names"]][cc] else ax[["row.names"]] <- .set_row_names(length(cc)) # best ??
        x <- .Call(C_subsetDT, x, cc, seq_along(unclass(x)), FALSE)
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) { # fl is a list !!
    fc <- .Call(C_vtypes, fl, 2L) # vapply(unattrib(fl), is.factor, TRUE)
    fcl <- any(fc)
    # if(!fcl && !missing(...)) unused_arg_action(match.call(), ...)
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demean(do.call(cbind, fl[!fc]), fl[fc], w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(alloc(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    # if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(Csv(fl, cc)) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    # if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl = Csv(fl, cc)) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }
  # Only this part of the code is different from fhdwithin !!
  if(variable.wise) {
    if(na.rm) { # this means there were mising values in fl, which were already removed!
      return(setAttributes(lapply(unattrib(x), function(y) {
        y[-cc] <- NA # which is not faster !!
        ycc <- whichv(y, NA, TRUE)
        y_cc <- Csv(y, cc)
        YC <- whichv(y_cc, NA, TRUE)
        y_cc <- Csv(y_cc, YC)
        wc <- w[YC]
        y[ycc] <- if(nallfc) y_cc %-=% flmres(demean(y_cc, subsetfl(fl, YC), wc, ...), xmat[YC, , drop = FALSE], wc, lm.method, ...) else if(fcl)
          demean(y_cc, subsetfl(fl, YC), wc, ..., means = TRUE) else flmres(y_cc, xmat[YC, , drop = FALSE], wc, lm.method, FALSE, ...)
        return(y)
      }), ax))
    }
    return(setAttributes(lapply(unattrib(x), function(y) {
      ycc <- whichv(y, NA, TRUE)
      y_cc <- Csv(y, ycc)
      wc <- w[ycc]
      y[ycc] <- if(nallfc) y_cc %-=% flmres(demean(y_cc, subsetfl(fl, ycc), wc, ...), xmat[ycc, , drop = FALSE], wc, lm.method, ...) else if(fcl)
        demean(y_cc, subsetfl(fl, ycc), wc, ..., means = TRUE) else flmres(y_cc, xmat[ycc, , drop = FALSE], wc, lm.method, FALSE, ...)
      return(y)
    }), ax)) # Rfast fastlm??
  } else { # at this point missing values are already removed from x and fl !!
    if(nallfc || !fcl) {
      Y <- if(nallfc) x %c-% flmres(demean(x, fl, w, ...), xmat, w, lm.method, ...) else flmres(x, xmat, w, lm.method, FALSE, ...)
    } else Y <- demean(x, fl, w, ..., means = TRUE)
    if(na.rm && fill)  # x[cc, ] <- Y; x[-cc, ] <- NA
      return(setAttributes(.Call(C_lassign, Y, nrx, cc, NA_real_), ax))
    return(setAttributes(Y, ax))
  }
}

fhdbetween.list <- function(x, ...) fhdbetween.data.frame(x, ...)


HDB <- function(x, ...) UseMethod("HDB") # , x

HDB.default <- function(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, lm.method = "qr", ...) {
  # if(is.matrix(x) && !inherits(x, "matrix")) return(HDB.matrix(x, fl, w, na.rm, fill, lm.method, ...))
  fhdbetween.default(x, fl, w, na.rm, fill, lm.method, ...)
}

HDB.pseries <- function(x, effect = "all", w = NULL, na.rm = .op[["na.rm"]], fill = TRUE, ...)
  fhdwithin.pseries(x, effect, w, na.rm, fill, ..., means = TRUE)


HDB.matrix <- function(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, stub = .op[["stub"]], lm.method = "qr", ...) {
  res <- fhdbetween.matrix(x, fl, w, na.rm, fill, lm.method, ...)
  if(isTRUE(stub) || is.character(stub)) return(add_stub(res, if(is.character(stub)) stub else "HDB."))
  res
}

HDB.zoo <- function(x, ...) if(is.matrix(x)) HDB.matrix(x, ...) else HDB.default(x, ...)
HDB.units <- HDB.zoo

HDB.data.frame <- function(x, fl, w = NULL, cols = is.numeric, na.rm = .op[["na.rm"]], fill = FALSE,
                           variable.wise = FALSE, stub = .op[["stub"]], lm.method = "qr", ...) {
  if(is.call(fl)) {
    ax <- attributes(x)
    nam <- ax[["names"]]
    if(length(fl) == 3L) {
      fvars <- ckmatch(all.vars(fl[[3L]]), nam)
      Xvars <- ckmatch(all.vars(fl[[2L]]), nam)
      fl[[2L]] <- NULL
    } else {
      fvars <- ckmatch(all.vars(fl), nam)
      Xvars <- cols2intrmgn(fvars, cols, x) #  if(length(cols)) fsetdiff(cols2int(cols, x, nam), fvars) else seq_along(unclass(x))[-fvars]
    }
    ax[["names"]] <- do_stub(stub, nam[Xvars], "HDB.")

    if(na.rm) {
      miss <- missDF(x, if(variable.wise) fvars else c(Xvars, fvars))
      if(missw <- length(w) && anyNA(w)) miss <- miss | is.na(w)
      if(missw || any(miss)) {
        ax[["na.rm"]] <- which(miss)
        cc <- whichv(miss, FALSE)
        w <- w[cc]
        if(!variable.wise) if(fill) nrx <- fnrow(x) else if(is.character(ax[["row.names"]]))
          ax[["row.names"]] <- ax[["row.names"]][cc] else ax[["row.names"]] <- .set_row_names(length(cc)) # best ??
      } else na.rm <- FALSE
    }

    xmat <- NULL
    list2env(getfl(myModFrame(fl, if(na.rm) .Call(C_subsetDT, x, cc, fvars, FALSE) else .subset(x, fvars))), envir = environment())
    fcl <- !is.null(fl)
    nallfc <- fcl && !is.null(xmat)
    if(nallfc) xmat <- demean(xmat, fl, w, ...)


    # Only this part of the code is different from fhdwithin !!
    if(variable.wise) {
      if(na.rm) { # this means there were mising values in fl, which were already removed!
        return(setAttributes(lapply(.subset(x, Xvars), function(y) {
          y[-cc] <- NA # which is not faster !!
          ycc <- whichv(y, NA, TRUE)
          y_cc <- Csv(y, cc)
          YC <- whichv(y_cc, NA, TRUE)
          y_cc <- Csv(y_cc, YC)
          wc <- w[YC]
          y[ycc] <- if(nallfc) y_cc %-=% flmres(demean(y_cc, subsetfl(fl, YC), wc, ...), xmat[YC, , drop = FALSE], wc, lm.method, ...) else if(fcl)
            demean(y_cc, subsetfl(fl, YC), wc, ..., means = TRUE) else flmres(y_cc, xmat[YC, , drop = FALSE], wc, lm.method, FALSE, ...)
          return(y)
        }), ax))
      }
      return(setAttributes(lapply(.subset(x, Xvars), function(y) {
        ycc <- whichv(y, NA, TRUE)
        y_cc <- Csv(y, ycc)
        wc <- w[ycc]
        y[ycc] <- if(nallfc) y_cc %-=% flmres(demean(y_cc, subsetfl(fl, ycc), wc, ...), xmat[ycc, , drop = FALSE], wc, lm.method, ...) else if(fcl)
          demean(y_cc, subsetfl(fl, ycc), wc, ..., means = TRUE) else flmres(y_cc, xmat[ycc, , drop = FALSE], wc, lm.method, FALSE, ...)
        return(y)
      }), ax))
    } else { # at this point missing values are already removed from  fl !!
      x <- if(na.rm) .Call(C_subsetDT, x, cc, Xvars, FALSE) else .subset(x, Xvars)
      if(nallfc || !fcl) {
        Y <- if(nallfc) x %c-% flmres(demean(x, fl, w, ...), xmat, w, lm.method, ...) else flmres(x, xmat, w, lm.method, FALSE, ...)
      } else Y <- demean(x, fl, w, ..., means = TRUE)
      if(na.rm && fill)  # x[cc, ] <- Y; x[-cc, ] <- NA
        return(setAttributes(.Call(C_lassign, Y, nrx, cc, NA_real_), ax))
      return(setAttributes(Y, ax))
    }
  }  # fl is not a formula !!
  res <- fhdbetween.data.frame(if(is.null(cols)) x else colsubset(x, cols), fl, w, na.rm, fill, variable.wise, lm.method, ...)
  if(isTRUE(stub) || is.character(stub)) return(add_stub(res, if(is.character(stub)) stub else "HDB."))
  res
}


HDB.pdata.frame <- function(x, effect = "all", w = NULL, cols = is.numeric, na.rm = .op[["na.rm"]], fill = TRUE,
                            variable.wise = TRUE, stub = .op[["stub"]], ...) {
  res <- fhdwithin.pdata.frame(fcolsubset(x, cols2intrmgn(which(attr(x, "names") %in% attr(findex(x), "nam")), cols, x)), effect, w, na.rm, fill, variable.wise, ..., means = TRUE)
  if(isTRUE(stub) || is.character(stub)) return(add_stub(res, if(is.character(stub)) stub else "HDB."))
  res
}

HDB.list <- function(x, ...) HDB.data.frame(x, ...)

fHDbetween <- function(x, ...) {
  message("Note that 'fHDbetween' was renamed to 'fhdbetween'. The S3 generic will not be removed anytime soon, but please use updated function names in new code, see help('collapse-renamed')")
  UseMethod("fhdbetween")
}
fHDbetween.default <- function(x, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fhdbetween.matrix(x, ...))
  # .Deprecated(msg = "This method belongs to a renamed function and will be removed end of 2022, see help('collapse-renamed')")
  fhdbetween.default(x, ...)
}
fHDbetween.matrix <- function(x, ...) {
  # .Deprecated(msg = "This method belongs to a renamed function and will be removed end of 2022, see help('collapse-renamed')")
  fhdbetween.matrix(x, ...)
}
fHDbetween.data.frame <- function(x, ...) {
  # .Deprecated(msg = "This method belongs to a renamed function and will be removed end of 2022, see help('collapse-renamed')")
  fhdbetween.data.frame(x, ...)
}


fHDwithin <- function(x, ...) {
  message("Note that 'fHDwithin' was renamed to 'fhdwithin'. The S3 generic will not be removed anytime soon, but please use updated function names in new code, see help('collapse-renamed')")
  UseMethod("fhdwithin")
}
fHDwithin.default <- function(x, ...) {
  if(is.matrix(x) && !inherits(x, "matrix")) return(fhdwithin.matrix(x, ...))
  # .Deprecated(msg = "This method belongs to a renamed function and will be removed end of 2022, see help('collapse-renamed')")
  fhdwithin.default(x, ...)
}
fHDwithin.matrix <- function(x, ...) {
  # .Deprecated(msg = "This method belongs to a renamed function and will be removed end of 2022, see help('collapse-renamed')")
  fhdwithin.matrix(x, ...)
}
fHDwithin.data.frame <- function(x, ...) {
  # .Deprecated(msg = "This method belongs to a renamed function and will be removed end of 2022, see help('collapse-renamed')")
  fhdwithin.data.frame(x, ...)
}


#
# HDW(x = mtcars, fl = ~ factor(cyl)*carb)
#
# HDW(x = mtcars, fl = ~ factor(cyl):vs)
#
# lm(mpg ~ factor(cyl):factor(vs), data = mtcars)
#
# HDW(x = mtcars, fl = ~ factor(cyl)*carb + vs + wt:gear + wt:gear:carb)
#
# # Works!! although there is a further interaction with carb!!
# lm(mpg ~ hp, data = HDW(mtcars, ~ factor(cyl)*carb + vs + wt:gear + wt:gear:carb))
# lm(mpg ~ hp + factor(cyl)*carb + vs + wt:gear + wt:gear:carb, data = mtcars)
#
# lm(mpg ~ hp, data = HDW(mtcars, ~ factor(cyl)*carb + vs + wt:gear))
# lm(mpg ~ hp + factor(cyl)*carb + vs + wt:gear, data = mtcars)
#
# lm(mpg ~ hp, data = HDW(mtcars, ~ cyl*carb + vs + wt:gear))
# lm(mpg ~ hp + cyl*carb + vs + wt:gear, data = mtcars)
#
#
# lm(mpg ~ hp, data = HDW(mtcars, mpg + hp ~ cyl*carb + factor(cyl)*poly(drat,2)))
# lm(mpg ~ hp + cyl*carb + factor(cyl)*poly(drat,2), data = mtcars)
#
