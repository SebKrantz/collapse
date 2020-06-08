

myModFrame <- function(f, data) {
  t <- terms.formula(f)
  v <- attr(t, "variables")
  # res <- eval(substitute(with(data, e), list(e = v)))
  res <- eval(v, data, parent.frame()) # faster !
  attributes(res) <- list(names = as.character(v[-1]),
                          row.names = .set_row_names(fnrow2(data)),
                          class = "data.frame",
                          terms = t)
  return(res)
}
# Example:
# mf <- myModFrame( ~ factor(cyl)*carb + factor(cyl):factor(vs) + vs + carb:am, data = mtcars)

getfl <- function(mf) {

  facts <- vapply(unattrib(mf), is.factor, TRUE)

  if(any(facts)) {
    terms <- attributes(attr(mf, "terms"))
    clmf <- class(mf)
    class(mf) <- NULL # good ??
    tl <- terms[["term.labels"]]
    factors <- terms[[2L]]
    fctterms <- colSums(factors[facts, , drop = FALSE]) > 0
    fctinteract <- fctterms & colSums(factors) > 1

    if(any(fctinteract)) { # if any interactions involving factors
      singlefct <- match(tl[fctterms & !fctinteract], names(mf))
      intterms <- lapply(which(fctinteract), function(i) factors[,i] > 0) # names(which( # better way ?
      factors <- factors[!facts, fctinteract, drop = FALSE]
      fctfct <- colSums(factors) == 0
      globalslopes <- tl %in% names(which(rowSums(factors) > 0))
      tvec <- c(length(singlefct), sum(globalslopes), sum(fctfct), sum(fctinteract)-sum(fctfct))
      ctvec <- cumsum(tvec)
      fctterms[globalslopes] <- TRUE
      fctdat <- vector("list", ctvec[4])
      if(tvec[1] != 0) fctdat[seq_along(singlefct)] <- mf[singlefct]
      if(tvec[2] != 0) fctdat[(ctvec[1]+1):ctvec[2]] <- lapply(mf[globalslopes], function(x)
        setAttributes(rep(1L, NROW(x)), list(levels = "1", class = "factor", x = x)))
      if(tvec[3] != 0)
        fctdat[(ctvec[2]+1):ctvec[3]] <- lapply(intterms[fctfct], function(x) if(length(x) == 2L) do.call(`:`, mf[x]) else as.factor.GRP(GRP.default(mf[x]))) # interaction(mf[x])) # or as.factor.GRP(GRP(mf[x]))
      if(tvec[4] != 0)
        fctdat[(ctvec[3]+1):ctvec[4]] <- lapply(intterms[!fctfct], function(x) {
          f <- x & facts
          nf <- x & !f
          f <- if(sum(f) == 1) mf[[which(f)]] else if(sum(f) == 2) do.call(`:`, mf[f]) else as.factor.GRP(GRP.default(mf[f])) # interaction(mf[f])
          attr(f, "x") <- if(sum(nf) > 1) do.call("*", mf[nf]) else mf[[which(nf)]]
          return(f)
        })
      # for (i in intersect(intns,fct)) { # Drop single fctor levels..
      #   tab = table(mf[[i]])<=1
      #   rm[[i]] = which(mf[[i]] %in% names(tab)[tab])
      #   if (length(rm[[i]])) mf = mf[-rm[[i]], , drop = FALSE]
      # }
    } else fctdat <- mf[facts]

    modelterms <- tl[!fctterms]
    if(length(modelterms)) {
      # modelterms <- terms.formula(as.formula(paste0("~ ",paste(modelterms, collapse = " + "))))
      # return(modelterms)
      # moddat <-  model.matrix.default(modelterms, data = mf) #  eval(substitute(with(mf, e), list(e = attr(modelterms, "variables")))))
        # .External2(stats:::C_modelmatrix, modelterms, # removes NA's?? -> Nope !!
        #                   eval(substitute(with(mf, e),
        #                   list(e = attr(modelterms, "variables"))))) # model.matrix.default(as.formula(paste0("~ ",paste(modelterms, collapse = " + "))), data = mf) else moddat <- NULL
      moddat <- model.matrix.default(as.formula(paste0("~ ",paste(modelterms, collapse = " + "))), data = `oldClass<-`(mf, clmf))
    } else moddat <- NULL

  } else {
    fctdat <- NULL
    moddat <- model.matrix.default(attr(mf, "terms"), data = mf) # .External2(stats:::C_modelmatrix, attr(mf, "terms"), mf)
  }
  return(list(fl = fctdat, xmat = moddat))
}

subsetfl <- function(fl, cc) {
  lapply(fl, function(f) { # use CsubsetDT or CsubsetVector ?? also check NA in regressors ??
    x <- attr(f, "x")
    if(is.null(x)) return(.Call(C_subsetVector, f, cc)) else
      return(`attr<-`(.Call(C_subsetVector, f, cc), "x",
                      if(is.matrix(x)) x[cc, , drop = FALSE] else
                        .Call(C_subsetVector, x, cc)))
  })
}
# Examples:
# getfl(model.frame( ~ cyl + carb, data = mtcars))
# getfl(model.frame( ~ factor(cyl)*carb, data = mtcars))
# getfl(model.frame( ~ factor(cyl)*carb + vs + wt:gear + wt:gear:carb, data = mtcars))

fHDwithin <- function(x, ...) UseMethod("fHDwithin") # , x

fHDwithin.default <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  ax <- attributes(x)
  if(na.rm) {
    cc <- complete.cases(x, fl) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- which(!cc)
      cc <- which(cc)
      if(!fill) {
        if(!is.null(names(x))) ax[["names"]] <- names(x)[cc] # best ??
        x <- x[cc]
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) {
    fc <- vapply(unattrib(fl), is.factor, TRUE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) unused_arg_action(match.call(), ...)
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(unclass(fl))) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demeanlist(do.call(cbind, fl[!fc]), fl[fc], weights = w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(rep(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(fl[cc]) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(nallfc || !fcl) {
    xmat <- qr.default(xmat)
    if(na.rm && fill) {
      x[-cc] <- NA
      x[cc] <- if(nallfc) qr.resid(xmat, demeanlist(x[cc], fl, weights = w, ...)) else qr.resid(xmat, x[cc])
      return(setAttributes(x, ax))
    } else if(nallfc)
      return(setAttributes(qr.resid(xmat, demeanlist(x, fl, weights = w, ...)), ax)) else
        return(setAttributes(qr.resid(xmat, x), ax))
  } else if(na.rm && fill) {
    x[-cc] <- NA
    x[cc] <- demeanlist(x[cc], fl, weights = w, ...)
    return(setAttributes(x, ax))
  } else return(setAttributes(demeanlist(x, fl, weights = w, ...), ax))
}
fHDwithin.pseries <- function(x, w = NULL, na.rm = TRUE, fill = TRUE, ...) {
  if(na.rm && !all(xcc <- !is.na(x))) {
     g <- lapply(attr(x, "index"), function(y) y[xcc]) # good ! faster than subsetDT
    if(fill) {
      x[xcc] <- demeanlist(x[xcc], g, weights = w[xcc], ...) # keeps attributes ?? -> Yes !!
      return(x)
    } else return(addAttributes(demeanlist(x[xcc], g, weights = w[xcc], ...),
                                list(index = g, na.rm = which(!xcc)))) # keeps attributes ?? -> Nope !!
  } else {
    g <- attr(x, "index") # what about cases ?? -> nah, named !!
    return(`attr<-`(demeanlist(x, g, weights = w, ...), "index", g)) # keeps attributes ?? -> Nope !!
  }
}
fHDwithin.matrix <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  ax <- attributes(x)
  if(na.rm) {
    cc <- complete.cases(x, fl) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- which(!cc)
      cc <- which(cc)
      if(!fill) {
        if(length(dimnames(x)[[1L]])) ax[["dimnames"]][[1L]] <- dimnames(x)[[1L]][cc] # best ??
        ax[["dim"]][1L] <- length(cc)
        x <- x[cc, , drop = FALSE]
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) {
    fc <- vapply(unattrib(fl), is.factor, TRUE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) unused_arg_action(match.call(), ...)
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(unclass(fl))) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demeanlist(do.call(cbind, fl[!fc]), fl[fc], weights = w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(rep(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(fl[cc]) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(nallfc || !fcl) {
    xmat <- qr.default(xmat)
    if(na.rm && fill) {
      x[-cc, ] <- NA
      x[cc, ] <- if(nallfc) qr.resid(xmat, demeanlist(x[cc, ], fl, weights = w, ...)) else qr.resid(xmat, x[cc, ])
      return(setAttributes(x, ax))
    } else if(nallfc)
      return(setAttributes(qr.resid(xmat, demeanlist(x, fl, weights = w, ...)), ax)) else
        return(setAttributes(qr.resid(xmat, x), ax))
  } else if(na.rm && fill) {
    x[-cc, ] <- NA
    x[cc, ] <- demeanlist(x[cc, ], fl, weights = w, ...)
    return(setAttributes(x, ax))
  } else return(setAttributes(demeanlist(x, fl, weights = w, ...), ax))
}
fHDwithin.pdata.frame <- function(x, w = NULL, na.rm = TRUE, fill = TRUE, variable.wise = TRUE, ...) {
  ax <- attributes(x)
  if(na.rm && fill && variable.wise) {
    attributes(x) <- NULL
    varwisecomp <- function(x, fl, w, ...) lapply(x, function(y) {
      ycc <- which(!is.na(y))
      y[ycc] <- demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)
      return(y)
    })
    return(setAttributes(varwisecomp(x, ax[["index"]], w, ...), ax))
  } else if(na.rm && any(miss <- .Call(C_dt_na, x, seq_along(unclass(x))))) {
      ax[["na.rm"]] <- which(miss)
      cc <- which(!miss)
      Y <- demeanlist(.Call(C_subsetDT, x, cc, seq_along(unclass(x))),
                            lapply(ax[["index"]], function(y) y[cc]),
                            weights = w[cc], ...)
    if(fill) return(setAttributes(.Call(Cpp_lassign, Y, fnrow2(x), cc, NA), ax)) else {
      ax[["row.names"]] <- ax[["row.names"]][cc]
      return(setAttributes(Y, ax))
    }
  } else return(setAttributes(demeanlist(x, ax[["index"]], weights = w, ...), ax))
}
fHDwithin.data.frame <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
  ax <- attributes(x)

  if(na.rm) {
    cc <- if(variable.wise) complete.cases(fl) else complete.cases(x, fl) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- which(!cc)
      cc <- which(cc)
      if(!is.null(w)) w <- w[cc]
      if(!variable.wise) {
        if(fill) nrx <- fnrow2(x) else ax[["row.names"]] <- ax[["row.names"]][cc] # best ??
        x <- .Call(C_subsetDT, x, cc, seq_along(unclass(x)))
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) { # fl is a list !!
    fc <- vapply(unattrib(fl), is.factor, TRUE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) unused_arg_action(match.call(), ...)
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(unclass(fl))) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demeanlist(do.call(cbind, fl[!fc]), fl[fc], weights = w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(rep(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
      if(!missing(...)) unused_arg_action(match.call(), ...)
      xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
      nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(fl[cc]) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(variable.wise) {
    if(na.rm) { # this means there were mising values in fl, which were already removed!
      return(setAttributes(lapply(unattrib(x), function(y) {
        y[-cc] <- NA # which is not faster !!
        ycc <- !is.na(y)
        YC <- which(ycc[cc])
        y[ycc] <- if(nallfc) qr.resid(qr.default(xmat[YC, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], ...)) else if(fcl)
          demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], ...) else  qr.resid(qr.default(xmat[YC, , drop = FALSE]), y[ycc])
        return(y)
      }), ax)) # Rfast fastlm??
    } else {
      return(setAttributes(lapply(unattrib(x), function(y) {
        ycc <- which(!is.na(y))
        y[ycc] <- if(nallfc) qr.resid(qr.default(xmat[ycc, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)) else if(fcl)
          demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...) else qr.resid(qr.default(xmat[ycc, , drop = FALSE]), y[ycc])
        return(y)
      }), ax)) # Rfast fastlm??
    }
  } else { # at this point missing values are already removed from x and fl !!
    if(nallfc || !fcl) {
      xmat <- qr.default(xmat)
      Y <- if(nallfc) lapply(demeanlist(x, fl, weights = w, ...), function(y) qr.resid(xmat, y)) else
           lapply(x, function(y) qr.resid(xmat, y))
    } else Y <- demeanlist(x, fl, weights = w, ...)
    if(na.rm && fill) { # x[cc, ] <- Y; x[-cc, ] <- NA
      return(setAttributes(.Call(Cpp_lassign, Y, nrx, cc, NA), ax))
    } else return(setAttributes(Y, ax))
  }
}


# Note: could also do Mudlack and add means to second regression -> better than two-times centering ??
HDW <- function(x, ...) UseMethod("HDW") # , x

HDW.default <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...)
  fHDwithin.default(x, fl, w, na.rm, fill, ...)

HDW.pseries <- function(x, w = NULL, na.rm = TRUE, fill = TRUE, ...)
  fHDwithin.pseries(x, w, na.rm, fill, ...)

HDW.matrix <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, stub = "HDW.", ...)
  add_stub(fHDwithin.matrix(x, fl, w, na.rm, fill, ...), stub)



HDW.data.frame <- function(x, fl, w = NULL, cols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, stub = "HDW.", ...) {
  if(is.call(fl)) {
    ax <- attributes(x)
    nam <- ax[["names"]]
    if(length(fl) == 3L) {
      fvars <- ckmatch(all.vars(fl[[3L]]), nam)
      Xvars <- ckmatch(all.vars(fl[[2L]]), nam)
      fl[[2L]] <- NULL
    } else {
      fvars <- ckmatch(all.vars(fl), nam)
      Xvars <- if(!is.null(cols)) fsetdiff(cols2int(cols, x, nam), fvars) else seq_along(unclass(x))[-fvars]
    }
    ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[Xvars]) else nam[Xvars]

    if(na.rm) {
      miss <- if(variable.wise) .Call(C_dt_na, x, fvars) else .Call(C_dt_na, x, c(Xvars,fvars))
      if(any(miss)) {
        ax[["na.rm"]] <- which(miss)
        cc <- which(!miss)
        if(!is.null(w)) w <- w[cc]
        if(!variable.wise) if(fill) nrx <- fnrow2(x) else ax[["row.names"]] <- ax[["row.names"]][cc] # best ??
      } else na.rm <- FALSE
    }

    xmat <- NULL
    list2env(getfl(myModFrame(fl, if(na.rm) .Call(C_subsetDT, x, cc, fvars) else .subset(x, fvars))), envir = environment())
    fcl <- !is.null(fl)
    nallfc <- fcl && !is.null(xmat)
    if(nallfc) xmat <- demeanlist(xmat[, -1L, drop = FALSE], fl, weights = w, ...)

    if(variable.wise) {
      if(na.rm) {
        return(setAttributes(lapply(.subset(x, Xvars), function(y) {
          y[-cc] <- NA
          ycc <- !is.na(y)
          YC <- which(ycc[cc])
          y[ycc] <- if(nallfc) qr.resid(qr.default(xmat[YC, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], ...)) else if(fcl)
            demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], ...) else  qr.resid(qr.default(xmat[YC, , drop = FALSE]), y[ycc])
          return(y)
        }), ax)) # Rfast fastlm??
      } else {
        return(setAttributes(lapply(.subset(x, Xvars), function(y) {
          ycc <- which(!is.na(y))
          y[ycc] <- if(nallfc) qr.resid(qr.default(xmat[ycc, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)) else if(fcl)
            demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...) else qr.resid(qr.default(xmat[ycc, , drop = FALSE]), y[ycc])
          return(y)
        }), ax)) # Rfast fastlm??
      }
    } else { # at this point missing values are already removed from  fl !!
      # x <- if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else .subset(x, Xvars) # incorporated below -> more memory efficient
      if(nallfc || !fcl) {
        xmat <- qr.default(xmat)
        Y <- if(nallfc) lapply(demeanlist(if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else .subset(x, Xvars), fl, weights = w, ...), function(y) qr.resid(xmat, y)) else
          lapply(if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else .subset(x, Xvars), function(y) qr.resid(xmat, y))
      } else Y <- demeanlist(if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else .subset(x, Xvars), fl, weights = w, ...)
      if(na.rm && fill) { # x[cc, ] <- Y; x[-cc, ] <- NA
        return(setAttributes(.Call(Cpp_lassign, Y, nrx, cc, NA), ax))
      } else return(setAttributes(Y, ax))
    }
  } else # fl is not a formula !!
 return(add_stub(fHDwithin.data.frame(if(is.null(cols)) x else colsubset(x, cols), fl, w, na.rm, fill, variable.wise, ...), stub))
}

HDW.pdata.frame <- function(x, w = NULL, cols = is.numeric, na.rm = TRUE, fill = TRUE,
                            variable.wise = TRUE, stub = "HDW.", ...)
    return(add_stub(fHDwithin.pdata.frame(if(is.null(cols)) x else colsubset(x, cols), w, na.rm, fill, variable.wise, ...), stub))


# Theory: y = ?1 x1 + ?2 x2 + e
# FWT: M2 y = ?1 M2 x1 + e so residuals: e = M2 y - ?1 M2 x1 and fitted:
# Now M = I - x(x'x)-1x' = I - P.
# So (I-P2) y = ?1 (I-P2) x1 + e or y - P2 y = ?1 x1 - ?1 P2 x1 + e

# I want y - e = y^ = ?1 x1 + ?2 x2
# so
# P2 y = ?1 P2 x1 + ?2 x2
# Haven't quite figgured it out, but my solution is to just subtract the demeaned data !!

# Note: Only changes to fHDwithin is in the computation part: Perhaps you can combine the code in some better way to reduce code duplication ??

fHDbetween <- function(x, ...) UseMethod("fHDbetween") # , x

fHDbetween.default <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  ax <- attributes(x)
  if(na.rm) {
    cc <- complete.cases(x, fl) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- which(!cc)
      cc <- which(cc)
      if(!fill) {
        if(!is.null(names(x))) ax[["names"]] <- names(x)[cc] # best ??
        x <- x[cc]
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) {
    fc <- vapply(unattrib(fl), is.factor, TRUE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) unused_arg_action(match.call(), ...)
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(unclass(fl))) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demeanlist(do.call(cbind, fl[!fc]), fl[fc], weights = w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(rep(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(fl[cc]) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(nallfc || !fcl) {
    xmat <- qr.default(xmat)
    if(na.rm && fill) {
      x[-cc] <- NA
      x[cc] <- if(nallfc) x[cc] - qr.resid(xmat, demeanlist(x[cc], fl, weights = w, ...)) else qr.fitted(xmat, x[cc])
      return(setAttributes(x, ax))
    } else if(nallfc)
      return(setAttributes(x - qr.resid(xmat, demeanlist(x, fl, weights = w, ...)), ax)) else
      return(setAttributes(qr.fitted(xmat, x), ax))
  } else if(na.rm && fill) {
    x[-cc] <- NA
    x[cc] <- demeanlist(x[cc], fl, weights = w, means = TRUE, ...)
    return(setAttributes(x, ax))
  } else return(setAttributes(demeanlist(x, fl, weights = w, means = TRUE, ...), ax))
}
fHDbetween.pseries <- function(x, w = NULL, na.rm = TRUE, fill = TRUE, ...) {
  if(na.rm && !all(xcc <- !is.na(x))) {
    g <- lapply(attr(x, "index"), function(y) y[xcc]) # good ! faster than subsetDT
    if(fill) {
      x[xcc] <- demeanlist(x[xcc], g, weights = w[xcc], means = TRUE, ...) # keeps attributes ?? -> Yes !!
      return(x)
    } else return(addAttributes(demeanlist(x[xcc], g, weights = w[xcc], means = TRUE, ...),
                                list(index = g, na.rm = which(!xcc)))) # keeps attributes ?? -> Nope !!
  } else {
    g <- attr(x, "index") # what about cases ?? -> nah, named !!
    return(`attr<-`(demeanlist(x, g, weights = w, means = TRUE, ...), "index", g)) # keeps attributes ?? -> Nope !!
  }
}
fHDbetween.matrix <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  ax <- attributes(x)
  if(na.rm) {
    cc <- complete.cases(x, fl) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- which(!cc)
      cc <- which(cc)
      if(!fill) {
        if(length(dimnames(x)[[1L]])) ax[["dimnames"]][[1L]] <- dimnames(x)[[1L]][cc] # best ??
        ax[["dim"]][1L] <- length(cc)
        x <- x[cc, , drop = FALSE]
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) {
    fc <- vapply(unattrib(fl), is.factor, TRUE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) unused_arg_action(match.call(), ...)
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(unclass(fl))) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demeanlist(do.call(cbind, fl[!fc]), fl[fc], weights = w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(rep(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(fl[cc]) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(nallfc || !fcl) {
    xmat <- qr.default(xmat)
    if(na.rm && fill) {
      x[-cc, ] <- NA
      x[cc, ] <- if(nallfc) x[cc, ] - qr.resid(xmat, demeanlist(x[cc, ], fl, weights = w, ...)) else qr.fitted(xmat, x[cc, ])
      return(setAttributes(x, ax))
    } else if(nallfc)
      return(setAttributes(x - qr.resid(xmat, demeanlist(x, fl, weights = w, ...)), ax)) else
        return(setAttributes(qr.fitted(xmat, x), ax))
  } else if(na.rm && fill) {
    x[-cc, ] <- NA
    x[cc, ] <- demeanlist(x[cc, ], fl, weights = w, means = TRUE, ...)
    return(setAttributes(x, ax))
  } else return(setAttributes(demeanlist(x, fl, weights = w, means = TRUE, ...), ax))
}
fHDbetween.pdata.frame <- function(x, w = NULL, na.rm = TRUE, fill = TRUE, variable.wise = TRUE, ...) {
  ax <- attributes(x)
  if(na.rm && fill && variable.wise) {
    attributes(x) <- NULL
    varwisecomp <- function(x, fl, w, ...) lapply(x, function(y) {
      ycc <- which(!is.na(y))
      y[ycc] <- demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], means = TRUE, ...)
      return(y)
    })
    return(setAttributes(varwisecomp(x, ax[["index"]], w, ...), ax))
  } else if(na.rm && any(miss <- .Call(C_dt_na, x, seq_along(unclass(x))))) {
    ax[["na.rm"]] <- which(miss)
    cc <- which(!miss)
    Y <- demeanlist(.Call(C_subsetDT, x, cc, seq_along(unclass(x))),
                    lapply(ax[["index"]], function(y) y[cc]),
                    weights = w[cc], means = TRUE, ...)
    if(fill) return(setAttributes(.Call(Cpp_lassign, Y, fnrow2(x), cc, NA), ax)) else {
      ax[["row.names"]] <- ax[["row.names"]][cc]
      return(setAttributes(Y, ax))
    }
  } else return(setAttributes(demeanlist(x, ax[["index"]], weights = w, means = TRUE, ...), ax))
}
fHDbetween.data.frame <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
  ax <- attributes(x)

  if(na.rm) {
    cc <- if(variable.wise) complete.cases(fl) else complete.cases(x, fl) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- which(!cc)
      cc <- which(cc)
      if(!is.null(w)) w <- w[cc]
      if(!variable.wise) {
        if(fill) nrx <- fnrow2(x) else ax[["row.names"]] <- ax[["row.names"]][cc] # best ??
        x <- .Call(C_subsetDT, x, cc, seq_along(unclass(x)))
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) { # fl is a list !!
    fc <- vapply(unattrib(fl), is.factor, TRUE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) unused_arg_action(match.call(), ...)
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(unclass(fl))) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demeanlist(do.call(cbind, fl[!fc]), fl[fc], weights = w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(rep(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(fl[cc]) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    if(!missing(...)) unused_arg_action(match.call(), ...)
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(variable.wise) {
    if(na.rm) { # this means there were mising values in fl, which were already removed!
      return(setAttributes(lapply(unattrib(x), function(y) {
        y[-cc] <- NA # which is not faster !!
        ycc <- !is.na(y)
        YC <- which(ycc[cc])
        y[ycc] <- if(nallfc) y[ycc] - qr.resid(qr.default(xmat[YC, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], ...)) else if(fcl)
          demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], means = TRUE, ...) else  qr.fitted(qr.default(xmat[YC, , drop = FALSE]), y[ycc])
        return(y)
      }), ax)) # Rfast fastlm??
    } else {
      return(setAttributes(lapply(unattrib(x), function(y) {
        ycc <- which(!is.na(y))
        y[ycc] <- if(nallfc) y[ycc] - qr.resid(qr.default(xmat[ycc, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)) else if(fcl)
          demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], means = TRUE, ...) else qr.fitted(qr.default(xmat[ycc, , drop = FALSE]), y[ycc])
        return(y)
      }), ax)) # Rfast fastlm??
    }
  } else { # at this point missing values are already removed from x and fl !!
    if(nallfc || !fcl) {
      xmat <- qr.default(xmat)
      Y <- if(nallfc) lapply(x, function(y) y - qr.resid(xmat, demeanlist(y, fl, weights = w, ...))) else
           lapply(x, function(y) qr.fitted(xmat, y))
    } else Y <- demeanlist(x, fl, weights = w, means = TRUE, ...)
    if(na.rm && fill) { # x[cc, ] <- Y; x[-cc, ] <- NA
      return(setAttributes(.Call(Cpp_lassign, Y, nrx, cc, NA), ax))
    } else return(setAttributes(Y, ax))
  }
}


HDB <- function(x, ...) UseMethod("HDB") # , x

HDB.default <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...)
  fHDbetween.default(x, fl, w, na.rm, fill, ...)

HDB.pseries <- function(x, w = NULL, na.rm = TRUE, fill = TRUE, ...)
  fHDbetween.pseries(x, w, na.rm, fill, ...)

HDB.matrix <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, stub = "HDB.", ...)
  add_stub(fHDbetween.matrix(x, fl, w, na.rm, fill, ...), stub)

HDB.data.frame <- function(x, fl, w = NULL, cols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, stub = "HDB.", ...) {
  if(is.call(fl)) {
    ax <- attributes(x)
    nam <- ax[["names"]]
    if(length(fl) == 3L) {
      fvars <- ckmatch(all.vars(fl[[3L]]), nam)
      Xvars <- ckmatch(all.vars(fl[[2L]]), nam)
      fl[[2L]] <- NULL
    } else {
      fvars <- ckmatch(all.vars(fl), nam)
      Xvars <- if(!is.null(cols)) fsetdiff(cols2int(cols, x, nam), fvars) else seq_along(unclass(x))[-fvars]
    }
    ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[Xvars]) else nam[Xvars]

    if(na.rm) {
      miss <- if(variable.wise) .Call(C_dt_na, x, fvars) else .Call(C_dt_na, x, c(Xvars,fvars))
      if(any(miss)) {
        ax[["na.rm"]] <- which(miss)
        cc <- which(!miss)
        if(!is.null(w)) w <- w[cc]
        if(!variable.wise) if(fill) nrx <- fnrow2(x) else ax[["row.names"]] <- ax[["row.names"]][cc] # best ??
      } else na.rm <- FALSE
    }

    xmat <- NULL
    list2env(getfl(myModFrame(fl, if(na.rm) .Call(C_subsetDT, x, cc, fvars) else .subset(x, fvars))), envir = environment())
    fcl <- !is.null(fl)
    nallfc <- fcl && !is.null(xmat)
    if(nallfc) xmat <- demeanlist(xmat[, -1L, drop = FALSE], fl, weights = w, ...)

    if(variable.wise) {
      if(na.rm) {
        return(setAttributes(lapply(.subset(x, Xvars), function(y) {
          y[-cc] <- NA
          ycc <- !is.na(y)
          YC <- which(ycc[cc])
          y[ycc] <- if(nallfc) y[ycc] - qr.resid(qr.default(xmat[YC, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], ...)) else if(fcl)
            demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], means = TRUE, ...) else  qr.fitted(qr.default(xmat[YC, , drop = FALSE]), y[ycc])
          return(y)
        }), ax)) # Rfast fastlm??
      } else {
        return(setAttributes(lapply(.subset(x, Xvars), function(y) {
          ycc <- which(!is.na(y))
          y[ycc] <- if(nallfc) y[ycc] - qr.resid(qr.default(xmat[ycc, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)) else if(fcl)
            demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], means = TRUE, ...) else qr.fitted(qr.default(xmat[ycc, , drop = FALSE]), y[ycc])
          return(y)
        }), ax)) # Rfast fastlm??
      }
    } else { # at this point missing values are already removed from  fl !!
      # x <- if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else .subset(x, Xvars) # incorporated below -> more memory efficient
      if(nallfc || !fcl) {
        xmat <- qr.default(xmat)
        Y <- if(nallfc) lapply(if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else .subset(x, Xvars), function(y) y - qr.resid(xmat, demeanlist(y, fl, weights = w, ...))) else
          lapply(if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else .subset(x, Xvars), function(y) qr.fitted(xmat, y))
      } else Y <- demeanlist(if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else .subset(x, Xvars), fl, weights = w, means = TRUE, ...)
      if(na.rm && fill) { # x[cc, ] <- Y; x[-cc, ] <- NA
        return(setAttributes(.Call(Cpp_lassign, Y, nrx, cc, NA), ax))
      } else return(setAttributes(Y, ax))
    }
  } else # fl is not a formula !!
    return(add_stub(fHDbetween.data.frame(if(is.null(cols)) x else colsubset(x, cols), fl, w, na.rm, fill, variable.wise, ...), stub))
}
HDB.pdata.frame <- function(x, w = NULL, cols = is.numeric, na.rm = TRUE, fill = TRUE,
                            variable.wise = TRUE, stub = "HDB.", ...)
  return(add_stub(fHDbetween.pdata.frame(if(is.null(cols)) x else colsubset(x, cols), w, na.rm, fill, variable.wise, ...), stub))





# todo: good performance with missing values ??
# -> do pseries and pdata.frame method !!


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
