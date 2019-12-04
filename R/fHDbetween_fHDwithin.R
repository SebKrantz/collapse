# # HDB and HDW
# library(lfe)
# library(Rcpp)

myModFrame <- function(f, data) {
  t <- terms.formula(f)
  v <- attr(t, "variables")
  res <- eval(substitute(with(data, e), list(e = v)))
  attributes(res) <- list(names = as.character(v[-1]),
                          row.names = .set_row_names(nrow(data)),
                          class = "data.frame",
                          terms = t)
  return(res)
}
# mf <- myModFrame( ~ factor(cyl)*carb + factor(cyl):factor(vs) + vs + carb:am, data = mtcars)

getfl <- function(mf) {

  facts <- vapply(mf, is.factor, TRUE)

  if(any(facts)) {
    terms <- attributes(attr(mf, "terms"))
    clmf <- class(mf)
    class(mf) <- NULL # good ??
    tl <- terms[["term.labels"]]
    factors <- terms[[2L]]
    fctterms <- colSums(factors[facts, , drop = FALSE]) > 0
    fctinteract <- fctterms & colSums(factors) > 1

    if(any(fctinteract)) { # if any interactions involving factors
      singlefct <- match(tl[fctterms & !fctinteract], names(facts))
      intterms <- lapply(which(fctinteract), function(i) factors[,i] > 0) # names(which( # better way ??
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
        fctdat[(ctvec[2]+1):ctvec[3]] <- lapply(intterms[fctfct], function(x) if(length(x) == 2) do.call(`:`,mf[x]) else interaction(mf[x])) # or as.factor.GRP(GRP(mf[x]))
      if(tvec[4] != 0)
        fctdat[(ctvec[3]+1):ctvec[4]] <- lapply(intterms[!fctfct], function(x) {
          f <- x & facts
          nf <- x & !f
          f <- if(sum(f) == 1) mf[[which(f)]] else if(sum(f) == 2) do.call(`:`, mf[f]) else interaction(mf[f])
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

# getfl(model.frame( ~ cyl + carb, data = mtcars))
# getfl(model.frame( ~ factor(cyl)*carb, data = mtcars))
# mf =  model.frame( ~ factor(cyl)*carb + vs + wt:gear + wt:gear:carb, data = mtcars)
# getfl(model.frame( ~ factor(cyl)*carb + vs + wt:gear + wt:gear:carb, data = mtcars))

fHDwithin <- function(x, ...) {
  UseMethod("fHDwithin", x)
}
fHDwithin.default <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  ax <- attributes(x)
  if(na.rm) {
    cc <- complete.cases(x, fl) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- which(!cc)
      cc <- which(cc)
      if(!fill) {
        if(!is.null(ax[["names"]])) ax[["names"]] <- ax[["names"]][cc] # best ??
        x <- x[cc]
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) {
    fc <- vapply(fl, is.factor, TRUE, USE.NAMES = FALSE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) stop("Unknown argument ", dotstostr(...))
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demeanlist(do.call(cbind, fl[!fc]), fl[fc], weights = w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(rep(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(fl[cc]) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(nallfc || !fcl) {
    xmat <- qr.default(xmat)
    Y <- if(nallfc) qr.resid(xmat, demeanlist(x, fl, weights = w, ...)) else
      qr.resid(xmat, x)
  } else Y <- demeanlist(x, fl, weights = w, ...)
  if(na.rm && fill) {
    x[cc] <- Y
    x[-cc] <- NA
    return(setAttributes(x, ax))
  } else return(setAttributes(Y, ax))
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
        ax[["dimnames"]][[1L]] <- ax[["dimnames"]][[1L]][cc] # best ??
        ax[["dim"]][1L] <- length(cc)
        x <- x[cc, , drop = FALSE]
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) {
    fc <- vapply(fl, is.factor, TRUE, USE.NAMES = FALSE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) stop("Unknown argument ", dotstostr(...))
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demeanlist(do.call(cbind, fl[!fc]), fl[fc], weights = w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(rep(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(fl[cc]) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(nallfc || !fcl) {
    xmat <- qr.default(xmat)
    Y <- if(nallfc) qr.resid(xmat, demeanlist(x, fl, weights = w, ...)) else
         qr.resid(xmat, x)
  } else Y <- demeanlist(x, fl, weights = w, ...)
  if(na.rm && fill) {
    x[cc, ] <- Y
    x[-cc, ] <- NA
    return(setAttributes(x, ax))
  } else return(setAttributes(Y, ax))
}
fHDwithin.pdata.frame <- function(x, w = NULL, na.rm = TRUE, fill = TRUE, variable.wise = FALSE, ...) {
  ax <- attributes(x)
  if(variable.wise) {
    attributes(x) <- NULL
    varwisecomp <- function(x, fl, w, ...) lapply(x, function(y) {
      ycc <- which(!is.na(y))
      y[ycc] <- demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)
      return(y)
    })
    return(setAttributes(varwisecomp(x, ax[["index"]], w, ...), ax))
  } else if(na.rm && any(miss <- .Call(C_dt_na, x, seq_along(x)))) {
      ax[["na.rm"]] <- which(miss)
      cc <- which(!miss)
      Y <- demeanlist(.Call(C_subsetDT, x, cc, seq_along(x)),
                            lapply(ax[["index"]], function(y) y[cc]),
                            weights = w[cc], ...)
    if(fill) return(setAttributes(.Call(Cpp_lassign, Y, nrow(x), cc, NA), ax)) else {
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
        x <- .Call(C_subsetDT, x, cc, seq_along(x))
        if(fill) nrx <- nrow(x) else ax[["row.names"]] <- ax[["row.names"]][cc] # best ??
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) { # fl is a list !!
    fc <- vapply(fl, is.factor, TRUE, USE.NAMES = FALSE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) stop("Unknown argument ", dotstostr(...))
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demeanlist(do.call(cbind, fl[!fc]), fl[fc], weights = w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(rep(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
      if(!missing(...)) stop("Unknown argument ", dotstostr(...))
      xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
      nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(fl[cc]) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(variable.wise) {
    if(na.rm) { # this means there were mising values in fl, which were already removed!
      return(setAttributes(lapply(x, function(y) {
        y[-cc] <- NA # which is not faster !!
        ycc <- !is.na(y)
        YC <- which(ycc[cc])
        y[ycc] <- if(nallfc) qr.resid(qr.default(xmat[YC, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], ...)) else if(fcl)
          demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], ...) else  qr.resid(qr.default(xmat[YC, , drop = FALSE]), y[ycc])
        return(y)
      }), ax)) # Rfast fastlm??
    } else {
      return(setAttributes(lapply(x, function(y) {
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
HDW <- function(x, ...) { # fl, w = NULL, na.rm = TRUE, fill = FALSE,
  UseMethod("HDW", x)
}
HDW.default <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  fHDwithin.default(x, fl, w, na.rm, fill, ...)
}
HDW.pseries <- function(x, w = NULL, na.rm = TRUE, fill = TRUE, ...) {
  fHDwithin.pseries(x, w, na.rm, fill, ...)
}
HDW.matrix <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, stub = "HDW_", ...) {
  add_stub(fHDwithin.matrix(x, fl, w, na.rm, fill, ...), stub)
}
HDW.data.frame <- function(x, fl, w = NULL, cols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, stub = "HDW_", ...) {
  if(is.call(fl)) {
    ax <- attributes(x)
    nam <- names(x)
    if(length(fl) == 3L) {
      fvars <- anyNAerror(match(all.vars(fl[[3L]]), nam), "Unknown variables in formula!")
      Xvars <- anyNAerror(match(all.vars(fl[[2L]]), nam), "Unknown variables in formula!")
      fl[[2L]] <- NULL
    } else {
      fvars <- anyNAerror(match(all.vars(fl), nam), "Unknown variables in formula!")
      Xvars <- if(!is.null(cols)) setdiff(cols2int(cols, x, nam), fvars) else seq_along(x)[-fvars]
    }
    ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[Xvars]) else nam[Xvars]

    if(na.rm) {
      miss <- if(variable.wise) .Call(C_dt_na, x, fvars) else .Call(C_dt_na, x, c(Xvars,fvars))
      if(any(miss)) {
        ax[["na.rm"]] <- which(miss)
        cc <- which(!miss)
        if(!is.null(w)) w <- w[cc]
        if(!variable.wise) if(fill) nrx <- nrow(x) else ax[["row.names"]] <- ax[["row.names"]][cc] # best ??
      } else na.rm <- FALSE
    }

    xmat <- NULL
    list2env(getfl(myModFrame(fl, if(na.rm) .Call(C_subsetDT, x, cc, fvars) else colsubset(x, fvars))), envir = environment())
    fcl <- !is.null(fl)
    nallfc <- fcl && !is.null(xmat)
    if(nallfc) xmat <- demeanlist(xmat[, -1L, drop = FALSE], fl, weights = w, ...)

    if(variable.wise) {
      if(na.rm) {
        return(setAttributes(lapply(unclass(x)[Xvars], function(y) {
          y[-cc] <- NA
          ycc <- !is.na(y)
          YC <- which(ycc[cc])
          y[ycc] <- if(nallfc) qr.resid(qr.default(xmat[YC, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], ...)) else if(fcl)
            demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], ...) else  qr.resid(qr.default(xmat[YC, , drop = FALSE]), y[ycc])
          return(y)
        }), ax)) # Rfast fastlm??
      } else {
        return(setAttributes(lapply(unclass(x)[Xvars], function(y) {
          ycc <- which(!is.na(y))
          y[ycc] <- if(nallfc) qr.resid(qr.default(xmat[ycc, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)) else if(fcl)
            demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...) else qr.resid(qr.default(xmat[ycc, , drop = FALSE]), y[ycc])
          return(y)
        }), ax)) # Rfast fastlm??
      }
    } else { # at this point missing values are already removed from  fl !!
      # x <- if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else unclass(x)[Xvars] # incorporated below -> more memory efficient
      if(nallfc || !fcl) {
        xmat <- qr.default(xmat)
        Y <- if(nallfc) lapply(demeanlist(if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else unclass(x)[Xvars], fl, weights = w, ...), function(y) qr.resid(xmat, y)) else
          lapply(if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else unclass(x)[Xvars], function(y) qr.resid(xmat, y))
      } else Y <- demeanlist(if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else unclass(x)[Xvars], fl, weights = w, ...)
      if(na.rm && fill) { # x[cc, ] <- Y; x[-cc, ] <- NA
        return(setAttributes(.Call(Cpp_lassign, Y, nrx, cc, NA), ax))
      } else return(setAttributes(Y, ax))
    }
  } else # fl is not a formula !!
 return(add_stub(fHDwithin.data.frame(if(is.null(cols)) x else colsubset(x, cols), fl, w, na.rm, fill, variable.wise, ...), stub))
}

HDW.pdata.frame <- function(x, w = NULL, cols = is.numeric, na.rm = TRUE, fill = TRUE,
                            variable.wise = FALSE, stub = "HDW_", ...) {
  return(add_stub(fHDwithin.data.frame(if(is.null(cols)) x else colsubset(x, cols), w, na.rm, fill, variable.wise, ...), stub))
}

# Theory: y = ?1 x1 + ?2 x2 + e
# FWT: M2 y = ?1 M2 x1 + e so residuals: e = M2 y - ?1 M2 x1 and fitted:
# Now M = I - x(x'x)-1x' = I - P.
# So (I-P2) y = ?1 (I-P2) x1 + e or y - P2 y = ?1 x1 - ?1 P2 x1 + e

# I want y - e = y^ = ?1 x1 + ?2 x2
# so
# P2 y = ?1 P2 x1 + ?2 x2
# Haven't quite figgured it out, but my solution is to just subtract the demeaned data !!

# Note: Only changes to fHDwithin is in the computation part: Perhaps you can combine the code in some better way to reduce code duplication ??

fHDbetween <- function(x, ...) {
  UseMethod("fHDbetween", x)
}
fHDbetween.default <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  ax <- attributes(x)
  if(na.rm) {
    cc <- complete.cases(x, fl) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- which(!cc)
      cc <- which(cc)
      if(!fill) {
        if(!is.null(ax[["names"]])) ax[["names"]] <- ax[["names"]][cc] # best ??
        x <- x[cc]
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) {
    fc <- vapply(fl, is.factor, TRUE, USE.NAMES = FALSE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) stop("Unknown argument ", dotstostr(...))
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demeanlist(do.call(cbind, fl[!fc]), fl[fc], weights = w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(rep(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(fl[cc]) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(nallfc || !fcl) {
    xmat <- qr.default(xmat)
    Y <- if(nallfc) x - qr.resid(xmat, demeanlist(x, fl, weights = w, ...)) else
          qr.fitted(xmat, x)
  } else Y <- demeanlist(x, fl, weights = w, means = TRUE, ...)
  if(na.rm && fill) {
    x[cc] <- Y
    x[-cc] <- NA
    return(setAttributes(x, ax))
  } else return(setAttributes(Y, ax))
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
        ax[["dimnames"]][[1L]] <- ax[["dimnames"]][[1L]][cc] # best ??
        ax[["dim"]][1L] <- length(cc)
        x <- x[cc, , drop = FALSE]
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) {
    fc <- vapply(fl, is.factor, TRUE, USE.NAMES = FALSE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) stop("Unknown argument ", dotstostr(...))
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demeanlist(do.call(cbind, fl[!fc]), fl[fc], weights = w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(rep(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(fl[cc]) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(nallfc || !fcl) {
    xmat <- qr.default(xmat)
    Y <- if(nallfc) x - qr.resid(xmat, demeanlist(x, fl, weights = w, ...)) else
         qr.fitted(xmat, x)
  } else Y <- demeanlist(x, fl, weights = w, means = TRUE, ...)
  if(na.rm && fill) {
    x[cc, ] <- Y
    x[-cc, ] <- NA
    return(setAttributes(x, ax))
  } else return(setAttributes(Y, ax))
}
fHDbetween.pdata.frame <- function(x, w = NULL, na.rm = TRUE, fill = TRUE, variable.wise = FALSE, ...) {
  ax <- attributes(x)
  if(variable.wise) {
    attributes(x) <- NULL
    varwisecomp <- function(x, fl, w, ...) lapply(x, function(y) {
      ycc <- which(!is.na(y))
      y[ycc] <- demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], means = TRUE, ...)
      return(y)
    })
    return(setAttributes(varwisecomp(x, ax[["index"]], w, ...), ax))
  } else if(na.rm && any(miss <- .Call(C_dt_na, x, seq_along(x)))) {
    ax[["na.rm"]] <- which(miss)
    cc <- which(!miss)
    Y <- demeanlist(.Call(C_subsetDT, x, cc, seq_along(x)),
                    lapply(ax[["index"]], function(y) y[cc]),
                    weights = w[cc], means = TRUE, ...)
    if(fill) return(setAttributes(.Call(Cpp_lassign, Y, nrow(x), cc, NA), ax)) else {
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
        x <- .Call(C_subsetDT, x, cc, seq_along(x))
        if(fill) nrx <- nrow(x) else ax[["row.names"]] <- ax[["row.names"]][cc] # best ??
      }
    } else na.rm <- FALSE
  }

  if(is.list(fl)) { # fl is a list !!
    fc <- vapply(fl, is.factor, TRUE, USE.NAMES = FALSE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) stop("Unknown argument ", dotstostr(...))
    nallfc <- fcl && !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(nallfc) {
      xmat <- demeanlist(do.call(cbind, fl[!fc]), fl[fc], weights = w, ...)
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(rep(1L, length(fl[[1L]]))), fl))
  } else if(is.matrix(fl)) {
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc, , drop = FALSE]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  } else if(is.factor(fl)) {
    fl <- if(na.rm) list(fl[cc]) else list(fl)
    fcl <- TRUE
    nallfc <- FALSE
  } else {
    if(!is.numeric(fl)) stop("fl must be a list of vectors / factors, a numeric matrix or a numeric vector")
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- if(na.rm) cbind(Intercept = 1L, fl[cc]) else cbind(Intercept = 1L, fl)
    nallfc <- fcl <- FALSE
  }

  if(variable.wise) {
    if(na.rm) { # this means there were mising values in fl, which were already removed!
      return(setAttributes(lapply(x, function(y) {
        y[-cc] <- NA # which is not faster !!
        ycc <- !is.na(y)
        YC <- which(ycc[cc])
        y[ycc] <- if(nallfc) y[ycc] - qr.resid(qr.default(xmat[YC, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], ...)) else if(fcl)
          demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], means = TRUE, ...) else  qr.fitted(qr.default(xmat[YC, , drop = FALSE]), y[ycc])
        return(y)
      }), ax)) # Rfast fastlm??
    } else {
      return(setAttributes(lapply(x, function(y) {
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


HDB <- function(x, ...) { # fl, w = NULL, na.rm = TRUE, fill = FALSE,
  UseMethod("HDB", x)
}
HDB.default <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  fHDbetween.default(x, fl, w, na.rm, fill, ...)
}
HDB.pseries <- function(x, w = NULL, na.rm = TRUE, fill = TRUE, ...) {
  fHDbetween.pseries(x, w, na.rm, fill, ...)
}
HDB.matrix <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, stub = "HDB_", ...) {
  add_stub(fHDbetween.matrix(x, fl, w, na.rm, fill, ...), stub)
}
HDB.data.frame <- function(x, fl, w = NULL, cols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, stub = "HDB_", ...) {
  if(is.call(fl)) {
    ax <- attributes(x)
    nam <- names(x)
    if(length(fl) == 3L) {
      fvars <- anyNAerror(match(all.vars(fl[[3L]]), nam), "Unknown variables in formula!")
      Xvars <- anyNAerror(match(all.vars(fl[[2L]]), nam), "Unknown variables in formula!")
      fl[[2L]] <- NULL
    } else {
      fvars <- anyNAerror(match(all.vars(fl), nam), "Unknown variables in formula!")
      Xvars <- if(!is.null(cols)) setdiff(cols2int(cols, x, nam), fvars) else seq_along(x)[-fvars]
    }
    ax[["names"]] <- if(is.character(stub)) paste0(stub, nam[Xvars]) else nam[Xvars]

    if(na.rm) {
      miss <- if(variable.wise) .Call(C_dt_na, x, fvars) else .Call(C_dt_na, x, c(Xvars,fvars))
      if(any(miss)) {
        ax[["na.rm"]] <- which(miss)
        cc <- which(!miss)
        if(!is.null(w)) w <- w[cc]
        if(!variable.wise) if(fill) nrx <- nrow(x) else ax[["row.names"]] <- ax[["row.names"]][cc] # best ??
      } else na.rm <- FALSE
    }

    xmat <- NULL
    list2env(getfl(myModFrame(fl, if(na.rm) .Call(C_subsetDT, x, cc, fvars) else colsubset(x, fvars))), envir = environment())
    fcl <- !is.null(fl)
    nallfc <- fcl && !is.null(xmat)
    if(nallfc) xmat <- demeanlist(xmat[, -1L, drop = FALSE], fl, weights = w, ...)

    if(variable.wise) {
      if(na.rm) {
        return(setAttributes(lapply(unclass(x)[Xvars], function(y) {
          y[-cc] <- NA
          ycc <- !is.na(y)
          YC <- which(ycc[cc])
          y[ycc] <- if(nallfc) y[ycc] - qr.resid(qr.default(xmat[YC, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], ...)) else if(fcl)
            demeanlist(y[ycc], subsetfl(fl, YC), weights = w[YC], means = TRUE, ...) else  qr.fitted(qr.default(xmat[YC, , drop = FALSE]), y[ycc])
          return(y)
        }), ax)) # Rfast fastlm??
      } else {
        return(setAttributes(lapply(unclass(x)[Xvars], function(y) {
          ycc <- which(!is.na(y))
          y[ycc] <- if(nallfc) y[ycc] - qr.resid(qr.default(xmat[ycc, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)) else if(fcl)
            demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], means = TRUE, ...) else qr.fitted(qr.default(xmat[ycc, , drop = FALSE]), y[ycc])
          return(y)
        }), ax)) # Rfast fastlm??
      }
    } else { # at this point missing values are already removed from  fl !!
      # x <- if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else unclass(x)[Xvars] # incorporated below -> more memory efficient
      if(nallfc || !fcl) {
        xmat <- qr.default(xmat)
        Y <- if(nallfc) lapply(if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else unclass(x)[Xvars], function(y) y - qr.resid(xmat, demeanlist(y, fl, weights = w, ...))) else
          lapply(if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else unclass(x)[Xvars], function(y) qr.fitted(xmat, y))
      } else Y <- demeanlist(if(na.rm) .Call(C_subsetDT, x, cc, Xvars) else unclass(x)[Xvars], fl, weights = w, means = TRUE, ...)
      if(na.rm && fill) { # x[cc, ] <- Y; x[-cc, ] <- NA
        return(setAttributes(.Call(Cpp_lassign, Y, nrx, cc, NA), ax))
      } else return(setAttributes(Y, ax))
    }
  } else # fl is not a formula !!
    return(add_stub(fHDbetween.data.frame(if(is.null(cols)) x else colsubset(x, cols), fl, w, na.rm, fill, variable.wise, ...), stub))
}
HDB.pdata.frame <- function(x, w = NULL, cols = is.numeric, na.rm = TRUE, fill = TRUE,
                            variable.wise = FALSE, stub = "HDB_", ...) {
  return(add_stub(fHDbetween.data.frame(if(is.null(cols)) x else colsubset(x, cols), w, na.rm, fill, variable.wise, ...), stub))
}




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


# Previous Versions: (for certain tasks slightly faster !!), and also implement keep.xt and keep.w arguments + weight formula !!
# getPartData <- function(X, fl, cols, na.rm, variable.wise) {
#   res <- list(fl = NULL, dat = NULL, cc = NULL, Xvars = NULL, na.rm = na.rm)
#   nam <- names(X)
#   form <- all(class(fl) == "formula")
#   twosided <- form && length(fl) == 3L
#   colsl <- !is.null(cols)
#   ind <- if(!twosided && colsl) cols2int(cols, X, nam) else seq_along(X)
#   if(form) { # fl is formula
#     if(twosided) {
#       fvars <- anyNAerror(match(all.vars(fl[[3L]]), nam), "Unknown variables in formula!")
#       Xvars <- anyNAerror(match(all.vars(fl[[2L]]), nam), "Unknown variables in formula!")
#       fl[[2L]] <- NULL
#     } else {
#       fvars <- anyNAerror(match(all.vars(fl), nam), "Unknown variables in formula!")
#       Xvars <- setdiff(ind, fvars)
#     }
#     res[[4L]] <- Xvars
#     if(na.rm) {
#       miss <- if(variable.wise) .Call(C_dt_na, X, fvars) else .Call(C_dt_na, X, c(Xvars,fvars))
#       if(any(miss)) {
#         cc <- which(!miss)
#         res[1:2] <- getfl(myModFrame(fl, .Call(C_subsetDT, X, cc, fvars)))
#         res[[3L]] <- cc
#       } else {
#         res[1:2] <- getfl(myModFrame(fl, X[fvars]))
#         if(variable.wise) res[[3L]] <- "no missing cases in fl!" else {
#           res[[3L]] <- "no missing cases!"
#           res[[5L]] <- FALSE
#         }
#       }
#     } else res[1:2] <- getfl(myModFrame(fl, X))
#   } else if(is.list(fl)) { # fl is factor list
#     class(fl) <- NULL
#     res[[4L]] <- ind
#     fc <- vapply(fl, is.factor, TRUE)
#     fcl <- any(!fc)
#     if(na.rm) {
#       cc <- if(variable.wise) which(!.Call(C_dt_na, fl, seq_along(fl))) else
#         which(!.Call(C_dt_na, c(X,fl), seq_len(length(X)+length(fl))))
#       missc <- length(cc) != nrow(X)
#       if(fcl) {
#         if(missc) fl <- subsetfl(fl, cc)
#         res[[1L]] <- fl[fc]
#         res[[2L]] <- do.call(cbind, c(list(Intercept = rep(1, length(fl[[1L]]))), fl[!fc]))
#       } else res[[1L]] <- if(missc) subsetfl(fl, cc) else fl
#       if(missc) res[[3L]] <- cc else if(variable.wise) res[[3L]] <- "no missing cases in fl!" else {
#         res[[3L]] <- "no missing cases!"
#         res[[5L]] <- FALSE
#       }
#     } else {
#       if(fcl) {
#         res[[1L]] <- fl[fc]
#         res[[2L]] <- do.call(cbind, c(list(Intercept = rep(1, length(fl[[1L]]))), fl[!fc]))
#       } else res[[1L]] <- fl
#     }
#   } else { # fl is factor, vector or matrix !!
#     res[[4L]] <- ind
#     if(na.rm) {
#       cc <- if(variable.wise) which(complete.cases(fl)) else which(complete.cases(X, fl))
#       if(length(cc) != nrow(X)) {
#         if(is.factor(fl)) res[[1L]] <- list(.Call(C_subsetVector, fl, cc)) else res[[2L]] <- cbind(Intercept = 1L, .Call(C_subsetVector, fl, cc))
#         res[[3L]] <- cc
#       } else {
#         if(is.factor(fl)) res[[1L]] <- list(fl) else res[[2L]] <- cbind(Intercept = 1L, fl)
#         if(variable.wise) res[[3L]] <- "no missing cases in fl!" else {
#           res[[3L]] <- "no missing cases!"
#           res[[5L]] <- FALSE
#         }
#       }
#     } else if(is.factor(fl)) res[[1L]] <- list(fl) else res[[2L]] <- cbind(Intercept = 1L, fl)
#   }
#   return(res)
# }
# HDW.data.frame.old <- function(x, fl, w = NULL, cols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
#   # return(getPartData(x, fl, cols, na.rm, variable.wise))
#   list2env(getPartData(x, fl, cols, na.rm, variable.wise), envir = environment())
#
#   if(variable.wise) fill <- TRUE
#   # Best solution ????:
#   if(na.rm && !fill) x <- .Call(C_subsetDT, x, cc, Xvars) else if(na.rm && fill && !variable.wise)
#     x <- .Call(C_subsetDT, x, NULL, Xvars) else if(!identical(Xvars,seq_along(x))) x <- x[Xvars]
#     Xvars <- seq_along(x)
#
#     if(!is.null(dat)) {
#       if(na.rm && fill) {
#         if(variable.wise) {
#           ax <- attributes(x)
#           ax[["cases"]] <- cc
#           lfl <- length(fl) > 0
#           if(cc[1] == "no missing cases in fl!") {
#             if(lfl) dat <- demeanlist(dat[,-1, drop = FALSE], fl, weights = w, ...)
#             return(setAttributes(lapply(x, function(y) {
#               ycc <- which(!is.na(y))
#               if(lfl) y[ycc] <- qr.resid(qr.default(dat[ycc, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)) else
#                 y[ycc] <- qr.resid(qr.default(dat[ycc, , drop = FALSE]), y[ycc])
#               return(y)
#             }), ax)) # Rfast fastlm??
#           } else {
#             if(lfl) dat <- demeanlist(dat[,-1, drop = FALSE], fl, weights = w[cc], ...)
#             return(setAttributes(lapply(x, function(y) {
#               ycc <- !is.na(y)
#               ycc[-cc] <- FALSE
#               y[-cc] <- NA # which is not faster !!
#               if(lfl) {
#                 YC <- which(ycc[cc])
#                 y[ycc] <- qr.resid(qr.default(dat[YC, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, YC), weights = w[ycc], ...))
#               } else y[ycc] <- qr.resid(qr.default(dat[ycc[cc], , drop = FALSE]), y[ycc])
#               return(y)
#             }), ax)) # Rfast fastlm??
#           }
#         } else {
#           attr(x, "cases") <- cc
#           if(length(fl)) { # fastest so far !! (and memory efficient)
#             dat <- qr.default(demeanlist(dat[,-1], fl, weights = w[cc], ...)) # This is necessary!!!
#             Y <- lapply(demeanlist(.Call(C_subsetDT, x, cc, Xvars), fl, weights = w[cc], ...), function(y) qr.resid(dat, y))
#             return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
#                          seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
#           } else {
#             dat <- qr.default(dat)
#             Y <- lapply(.Call(C_subsetDT, x, cc, Xvars), function(y) qr.resid(dat, y))
#             return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
#                          seq_row(X)[-cc], Xvars, NULL, NA, FALSE))
#           }
#         }
#       } else {
#         ax <- attributes(x)
#         if(na.rm) ax[["cases"]] <- cc
#         if(length(fl)) {
#           if(!is.null(w) && na.rm) w <- w[cc]
#           dat <- qr.default(demeanlist(dat[,-1], fl, weights = w, ...)) # This is necessary!!!
#           return(setAttributes(lapply(demeanlist(x, fl, weights = w, ...),
#                                       function(y) qr.resid(dat, y)), ax))
#         } else {
#           dat <- qr.default(dat)
#           return(setAttributes(lapply(x, function(y) qr.resid(dat, y)), ax))
#         }
#       }
#     } else {
#       if(na.rm && fill) {
#         if(variable.wise) {
#           ax <- attributes(x)
#           ax[["cases"]] <- cc
#           if(cc[1] == "no missing cases in fl!") {
#             return(setAttributes(lapply(x, function(y) {
#               ycc <- which(!is.na(y))
#               y[ycc] <- demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)
#               return(y)
#             }), ax))
#           } else {
#             return(setAttributes(lapply(x, function(y) {
#               ycc <- !is.na(y)
#               ycc[-cc] <- FALSE
#               y[-cc] <- NA # which is not faster !!
#               y[xcc] <- demeanlist(y[ycc], subsetfl(fl, which(ycc[cc])), weights = w[ycc], ...)
#               return(y)
#             }), ax))
#           }
#         } else {
#           attr(x, "cases") <- cc
#           Y <- demeanlist(.Call(C_subsetDT, x, cc, Xvars), fl, weights = w[cc], ...)
#           # NUMlassignCpp(Y, cc, Xvars)
#           return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
#                        seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
#         }
#       } else {
#         if(na.rm) return(`attr<-`(demeanlist(x, fl, weights = w, ...), "cases", cc)) else
#           return(demeanlist(x, fl, weights = w, ...))
#       }
#     }
# }
# HDW.pdata.frame <- function(x, w = NULL, cols = is.numeric, na.rm = TRUE, fill = TRUE,
#                             variable.wise = FALSE, drop.xt = FALSE, drop.w = TRUE, ...) {
#   ax <- attributes(x)
#   nam <- ax[["names"]]
#   gn <- match(names(ax[["index"]]), nam)
#   gn2 = gn <- gn[!is.na(gn)]
#   if(!is.null(cols)) {
#     if(is.function(cols)) cols <- seq_along(x)[!vapply(x, cols, TRUE)] else if(is.character(cols))
#       cols <- seq_along(x)[-match(cols, nam)] else if(is.logical(cols))
#         cols <- seq_along(x)[!cols] else if(is.numeric(cols))
#           cols <- seq_along(x)[-cols] else stop("cols needs to be a function, column names, indices or a logical vector")
#         if(drop.xt) gn <- unique.default(c(gn, cols)) else if(length(gn)) gn2 <- unique.default(c(gn2, cols)) else gn <- cols
#   }
#   if(!is.null(w)) {
#     if(is.call(w)) w <- all.vars(w)
#     if(length(w) == 1) {
#       wn <- match(w, nam)
#       if(any(gn == wn)) stop("Weights coincide with grouping variables!")
#       w <- x[[wn]]
#       if(drop.w) if(!length(gn)) x[[wn]] <- NULL else if(drop.xt) gn <- c(gn,wn) else gn2 <- c(gn2,wn)
#     }
#   }
#   varwisecomp <- function(x, fl, w, ...) lapply(x, function(y) {
#     ycc <- which(!is.na(y))
#     y[ycc] <- demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)
#     return(y)
#   })
#   fillcomp <- function(x, fl, w, cc, Xvars, ...) {
#     if(length(cc) != length(x[[1L]])) {
#       Y <- demeanlist(.Call(C_subsetDT, x, cc, Xvars), lapply(fl, function(y) y[cc]), weights = w[cc], ...)
#       return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
#                    seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
#     } else return(demeanlist(x[Xvars], fl, weights = w, ...))
#   }
#   nonfillcomp <- function(x, fl, w, cc, Xvars, ...) {
#     if(length(cc) != length(x[[1L]]))
#       return(demeanlist(.Call(C_subsetDT, x, cc, Xvars), lapply(fl, function(y) y[cc]), weights = w[cc], ...)) else
#         return(demeanlist(x[Xvars], fl, weights = w, ...))
#   }
#   if(na.rm && fill) {
#     if(variable.wise) {
#       attributes(x) <- NULL
#       if(length(gn)) {
#         if(drop.xt) {
#           ax[["names"]] <- if(give.names) paste0("HDW.",nam[-gn]) else nam[-gn]
#           return(setAttributes(varwisecomp(x[-gn], ax[["index"]], w, ...), ax))
#         } else {
#           ax[["names"]] <- c(nam[gn], if(give.names) paste0("HDW.",nam[-gn2]) else nam[-gn2])
#           return(setAttributes(c(x[gn],varwisecomp(x[-gn2], ax[["index"]], w, ...)), ax))
#         }
#       } else {
#         if(give.names) paste0("HDW.",nam)
#         return(setAttributes(varwisecomp(x, ax[["index"]], w, ...), ax))
#       }
#     } else {
#       cc <- which(!.Call(C_dt_na, x, seq_along(x)[-c(gn,gn2)])) # good??
#       ax[["cases"]] <- if(length(cc) == length(x[[1L]])) "no missing cases!" else cc
#       if(length(gn)) {
#         if(drop.xt) {
#           ax[["names"]] <- if(give.names) paste0("HDW.",nam[-gn]) else nam[-gn]
#           return(setAttributes(fillcomp(x, ax[["index"]], w, cc, seq_along(x)[-gn], ...), ax))
#         } else {
#           ax[["names"]] <- c(nam[gn], if(give.names) paste0("HDW.",nam[-gn2]) else nam[-gn2])
#           return(setAttributes(c(x[gn], fillcomp(x, ax[["index"]], w, cc, seq_along(x)[-gn2], ...)), ax))
#         }
#       } else {
#         if(give.names) paste0("HDW.",nam)
#         return(setAttributes(fillcomp(x, ax[["index"]], w, cc, seq_along(x), ...), ax))
#       }
#     }
#   } else {
#     if(na.rm) {
#       cc <- which(!.Call(C_dt_na, x, seq_along(x)[-c(gn,gn2)])) # good??
#       if(length(cc) == length(x[[1L]])) ax[["cases"]] <- "no missing cases!" else {
#         ax[["row.names"]] <- ax[["row.names"]][cc]
#         ax[["cases"]] <- cc
#       }
#       if(length(gn)) {
#         if(drop.xt) {
#           ax[["names"]] <- if(give.names) paste0("HDW.",nam[-gn]) else nam[-gn]
#           return(setAttributes(nonfillcomp(x, ax[["index"]], w, cc, seq_along(x)[-gn], ...), ax))
#         } else {
#           ax[["names"]] <- c(nam[gn], if(give.names) paste0("HDW.",nam[-gn2]) else nam[-gn2])
#           return(setAttributes(c(.Call(C_subsetDT, x, cc, gn), nonfillcomp(x, ax[["index"]], w, cc, seq_along(x)[-gn2], ...)), ax))
#         }
#       } else {
#         if(give.names) paste0("HDW.",nam)
#         return(setAttributes(nonfillcomp(x, ax[["index"]], w, cc, seq_along(x), ...), ax))
#       }
#     } else {
#       attributes(x) <- NULL
#       if(length(gn)) {
#         if(drop.xt) {
#           ax[["names"]] <- if(give.names) paste0("HDW.",nam[-gn]) else nam[-gn]
#           return(setAttributes(demeanlist(x[-gn], ax[["index"]], weights = w, ...), ax))
#         } else {
#           ax[["names"]] <- c(nam[gn], if(give.names) paste0("HDW.",nam[-gn2]) else nam[-gn2])
#           return(setAttributes(c(x[gn],demeanlist(x[-gn2], ax[["index"]], weights = w, ...)), ax))
#         }
#       } else {
#         if(give.names) paste0("HDW.",nam)
#         return(setAttributes(demeanlist(x, ax[["index"]], weights = w, ...), ax))
#       }
#     }
#   }
# }

# HDB.default.old <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
#   ax <- attributes(x)
#   if(na.rm) {
#     cc <- which(complete.cases(x, fl)) # gives error if lengths don't match, otherwise demeanlist and qr.fitted give errors !!
#     if(length(cc) != length(x)) {
#       if(!fill) ax[["names"]] <- NULL
#       ax[["cases"]] <- cc
#     } else {
#       ax[["cases"]] <- "no missing cases!"
#       na.rm <- FALSE
#     }
#   }
#   if(is.list(fl)) {
#     fc <- vapply(fl, is.factor, TRUE)
#     fcl <- any(fc)
#     nallfc <- !all(fc)
#     if(na.rm) fl <- if(is.data.frame(fl)) .Call(C_subsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
#     if(fcl && nallfc) {
#       fl <- fl[fc]
#       xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl[!fc]))
#     } else if(!fcl) xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl))
#     if(na.rm) {
#       Y <- if(fcl && nallfc) x[cc] - qr.resid(qr.default(demeanlist(xmat, fl, weights = w[cc], ...)),
#                                               demeanlist(x[cc], fl, weights = w[cc], ...)) else if(fcl)
#                                                 demeanlist(x[cc], fl, weights = w[cc], means = TRUE, ...) else qr.fitted(qr.default(xmat), x[cc])
#       if(fill) {
#         x[cc] <- Y
#         x[-cc] <- NA
#         return(setAttributes(x, ax))
#       } else return(setAttributes(Y, ax))
#     } else {
#       if(fcl && nallfc)
#         return(setAttributes(x - qr.resid(qr.default(demeanlist(xmat, fl, weights = w, ...)),
#                                           demeanlist(x, fl, weights = w, ...)), ax)) else if(fcl)
#                                             return(setAttributes(demeanlist(x, fl, weights = w, means = TRUE, ...), ax)) else
#                                               return(setAttributes(qr.fitted(qr.default(xmat), x), ax))
#     }
#   } else if (is.matrix(fl)) {
#     xmat <- cbind(Intercept = 1L, fl)
#     if(na.rm) {
#       if(fill) {
#         x[cc] <- qr.fitted(qr.default(xmat[cc, , drop = FALSE]), x[cc])
#         x[-cc] <- NA
#         return(setAttributes(x, ax))
#       } else return(setAttributes(qr.fitted(qr.default(xmat[cc, , drop = FALSE]), x[cc]), ax)) # best ??
#     } else return(setAttributes(qr.fitted(qr.default(xmat), x), ax))
#   } else {
#     if(na.rm) {
#       Y <- if(is.factor(fl)) demeanlist(x[cc], list(fl[cc]), weights = w[cc], means = TRUE, ...) else
#         qr.fitted(qr.default(cbind(1L,fl[cc])), x[cc])
#       if(fill) {
#         x[cc] <- Y
#         x[-cc] <- NA
#         return(setAttributes(x, ax))
#       } else return(setAttributes(Y, ax))
#     } else {
#       if(is.factor(fl)) return(setAttributes(demeanlist(x, list(fl), weights = w, means = TRUE, ...), ax)) else
#         return(setAttributes(qr.fitted(qr.default(cbind(1L,fl)), x), ax))
#     }
#   }
# }
# HDB.matrix.old <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
#   ax <- attributes(x)
#   if(na.rm) {
#     cc <- which(complete.cases(x, fl)) # gives error if lengths don't match, otherwise demeanlist and qr.fitted give errors !!
#     if(length(cc) != nrow(x)) {
#       if(!fill) {
#         ax[["dimnames"]][[1L]] <- ax[["dimnames"]][[1L]][cc] # best ??
#         ax[["dim"]][1] <- length(cc)
#       }
#       ax[["cases"]] <- cc
#     } else {
#       na.rm <- FALSE
#       ax[["cases"]] <- "no missing cases!"
#     }
#   }
#   if(is.list(fl)) {
#     fc <- vapply(fl, is.factor, TRUE)
#     fcl <- any(fc)
#     nallfc <- !all(fc)
#     if(na.rm) fl <- if(is.data.frame(fl)) .Call(C_subsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
#     if(fcl && nallfc) {
#       fl <- fl[fc]
#       xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl[!fc]))
#     } else if(!fcl) xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl))
#     if(na.rm) {
#       Y <- if(fcl && nallfc) x[cc, , drop = FALSE] - qr.resid(qr.default(demeanlist(xmat, fl, weights = w[cc], ...)),
#                                                               demeanlist(x[cc, , drop = FALSE], fl, weights = w[cc], ...)) else if(fcl)
#                                                                 demeanlist(x[cc, , drop = FALSE], fl, weights = w[cc], means = TRUE, ...) else qr.fitted(qr.default(xmat), x[cc, , drop = FALSE])
#       if(fill) {
#         x[cc, ] <- Y
#         x[-cc, ] <- NA
#         return(setAttributes(x, ax))
#       } else return(setAttributes(Y, ax))
#     } else {
#       if(fcl && nallfc)
#         return(setAttributes(x - qr.resid(qr.default(demeanlist(xmat, fl, weights = w, ...)),
#                                           demeanlist(x, fl, weights = w, ...)), ax)) else if(fcl)
#                                             return(setAttributes(demeanlist(x, fl, weights = w, means = TRUE, ...), ax)) else
#                                               return(setAttributes(qr.fitted(qr.default(xmat), x), ax))
#     }
#   } else if (is.matrix(fl)) {
#     xmat <- cbind(Intercept = 1L, fl)
#     if(na.rm) {
#       if(fill) {
#         x[cc, ] <- qr.fitted(qr.default(xmat[cc, , drop = FALSE]), x[cc, , drop = FALSE])
#         x[-cc, ] <- NA
#         return(setAttributes(x, ax))
#       } else return(setAttributes(qr.fitted(qr.default(xmat[cc, , drop = FALSE]), x[cc, , drop = FALSE]), ax)) # best ??
#     } else return(setAttributes(qr.fitted(qr.default(xmat), x), ax))
#   } else {
#     if(na.rm) {
#       Y <- if(is.factor(fl)) demeanlist(x[cc, , drop = FALSE], list(fl[cc]), weights = w[cc], means = TRUE, ...) else
#         qr.fitted(qr.default(cbind(1L,fl[cc])), x[cc, , drop = FALSE])
#       if(fill) {
#         x[cc, ] <- Y
#         x[-cc, ] <- NA
#         return(setAttributes(x, ax))
#       } else return(setAttributes(Y, ax))
#     } else {
#       if(is.factor(fl)) return(setAttributes(demeanlist(x, list(fl), weights = w, means = TRUE, ...), ax)) else
#         return(setAttributes(qr.fitted(qr.default(cbind(1L,fl)), x), ax))
#     }
#   }
# }
# HDB.data.frame.old <- function(x, fl, w = NULL, cols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
#   list2env(getPartData(x, fl, cols, na.rm, variable.wise), envir = environment())
#
#   if(variable.wise) fill <- TRUE
#   # Best solution ????:
#   if(na.rm && !fill) x <- .Call(C_subsetDT, x, cc, Xvars) else if(na.rm && fill && !variable.wise)
#     x <- .Call(C_subsetDT, x, NULL, Xvars) else if(!identical(Xvars,seq_along(x))) x <- x[Xvars]
#     Xvars <- seq_along(x)
#
#     if(!is.null(dat)) {
#       if(na.rm && fill) {
#         if(variable.wise) {
#           ax <- attributes(x)
#           ax[["cases"]] <- cc
#           lfl <- length(fl) > 0
#           if(cc[1] == "no missing cases in fl!") {
#             if(lfl) dat <- demeanlist(dat[,-1, drop = FALSE], fl, weights = w, ...)
#             return(setAttributes(lapply(x, function(y) {
#               ycc <- which(!is.na(y))
#               if(lfl) y[ycc] <- y[ycc] - qr.resid(qr.default(dat[ycc, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)) else
#                 y[ycc] <- qr.fitted(qr.default(dat[ycc, , drop = FALSE]), y[ycc])
#               return(y)
#             }), ax)) # Rfast fastlm??
#           } else {
#             if(lfl) dat <- demeanlist(dat[,-1, drop = FALSE], fl, weights = w[cc], ...)
#             return(setAttributes(lapply(x, function(y) { # change !!
#               ycc <- !is.na(y)
#               ycc[-cc] <- FALSE
#               y[-cc] <- NA # which is not faster !!
#               if(lfl) {
#                 YC <- which(ycc[cc])
#                 y[ycc] <- y[ycc] - qr.resid(qr.default(dat[YC, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, YC), weights = w[ycc], ...))
#               } else y[ycc] <- qr.fitted(qr.default(dat[ycc[cc], , drop = FALSE]), y[ycc])
#               return(y)
#             }), ax)) # Rfast fastlm??
#           }
#         } else {
#           attr(x, "cases") <- cc
#           if(length(fl)) { # fastest so far !! (and memory efficient)
#             dat <- qr.default(demeanlist(dat[,-1], fl, weights = w[cc], ...)) # This is necessary!!!
#             Y <- lapply(.Call(C_subsetDT, x, cc, Xvars), function(y) y - qr.resid(dat, demeanlist(y, fl, weights = w[cc], ...))) # good ??
#             return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
#                          seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
#           } else {
#             dat <- qr.default(dat)
#             Y <- lapply(.Call(C_subsetDT, x, cc, Xvars), function(y) qr.fitted(dat, y))
#             return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
#                          seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
#           }
#         }
#       } else {
#         ax <- attributes(x)
#         if(na.rm) ax[["cases"]] <- cc
#         if(length(fl)) {
#           if(!is.null(w) && na.rm) w <- w[cc]
#           dat <- qr.default(demeanlist(dat[,-1], fl, weights = w, ...)) # This is necessary!!!
#           return(setAttributes(Map(function(x, y) x - qr.resid(dat, y), x, demeanlist(x, fl, weights = w, ...)), ax))
#         } else {
#           dat <- qr.default(dat)
#           return(setAttributes(lapply(x, function(y) qr.fitted(dat, y)), ax))
#         }
#       }
#     } else {
#       if(na.rm && fill) {
#         if(variable.wise) {
#           ax <- attributes(x)
#           ax[["cases"]] <- cc
#           if(cc[1] == "no missing cases in fl!") {
#             return(setAttributes(lapply(x, function(y) {
#               ycc <- which(!is.na(y))
#               y[ycc] <- demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], means = TRUE, ...)
#               return(y)
#             }), ax))
#           } else {
#             return(setAttributes(lapply(x, function(y) {
#               ycc <- !is.na(y)
#               ycc[-cc] <- FALSE
#               y[-cc] <- NA # which is not faster !!
#               y[ycc] <- demeanlist(y[ycc], subsetfl(fl, which(ycc[cc])), weights = w[ycc], means = TRUE, ...)
#               return(y)
#             }), ax))
#           }
#         } else {
#           attr(x, "cases") <- cc
#           Y <- demeanlist(.Call(C_subsetDT, x, cc, Xvars), fl, weights = w[cc], means = TRUE, ...)
#           return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
#                        seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
#         }
#       } else {
#         if(na.rm) return(`attr<-`(demeanlist(x, fl, weights = w, means = TRUE, ...), "cases", cc)) else
#           return(demeanlist(x, fl, weights = w, means = TRUE, ...))
#       }
#     }
# }
