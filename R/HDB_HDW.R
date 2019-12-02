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

# # C functions used:
# CSV <- data.table:::CsubsetVector
# CSDT <- data.table:::CsubsetDT
# NADT <- data.table:::Cdt_na
# AS <- data.table:::Cassign
# MM <- stats:::C_modelmatrix

getfl <- function(mf) {

  facts <- vapply(mf, is.factor, TRUE)

  if(any(facts)) {
    terms <- attributes(attr(mf, "terms"))
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
      modelterms <- terms.formula(as.formula(paste0("~ ",paste(modelterms, collapse = " + "))))
      moddat <- .External2(stats:::C_modelmatrix, modelterms, # removes NA's?? -> Nope !!
                           eval(substitute(with(mf, e),
                           list(e = attr(modelterms, "variables"))))) # model.matrix.default(as.formula(paste0("~ ",paste(modelterms, collapse = " + "))), data = mf) else moddat <- NULL
    } else moddat <- NULL

  } else {
    fctdat <- NULL
    moddat <- .External2(stats:::C_modelmatrix, attr(mf, "terms"), mf) # model.matrix(attr(mf, "terms"), data = mf)
  }
  return(list(fl = fctdat, data = moddat))
}


# getfl(model.frame( ~ factor(cyl)*carb, data = mtcars))
# mf =  model.frame( ~ factor(cyl)*carb + vs + wt:gear + wt:gear:carb, data = mtcars)
# getfl(model.frame( ~ factor(cyl)*carb + vs + wt:gear + wt:gear:carb, data = mtcars))

subsetfl <- function(fl, cc) {
  lapply(fl, function(f) { # use CsubsetDT or CsubsetVector ?? also check NA in regressors ??
    x <- attr(f, "x")
    if(is.null(x)) return(.Call(C_subsetVector, f, cc)) else
      return(`attr<-`(.Call(C_subsetVector, f, cc), "x",
                      if(is.matrix(x)) x[cc, , drop = FALSE] else
                        .Call(C_subsetVector, x, cc)))
  })
}

getPartData <- function(X, fl, cols, na.rm, variable.wise) {
  res <- list(fl = NULL, dat = NULL, cc = NULL, Xvars = NULL, na.rm = na.rm)
  nam <- names(X)
  form <- all(class(fl) == "formula")
  twosided <- form && length(fl) == 3
  colsl <- !is.null(cols)
  if(!twosided && colsl) {
    if(is.function(cols))
    ind <- which(vapply(X, cols, TRUE)) else if(is.numeric(cols))
    ind <- cols else if(is.character(cols))
    ind <- match(cols, nam) else if(is.logical(cols))
    ind <- which(cols) else
    stop("cols needs to be a function, column names, indices or a logical vector")
  } else ind <- seq_along(X)

    if(form) { # fl is formula
      if(twosided) {
        fvars <- match(all.vars(fl[[3L]]), nam)
        Xvars <- match(all.vars(fl[[2L]]), nam)
        fl[[2L]] <- NULL
      } else {
        fvars <- match(all.vars(fl), nam)
        Xvars <- setdiff(ind, fvars)
      }
      res[[4L]] <- Xvars
      if(na.rm) {
        cc <- if(variable.wise) which(!.Call(C_dt_na, X, fvars)) else which(!.Call(C_dt_na, X, c(Xvars,fvars)))
        if(length(cc) != nrow(X)) {
          res[1:2] <- getfl(myModFrame(fl, .Call(C_subsetDT, X, cc, fvars)))
          res[[3L]] <- cc
        } else {
          res[1:2] <- getfl(myModFrame(fl, X[fvars]))
          if(variable.wise) res[[3L]] <- "no missing cases in fl!" else {
            res[[3L]] <- "no missing cases!"
            res[[5L]] <- FALSE
          }
        }
      } else res[1:2] <- getfl(myModFrame(fl, X))
  } else if(is.list(fl)) { # fl is factor list
    class(fl) <- NULL
    res[[4L]] <- ind
    fc <- vapply(fl, is.factor, TRUE)
    fcl <- any(!fc)
    if(na.rm) {
      cc <- if(variable.wise) which(!.Call(C_dt_na, fl, seq_along(fl))) else
                              which(!.Call(C_dt_na, c(X,fl), seq_len(length(X)+length(fl))))
      missc <- length(cc) != nrow(X)
        if(fcl) {
          if(missc) fl <- subsetfl(fl, cc)
          res[[1L]] <- fl[fc]
          res[[2L]] <- do.call(cbind, c(list(Intercept = rep(1, length(fl[[1L]]))), fl[!fc]))
        } else res[[1L]] <- if(missc) subsetfl(fl, cc) else fl
      if(missc) res[[3L]] <- cc else if(variable.wise) res[[3L]] <- "no missing cases in fl!" else {
          res[[3L]] <- "no missing cases!"
          res[[5L]] <- FALSE
      }
    } else {
      if(fcl) {
        res[[1L]] <- fl[fc]
        res[[2L]] <- do.call(cbind, c(list(Intercept = rep(1, length(fl[[1L]]))), fl[!fc]))
      } else res[[1L]] <- fl
    }
  } else { # fl is factor, vector or matrix !!
    res[[4L]] <- ind
    if(na.rm) {
      cc <- if(variable.wise) which(complete.cases(fl)) else which(complete.cases(X, fl))
      if(length(cc) != nrow(X)) {
        if(is.factor(fl)) res[[1L]] <- list(.Call(C_subsetVector, fl, cc)) else res[[2L]] <- cbind(Intercept = 1L, .Call(C_subsetVector, fl, cc))
        res[[3L]] <- cc
      } else {
        if(is.factor(fl)) res[[1L]] <- list(fl) else res[[2L]] <- cbind(Intercept = 1L, fl)
        if(variable.wise) res[[3L]] <- "no missing cases in fl!" else {
          res[[3L]] <- "no missing cases!"
          res[[5L]] <- FALSE
        }
      }
    } else if(is.factor(fl)) res[[1L]] <- list(fl) else res[[2L]] <- cbind(Intercept = 1L, fl)
  }
  return(res)
}


fHDwithin <- function(x, ...) {
  UseMethod("fHDwithin", x)
}
fHDwithin.default <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  ax <- attributes(x)
  if(na.rm) {
    cc <- complete.cases(x, fl) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      if(!fill) ax[["names"]] <- NULL
      ax[["na.rm"]] <- which(!cc)
      cc <- which(cc)
    } else na.rm <- FALSE
  }
  if(is.list(fl)) {
    fc <- vapply(fl, is.factor, TRUE, USE.NAMES = FALSE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) stop("Unknown argument ", dotstostr(...))
    nallfc <- !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(fcl && nallfc) {
      xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl[!fc]))
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl))
    if(na.rm) {
      Y <- if(fcl && nallfc) qr.resid(qr(demeanlist(xmat, fl, weights = w[cc], ...)),
                                      demeanlist(x[cc], fl, weights = w[cc], ...)) else if(fcl)
                                        demeanlist(x[cc], fl, weights = w[cc], ...) else qr.resid(qr(xmat), x[cc])
      if(fill) {
        x[cc] <- Y
        x[-cc] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(Y, ax))
    } else {
      if(fcl && nallfc)
        return(setAttributes(qr.resid(qr(demeanlist(xmat, fl, weights = w, ...)),
                                      demeanlist(x, fl, weights = w, ...)), ax)) else if(fcl)
                                        return(setAttributes(demeanlist(x, fl, weights = w, ...), ax)) else
                                          return(setAttributes(qr.resid(qr(xmat), x), ax))
    }
  } else if (is.matrix(fl)) {
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- cbind(Intercept = 1L, fl)
    if(na.rm) {
      if(fill) {
        x[cc] <- qr.resid(qr(xmat[cc, , drop = FALSE]), x[cc])
        x[-cc] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(qr.resid(qr(xmat[cc, , drop = FALSE]), x[cc]), ax)) # best ??
    } else return(setAttributes(qr.resid(qr(xmat), x), ax))
  } else {
    ffll <- is.factor(fl)
    if(!ffll && !missing(...)) stop("Unknown argument ", dotstostr(...))
    if(na.rm) {
      Y <- if(ffll) demeanlist(x[cc], list(fl[cc]), weights = w[cc], ...) else
        qr.resid(qr(cbind(1L,fl[cc])), x[cc])
      if(fill) {
        x[cc] <- Y
        x[-cc] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(Y, ax))
    } else {
      if(ffll) return(setAttributes(demeanlist(x, list(fl), weights = w, ...), ax)) else
        return(setAttributes(qr.resid(qr(cbind(1L,fl)), x), ax))
    }
  }
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
      }
    } else na.rm <- FALSE
  }
  if(is.list(fl)) {
    fc <- vapply(fl, is.factor, TRUE, USE.NAMES = FALSE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) stop("Unknown argument ", dotstostr(...))
    nallfc <- !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(fcl && nallfc) {
      xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl[!fc]))
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl))
    if(na.rm) {
      Y <- if(fcl && nallfc) qr.resid(qr(demeanlist(xmat, fl, weights = w[cc], ...)),
                                      demeanlist(x[cc, , drop = FALSE], fl, weights = w[cc], ...)) else if(fcl)
                                        demeanlist(x[cc, , drop = FALSE], fl, weights = w[cc], ...) else qr.resid(qr(xmat), x[cc, , drop = FALSE])
      if(fill) {
        x[cc, ] <- Y
        x[-cc, ] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(Y, ax))
    } else {
      if(fcl && nallfc)
        return(setAttributes(qr.resid(qr(demeanlist(xmat, fl, weights = w, ...)),
                                      demeanlist(x, fl, weights = w, ...)), ax)) else if(fcl)
                                        return(setAttributes(demeanlist(x, fl, weights = w, ...), ax)) else
                                          return(setAttributes(qr.resid(qr(xmat), x), ax))
    }
  } else if (is.matrix(fl)) {
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- cbind(Intercept = 1L, fl)
    if(na.rm) {
      if(fill) {
        x[cc, ] <- qr.resid(qr(xmat[cc, , drop = FALSE]), x[cc, , drop = FALSE])
        x[-cc, ] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(qr.resid(qr(xmat[cc, , drop = FALSE]), x[cc, , drop = FALSE]), ax)) # best ??
    } else return(setAttributes(qr.resid(qr(xmat), x), ax))
  } else {
    ffll <- is.factor(fl)
    if(!ffll && !missing(...)) stop("Unknown argument ", dotstostr(...))
    if(na.rm) {
      Y <- if(ffll) demeanlist(x[cc, , drop = FALSE], list(fl[cc]), weights = w[cc], ...) else
        qr.resid(qr(cbind(1L,fl[cc])), x[cc, , drop = FALSE])
      if(fill) {
        x[cc, ] <- Y
        x[-cc, ] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(Y, ax))
    } else {
      if(ffll) return(setAttributes(demeanlist(x, list(fl), weights = w, ...), ax)) else
        return(setAttributes(qr.resid(qr(cbind(1L,fl)), x), ax))
    }
  }
}

fHDwithin.pdata.frame <- function(x, w = NULL, na.rm = TRUE, fill = TRUE, variable.wise = FALSE, ...) {
  ax <- attributes(x)
  varwisecomp <- function(x, fl, w, ...) lapply(x, function(y) {
    ycc <- which(!is.na(y))
    y[ycc] <- demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)
    return(y)
  })
  if(na.rm && fill) {
    if(variable.wise) {
      attributes(x) <- NULL
      return(setAttributes(varwisecomp(x, ax[["index"]], w, ...), ax))
    } else {
      miss <- .Call(C_dt_na, x, seq_along(x)) # good??
      if(any(miss)) {
        ax[["na.rm"]] <- which(miss)
        cc <- which(!miss)
        x[cc, ] <- demeanlist(.Call(C_subsetDT, x, cc, seq_along(x)),
                              lapply(ax[["index"]], function(y) y[cc]),
                              weights = w[cc], ...)
              # return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
              #             seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
        x[miss, ] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(demeanlist(x, ax[["index"]], weights = w, ...), ax))
    }
  } else {
    if(na.rm) {
      miss <- .Call(C_dt_na, x, seq_along(x)) # good??
      if(any(miss)) {
        ax[["na.rm"]] <- which(miss)
        cc <- which(!miss)
        ax[["row.names"]] <- ax[["row.names"]][cc]
        return(setAttributes(demeanlist(.Call(C_subsetDT, x, cc, seq_along(x)),
                                        lapply(ax[["index"]], function(y) y[cc]),
                                        weights = w[cc], ...), ax))
      } else return(setAttributes(demeanlist(x, ax[["index"]], weights = w, ...), ax))
    } else {
      attributes(x) <- NULL # good ??
      return(setAttributes(demeanlist(x, ax[["index"]], weights = w, ...), ax))
    }
  }
}

fHDwithin.df <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
  ax <- attributes(x)
  if(variable.wise) fill <- TRUE

  if(na.rm) {
    cc <- if(variable.wise) complete.cases(fl) else complete.cases(x, fl) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(!all(cc)) {
      ax[["na.rm"]] <- which(!cc)
      cc <- which(cc)
      if(!fill) ax[["row.names"]] <- ax[["row.names"]][cc] # best ??
      x <- .Call(C_subsetDT, x, cc, seq_along(x))
    } else na.rm <- FALSE
  }
  if(is.list(fl)) {
    fc <- vapply(fl, is.factor, TRUE, USE.NAMES = FALSE)
    fcl <- any(fc)
    if(!fcl && !missing(...)) stop("Unknown argument ", dotstostr(...))
    nallfc <- !all(fc)
    if(na.rm) fl <- if(inherits(fl, "data.frame")) .Call(C_subsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    attributes(fl) <- NULL # good here ??
    if(fcl && nallfc) {
      xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl[!fc]))
      fl <- fl[fc]
    } else if(!fcl) xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl))
    if(na.rm) {
      Y <- if(fcl && nallfc) qr.resid(qr(demeanlist(xmat, fl, weights = w[cc], ...)),
                                      demeanlist(x[cc, , drop = FALSE], fl, weights = w[cc], ...)) else if(fcl)
                                        demeanlist(x[cc, , drop = FALSE], fl, weights = w[cc], ...) else qr.resid(qr(xmat), x[cc, , drop = FALSE])
      if(fill) {
        x[cc, ] <- Y
        x[-cc, ] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(Y, ax))
    } else {
      if(fcl && nallfc)
        return(setAttributes(qr.resid(qr(demeanlist(xmat, fl, weights = w, ...)),
                                      demeanlist(x, fl, weights = w, ...)), ax)) else if(fcl)
                                        return(setAttributes(demeanlist(x, fl, weights = w, ...), ax)) else
                                          return(setAttributes(qr.resid(qr(xmat), x), ax))
    }
  } else if (is.matrix(fl)) {
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    xmat <- cbind(Intercept = 1L, fl)
    if(na.rm) {
      if(fill) {
        x[cc, ] <- qr.resid(qr(xmat[cc, , drop = FALSE]), x[cc, , drop = FALSE])
        x[-cc, ] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(qr.resid(qr(xmat[cc, , drop = FALSE]), x[cc, , drop = FALSE]), ax)) # best ??
    } else return(setAttributes(qr.resid(qr(xmat), x), ax))
  } else {
    ffll <- is.factor(fl)
    if(!ffll && !missing(...)) stop("Unknown argument ", dotstostr(...))
    if(na.rm) {
      Y <- if(ffll) demeanlist(x[cc, , drop = FALSE], list(fl[cc]), weights = w[cc], ...) else
        qr.resid(qr(cbind(1L,fl[cc])), x[cc, , drop = FALSE])
      if(fill) {
        x[cc, ] <- Y
        x[-cc, ] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(Y, ax))
    } else {
      if(ffll) return(setAttributes(demeanlist(x, list(fl), weights = w, ...), ax)) else
        return(setAttributes(qr.resid(qr(cbind(1L,fl)), x), ax))
    }
  }
}


fHDwithin.data.frame <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
  if(variable.wise) fill <- TRUE
  # Best solution ????:
  if(na.rm && !fill) x <- .Call(C_subsetDT, x, cc, Xvars) else if(na.rm && fill && !variable.wise)
    x <- .Call(C_subsetDT, x, NULL, Xvars) else if(!identical(Xvars,seq_along(x))) x <- x[Xvars]
    Xvars <- seq_along(x)

    if(!is.null(dat)) {
      if(na.rm && fill) {
        if(variable.wise) {
          ax <- attributes(x)
          ax[["cases"]] <- cc
          lfl <- length(fl) > 0
          if(cc[1] == "no missing cases in fl!") {
            if(lfl) dat <- demeanlist(dat[,-1, drop = FALSE], fl, weights = w, ...)
            return(setAttributes(lapply(x, function(y) {
              ycc <- which(!is.na(y))
              if(lfl) y[ycc] <- qr.resid(qr(dat[ycc, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)) else
                y[ycc] <- qr.resid(qr(dat[ycc, , drop = FALSE]), y[ycc])
              return(y)
            }), ax)) # Rfast fastlm??
          } else {
            if(lfl) dat <- demeanlist(dat[,-1, drop = FALSE], fl, weights = w[cc], ...)
            return(setAttributes(lapply(x, function(y) {
              ycc <- !is.na(y)
              ycc[-cc] <- FALSE
              y[-cc] <- NA # which is not faster !!
              if(lfl) {
                YC <- which(ycc[cc])
                y[ycc] <- qr.resid(qr(dat[YC, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, YC), weights = w[ycc], ...))
              } else y[ycc] <- qr.resid(qr(dat[ycc[cc], , drop = FALSE]), y[ycc])
              return(y)
            }), ax)) # Rfast fastlm??
          }
        } else {
          attr(x, "cases") <- cc
          if(length(fl)) { # fastest so far !! (and memory efficient)
            dat <- qr(demeanlist(dat[,-1], fl, weights = w[cc], ...)) # This is necessary!!!
            Y <- lapply(demeanlist(.Call(C_subsetDT, x, cc, Xvars), fl, weights = w[cc], ...), function(y) qr.resid(dat, y))
            return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
                         seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
          } else {
            dat <- qr(dat)
            Y <- lapply(.Call(C_subsetDT, x, cc, Xvars), function(y) qr.resid(dat, y))
            return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
                         seq_row(X)[-cc], Xvars, NULL, NA, FALSE))
          }
        }
      } else {
        ax <- attributes(x)
        if(na.rm) ax[["cases"]] <- cc
        if(length(fl)) {
          if(!is.null(w) && na.rm) w <- w[cc]
          dat <- qr(demeanlist(dat[,-1], fl, weights = w, ...)) # This is necessary!!!
          return(setAttributes(lapply(demeanlist(x, fl, weights = w, ...),
                                      function(y) qr.resid(dat, y)), ax))
        } else {
          dat <- qr(dat)
          return(setAttributes(lapply(x, function(y) qr.resid(dat, y)), ax))
        }
      }
    } else {
      if(na.rm && fill) {
        if(variable.wise) {
          ax <- attributes(x)
          ax[["cases"]] <- cc
          if(cc[1] == "no missing cases in fl!") {
            return(setAttributes(lapply(x, function(y) {
              ycc <- which(!is.na(y))
              y[ycc] <- demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)
              return(y)
            }), ax))
          } else {
            return(setAttributes(lapply(x, function(y) {
              ycc <- !is.na(y)
              ycc[-cc] <- FALSE
              y[-cc] <- NA # which is not faster !!
              y[xcc] <- demeanlist(y[ycc], subsetfl(fl, which(ycc[cc])), weights = w[ycc], ...)
              return(y)
            }), ax))
          }
        } else {
          attr(x, "cases") <- cc
          Y <- demeanlist(.Call(C_subsetDT, x, cc, Xvars), fl, weights = w[cc], ...)
          # NUMlassignCpp(Y, cc, Xvars)
          return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
                       seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
        }
      } else {
        if(na.rm) return(`attr<-`(demeanlist(x, fl, weights = w, ...), "cases", cc)) else
          return(demeanlist(x, fl, weights = w, ...))
      }
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
HDW.matrix <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  fHDwithin.matrix(x, fl, w, na.rm, fill, ...)
}


HDW.data.frame <- function(x, fl, w = NULL, cols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
  list2env(getPartData(x, fl, cols, na.rm, variable.wise), envir = environment())

  if(variable.wise) fill <- TRUE
  # Best solution ????:
  if(na.rm && !fill) x <- .Call(C_subsetDT, x, cc, Xvars) else if(na.rm && fill && !variable.wise)
    x <- .Call(C_subsetDT, x, NULL, Xvars) else if(!identical(Xvars,seq_along(x))) x <- x[Xvars]
    Xvars <- seq_along(x)

    if(!is.null(dat)) {
      if(na.rm && fill) {
        if(variable.wise) {
          ax <- attributes(x)
          ax[["cases"]] <- cc
          lfl <- length(fl) > 0
          if(cc[1] == "no missing cases in fl!") {
            if(lfl) dat <- demeanlist(dat[,-1, drop = FALSE], fl, weights = w, ...)
            return(setAttributes(lapply(x, function(y) {
              ycc <- which(!is.na(y))
              if(lfl) y[ycc] <- qr.resid(qr(dat[ycc, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)) else
                      y[ycc] <- qr.resid(qr(dat[ycc, , drop = FALSE]), y[ycc])
              return(y)
            }), ax)) # Rfast fastlm??
          } else {
            if(lfl) dat <- demeanlist(dat[,-1, drop = FALSE], fl, weights = w[cc], ...)
            return(setAttributes(lapply(x, function(y) {
              ycc <- !is.na(y)
              ycc[-cc] <- FALSE
              y[-cc] <- NA # which is not faster !!
              if(lfl) {
                YC <- which(ycc[cc])
                y[ycc] <- qr.resid(qr(dat[YC, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, YC), weights = w[ycc], ...))
              } else y[ycc] <- qr.resid(qr(dat[ycc[cc], , drop = FALSE]), y[ycc])
              return(y)
            }), ax)) # Rfast fastlm??
          }
        } else {
          attr(x, "cases") <- cc
          if(length(fl)) { # fastest so far !! (and memory efficient)
            dat <- qr(demeanlist(dat[,-1], fl, weights = w[cc], ...)) # This is necessary!!!
            Y <- lapply(demeanlist(.Call(C_subsetDT, x, cc, Xvars), fl, weights = w[cc], ...), function(y) qr.resid(dat, y))
            return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
                         seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
          } else {
            dat <- qr(dat)
            Y <- lapply(.Call(C_subsetDT, x, cc, Xvars), function(y) qr.resid(dat, y))
            return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
                         seq_row(X)[-cc], Xvars, NULL, NA, FALSE))
          }
        }
      } else {
        ax <- attributes(x)
        if(na.rm) ax[["cases"]] <- cc
        if(length(fl)) {
          if(!is.null(w) && na.rm) w <- w[cc]
          dat <- qr(demeanlist(dat[,-1], fl, weights = w, ...)) # This is necessary!!!
          return(setAttributes(lapply(demeanlist(x, fl, weights = w, ...),
                                       function(y) qr.resid(dat, y)), ax))
        } else {
          dat <- qr(dat)
          return(setAttributes(lapply(x, function(y) qr.resid(dat, y)), ax))
        }
      }
    } else {
      if(na.rm && fill) {
        if(variable.wise) {
          ax <- attributes(x)
          ax[["cases"]] <- cc
          if(cc[1] == "no missing cases in fl!") {
            return(setAttributes(lapply(x, function(y) {
              ycc <- which(!is.na(y))
              y[ycc] <- demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)
              return(y)
            }), ax))
          } else {
            return(setAttributes(lapply(x, function(y) {
              ycc <- !is.na(y)
              ycc[-cc] <- FALSE
              y[-cc] <- NA # which is not faster !!
              y[xcc] <- demeanlist(y[ycc], subsetfl(fl, which(ycc[cc])), weights = w[ycc], ...)
              return(y)
            }), ax))
          }
        } else {
          attr(x, "cases") <- cc
          Y <- demeanlist(.Call(C_subsetDT, x, cc, Xvars), fl, weights = w[cc], ...)
          # NUMlassignCpp(Y, cc, Xvars)
          return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
                       seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
        }
      } else {
        if(na.rm) return(`attr<-`(demeanlist(x, fl, weights = w, ...), "cases", cc)) else
        return(demeanlist(x, fl, weights = w, ...))
      }
    }
}
HDW.pdata.frame <- function(x, w = NULL, cols = is.numeric, na.rm = TRUE, fill = TRUE,
                            variable.wise = FALSE, drop.xt = FALSE, drop.w = TRUE, ...) {
  ax <- attributes(x)
  nam <- ax[["names"]]
  gn <- match(names(ax[["index"]]), nam)
  gn2 = gn <- gn[!is.na(gn)]
  if(!is.null(cols)) {
    if(is.function(cols)) cols <- seq_along(x)[!vapply(x, cols, TRUE)] else if(is.character(cols))
      cols <- seq_along(x)[-match(cols, nam)] else if(is.logical(cols))
        cols <- seq_along(x)[!cols] else if(is.numeric(cols))
          cols <- seq_along(x)[-cols] else stop("cols needs to be a function, column names, indices or a logical vector")
        if(drop.xt) gn <- unique.default(c(gn, cols)) else if(length(gn)) gn2 <- unique.default(c(gn2, cols)) else gn <- cols
  }
  if(!is.null(w)) {
    if(is.call(w)) w <- all.vars(w)
    if(length(w) == 1) {
      wn <- match(w, nam)
      if(any(gn == wn)) stop("Weights coincide with grouping variables!")
      w <- x[[wn]]
      if(drop.w) if(!length(gn)) x[[wn]] <- NULL else if(drop.xt) gn <- c(gn,wn) else gn2 <- c(gn2,wn)
    }
  }
  varwisecomp <- function(x, fl, w, ...) lapply(x, function(y) {
    ycc <- which(!is.na(y))
    y[ycc] <- demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)
    return(y)
  })
  fillcomp <- function(x, fl, w, cc, Xvars, ...) {
    if(length(cc) != length(x[[1L]])) {
      Y <- demeanlist(.Call(C_subsetDT, x, cc, Xvars), lapply(fl, function(y) y[cc]), weights = w[cc], ...)
      return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
                   seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
    } else return(demeanlist(x[Xvars], fl, weights = w, ...))
  }
  nonfillcomp <- function(x, fl, w, cc, Xvars, ...) {
    if(length(cc) != length(x[[1L]]))
      return(demeanlist(.Call(C_subsetDT, x, cc, Xvars), lapply(fl, function(y) y[cc]), weights = w[cc], ...)) else
      return(demeanlist(x[Xvars], fl, weights = w, ...))
  }
  if(na.rm && fill) {
    if(variable.wise) {
      attributes(x) <- NULL
      if(length(gn)) {
        if(drop.xt) {
          ax[["names"]] <- if(give.names) paste0("HDW.",nam[-gn]) else nam[-gn]
          return(setAttributes(varwisecomp(x[-gn], ax[["index"]], w, ...), ax))
        } else {
            ax[["names"]] <- c(nam[gn], if(give.names) paste0("HDW.",nam[-gn2]) else nam[-gn2])
            return(setAttributes(c(x[gn],varwisecomp(x[-gn2], ax[["index"]], w, ...)), ax))
        }
      } else {
        if(give.names) paste0("HDW.",nam)
        return(setAttributes(varwisecomp(x, ax[["index"]], w, ...), ax))
      }
    } else {
      cc <- which(!.Call(C_dt_na, x, seq_along(x)[-c(gn,gn2)])) # good??
      ax[["cases"]] <- if(length(cc) == length(x[[1L]])) "no missing cases!" else cc
      if(length(gn)) {
        if(drop.xt) {
          ax[["names"]] <- if(give.names) paste0("HDW.",nam[-gn]) else nam[-gn]
          return(setAttributes(fillcomp(x, ax[["index"]], w, cc, seq_along(x)[-gn], ...), ax))
        } else {
          ax[["names"]] <- c(nam[gn], if(give.names) paste0("HDW.",nam[-gn2]) else nam[-gn2])
          return(setAttributes(c(x[gn], fillcomp(x, ax[["index"]], w, cc, seq_along(x)[-gn2], ...)), ax))
        }
      } else {
        if(give.names) paste0("HDW.",nam)
        return(setAttributes(fillcomp(x, ax[["index"]], w, cc, seq_along(x), ...), ax))
      }
    }
  } else {
    if(na.rm) {
      cc <- which(!.Call(C_dt_na, x, seq_along(x)[-c(gn,gn2)])) # good??
      if(length(cc) == length(x[[1L]])) ax[["cases"]] <- "no missing cases!" else {
        ax[["row.names"]] <- ax[["row.names"]][cc]
        ax[["cases"]] <- cc
      }
      if(length(gn)) {
        if(drop.xt) {
          ax[["names"]] <- if(give.names) paste0("HDW.",nam[-gn]) else nam[-gn]
          return(setAttributes(nonfillcomp(x, ax[["index"]], w, cc, seq_along(x)[-gn], ...), ax))
        } else {
          ax[["names"]] <- c(nam[gn], if(give.names) paste0("HDW.",nam[-gn2]) else nam[-gn2])
          return(setAttributes(c(.Call(C_subsetDT, x, cc, gn), nonfillcomp(x, ax[["index"]], w, cc, seq_along(x)[-gn2], ...)), ax))
        }
      } else {
        if(give.names) paste0("HDW.",nam)
        return(setAttributes(nonfillcomp(x, ax[["index"]], w, cc, seq_along(x), ...), ax))
      }
    } else {
      attributes(x) <- NULL
      if(length(gn)) {
        if(drop.xt) {
          ax[["names"]] <- if(give.names) paste0("HDW.",nam[-gn]) else nam[-gn]
          return(setAttributes(demeanlist(x[-gn], ax[["index"]], weights = w, ...), ax))
        } else {
          ax[["names"]] <- c(nam[gn], if(give.names) paste0("HDW.",nam[-gn2]) else nam[-gn2])
          return(setAttributes(c(x[gn],demeanlist(x[-gn2], ax[["index"]], weights = w, ...)), ax))
        }
      } else {
        if(give.names) paste0("HDW.",nam)
        return(setAttributes(demeanlist(x, ax[["index"]], weights = w, ...), ax))
      }
    }
  }
}

# Theory: y = ?1 x1 + ?2 x2 + e
# FWT: M2 y = ?1 M2 x1 + e so residuals: e = M2 y - ?1 M2 x1 and fitted:
# Now M = I - x(x'x)-1x' = I - P.
# So (I-P2) y = ?1 (I-P2) x1 + e or y - P2 y = ?1 x1 - ?1 P2 x1 + e

# I want y - e = y^ = ?1 x1 + ?2 x2
# so
# P2 y = ?1 P2 x1 + ?2 x2
# Haven't quite figgured it out, but my solution is to just subtract the demeaned data !!

HDB <- function(x, ...) { # fl, w = NULL, na.rm = TRUE, fill = FALSE,
  UseMethod("HDB", x)
}
HDB.default <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  ax <- attributes(x)
  if(na.rm) {
    cc <- which(complete.cases(x, fl)) # gives error if lengths don't match, otherwise demeanlist and qr.fitted give errors !!
    if(length(cc) != length(x)) {
      if(!fill) ax[["names"]] <- NULL
      ax[["cases"]] <- cc
    } else {
      ax[["cases"]] <- "no missing cases!"
      na.rm <- FALSE
    }
  }
  if(is.list(fl)) {
    fc <- vapply(fl, is.factor, TRUE)
    fcl <- any(fc)
    nallfc <- !all(fc)
    if(na.rm) fl <- if(is.data.frame(fl)) .Call(C_subsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    if(fcl && nallfc) {
      fl <- fl[fc]
      xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl[!fc]))
    } else if(!fcl) xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl))
    if(na.rm) {
      Y <- if(fcl && nallfc) x[cc] - qr.resid(qr(demeanlist(xmat, fl, weights = w[cc], ...)),
                                          demeanlist(x[cc], fl, weights = w[cc], ...)) else if(fcl)
                                          demeanlist(x[cc], fl, weights = w[cc], means = TRUE, ...) else qr.fitted(qr(xmat), x[cc])
      if(fill) {
        x[cc] <- Y
        x[-cc] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(Y, ax))
    } else {
      if(fcl && nallfc)
        return(setAttributes(x - qr.resid(qr(demeanlist(xmat, fl, weights = w, ...)),
                                       demeanlist(x, fl, weights = w, ...)), ax)) else if(fcl)
                                         return(setAttributes(demeanlist(x, fl, weights = w, means = TRUE, ...), ax)) else
                                           return(setAttributes(qr.fitted(qr(xmat), x), ax))
    }
  } else if (is.matrix(fl)) {
    xmat <- cbind(Intercept = 1L, fl)
    if(na.rm) {
      if(fill) {
        x[cc] <- qr.fitted(qr(xmat[cc, , drop = FALSE]), x[cc])
        x[-cc] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(qr.fitted(qr(xmat[cc, , drop = FALSE]), x[cc]), ax)) # best ??
    } else return(setAttributes(qr.fitted(qr(xmat), x), ax))
  } else {
    if(na.rm) {
      Y <- if(is.factor(fl)) demeanlist(x[cc], list(fl[cc]), weights = w[cc], means = TRUE, ...) else
        qr.fitted(qr(cbind(1L,fl[cc])), x[cc])
      if(fill) {
        x[cc] <- Y
        x[-cc] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(Y, ax))
    } else {
      if(is.factor(fl)) return(setAttributes(demeanlist(x, list(fl), weights = w, means = TRUE, ...), ax)) else
        return(setAttributes(qr.fitted(qr(cbind(1L,fl)), x), ax))
    }
  }
}
HDB.matrix <- function(x, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  ax <- attributes(x)
  if(na.rm) {
    cc <- which(complete.cases(x, fl)) # gives error if lengths don't match, otherwise demeanlist and qr.fitted give errors !!
    if(length(cc) != nrow(x)) {
      if(!fill) {
        ax[["dimnames"]][[1L]] <- ax[["dimnames"]][[1L]][cc] # best ??
        ax[["dim"]][1] <- length(cc)
      }
      ax[["cases"]] <- cc
    } else {
      na.rm <- FALSE
      ax[["cases"]] <- "no missing cases!"
    }
  }
  if(is.list(fl)) {
    fc <- vapply(fl, is.factor, TRUE)
    fcl <- any(fc)
    nallfc <- !all(fc)
    if(na.rm) fl <- if(is.data.frame(fl)) .Call(C_subsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    if(fcl && nallfc) {
      fl <- fl[fc]
      xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl[!fc]))
    } else if(!fcl) xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1L]]))), fl))
    if(na.rm) {
      Y <- if(fcl && nallfc) x[cc, , drop = FALSE] - qr.resid(qr(demeanlist(xmat, fl, weights = w[cc], ...)),
                                        demeanlist(x[cc, , drop = FALSE], fl, weights = w[cc], ...)) else if(fcl)
                                        demeanlist(x[cc, , drop = FALSE], fl, weights = w[cc], means = TRUE, ...) else qr.fitted(qr(xmat), x[cc, , drop = FALSE])
      if(fill) {
        x[cc, ] <- Y
        x[-cc, ] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(Y, ax))
    } else {
      if(fcl && nallfc)
        return(setAttributes(x - qr.resid(qr(demeanlist(xmat, fl, weights = w, ...)),
                                       demeanlist(x, fl, weights = w, ...)), ax)) else if(fcl)
                                         return(setAttributes(demeanlist(x, fl, weights = w, means = TRUE, ...), ax)) else
                                           return(setAttributes(qr.fitted(qr(xmat), x), ax))
    }
  } else if (is.matrix(fl)) {
    xmat <- cbind(Intercept = 1L, fl)
    if(na.rm) {
      if(fill) {
        x[cc, ] <- qr.fitted(qr(xmat[cc, , drop = FALSE]), x[cc, , drop = FALSE])
        x[-cc, ] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(qr.fitted(qr(xmat[cc, , drop = FALSE]), x[cc, , drop = FALSE]), ax)) # best ??
    } else return(setAttributes(qr.fitted(qr(xmat), x), ax))
  } else {
    if(na.rm) {
      Y <- if(is.factor(fl)) demeanlist(x[cc, , drop = FALSE], list(fl[cc]), weights = w[cc], means = TRUE, ...) else
           qr.fitted(qr(cbind(1L,fl[cc])), x[cc, , drop = FALSE])
      if(fill) {
        x[cc, ] <- Y
        x[-cc, ] <- NA
        return(setAttributes(x, ax))
      } else return(setAttributes(Y, ax))
    } else {
      if(is.factor(fl)) return(setAttributes(demeanlist(x, list(fl), weights = w, means = TRUE, ...), ax)) else
        return(setAttributes(qr.fitted(qr(cbind(1L,fl)), x), ax))
    }
  }
}
HDB.data.frame <- function(x, fl, w = NULL, cols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
  list2env(getPartData(x, fl, cols, na.rm, variable.wise), envir = environment())

  if(variable.wise) fill <- TRUE
  # Best solution ????:
  if(na.rm && !fill) x <- .Call(C_subsetDT, x, cc, Xvars) else if(na.rm && fill && !variable.wise)
    x <- .Call(C_subsetDT, x, NULL, Xvars) else if(!identical(Xvars,seq_along(x))) x <- x[Xvars]
    Xvars <- seq_along(x)

    if(!is.null(dat)) {
      if(na.rm && fill) {
        if(variable.wise) {
          ax <- attributes(x)
          ax[["cases"]] <- cc
          lfl <- length(fl) > 0
          if(cc[1] == "no missing cases in fl!") {
            if(lfl) dat <- demeanlist(dat[,-1, drop = FALSE], fl, weights = w, ...)
            return(setAttributes(lapply(x, function(y) {
              ycc <- which(!is.na(y))
              if(lfl) y[ycc] <- y[ycc] - qr.resid(qr(dat[ycc, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], ...)) else
                y[ycc] <- qr.fitted(qr(dat[ycc, , drop = FALSE]), y[ycc])
              return(y)
            }), ax)) # Rfast fastlm??
          } else {
            if(lfl) dat <- demeanlist(dat[,-1, drop = FALSE], fl, weights = w[cc], ...)
            return(setAttributes(lapply(x, function(y) { # change !!
              ycc <- !is.na(y)
              ycc[-cc] <- FALSE
              y[-cc] <- NA # which is not faster !!
              if(lfl) {
                YC <- which(ycc[cc])
                y[ycc] <- y[ycc] - qr.resid(qr(dat[YC, , drop = FALSE]), demeanlist(y[ycc], subsetfl(fl, YC), weights = w[ycc], ...))
              } else y[ycc] <- qr.fitted(qr(dat[ycc[cc], , drop = FALSE]), y[ycc])
              return(y)
            }), ax)) # Rfast fastlm??
          }
        } else {
          attr(x, "cases") <- cc
          if(length(fl)) { # fastest so far !! (and memory efficient)
            dat <- qr(demeanlist(dat[,-1], fl, weights = w[cc], ...)) # This is necessary!!!
            Y <- lapply(.Call(C_subsetDT, x, cc, Xvars), function(y) y - qr.resid(dat, demeanlist(y, fl, weights = w[cc], ...))) # good ??
            return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
                         seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
          } else {
            dat <- qr(dat)
            Y <- lapply(.Call(C_subsetDT, x, cc, Xvars), function(y) qr.fitted(dat, y))
            return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
                         seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
          }
        }
      } else {
        ax <- attributes(x)
        if(na.rm) ax[["cases"]] <- cc
        if(length(fl)) {
          if(!is.null(w) && na.rm) w <- w[cc]
          dat <- qr(demeanlist(dat[,-1], fl, weights = w, ...)) # This is necessary!!!
          return(setAttributes(Map(function(x, y) x - qr.resid(dat, y), x, demeanlist(x, fl, weights = w, ...)), ax))
        } else {
          dat <- qr(dat)
          return(setAttributes(lapply(x, function(y) qr.fitted(dat, y)), ax))
        }
      }
    } else {
      if(na.rm && fill) {
        if(variable.wise) {
          ax <- attributes(x)
          ax[["cases"]] <- cc
          if(cc[1] == "no missing cases in fl!") {
            return(setAttributes(lapply(x, function(y) {
              ycc <- which(!is.na(y))
              y[ycc] <- demeanlist(y[ycc], subsetfl(fl, ycc), weights = w[ycc], means = TRUE, ...)
              return(y)
            }), ax))
          } else {
            return(setAttributes(lapply(x, function(y) {
              ycc <- !is.na(y)
              ycc[-cc] <- FALSE
              y[-cc] <- NA # which is not faster !!
              y[ycc] <- demeanlist(y[ycc], subsetfl(fl, which(ycc[cc])), weights = w[ycc], means = TRUE, ...)
              return(y)
            }), ax))
          }
        } else {
          attr(x, "cases") <- cc
          Y <- demeanlist(.Call(C_subsetDT, x, cc, Xvars), fl, weights = w[cc], means = TRUE, ...)
          return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, x, cc, Xvars, NULL, Y, FALSE),
                       seq_row(x)[-cc], Xvars, NULL, NA, FALSE))
        }
      } else {
        if(na.rm) return(`attr<-`(demeanlist(x, fl, weights = w, means = TRUE, ...), "cases", cc)) else
          return(demeanlist(x, fl, weights = w, means = TRUE, ...))
      }
    }
}


# todo: good performance with missing values ??
# -> do pseries and pdata.frame method !!


#
# HDW(X = mtcars, fl = ~ factor(cyl)*carb)
#
# HDW(X = mtcars, fl = ~ factor(cyl):vs)
#
# lm(mpg ~ factor(cyl):factor(vs), data = mtcars)
#
# HDW(X = mtcars, fl = ~ factor(cyl)*carb + vs + wt:gear + wt:gear:carb)
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



# Revious Version: Not subsetting X at the top !!
# HDW.data.frame.old <- function(X, fl, w = NULL, cols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
#   list2env(getPartData(X, fl, cols, na.rm, variable.wise), envir = environment())
#   ay <- attributes(Y)
#   if(na.rm) ay[["cases"]] <- cc # good !!
#   if(!is.null(dat)) {
#     if(na.rm && fill) {
#       if(variable.wise) {
#         lfl <- length(fl) > 0
#         if(lfl) dat <- demeanlist(dat[,-1], fl, weights = w[cc], ...)
#         return(setAttributes(lapply(Y, function(x) {
#           xcc <- !is.na(x)
#           xcc[-cc] <- FALSE
#           x[-cc] <- NA # which is not faster !!
#           if(lfl) {
#             XC <- which(xcc[cc])
#             x[xcc] <- qr.resid(qr(dat[XC, , drop = FALSE]), demeanlist(x[xcc], subsetfl(fl, XC), weights = w[xcc], ...))
#           } else x[xcc] <- qr.resid(qr(dat[xcc[cc], , drop = FALSE]), x[xcc])
#           return(x)
#         }), ay)) # Rfast fastlm??
#       } else {
#         if(!is.null(fl)) { # fastest so far !! (and memory efficient)
#           dat <- qr(demeanlist(dat[,-1], fl, weights = w[cc], ...)) # This is necessary!!!
#           Y <- lapply(demeanlist(.Call(CSDT, X, cc, Xvars), fl, weights = w[cc], ...), function(x) qr.resid(dat,x))
#           return(.Call(AS, .Call(AS, .Call(CSDT, X, NULL, Xvars), cc, seq_along(Xvars), NULL, Y, FALSE),
#                        seq_row(X)[-cc], seq_along(Xvars), NULL, NA, FALSE))
#         } else {
#           dat <- qr(dat)
#           Y <- lapply(.Call(CSDT, X, cc, Xvars), function(x) qr.resid(dat,x))
#           return(.Call(AS, .Call(AS, .Call(CSDT, X, NULL, Xvars), cc, seq_along(Xvars), NULL, Y, FALSE),
#                        seq_row(X)[-cc], seq_along(Xvars), NULL, NA, FALSE))
#         }
#       }
#     } else {
#       if(length(fl)) {
#         if(!is.null(w) && na.rm) w <- w[cc]
#         dat <- qr(demeanlist(dat[,-1], fl, weights = w, ...)) # This is necessary!!!
#         if(na.rm) X <- .Call(CSDT, X, cc, Xvars) else if(!identical(Xvars,seq_along(X))) X <- X[Xvars]
#         return(setAttributes(lapply(demeanlist(X, fl, weights = w, ...),
#                                      function(x) qr.resid(dat, x)), ay))
#       } else {
#         dat <- qr(dat)
#         if(na.rm) X <- .Call(CSDT, X, cc, Xvars) else if(!identical(Xvars,seq_along(X))) X <- X[Xvars]
#         return(setAttributes(lapply(X, function(x) qr.resid(dat, x)), ay))
#       }
#     }
#   } else {
#     if(na.rm && fill) {
#       if(variable.wise) {
#         return(setAttributes(lapply(Y, function(x) {
#           xcc <- !is.na(x)
#           xcc[-cc] <- FALSE
#           x[-cc] <- NA # which is not faster !!
#           x[xcc] <- demeanlist(x[xcc], subsetfl(fl, which(xcc[cc])), weights = w[xcc], ...)
#           return(x)
#         }), ay))
#       } else return(.Call(AS, .Call(AS, .Call(CSDT, X, NULL, Xvars), cc, seq_along(Xvars), NULL,
#                                     demeanlist(Y, fl, weights = w[cc], ...), FALSE),
#                           seq_len(nrow(X))[-cc], seq_along(Xvars),
#                           NULL, NA, FALSE))
#
#     } else return(demeanlist(Y, fl, weights = w, ...))
#   }
# }


# HDW(X = na.omit(GGDC),fl = ~factor(Country))

# Previous Getpartdata -> generating Y before passing to HDW.data.frame

# getPartData <- function(X, fl, cols, na.rm, variable.wise) {
#   res <- list(fl = NULL, dat = NULL, Y = NULL, cc = NULL, Xvars = NULL)
#   nam <- names(X)
#   form <- class(fl) == "formula"
#   twosided <- form && length(fl) == 3
#   colsl <- !is.null(cols)
#   if(!twosided && colsl) {
#     if(is.function(cols))
#       ind <- which(vapply(X, cols, TRUE)) else if(is.numeric(cols))
#         ind <- cols else if(is.character(cols))
#           ind <- match(cols, nam) else if(is.logical(cols))
#             ind <- which(cols) else
#               stop("cols needs to be a function, column names, indices or a logical vector")
#   } else ind <- seq_along(nam)
#
#   if(form) { # fl is formula
#     if(twosided) {
#       fvars <- match(all.vars(fl[[3L]]), nam)
#       Xvars <- match(all.vars(fl[[2L]]), nam)
#       fl[[2L]] <- NULL
#     } else {
#       fvars <- match(all.vars(fl), nam)
#       Xvars <- setdiff(ind, fvars)
#     }
#     res[[5L]] <- Xvars
#     if(na.rm) {
#       if(variable.wise) {
#         cc <- which(!.Call(NADT, X, fvars)) # best ??
#         res[1:2] <- getfl(myModFrame(fl, .Call(CSDT, X, cc, fvars)))
#         res[[3L]] <- .Call(CSDT, X, NULL, Xvars)
#       } else {
#         cc <- which(!.Call(NADT, X, c(Xvars,fvars))) # best ??
#         res[1:2] <- getfl(myModFrame(fl, .Call(CSDT, X, cc, fvars)))
#         res[[3L]] <- .Call(CSDT, X, cc, Xvars)
#       }
#       res[[4L]] <- cc
#     } else {
#       res[1:2] <- getfl(myModFrame(fl, X))
#       res[[3L]] <- .Call(CSDT, X, NULL, Xvars)
#     }
#   } else if(is.list(fl)) { # fl is factor list
#     class(fl) <- NULL
#     res[[5L]] <- ind
#     fc <- vapply(fl, is.factor, TRUE)
#     fcl <- any(!fc)
#     if(na.rm) {
#       if(variable.wise) {
#         cc <- which(!.Call(NADT, fl, seq_along(fl))) # best ??
#         if(fcl) {
#           fl <- subsetfl(fl, cc)
#           res[[1L]] <- fl[fc]
#           res[[2L]] <- do.call(cbind, c(Intercept = rep(1, length(fl[[1L]])), fl[!fc]))
#         } else res[[1L]] <- subsetfl(fl, cc)
#         res[[3L]] <- if(colsl) .Call(CSDT, X, NULL, ind) else X
#       } else {
#         cc <- which(!.Call(NADT, c(X,fl), seq_len(length(X)+length(fl)))) # best ??
#         if(fcl) {
#           fl <- subsetfl(fl, cc)
#           res[[1L]] <- fl[fc]
#           res[[2L]] <- do.call(cbind, c(Intercept = rep(1, length(fl[[1L]])), fl[!fc]))
#         } else res[[1L]] <- subsetfl(fl, cc)
#         res[[3L]] <- .Call(CSDT, X, cc, ind)
#       }
#       res[[4L]] <- cc
#     } else {
#       if(fcl) {
#         res[[1L]] <- fl[fc]
#         res[[2L]] <- do.call(cbind, c(Intercept = rep(1, length(fl[[1L]])), fl[!fc]))
#       } else res[[1L]] <- fl
#       res[[3L]] <- if(colsl) .Call(CSDT, X, NULL, ind) else X
#     }
#   } else { # fl is vector or matrix !!
#     res[[5L]] <- ind
#     if(is.factor(fl)) stop("For single factors please use W(), not HDW(), or pass the factor inside a list.")
#     xmat <- cbind(Intercept = 1L, fl)
#     if(na.rm) {
#       if(variable.wise) {
#         cc <- which(complete.cases(xmat))
#         res[[2L]] <- xmat[cc, , drop = FALSE]
#         res[[3L]] <- if(colsl) .Call(CSDT, X, NULL, ind) else X
#       } else {
#         cc <- which(complete.cases(X, xmat))
#         res[[2L]] <- xmat[cc, , drop = FALSE]
#         res[[3L]] <- .Call(CSDT, X, cc, ind)
#       }
#       res[[4L]] <- cc
#     } else {
#       res[[2L]] <- xmat
#       res[[3L]] <- if(colsl) .Call(CSDT, X, NULL, ind) else X
#     }
#   }
#   return(res)
# }

# regression functions:
# regres <- function(y, x, residuals = TRUE) { # This is the problem, in doing this projection, you are already assuming intercept removal!!
#   dmx = x-ave(x) # NOT demeaning -> maybe that solves the interaction problem!!-> Nope, does not converge!!
#   ay = ave(y)
#   dmy = y-ay # -> maybe you need to add the group average for y, or for the interaction variable??
#   beta = dmx%*%dmy/dmx%*%dmx # much faster than cov(dmx,dmy)/var(dmx) !!
#   if (residuals) {
#     dmy - dmx%*%beta # residuals -> add ave(y)?? -> nope!!
#   } else
#     ay + dmx%*%beta # fitted values
# }

# regresmat <- function(y, x) { # Write in C++!!, because the solve option is still faster!!
#   dmy = y-ave(y) # also does not give the same result!!
#   xlist = mctl(x - colMeans(x))
#   beta = numeric(1)
#   for (j in xlist) {
#     beta = j%*%dmy/j%*%j
#     dmy = dmy - j%*%beta
#   }
#   dmy
# }


# mydml without na.rm -> same speed as demeanlist
#
# mydml <- function(mtx, fl, icpt = 0L, eps = getOption("lfe.eps"), threads = getOption("lfe.threads"),
#                   progress = getOption("lfe.pint"), accel = getOption("lfe.accel"),
#                   randfact = TRUE, means = FALSE, weights = NULL, scale = TRUE, attrs = NULL) {
#   if (length(fl) == 0) {
#     if (means) {
#       foo <- unlist(utils::as.relistable(mtx))
#       foo[] <- 0
#       return(utils::relist(foo))
#     }
#     return(eval.parent(substitute(mtx)))
#   }
#   mf <- match.call()
#   ff <- formals(sys.function())
#   ff <- ff[match(c("mtx", "fl", "icpt", "eps", "threads",
#                    "progress", "accel", "means", "weights", "scale", "attrs"),
#                  names(ff), 0L)]
#   m <- names(ff)[names(ff) %in% names(mf)]
#   ff[m] <- mf[m]
#   env <- new.env(parent = parent.frame())
#   assign("C_demeanlist", lfe:::C_demeanlist, envir = env)
#   assign("unnamed", lfe:::unnamed, envir = env)
#   assign(".fl", eval.parent(ff[["fl"]]), envir = env)
#   .fl <- NULL
#   ff[["fl"]] <- quote(.fl)
#   if(randfact && length(get(".fl", env)) > 2) {
#     delayedAssign("..fl", .fl[order(runif(length(.fl)))],
#                   env, env)
#     ff[["fl"]] <- quote(..fl)
#   }
#   eval(as.call(c(list(quote(.Call), quote(C_demeanlist)), ff)), env)
# }


# Previous mydml and HDW: na.rm inside mydml:
# mydml <- function(mtx, fl, icpt = 0L, eps = getOption("lfe.eps"), threads = getOption("lfe.threads"),
#                   progress = getOption("lfe.pint"), accel = getOption("lfe.accel"),
#                   randfact = TRUE, means = FALSE, weights = NULL, scale = TRUE,
#                   na.rm = FALSE, fill = TRUE, attrs = NULL) {
#   if (length(fl) == 0) {
#     if (means) {
#       foo <- unlist(utils::as.relistable(mtx))
#       foo[] <- 0
#       return(utils::relist(foo))
#     }
#     return(eval.parent(substitute(mtx)))
#   }
#   mf <- match.call()
#   ff <- formals(sys.function())
#   ff <- ff[match(c("mtx", "fl", "icpt", "eps", "threads",
#                    "progress", "accel", "means", "weights", "scale", "attrs"),
#                  names(ff), 0L)]
#   m <- names(ff)[names(ff) %in% names(mf)]
#   ff[m] <- mf[m]
#   env <- new.env(parent = parent.frame())
#   assign("C_demeanlist", lfe:::C_demeanlist, envir = env)
#   assign("unnamed", lfe:::unnamed, envir = env)
#   assign(".fl", eval.parent(ff[["fl"]]), envir = env)
#   .fl <- NULL
#   ff[["fl"]] <- quote(.fl)
#   if(randfact && length(get(".fl", env)) > 2) {
#     delayedAssign("..fl", .fl[order(runif(length(.fl)))],
#                   env, env)
#     ff[["fl"]] <- quote(..fl)
#   }
#   if(na.rm) {
#     mtx <- eval.parent(ff[["mtx"]])
#     imtx <- seq_along(mtx)
#     badrows <- which(.Call("Cdt_na", mtx, imtx)) # !complete.cases(mtx) # faster !!
#     if(length(badrows) > 0) {
#       newmtx <- .Call("CsubsetDT", mtx, -badrows, imtx)
#       if(!fill) rm(mtx)
#       fl <- eval(ff[["fl"]], env)
#       N <- length(fl[[1L]])
#       newfl <- lapply(fl, function(f) { # use CsubsetDT or CsubsetVector ?? also check NA in regressors ??
#         x <- attr(f, "x")
#         if(is.null(x)) return(f[-badrows]) else
#           return(`attr<-`(f[-badrows], "x",
#                           if(is.matrix(x)) x[-badrows, , drop = FALSE] else
#                             x[-badrows]))
#       })
#       assign(".mtx", newmtx, envir = env)
#       assign(".fl", newfl, envir = env)
#       ff[["fl"]] <- quote(.fl)
#       rm(newfl)
#       ff[["attrs"]] <- c(ff[["attrs"]], list(na.rm = badrows))
#       ff[["mtx"]] <- quote(unnamed(.mtx))
#     } else {
#       assign(".mtx", mtx, envir = env)
#       ff[["mtx"]] <- quote(.mtx)
#     }
#   }
#   res <- eval(as.call(c(list(quote(.Call), quote(C_demeanlist)), ff)), env)
#   if(na.rm && fill && length(badrows) > 0) {
#     # mtx <- .Call("Cassign", mtx, badrows, imtx, NULL, NA, FALSE) #     mtx[badrows, ] <- NA #  .Call(Cassign,x,irows,cols,newnames,jval,verbose)
#     return(.Call("Cassign", .Call("Cassign", mtx, badrows, imtx, NULL, NA, FALSE),
#                  seq_len(nrow(mtx))[-badrows], imtx, NULL, res, FALSE)) # mtx[-badrows, ] <- res
#   } else return(res)
# }

# HDW.data.frame <- function(X, fl, w = NULL, cols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
#   list2env(getPartData(X, fl, cols, na.rm, variable.wise), envir = environment())
#
#   if(length(dat)) {
#     if(length(fl)) {
#       r <- mydml(Y, fl, weights = w, na.rm = na.rm, fill = fill, ...)
#       qrdat <- qr(mydml(dat[,-1], fl, weights = w, na.rm = na.rm, fill = fill, ...)) # This is necessary!!!
#       if(is.atomic(r)) qr.resid(qrdat,r) else
#         duplicate_attributes(lapply(r, function(x) qr.resid(qrdat,x)),r) # Rfast fastlm??
#     } else {
#       if(na.rm) {
#         if(fill && variable.wise) {
#           cc <- complete.cases(dat) # fastest ?? which ??
#           qrdat <- qr(dat[cc, , drop = FALSE])
#           res <- if(is.atomic(Y)) qr.resid(qrdat,Y[cc & complete.cases(Y)]) else
#             duplicate_attributes(lapply(Y, function(x) {
#               xNA <- is.na(x)
#               ccc <- cc & !xNA # fastest ?? which ??
#               x[ccc] <- qr.resid(qrdat[!xNA[cc],], x[ccc])
#               x[!ccc] <- NA
#               return(x)
#             }), Y) # Rfast fastlm??
#         } else {
#           cc <- which(complete.cases(dat, Y)) # fastest ?? what abour complete.cases(c(Y,dat))
#           qrdat <- qr(dat[cc, , drop = FALSE])
#           if(is.atomic(Y)) { # can it be a matrix ??
#             if(fill) {
#               Y[cc] <- qr.resid(qrdat, Y[cc])
#               Y[-cc] <- NA
#               return(Y)
#             } else return(`attr<-`(qr.resid(qrdat, Y[cc]), "cases", cc))
#           } else {
#             if(fill) {
#               Y[cc,] <- lapply(res, function(x) qr.resid(qrdat,x)) # Rfast fastlm??
#               Y[-cc,] <- NA
#               return(Y)
#             } else {
#               res <- Y[cc,]
#               ar <- attributes(res)
#               ar[["cases"]] <- cc
#               return(setAttributes(lapply(res, function(x) qr.resid(qrdat,x)), ar)) # Rfast fastlm??
#             }
#           }
#         }
#       } else {
#         qrdat <- qr(dat)
#         if(is.atomic(Y)) qr.resid(qrdat,Y) else
#           duplicate_attributes(lapply(Y,function(x)qr.resid(qrdat,x)),Y)
#       }
#     }
#   } else
#     mydml(Y, fl, weights = w, na.rm = na.rm, fill = fill, ...)
# }


# Old datapart : Also meant for vectors and matrices !!

# getPartData <- function(X, fl, cols) {
#   res <- list(fl = NULL, dat = NULL, Y = NULL)
#   if(is.list(X)) {
#     nam <- names(X)
#     ar <- FALSE
#   } else if (is.array(X)) {
#     nam <- colnames(X)
#     ar <- TRUE
#   } else {
#     nam <- NULL
#     ar <- FALSE
#   }
#   if(class(fl) == "formula") {
#     if(length(fl)<3) { # one -sided formula
#       if(!is.null(cols)) {
#         if(is.function(cols) && !ar)
#           ind <- which(vapply(X, cols, TRUE)) else if(is.numeric(cols))
#             ind <- cols else if(is.character(cols))
#               ind <- match(cols,nam) else if(is.logical(cols))
#                 ind <- which(cols) else
#                   stop("cols needs to be a function, column names, indices or a logical vector")
#       } else ind <- seq_along(nam)
#       fvars <- all.vars(fl)
#       res[1:2] <- getfl(myModFrame(fl, X))
#       if (!is.null(nam)) {
#         ind <- setdiff(ind, match(fvars, nam))
#         res[[3L]] <- if(ar) X[,ind] else X[ind]
#       } else res[[3L]] <- X
#     } else { # two-sided formula
#       fvars <- all.vars(fl[[3L]])
#       Xvars <- all.vars(fl[[2L]])
#       fl[[2L]] <- NULL
#       res[1:2] <- getfl(myModFrame(fl, X))
#       if(!is.null(nam)) {
#         res[[3L]] <- if(ar) X[,Xvars] else X[Xvars]
#       } else res[[3L]] <- X
#     }
#   } else { # fl is factor list !!
#     if(!is.null(cols)) {
#       if(is.function(cols) && !ar)
#         ind <- which(vapply(X, cols, TRUE)) else if(is.numeric(cols))
#           ind <- cols else if(is.character(cols))
#             ind <- match(cols, nam) else if (is.logical(cols))
#               ind <- which(cols) else
#                 stop("cols needs to be a function, column names, indices or a logical vector")
#             res[[3L]] <- if(ar) X[,ind] else X[ind]
#     } else res[[3L]] <- X
#     if(is.list(fl)) {
#       fc <- vapply(fl, is.factor, TRUE)
#       res[[1L]] <- c(Intercept = rep(1, length(fl[[1L]])), fl[!fc])
#       res[[2L]] <- fl[fc]
#     } else {
#       res[[1L]] <- c(Intercept = rep(1, length(fl[[1L]])), fl) # needed??
#       res[[3L]] <- NULL
#     }
#   } # still do else list of factors or data.frame or model matrix!!
#   return(res)
# }
