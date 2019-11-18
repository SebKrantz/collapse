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
    factors <- terms[[2]]
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
        `attributes<-`(rep(1L, NROW(x)), list(levels = "1", class = "factor", x = x)))
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
    if(is.null(x)) return(.Call(data.table:::CsubsetVector, f, cc)) else
      return(`attr<-`(.Call(data.table:::CsubsetVector, f, cc), "x",
                      if(is.matrix(x)) x[cc, , drop = FALSE] else
                        .Call(data.table:::CsubsetVector, x, cc)))
  })
}

getPartData <- function(X, fl, Xcols, na.rm, variable.wise) {
  res <- list(fl = NULL, dat = NULL, cc = NULL, Xvars = NULL, na.rm = na.rm)
  nam <- names(X)
  form <- all(class(fl) == "formula")
  twosided <- form && length(fl) == 3
  Xcolsl <- !is.null(Xcols)
  if(!twosided && Xcolsl) {
    if(is.function(Xcols))
    ind <- which(vapply(X, Xcols, TRUE)) else if(is.numeric(Xcols))
    ind <- Xcols else if(is.character(Xcols))
    ind <- match(Xcols, nam) else if(is.logical(Xcols))
    ind <- which(Xcols) else
    stop("Xcols needs to be a function, column names, indices or a logical vector")
  } else ind <- seq_along(X)

    if(form) { # fl is formula
      if(twosided) {
        fvars <- match(all.vars(fl[[3]]), nam)
        Xvars <- match(all.vars(fl[[2]]), nam)
        fl[[2]] <- NULL
      } else {
        fvars <- match(all.vars(fl), nam)
        Xvars <- setdiff(ind, fvars)
      }
      res[[4]] <- Xvars
      if(na.rm) {
        cc <- if(variable.wise) which(!.Call(data.table:::Cdt_na, X, fvars)) else which(!.Call(data.table:::Cdt_na, X, c(Xvars,fvars)))
        if(length(cc) != nrow(X)) {
          res[1:2] <- getfl(myModFrame(fl, .Call(data.table:::CsubsetDT, X, cc, fvars)))
          res[[3]] <- cc
        } else {
          res[1:2] <- getfl(myModFrame(fl, X[fvars]))
          if(variable.wise) res[[3]] <- "no missing cases in fl!" else {
            res[[3]] <- "no missing cases!"
            res[[5]] <- FALSE
          }
        }
      } else res[1:2] <- getfl(myModFrame(fl, X))
  } else if(is.list(fl)) { # fl is factor list
    class(fl) <- NULL
    res[[4]] <- ind
    fc <- vapply(fl, is.factor, TRUE)
    fcl <- any(!fc)
    if(na.rm) {
      cc <- if(variable.wise) which(!.Call(data.table:::Cdt_na, fl, seq_along(fl))) else
                              which(!.Call(data.table:::Cdt_na, c(X,fl), seq_len(length(X)+length(fl))))
      missc <- length(cc) != nrow(X)
        if(fcl) {
          if(missc) fl <- subsetfl(fl, cc)
          res[[1]] <- fl[fc]
          res[[2]] <- do.call(cbind, c(list(Intercept = rep(1, length(fl[[1]]))), fl[!fc]))
        } else res[[1]] <- if(missc) subsetfl(fl, cc) else fl
      if(missc) res[[3]] <- cc else if(variable.wise) res[[3]] <- "no missing cases in fl!" else {
          res[[3]] <- "no missing cases!"
          res[[5]] <- FALSE
      }
    } else {
      if(fcl) {
        res[[1]] <- fl[fc]
        res[[2]] <- do.call(cbind, c(list(Intercept = rep(1, length(fl[[1]]))), fl[!fc]))
      } else res[[1]] <- fl
    }
  } else { # fl is factor, vector or matrix !!
    res[[4]] <- ind
    if(na.rm) {
      cc <- if(variable.wise) which(complete.cases(fl)) else which(complete.cases(X, fl))
      if(length(cc) != nrow(X)) {
        if(is.factor(fl)) res[[1]] <- list(.Call(data.table:::CsubsetVector, fl, cc)) else res[[2]] <- cbind(Intercept = 1L, .Call(data.table:::CsubsetVector, fl, cc))
        res[[3]] <- cc
      } else {
        if(is.factor(fl)) res[[1]] <- list(fl) else res[[2]] <- cbind(Intercept = 1L, fl)
        if(variable.wise) res[[3]] <- "no missing cases in fl!" else {
          res[[3]] <- "no missing cases!"
          res[[5]] <- FALSE
        }
      }
    } else if(is.factor(fl)) res[[1]] <- list(fl) else res[[2]] <- cbind(Intercept = 1L, fl)
  }
  return(res)
}

# Note: could also do Mudlack and add means to second regression -> better than two-times centering ??
HDW <- function(X, ...) { # fl, w = NULL, na.rm = TRUE, fill = FALSE,
  UseMethod("HDW", X)
}
HDW.default <- function(X, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  ax <- attributes(X)
  if(na.rm) {
    cc <- which(complete.cases(X, fl)) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(length(cc) != length(X)) {
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
    if(na.rm) fl <- if(is.data.frame(fl)) .Call(data.table:::CsubsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    if(fcl && nallfc) {
      fl <- fl[fc]
      xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1]]))), fl[!fc]))
    } else if(!fcl) xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1]]))), fl))
    if(na.rm) {
      Y <- if(fcl && nallfc) qr.resid(qr(lfe::demeanlist(xmat, fl, weights = w[cc], ...)),
            lfe::demeanlist(X[cc], fl, weights = w[cc], ...)) else if(fcl)
            lfe::demeanlist(X[cc], fl, weights = w[cc], ...) else qr.resid(qr(xmat), X[cc])
      if(fill) {
        X[cc] <- Y
        X[-cc] <- NA
        attributes(X) <- ax
        return(X)
      } else return(`attributes<-`(Y, ax))
    } else {
      if(fcl && nallfc)
      return(`attributes<-`(qr.resid(qr(lfe::demeanlist(xmat, fl, weights = w, ...)),
      lfe::demeanlist(X, fl, weights = w, ...)), ax)) else if(fcl)
      return(`attributes<-`(lfe::demeanlist(X, fl, weights = w, ...), ax)) else
      return(`attributes<-`(qr.resid(qr(xmat), X), ax))
    }
  } else if (is.matrix(fl)) {
    xmat <- cbind(Intercept = 1L, fl)
    if(na.rm) {
      if(fill) {
        X[cc] <- qr.resid(qr(xmat[cc, , drop = FALSE]), X[cc])
        X[-cc] <- NA
        attributes(X) <- ax
        return(X)
      } else return(`attributes<-`(qr.resid(qr(xmat[cc, , drop = FALSE]), X[cc]), ax)) # best ??
    } else return(`attributes<-`(qr.resid(qr(xmat), X), ax))
  } else {
    if(na.rm) {
      Y <- if(is.factor(fl)) lfe::demeanlist(X[cc], list(fl[cc]), weights = w[cc], ...) else
            qr.resid(qr(cbind(1L,fl[cc])), X[cc])
      if(fill) {
        X[cc] <- Y
        X[-cc] <- NA
        attributes(X) <- ax
        return(X)
      } else return(`attributes<-`(Y, ax))
    } else {
      if(is.factor(fl)) return(`attributes<-`(lfe::demeanlist(X, list(fl), weights = w, ...), ax)) else
      return(`attributes<-`(qr.resid(qr(cbind(1L,fl)), X), ax))
    }
  }
}
HDW.pseries <- function(X, w = NULL, na.rm = TRUE, fill = TRUE, ...) {
  if(na.rm) {
    xcc <- which(!is.na(X))
    g <- lapply(attr(X, "index"), function(x) x[xcc])
    if(fill) {
      X[xcc] <- lfe::demeanlist(X[xcc], g, weights = w[xcc], ...) # keeps attributes ?? -> Yes !!
      return(X)
    } else return(`attr<-`(lfe::demeanlist(X[xcc], g, weights = w[xcc], ...), "index", g)) # keeps attributes ?? -> Nope !!
  } else {
    g <- attr(X, "index") # what about cases ?? -> nah, named !!
    return(`attr<-`(lfe::demeanlist(X, g, weights = w, ...), "index", g)) # keeps attributes ?? -> Nope !!
  }
}
HDW.matrix <- function(X, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  ax <- attributes(X)
  if(na.rm) {
    cc <- which(complete.cases(X, fl)) # gives error if lengths don't match, otherwise demeanlist and qr.resid give errors !!
    if(length(cc) != nrow(X)) {
      if(!fill) {
        ax[["dimnames"]][[1]] <- ax[["dimnames"]][[1]][cc] # best ??
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
    if(na.rm) fl <- if(is.data.frame(fl)) .Call(data.table:::CsubsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    if(fcl && nallfc) {
      fl <- fl[fc]
      xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1]]))), fl[!fc]))
    } else if(!fcl) xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1]]))), fl))
    if(na.rm) {
      Y <- if(fcl && nallfc) qr.resid(qr(lfe::demeanlist(xmat, fl, weights = w[cc], ...)),
            lfe::demeanlist(X[cc, , drop = FALSE], fl, weights = w[cc], ...)) else if(fcl)
            lfe::demeanlist(X[cc, , drop = FALSE], fl, weights = w[cc], ...) else qr.resid(qr(xmat), X[cc, , drop = FALSE])
      if(fill) {
        X[cc, ] <- Y
        X[-cc, ] <- NA
        attributes(X) <- ax
        return(X)
      } else return(`attributes<-`(Y, ax))
    } else {
      if(fcl && nallfc)
        return(`attributes<-`(qr.resid(qr(lfe::demeanlist(xmat, fl, weights = w, ...)),
        lfe::demeanlist(X, fl, weights = w, ...)), ax)) else if(fcl)
        return(`attributes<-`(lfe::demeanlist(X, fl, weights = w, ...), ax)) else
        return(`attributes<-`(qr.resid(qr(xmat), X), ax))
    }
  } else if (is.matrix(fl)) {
    xmat <- cbind(Intercept = 1L, fl)
    if(na.rm) {
      if(fill) {
        X[cc, ] <- qr.resid(qr(xmat[cc, , drop = FALSE]), X[cc, , drop = FALSE])
        X[-cc, ] <- NA
        attributes(X) <- ax
        return(X)
      } else return(`attributes<-`(qr.resid(qr(xmat[cc, , drop = FALSE]), X[cc, , drop = FALSE]), ax)) # best ??
    } else return(`attributes<-`(qr.resid(qr(xmat), X), ax))
  } else {
    if(na.rm) {
      Y <- if(is.factor(fl)) lfe::demeanlist(X[cc, , drop = FALSE], list(fl[cc]), weights = w[cc], ...) else
        qr.resid(qr(cbind(1L,fl[cc])), X[cc, , drop = FALSE])
      if(fill) {
        X[cc, ] <- Y
        X[-cc, ] <- NA
        attributes(X) <- ax
        return(X)
      } else return(`attributes<-`(Y, ax))
    } else {
      if(is.factor(fl)) return(`attributes<-`(lfe::demeanlist(X, list(fl), weights = w, ...), ax)) else
        return(`attributes<-`(qr.resid(qr(cbind(1L,fl)), X), ax))
    }
  }
}
HDW.data.frame <- function(X, fl, w = NULL, Xcols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
  list2env(getPartData(X, fl, Xcols, na.rm, variable.wise), envir = environment())

  if(variable.wise) fill <- TRUE
  # Best solution ????:
  if(na.rm && !fill) X <- .Call(data.table:::CsubsetDT, X, cc, Xvars) else if(na.rm && fill && !variable.wise)
    X <- .Call(data.table:::CsubsetDT, X, NULL, Xvars) else if(!identical(Xvars,seq_along(X))) X <- X[Xvars]
    Xvars <- seq_along(X)

    if(!is.null(dat)) {
      if(na.rm && fill) {
        if(variable.wise) {
          ax <- attributes(X)
          ax[["cases"]] <- cc
          lfl <- length(fl) > 0
          if(cc[1] == "no missing cases in fl!") {
            if(lfl) dat <- lfe::demeanlist(dat[,-1, drop = FALSE], fl, weights = w, ...)
            return(`attributes<-`(lapply(X, function(x) {
              xcc <- which(!is.na(x))
              if(lfl) x[xcc] <- qr.resid(qr(dat[xcc, , drop = FALSE]), lfe::demeanlist(x[xcc], subsetfl(fl, xcc), weights = w[xcc], ...)) else
                      x[xcc] <- qr.resid(qr(dat[xcc, , drop = FALSE]), x[xcc])
              return(x)
            }), ax)) # Rfast fastlm??
          } else {
            if(lfl) dat <- lfe::demeanlist(dat[,-1, drop = FALSE], fl, weights = w[cc], ...)
            return(`attributes<-`(lapply(X, function(x) {
              xcc <- !is.na(x)
              xcc[-cc] <- FALSE
              x[-cc] <- NA # which is not faster !!
              if(lfl) {
                XC <- which(xcc[cc])
                x[xcc] <- qr.resid(qr(dat[XC, , drop = FALSE]), lfe::demeanlist(x[xcc], subsetfl(fl, XC), weights = w[xcc], ...))
              } else x[xcc] <- qr.resid(qr(dat[xcc[cc], , drop = FALSE]), x[xcc])
              return(x)
            }), ax)) # Rfast fastlm??
          }
        } else {
          attr(X, "cases") <- cc
          if(length(fl)) { # fastest so far !! (and memory efficient)
            dat <- qr(lfe::demeanlist(dat[,-1], fl, weights = w[cc], ...)) # This is necessary!!!
            Y <- lapply(lfe::demeanlist(.Call(data.table:::CsubsetDT, X, cc, Xvars), fl, weights = w[cc], ...), function(x) qr.resid(dat, x))
            return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, X, cc, Xvars, NULL, Y, FALSE),
                         seq_row(X)[-cc], Xvars, NULL, NA, FALSE))
          } else {
            dat <- qr(dat)
            Y <- lapply(.Call(data.table:::CsubsetDT, X, cc, Xvars), function(x) qr.resid(dat,x))
            return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, X, cc, Xvars, NULL, Y, FALSE),
                         seq_row(X)[-cc], Xvars, NULL, NA, FALSE))
          }
        }
      } else {
        ax <- attributes(X)
        if(na.rm) ax[["cases"]] <- cc
        if(length(fl)) {
          if(!is.null(w) && na.rm) w <- w[cc]
          dat <- qr(lfe::demeanlist(dat[,-1], fl, weights = w, ...)) # This is necessary!!!
          return(`attributes<-`(lapply(lfe::demeanlist(X, fl, weights = w, ...),
                                       function(x) qr.resid(dat, x)), ax))
        } else {
          dat <- qr(dat)
          return(`attributes<-`(lapply(X, function(x) qr.resid(dat, x)), ax))
        }
      }
    } else {
      if(na.rm && fill) {
        if(variable.wise) {
          ax <- attributes(X)
          ax[["cases"]] <- cc
          if(cc[1] == "no missing cases in fl!") {
            return(`attributes<-`(lapply(X, function(x) {
              xcc <- which(!is.na(x))
              x[xcc] <- lfe::demeanlist(x[xcc], subsetfl(fl, xcc), weights = w[xcc], ...)
              return(x)
            }), ax))
          } else {
            return(`attributes<-`(lapply(X, function(x) {
              xcc <- !is.na(x)
              xcc[-cc] <- FALSE
              x[-cc] <- NA # which is not faster !!
              x[xcc] <- lfe::demeanlist(x[xcc], subsetfl(fl, which(xcc[cc])), weights = w[xcc], ...)
              return(x)
            }), ax))
          }
        } else {
          attr(X, "cases") <- cc
          Y <- lfe::demeanlist(.Call(data.table:::CsubsetDT, X, cc, Xvars), fl, weights = w[cc], ...)
          return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, X, cc, Xvars, NULL, Y, FALSE),
                       seq_row(X)[-cc], Xvars, NULL, NA, FALSE))
        }
      } else {
        if(na.rm) return(`attr<-`(lfe::demeanlist(X, fl, weights = w, ...), "cases", cc)) else
        return(lfe::demeanlist(X, fl, weights = w, ...))
      }
    }
}
HDW.pdata.frame <- function(X, w = NULL, Xcols = is.numeric, na.rm = TRUE, fill = TRUE,
                            variable.wise = FALSE, drop.xt = FALSE, drop.w = TRUE, ...) {
  ax <- attributes(X)
  nam <- ax[["names"]]
  gn <- match(names(ax[["index"]]), nam)
  gn2 = gn <- gn[!is.na(gn)]
  if(!is.null(Xcols)) {
    if(is.function(Xcols)) Xcols <- seq_along(X)[!vapply(X, Xcols, TRUE)] else if(is.character(Xcols))
      Xcols <- seq_along(X)[-match(Xcols, nam)] else if(is.logical(Xcols))
        Xcols <- seq_along(X)[!Xcols] else if(is.numeric(Xcols))
          Xcols <- seq_along(X)[-Xcols] else stop("Xcols needs to be a function, column names, indices or a logical vector")
        if(drop.xt) gn <- unique.default(c(gn, Xcols)) else if(length(gn)) gn2 <- unique.default(c(gn2, Xcols)) else gn <- Xcols
  }
  if(!is.null(w)) {
    if(is.call(w)) w <- all.vars(w)
    if(length(w) == 1) {
      wn <- match(w, nam)
      if(any(gn == wn)) stop("Weights coincide with grouping variables!")
      w <- X[[wn]]
      if(drop.w) if(!length(gn)) X[[wn]] <- NULL else if(drop.xt) gn <- c(gn,wn) else gn2 <- c(gn2,wn)
    }
  }
  varwisecomp <- function(X, fl, w, ...) lapply(X, function(x) {
    xcc <- which(!is.na(x))
    x[xcc] <- lfe::demeanlist(x[xcc], subsetfl(fl, xcc), weights = w[xcc], ...)
    return(x)
  })
  fillcomp <- function(X, fl, w, cc, Xvars, ...) {
    if(length(cc) != length(X[[1]])) {
      Y <- lfe::demeanlist(.Call(data.table:::CsubsetDT, X, cc, Xvars), lapply(fl, function(x) x[cc]), weights = w[cc], ...)
      return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, X, cc, Xvars, NULL, Y, FALSE),
                   seq_row(X)[-cc], Xvars, NULL, NA, FALSE))
    } else return(lfe::demeanlist(X[Xvars], fl, weights = w, ...))
  }
  nonfillcomp <- function(X, fl, w, cc, Xvars, ...) {
    if(length(cc) != length(X[[1]]))
      return(lfe::demeanlist(.Call(data.table:::CsubsetDT, X, cc, Xvars), lapply(fl, function(x) x[cc]), weights = w[cc], ...)) else
      return(lfe::demeanlist(X[Xvars], fl, weights = w, ...))
  }
  if(na.rm && fill) {
    if(variable.wise) {
      attributes(X) <- NULL
      if(length(gn)) {
        if(drop.xt) {
          ax[["names"]] <- if(give.names) paste0("HDW.",nam[-gn]) else nam[-gn]
          return(`attributes<-`(varwisecomp(X[-gn], ax[["index"]], w, ...), ax))
        } else {
            ax[["names"]] <- c(nam[gn], if(give.names) paste0("HDW.",nam[-gn2]) else nam[-gn2])
            return(`attributes<-`(c(X[gn],varwisecomp(X[-gn2], ax[["index"]], w, ...)), ax))
        }
      } else {
        if(give.names) paste0("HDW.",nam)
        return(`attributes<-`(varwisecomp(X, ax[["index"]], w, ...), ax))
      }
    } else {
      cc <- which(!.Call(data.table:::Cdt_na, X, seq_along(X)[-c(gn,gn2)])) # good??
      ax[["cases"]] <- if(length(cc) == length(X[[1]])) "no missing cases!" else cc
      if(length(gn)) {
        if(drop.xt) {
          ax[["names"]] <- if(give.names) paste0("HDW.",nam[-gn]) else nam[-gn]
          return(`attributes<-`(fillcomp(X, ax[["index"]], w, cc, seq_along(X)[-gn], ...), ax))
        } else {
          ax[["names"]] <- c(nam[gn], if(give.names) paste0("HDW.",nam[-gn2]) else nam[-gn2])
          return(`attributes<-`(c(X[gn], fillcomp(X, ax[["index"]], w, cc, seq_along(X)[-gn2], ...)), ax))
        }
      } else {
        if(give.names) paste0("HDW.",nam)
        return(`attributes<-`(fillcomp(X, ax[["index"]], w, cc, seq_along(X), ...), ax))
      }
    }
  } else {
    if(na.rm) {
      cc <- which(!.Call(data.table:::Cdt_na, X, seq_along(X)[-c(gn,gn2)])) # good??
      if(length(cc) == length(X[[1]])) ax[["cases"]] <- "no missing cases!" else {
        ax[["row.names"]] <- ax[["row.names"]][cc]
        ax[["cases"]] <- cc
      }
      if(length(gn)) {
        if(drop.xt) {
          ax[["names"]] <- if(give.names) paste0("HDW.",nam[-gn]) else nam[-gn]
          return(`attributes<-`(nonfillcomp(X, ax[["index"]], w, cc, seq_along(X)[-gn], ...), ax))
        } else {
          ax[["names"]] <- c(nam[gn], if(give.names) paste0("HDW.",nam[-gn2]) else nam[-gn2])
          return(`attributes<-`(c(.Call(data.table:::CsubsetDT, X, cc, gn), nonfillcomp(X, ax[["index"]], w, cc, seq_along(X)[-gn2], ...)), ax))
        }
      } else {
        if(give.names) paste0("HDW.",nam)
        return(`attributes<-`(nonfillcomp(X, ax[["index"]], w, cc, seq_along(X), ...), ax))
      }
    } else {
      attributes(X) <- NULL
      if(length(gn)) {
        if(drop.xt) {
          ax[["names"]] <- if(give.names) paste0("HDW.",nam[-gn]) else nam[-gn]
          return(`attributes<-`(lfe::demeanlist(X[-gn], ax[["index"]], weights = w, ...), ax))
        } else {
          ax[["names"]] <- c(nam[gn], if(give.names) paste0("HDW.",nam[-gn2]) else nam[-gn2])
          return(`attributes<-`(c(X[gn],lfe::demeanlist(X[-gn2], ax[["index"]], weights = w, ...)), ax))
        }
      } else {
        if(give.names) paste0("HDW.",nam)
        return(`attributes<-`(lfe::demeanlist(X, ax[["index"]], weights = w, ...), ax))
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

HDB <- function(X, ...) { # fl, w = NULL, na.rm = TRUE, fill = FALSE,
  UseMethod("HDB", X)
}
HDB.default <- function(X, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  ax <- attributes(X)
  if(na.rm) {
    cc <- which(complete.cases(X, fl)) # gives error if lengths don't match, otherwise demeanlist and qr.fitted give errors !!
    if(length(cc) != length(X)) {
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
    if(na.rm) fl <- if(is.data.frame(fl)) .Call(data.table:::CsubsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    if(fcl && nallfc) {
      fl <- fl[fc]
      xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1]]))), fl[!fc]))
    } else if(!fcl) xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1]]))), fl))
    if(na.rm) {
      Y <- if(fcl && nallfc) X[cc] - qr.resid(qr(lfe::demeanlist(xmat, fl, weights = w[cc], ...)),
                                          lfe::demeanlist(X[cc], fl, weights = w[cc], ...)) else if(fcl)
                                          lfe::demeanlist(X[cc], fl, weights = w[cc], means = TRUE, ...) else qr.fitted(qr(xmat), X[cc])
      if(fill) {
        X[cc] <- Y
        X[-cc] <- NA
        attributes(X) <- ax
        return(X)
      } else return(`attributes<-`(Y, ax))
    } else {
      if(fcl && nallfc)
        return(`attributes<-`(X - qr.resid(qr(lfe::demeanlist(xmat, fl, weights = w, ...)),
                                       lfe::demeanlist(X, fl, weights = w, ...)), ax)) else if(fcl)
                                         return(`attributes<-`(lfe::demeanlist(X, fl, weights = w, means = TRUE, ...), ax)) else
                                           return(`attributes<-`(qr.fitted(qr(xmat), X), ax))
    }
  } else if (is.matrix(fl)) {
    xmat <- cbind(Intercept = 1L, fl)
    if(na.rm) {
      if(fill) {
        X[cc] <- qr.fitted(qr(xmat[cc, , drop = FALSE]), X[cc])
        X[-cc] <- NA
        attributes(X) <- ax
        return(X)
      } else return(`attributes<-`(qr.fitted(qr(xmat[cc, , drop = FALSE]), X[cc]), ax)) # best ??
    } else return(`attributes<-`(qr.fitted(qr(xmat), X), ax))
  } else {
    if(na.rm) {
      Y <- if(is.factor(fl)) lfe::demeanlist(X[cc], list(fl[cc]), weights = w[cc], means = TRUE, ...) else
        qr.fitted(qr(cbind(1L,fl[cc])), X[cc])
      if(fill) {
        X[cc] <- Y
        X[-cc] <- NA
        attributes(X) <- ax
        return(X)
      } else return(`attributes<-`(Y, ax))
    } else {
      if(is.factor(fl)) return(`attributes<-`(lfe::demeanlist(X, list(fl), weights = w, means = TRUE, ...), ax)) else
        return(`attributes<-`(qr.fitted(qr(cbind(1L,fl)), X), ax))
    }
  }
}
HDB.matrix <- function(X, fl, w = NULL, na.rm = TRUE, fill = FALSE, ...) {
  ax <- attributes(X)
  if(na.rm) {
    cc <- which(complete.cases(X, fl)) # gives error if lengths don't match, otherwise demeanlist and qr.fitted give errors !!
    if(length(cc) != nrow(X)) {
      if(!fill) {
        ax[["dimnames"]][[1]] <- ax[["dimnames"]][[1]][cc] # best ??
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
    if(na.rm) fl <- if(is.data.frame(fl)) .Call(data.table:::CsubsetDT, fl, cc, seq_along(fl)) else subsetfl(fl, cc)
    if(fcl && nallfc) {
      fl <- fl[fc]
      xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1]]))), fl[!fc]))
    } else if(!fcl) xmat <- do.call(cbind, c(list(Intercept = rep(1L, length(fl[[1]]))), fl))
    if(na.rm) {
      Y <- if(fcl && nallfc) X[cc, , drop = FALSE] - qr.resid(qr(lfe::demeanlist(xmat, fl, weights = w[cc], ...)),
                                        lfe::demeanlist(X[cc, , drop = FALSE], fl, weights = w[cc], ...)) else if(fcl)
                                        lfe::demeanlist(X[cc, , drop = FALSE], fl, weights = w[cc], means = TRUE, ...) else qr.fitted(qr(xmat), X[cc, , drop = FALSE])
      if(fill) {
        X[cc, ] <- Y
        X[-cc, ] <- NA
        attributes(X) <- ax
        return(X)
      } else return(`attributes<-`(Y, ax))
    } else {
      if(fcl && nallfc)
        return(`attributes<-`(X - qr.resid(qr(lfe::demeanlist(xmat, fl, weights = w, ...)),
                                       lfe::demeanlist(X, fl, weights = w, ...)), ax)) else if(fcl)
                                         return(`attributes<-`(lfe::demeanlist(X, fl, weights = w, means = TRUE, ...), ax)) else
                                           return(`attributes<-`(qr.fitted(qr(xmat), X), ax))
    }
  } else if (is.matrix(fl)) {
    xmat <- cbind(Intercept = 1L, fl)
    if(na.rm) {
      if(fill) {
        X[cc, ] <- qr.fitted(qr(xmat[cc, , drop = FALSE]), X[cc, , drop = FALSE])
        X[-cc, ] <- NA
        attributes(X) <- ax
        return(X)
      } else return(`attributes<-`(qr.fitted(qr(xmat[cc, , drop = FALSE]), X[cc, , drop = FALSE]), ax)) # best ??
    } else return(`attributes<-`(qr.fitted(qr(xmat), X), ax))
  } else {
    if(na.rm) {
      Y <- if(is.factor(fl)) lfe::demeanlist(X[cc, , drop = FALSE], list(fl[cc]), weights = w[cc], means = TRUE, ...) else
           qr.fitted(qr(cbind(1L,fl[cc])), X[cc, , drop = FALSE])
      if(fill) {
        X[cc, ] <- Y
        X[-cc, ] <- NA
        attributes(X) <- ax
        return(X)
      } else return(`attributes<-`(Y, ax))
    } else {
      if(is.factor(fl)) return(`attributes<-`(lfe::demeanlist(X, list(fl), weights = w, means = TRUE, ...), ax)) else
        return(`attributes<-`(qr.fitted(qr(cbind(1L,fl)), X), ax))
    }
  }
}
HDB.data.frame <- function(X, fl, w = NULL, Xcols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
  list2env(getPartData(X, fl, Xcols, na.rm, variable.wise), envir = environment())

  if(variable.wise) fill <- TRUE
  # Best solution ????:
  if(na.rm && !fill) X <- .Call(data.table:::CsubsetDT, X, cc, Xvars) else if(na.rm && fill && !variable.wise)
    X <- .Call(data.table:::CsubsetDT, X, NULL, Xvars) else if(!identical(Xvars,seq_along(X))) X <- X[Xvars]
    Xvars <- seq_along(X)

    if(!is.null(dat)) {
      if(na.rm && fill) {
        if(variable.wise) {
          ax <- attributes(X)
          ax[["cases"]] <- cc
          lfl <- length(fl) > 0
          if(cc[1] == "no missing cases in fl!") {
            if(lfl) dat <- lfe::demeanlist(dat[,-1, drop = FALSE], fl, weights = w, ...)
            return(`attributes<-`(lapply(X, function(x) {
              xcc <- which(!is.na(x))
              if(lfl) x[xcc] <- x[xcc] - qr.resid(qr(dat[xcc, , drop = FALSE]), lfe::demeanlist(x[xcc], subsetfl(fl, xcc), weights = w[xcc], ...)) else
                x[xcc] <- qr.fitted(qr(dat[xcc, , drop = FALSE]), x[xcc])
              return(x)
            }), ax)) # Rfast fastlm??
          } else {
            if(lfl) dat <- lfe::demeanlist(dat[,-1, drop = FALSE], fl, weights = w[cc], ...)
            return(`attributes<-`(lapply(X, function(x) { # change !!
              xcc <- !is.na(x)
              xcc[-cc] <- FALSE
              x[-cc] <- NA # which is not faster !!
              if(lfl) {
                XC <- which(xcc[cc])
                x[xcc] <- x[xcc] - qr.resid(qr(dat[XC, , drop = FALSE]), lfe::demeanlist(x[xcc], subsetfl(fl, XC), weights = w[xcc], ...))
              } else x[xcc] <- qr.fitted(qr(dat[xcc[cc], , drop = FALSE]), x[xcc])
              return(x)
            }), ax)) # Rfast fastlm??
          }
        } else {
          attr(X, "cases") <- cc
          if(length(fl)) { # fastest so far !! (and memory efficient)
            dat <- qr(lfe::demeanlist(dat[,-1], fl, weights = w[cc], ...)) # This is necessary!!!
            Y <- lapply(.Call(data.table:::CsubsetDT, X, cc, Xvars), function(x) x - qr.resid(dat, lfe::demeanlist(x, fl, weights = w[cc], ...))) # good ??
            return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, X, cc, Xvars, NULL, Y, FALSE),
                         seq_row(X)[-cc], Xvars, NULL, NA, FALSE))
          } else {
            dat <- qr(dat)
            Y <- lapply(.Call(data.table:::CsubsetDT, X, cc, Xvars), function(x) qr.fitted(dat,x))
            return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, X, cc, Xvars, NULL, Y, FALSE),
                         seq_row(X)[-cc], Xvars, NULL, NA, FALSE))
          }
        }
      } else {
        ax <- attributes(X)
        if(na.rm) ax[["cases"]] <- cc
        if(length(fl)) {
          if(!is.null(w) && na.rm) w <- w[cc]
          dat <- qr(lfe::demeanlist(dat[,-1], fl, weights = w, ...)) # This is necessary!!!
          return(`attributes<-`(Map(function(x, y) x - qr.resid(dat, y), X, lfe::demeanlist(X, fl, weights = w, ...)), ax))
        } else {
          dat <- qr(dat)
          return(`attributes<-`(lapply(X, function(x) qr.fitted(dat, x)), ax))
        }
      }
    } else {
      if(na.rm && fill) {
        if(variable.wise) {
          ax <- attributes(X)
          ax[["cases"]] <- cc
          if(cc[1] == "no missing cases in fl!") {
            return(`attributes<-`(lapply(X, function(x) {
              xcc <- which(!is.na(x))
              x[xcc] <- lfe::demeanlist(x[xcc], subsetfl(fl, xcc), weights = w[xcc], means = TRUE, ...)
              return(x)
            }), ax))
          } else {
            return(`attributes<-`(lapply(X, function(x) {
              xcc <- !is.na(x)
              xcc[-cc] <- FALSE
              x[-cc] <- NA # which is not faster !!
              x[xcc] <- lfe::demeanlist(x[xcc], subsetfl(fl, which(xcc[cc])), weights = w[xcc], means = TRUE, ...)
              return(x)
            }), ax))
          }
        } else {
          attr(X, "cases") <- cc
          Y <- lfe::demeanlist(.Call(data.table:::CsubsetDT, X, cc, Xvars), fl, weights = w[cc], means = TRUE, ...)
          return(.Call(data.table:::Cassign, .Call(data.table:::Cassign, X, cc, Xvars, NULL, Y, FALSE),
                       seq_row(X)[-cc], Xvars, NULL, NA, FALSE))
        }
      } else {
        if(na.rm) return(`attr<-`(lfe::demeanlist(X, fl, weights = w, means = TRUE, ...), "cases", cc)) else
          return(lfe::demeanlist(X, fl, weights = w, means = TRUE, ...))
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
# HDW.data.frame.old <- function(X, fl, w = NULL, Xcols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
#   list2env(getPartData(X, fl, Xcols, na.rm, variable.wise), envir = environment())
#   ay <- attributes(Y)
#   if(na.rm) ay[["cases"]] <- cc # good !!
#   if(!is.null(dat)) {
#     if(na.rm && fill) {
#       if(variable.wise) {
#         lfl <- length(fl) > 0
#         if(lfl) dat <- lfe::demeanlist(dat[,-1], fl, weights = w[cc], ...)
#         return(`attributes<-`(lapply(Y, function(x) {
#           xcc <- !is.na(x)
#           xcc[-cc] <- FALSE
#           x[-cc] <- NA # which is not faster !!
#           if(lfl) {
#             XC <- which(xcc[cc])
#             x[xcc] <- qr.resid(qr(dat[XC, , drop = FALSE]), lfe::demeanlist(x[xcc], subsetfl(fl, XC), weights = w[xcc], ...))
#           } else x[xcc] <- qr.resid(qr(dat[xcc[cc], , drop = FALSE]), x[xcc])
#           return(x)
#         }), ay)) # Rfast fastlm??
#       } else {
#         if(!is.null(fl)) { # fastest so far !! (and memory efficient)
#           dat <- qr(lfe::demeanlist(dat[,-1], fl, weights = w[cc], ...)) # This is necessary!!!
#           Y <- lapply(lfe::demeanlist(.Call(CSDT, X, cc, Xvars), fl, weights = w[cc], ...), function(x) qr.resid(dat,x))
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
#         dat <- qr(lfe::demeanlist(dat[,-1], fl, weights = w, ...)) # This is necessary!!!
#         if(na.rm) X <- .Call(CSDT, X, cc, Xvars) else if(!identical(Xvars,seq_along(X))) X <- X[Xvars]
#         return(`attributes<-`(lapply(lfe::demeanlist(X, fl, weights = w, ...),
#                                      function(x) qr.resid(dat, x)), ay))
#       } else {
#         dat <- qr(dat)
#         if(na.rm) X <- .Call(CSDT, X, cc, Xvars) else if(!identical(Xvars,seq_along(X))) X <- X[Xvars]
#         return(`attributes<-`(lapply(X, function(x) qr.resid(dat, x)), ay))
#       }
#     }
#   } else {
#     if(na.rm && fill) {
#       if(variable.wise) {
#         return(`attributes<-`(lapply(Y, function(x) {
#           xcc <- !is.na(x)
#           xcc[-cc] <- FALSE
#           x[-cc] <- NA # which is not faster !!
#           x[xcc] <- lfe::demeanlist(x[xcc], subsetfl(fl, which(xcc[cc])), weights = w[xcc], ...)
#           return(x)
#         }), ay))
#       } else return(.Call(AS, .Call(AS, .Call(CSDT, X, NULL, Xvars), cc, seq_along(Xvars), NULL,
#                                     lfe::demeanlist(Y, fl, weights = w[cc], ...), FALSE),
#                           seq_len(nrow(X))[-cc], seq_along(Xvars),
#                           NULL, NA, FALSE))
#
#     } else return(lfe::demeanlist(Y, fl, weights = w, ...))
#   }
# }


# HDW(X = na.omit(GGDC),fl = ~factor(Country))

# Previous Getpartdata -> generating Y before passing to HDW.data.frame

# getPartData <- function(X, fl, Xcols, na.rm, variable.wise) {
#   res <- list(fl = NULL, dat = NULL, Y = NULL, cc = NULL, Xvars = NULL)
#   nam <- names(X)
#   form <- class(fl) == "formula"
#   twosided <- form && length(fl) == 3
#   Xcolsl <- !is.null(Xcols)
#   if(!twosided && Xcolsl) {
#     if(is.function(Xcols))
#       ind <- which(vapply(X, Xcols, TRUE)) else if(is.numeric(Xcols))
#         ind <- Xcols else if(is.character(Xcols))
#           ind <- match(Xcols, nam) else if(is.logical(Xcols))
#             ind <- which(Xcols) else
#               stop("Xcols needs to be a function, column names, indices or a logical vector")
#   } else ind <- seq_along(nam)
#
#   if(form) { # fl is formula
#     if(twosided) {
#       fvars <- match(all.vars(fl[[3]]), nam)
#       Xvars <- match(all.vars(fl[[2]]), nam)
#       fl[[2]] <- NULL
#     } else {
#       fvars <- match(all.vars(fl), nam)
#       Xvars <- setdiff(ind, fvars)
#     }
#     res[[5]] <- Xvars
#     if(na.rm) {
#       if(variable.wise) {
#         cc <- which(!.Call(NADT, X, fvars)) # best ??
#         res[1:2] <- getfl(myModFrame(fl, .Call(CSDT, X, cc, fvars)))
#         res[[3]] <- .Call(CSDT, X, NULL, Xvars)
#       } else {
#         cc <- which(!.Call(NADT, X, c(Xvars,fvars))) # best ??
#         res[1:2] <- getfl(myModFrame(fl, .Call(CSDT, X, cc, fvars)))
#         res[[3]] <- .Call(CSDT, X, cc, Xvars)
#       }
#       res[[4]] <- cc
#     } else {
#       res[1:2] <- getfl(myModFrame(fl, X))
#       res[[3]] <- .Call(CSDT, X, NULL, Xvars)
#     }
#   } else if(is.list(fl)) { # fl is factor list
#     class(fl) <- NULL
#     res[[5]] <- ind
#     fc <- vapply(fl, is.factor, TRUE)
#     fcl <- any(!fc)
#     if(na.rm) {
#       if(variable.wise) {
#         cc <- which(!.Call(NADT, fl, seq_along(fl))) # best ??
#         if(fcl) {
#           fl <- subsetfl(fl, cc)
#           res[[1]] <- fl[fc]
#           res[[2]] <- do.call(cbind, c(Intercept = rep(1, length(fl[[1]])), fl[!fc]))
#         } else res[[1]] <- subsetfl(fl, cc)
#         res[[3]] <- if(Xcolsl) .Call(CSDT, X, NULL, ind) else X
#       } else {
#         cc <- which(!.Call(NADT, c(X,fl), seq_len(length(X)+length(fl)))) # best ??
#         if(fcl) {
#           fl <- subsetfl(fl, cc)
#           res[[1]] <- fl[fc]
#           res[[2]] <- do.call(cbind, c(Intercept = rep(1, length(fl[[1]])), fl[!fc]))
#         } else res[[1]] <- subsetfl(fl, cc)
#         res[[3]] <- .Call(CSDT, X, cc, ind)
#       }
#       res[[4]] <- cc
#     } else {
#       if(fcl) {
#         res[[1]] <- fl[fc]
#         res[[2]] <- do.call(cbind, c(Intercept = rep(1, length(fl[[1]])), fl[!fc]))
#       } else res[[1]] <- fl
#       res[[3]] <- if(Xcolsl) .Call(CSDT, X, NULL, ind) else X
#     }
#   } else { # fl is vector or matrix !!
#     res[[5]] <- ind
#     if(is.factor(fl)) stop("For single factors please use W(), not HDW(), or pass the factor inside a list.")
#     xmat <- cbind(Intercept = 1L, fl)
#     if(na.rm) {
#       if(variable.wise) {
#         cc <- which(complete.cases(xmat))
#         res[[2]] <- xmat[cc, , drop = FALSE]
#         res[[3]] <- if(Xcolsl) .Call(CSDT, X, NULL, ind) else X
#       } else {
#         cc <- which(complete.cases(X, xmat))
#         res[[2]] <- xmat[cc, , drop = FALSE]
#         res[[3]] <- .Call(CSDT, X, cc, ind)
#       }
#       res[[4]] <- cc
#     } else {
#       res[[2]] <- xmat
#       res[[3]] <- if(Xcolsl) .Call(CSDT, X, NULL, ind) else X
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


# mydml without na.rm -> same speed as lfe::demeanlist
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
#       N <- length(fl[[1]])
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

# HDW.data.frame <- function(X, fl, w = NULL, Xcols = is.numeric, na.rm = TRUE, fill = FALSE, variable.wise = FALSE, ...) {
#   list2env(getPartData(X, fl, Xcols, na.rm, variable.wise), envir = environment())
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
#               return(`attributes<-`(lapply(res, function(x) qr.resid(qrdat,x)), ar)) # Rfast fastlm??
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

# getPartData <- function(X, fl, Xcols) {
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
#       if(!is.null(Xcols)) {
#         if(is.function(Xcols) && !ar)
#           ind <- which(vapply(X, Xcols, TRUE)) else if(is.numeric(Xcols))
#             ind <- Xcols else if(is.character(Xcols))
#               ind <- match(Xcols,nam) else if(is.logical(Xcols))
#                 ind <- which(Xcols) else
#                   stop("Xcols needs to be a function, column names, indices or a logical vector")
#       } else ind <- seq_along(nam)
#       fvars <- all.vars(fl)
#       res[1:2] <- getfl(myModFrame(fl, X))
#       if (!is.null(nam)) {
#         ind <- setdiff(ind, match(fvars, nam))
#         res[[3]] <- if(ar) X[,ind] else X[ind]
#       } else res[[3]] <- X
#     } else { # two-sided formula
#       fvars <- all.vars(fl[[3]])
#       Xvars <- all.vars(fl[[2]])
#       fl[[2]] <- NULL
#       res[1:2] <- getfl(myModFrame(fl, X))
#       if(!is.null(nam)) {
#         res[[3]] <- if(ar) X[,Xvars] else X[Xvars]
#       } else res[[3]] <- X
#     }
#   } else { # fl is factor list !!
#     if(!is.null(Xcols)) {
#       if(is.function(Xcols) && !ar)
#         ind <- which(vapply(X, Xcols, TRUE)) else if(is.numeric(Xcols))
#           ind <- Xcols else if(is.character(Xcols))
#             ind <- match(Xcols, nam) else if (is.logical(Xcols))
#               ind <- which(Xcols) else
#                 stop("Xcols needs to be a function, column names, indices or a logical vector")
#             res[[3]] <- if(ar) X[,ind] else X[ind]
#     } else res[[3]] <- X
#     if(is.list(fl)) {
#       fc <- vapply(fl, is.factor, TRUE)
#       res[[1]] <- c(Intercept = rep(1, length(fl[[1]])), fl[!fc])
#       res[[2]] <- fl[fc]
#     } else {
#       res[[1]] <- c(Intercept = rep(1, length(fl[[1]])), fl) # needed??
#       res[[3]] <- NULL
#     }
#   } # still do else list of factors or data.frame or model matrix!!
#   return(res)
# }
