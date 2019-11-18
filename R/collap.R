# Helper functions:
# Crbindlist <- data.table:::Crbindlist
# Csetcolorder <- data.table:::Csetcolorder

# global macro
.FAST_STAT_FUN <- c("fmean","fmedian","fmode","fsum","fprod","fsd","fvar",
                    "fmin","fmax","ffirst","flast","fNobs","fNunique",
                    "flag","flead","fdiff","fgrowth","fscale")

# todo speed up: use internals of order ??
collap <- function(X, by, FUN = fmean, catFUN = fmode, cols = NULL, custom = NULL,
                   keep.by = TRUE, sort.col = TRUE, sort.row = TRUE, parallel = FALSE, mc.cores = 1L,
                   multi.FUN.out = c("wide","list","long","long_dupl"), give.names = "auto", ...) {

  multi.FUN.out <- switch(multi.FUN.out[1L], wide = 1L, list = 2L, long = 3L, long_dupl = 4L, stop("Unknown multi.FUN output option"))
  widel <- multi.FUN.out == 1L
  customl <- !is.null(custom)
  if(!is.data.frame(X)) X <- qDF(X)
  ax <- attributes(X)
  nam <- ax[["names"]]
  # attributes(X) <- NULL
  class(X) <- NULL
  # attr(X, "class") <- "data.frame" # class needed for method dispatch of fast functions, not for BY !!

  aplyfun <- if(parallel) function(...) parallel::mclapply(..., mc.cores = mc.cores) else lapply

  # identifying by and cols
  vl <- TRUE
  if(is.call(by)) {
      if(length(by) == 3L) {
        v <- nam %in% all.vars(by[[2L]])
        namby <- all.vars(by[[3L]])
        numby <- anyNAerror(match(namby, nam), "Unknown 'by' columns selected!")
      } else {
        namby <- all.vars(by)
        numby <- anyNAerror(match(namby, nam), "Unknown 'by' columns selected!")
        if(!customl) {
          v <- if(is.null(cols)) !logical(length(X)) else cols2log(X, nam, cols)
          v[numby] <- FALSE
        }
      }
      if(length(numby) == 1L) by <- qF(X[[numby]], ordered = sort.row) else
                              by <- GRP(X, numby, sort = sort.row, return.groups = keep.by)
  } else if(is.atomic(by)) {
    namby <- deparse(substitute(by))
    numby <- 1L
    if(!customl) if(is.null(cols)) vl <- FALSE else v <- cols2log(X, nam, cols)
    if(!is.factor(by)) by <- qF(by, ordered = sort.row)
  } else {
    if(!customl) if(is.null(cols)) vl <- FALSE else v <- cols2log(X, nam, cols)
    if(!is.GRP(by)) {
      numby <- seq_along(by)
      namby <- names(by)
      if(is.null(namby)) namby <- paste0("GRP.", numby)
      by <- GRP(by, numby, sort = sort.row, return.groups = keep.by)
    } else {
      namby <- by[[5L]]
      if(is.null(namby)) namby <- paste0("GRP.", 1L:length(by[[4L]])) # necessary ??
      numby <- seq_along(namby)
    }
  }

  if(!customl) {

    # identifying data
    nu <- vapply(X, is.numeric, TRUE, USE.NAMES = FALSE)
    if(vl) {
      nnu <- which(!nu & v) # faster way ??
      nu <- which(nu & v)
    } else {
      nnu <- which(!nu)
      nu <- which(nu)
    }
    nul <- length(nu) > 0L
    nnul <- length(nnu) > 0L

    # Identifying FUN and catFUN:
    if(nul) if(is.character(FUN)) {
      # FUN <- unlist(strsplit(FUN,",",fixed = TRUE), use.names = FALSE)
      namFUN <- FUN
      FUN <- if(length(FUN) > 1L) lapply(FUN, match.fun, descend = FALSE) else
                                  match.fun(FUN, descend = FALSE)
    } else if(is.list(FUN)) {
      namFUN <- names(FUN)
      if(is.null(namFUN)) namFUN <- all.vars(substitute(FUN))
    } else namFUN <- deparse(substitute(FUN))

    if(nnul) if(is.character(catFUN)) {
      # catFUN <- unlist(strsplit(catFUN,",",fixed = TRUE), use.names = FALSE)
      namcatFUN <- catFUN
      catFUN <- if(length(catFUN) > 1L) lapply(catFUN, match.fun, descend = FALSE) else
                                        match.fun(catFUN, descend = FALSE)
    } else if(is.list(catFUN)) {
      namcatFUN <- names(catFUN)
      if(is.null(namcatFUN)) namcatFUN <- all.vars(substitute(catFUN))
    } else namcatFUN <- deparse(substitute(catFUN))

    if(give.names == "auto") give.names <- !widel || length(FUN) > 1L || length(catFUN) > 1L

    # Aggregator function ! # drop level of nesting i.e. make rest length(by)+length(FUN)+length(catFUN)  ??
    agg <- function(xnu, xnnu) { #by, FUN, namFUN, catFUN, namcatFUN, drop.by
      lr <- nul + nnul + keep.by
      res <- vector("list", lr)
        if(keep.by) {
          res[[1L]] <- if(is.atomic(by)) list(`names<-`(list(attr(by, "levels")), namby)) else list(by[[4L]]) # nah... could add later using "c" ??
          ind <- 2L
        } else ind <- 1L
        if(nul) {
          fFUN <- namFUN %in% .FAST_STAT_FUN
          if(is.list(FUN))
            res[[ind]] <- condsetn(aplyfun(seq_along(namFUN), function(i)
                          if(fFUN[i]) FUN[[i]](xnu, by, ..., use.g.names = FALSE) else
                          BY.data.frame(xnu, by, FUN[[i]], ..., use.g.names = FALSE)), namFUN, give.names) else
            res[[ind]] <- if(fFUN) condsetn(list(FUN(xnu, by, ..., use.g.names = FALSE)), namFUN, give.names) else # give.names || !widel
                          condsetn(list(BY.data.frame(xnu, by, FUN, ..., use.g.names = FALSE, parallel = parallel, mc.cores = mc.cores)), namFUN, give.names) # give.names || !widel
        }
        if(nnul) {
          fcatFUN <- namcatFUN %in% .FAST_STAT_FUN
          if(is.list(catFUN))
            res[[lr]] <- condsetn(aplyfun(seq_along(namcatFUN), function(i)
                         if(fcatFUN[i]) catFUN[[i]](xnnu, by, ..., use.g.names = FALSE) else
                         BY.data.frame(xnnu, by, catFUN[[i]], ..., use.g.names = FALSE)), namcatFUN, give.names) else
            res[[lr]] <- if(fcatFUN) condsetn(list(catFUN(xnnu, by, ..., use.g.names = FALSE)), namcatFUN, give.names) else # give.names || !widel
                         condsetn(list(BY.data.frame(xnnu, by, catFUN, ..., use.g.names = FALSE, parallel = parallel, mc.cores = mc.cores)), namcatFUN, give.names) # give.names || !widel
        }
      return(res)
    } # fastest isung res list ?? or better combine at the end ??
    res <- agg(if(nul) `oldClass<-`(X[nu], "data.frame") else NULL, if(nnul) `oldClass<-`(X[nnu], "data.frame") else NULL) #colsubset(X, nu)

    if(sort.col && widel) o <- sort.list(c(if(!keep.by) NULL else numby,
                                           if(nul) rep(nu,length(FUN)) else NULL,
                                           if(nnul) rep(nnu,length(catFUN)) else NULL), method = "radix")

  } else { # custom aggregation:
    if(give.names == "auto") give.names <- TRUE
    namFUN <- names(custom)
    if(is.null(namFUN)) stop("custom needs to be a named list, see ?collap")
    fFUN <- namFUN %in% .FAST_STAT_FUN
    if(!keep.by) {
      res <- vector("list", 1L)
      ind <- 1L
    } else {
      res <- vector("list", 2L)
      res[[1L]] <- if(is.atomic(by)) list(`names<-`(list(attr(by, "levels")), namby)) else list(by[[4L]]) # nah... could add later using "c" ??
      ind <- 2L
    }
    res[[ind]] <- condsetn(aplyfun(seq_along(namFUN), function(i)
            if(fFUN[i]) match.fun(namFUN[i])(`oldClass<-`(X[custom[[i]]], "data.frame"), by, ..., use.g.names = FALSE) else
            BY.data.frame(X[custom[[i]]], by, namFUN[i], ..., use.g.names = FALSE)), namFUN, give.names)
    if(sort.col && widel) {
      o <- unlist(custom, use.names = FALSE)
      o <- sort.list(c(if(!keep.by) NULL else if(is.character(o)) namby else numby, o), method = "radix")
    }
  }
  if(widel) res <- unlist(unlist(res, FALSE), FALSE) else {
    if(length(FUN) > 1L || length(catFUN) > 1L || length(custom) > 1L) {
      res <- unlist(res, FALSE)
      if(multi.FUN.out == 2L) {
        ax[["row.names"]] <- if(is.list(by)) .set_row_names(by[[1L]]) else .set_row_names(length(res[[1L]]))  # fnlevels(by) best ??
        if(!keep.by) return(lapply(res, function(e) {
          ax[["names"]] <- names(e)
          return(setAttributes(e, ax)) })) else
          return(lapply(res[-1L], function(e) {
              ax[["names"]] <- c(namby, names(e))
              setAttributes(c(res[[1L]], e), ax) }))
      } else {
        if(multi.FUN.out != 4L) {
          res <- if(!keep.by) .Call(rbindlist, res, TRUE, TRUE, "Function") else # data.table:::Crbindlist
                 .Call(rbindlist, lapply(res[-1L], function(e) c(res[[1L]], e)), TRUE, TRUE, "Function")
        } else {
          if(!(nul && nnul) || customl) stop("long_dupl is only meaningful for aggregations with both numeric and categorical data, and multiple functions used for only one of the two data types!")
          mFUN <- length(FUN) > 1L
          nid <- if(mFUN) length(res) else 2L-!keep.by
          if(!keep.by) {
            res <- if(mFUN) lapply(res[-nid], function(e) c(e, res[[nid]])) else
                            lapply(res[-nid], function(e) c(res[[nid]], e))
          } else res <- if(mFUN) lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], e, res[[nid]])) else
                                 lapply(res[-c(nid, 1L)], function(e) c(res[[1L]], res[[nid]], e))
          res <- .Call(rbindlist, res, FALSE, FALSE, "Function")
        }
        if(sort.col)  o <- sort.list(c(0L, if(!keep.by) NULL else numby, nu, nnu), method = "radix")
      }
    } else message("multi.FUN.out options are only meaningful if multiple functions are used!")
  }

  if(sort.col) .Call(setcolorder, res, o) # data.table:::Csetcolorder
  ax[["names"]] <- names(res)
  ax[["row.names"]] <- .set_row_names(length(res[[1L]]))
  return(setAttributes(res, ax))
}




# Old Helpers:
# FUNchar <- function(FUN) {
#   if (is.function(FUN)) list(deparse(substitute(FUN)), FALSE) else if (is.character(FUN))
#     list(unlist(strsplit(FUN, ",", fixed = TRUE), use.names = FALSE), FALSE) else if (is.list(FUN)) {
#       if(names(FUN)) stop("If a list of functions is passed to FUN/catFUN, it needs to be named!")
#       list(names(FUN), TRUE)
#     } else stop("FUN/catFUN needs to be a function, string of functions, or named list of functions")
# }
# getdots <- function() { # data.tables getdots.R
#   # return a string vector of the arguments in '...'
#   # My long winded way: gsub(" ","",unlist(strsplit(deparse(substitute(list(...))),"[(,)]")))[-1]
#   # Peter Dalgaard's & Brian Ripley helped out and ended up with :
#
#   # as.character(match.call(sys.function(-1L), call=sys.call(-1L), expand.dots=FALSE)$...)
#   match.call(sys.function(-1L), call=sys.call(-1L), expand.dots=FALSE)$...
# }


# # Previous Attempt:
# collap.data.frame <- function(X, by = NULL, FUN = fmean, catFUN = fmode, Xcols = NULL, custom = NULL,
#                               drop.by = FALSE, sort.col = TRUE, sort.row = TRUE, reshape.long = FALSE,
#                               list.out = FALSE, ...) {
#   nam <- attr(X, "names")
#   # ax <- attributes(X)
#   clx <- class(X)
#   #   isDT <- inherits(X, "data.table")
#   isDT <- any(clx == "data.table") # ax[["class"]]
#   byl <- !is.null(by)
#   Xcolsl <- !is.null(Xcols)
#
#   if (byl) { # Aggregation by Groups
#     listby = is.list(by)
#     if (!listby) { # necessary??
#       if (is.numeric(by)) {
#         num = by
#         by = nam[by]
#       } else if (is.character(by)) {
#         by = unlist(strsplit(by, ",", fixed = TRUE), use.names = FALSE)
#         num = match(by, nam)
#       } else {
#         if (length(by)>2) {
#           nam <- all.vars(by)
#           if(isDT) X <- X[, nam, with = FALSE] else X <- X[nam]
#           by <- all.vars(by[[3]])
#         } else by <- all.vars(by) # or  all.vars(by[[3]]) ??
#         num <- match(by, nam) # use setdiff??
#       }
#     } else {
#       if (any(lengths(by) != NROW(X))) stop("arguments must have same length")
#       indl <- which(!nzchar(names(by)))
#       names(by)[indl] <- paste0("Group.", indl)
#     }
#     lby <- length(by)
#     # iby <- seq_len(lby)
#
#     if(is.null(custom)) { # Default mode !!
#
#       # Characterizing the variables
#       cols <- setdiff(seq_along(X), num) # better option than stediff ??
#       nu <- setdiff(which(vapply(X, is.numeric, FUN.VALUE = logical(length(X)))), num) # speed gain ??
#       lnu <- length(nu) # needed ???
#       nul <- lnu > 0
#       nnu <- setdiff(cols, nu)
#       nnul <- length(nuu) > 0
#
#       if(nul && nnul) { # Both numeric and categorical variables !!
#         FUNc <- FUNchar(FUN)
#         tFUNc <- FUNc[[2]]
#         FUNc <- FUNc[[1]]
#         lFUNc <- length(FUNc)
#         catFUNc <- FUNchar(catFUN) # additional security checks ??
#         tcatFUNc <- catFUNc[[2]]
#         catFUNc <- catFUNc[[1]]
#         lcatFUNc <- length(catFUNc)
#
#         if(lFUNc + lcatFUNc > 2) {  # multi-function aggregation
#           clpf <- FUNc %in% clpfun # good ?? ????
#           aclpf <- any(clpf)
#           clpcf <- catFUNc %in% clpfun
#           aclpcf <- any(clpcf)
#
#           if(aclpf || aclpcf)  # some fast functions
#             g <- if(listby) GRP(by, seq_along(by), keep.groups = !drop.by, order = sort.row) else GRP(X, num, keep.groups = !drop.by, order = sort.row)
#           res <- list(NULL,NULL) # needed ?? or res <- list(list(),list())
#           if(aclpf) {
#             res[[1]] <- if(all(clpf)) lapply(FUNc, function(x) match.fun(x)(X[nu], g, ...)) else if(tFUNc)
#               lapply(seq_along(FUNc), function(i) if(clpf[i]) match.fun(FUNc[i])(X[nu], g, ...) else DTagg(FUN[[i]], nu, ...)) else
#                 lapply(seq_along(FUNc), function(i) if(clpf[i]) match.fun(FUNc[i])(X[nu], g, ...) else DTagg(as.name(FUNc[i]), nu, ...))
#           } else
#             res[[1]] <- if(tFUNc) lapply(FUN, DTagg, nu, ...) else lapply(FUNc, function(x) DTagg(as.name(x), nu, ...))
#           if(aclpcf) {
#             res[[2]] <- if(all(clpcf)) lapply(catFUNc, function(x) match.fun(x)(X[nnu], g, ...)) else if(tcatFUNc)
#               lapply(seq_along(catFUNc), function(i) if(clpcf[i]) match.fun(catFUNc[i])(X[nnu], g, ...) else DTagg(catFUN[[i]], nnu, ...)) else
#                 lapply(seq_along(catFUNc), function(i) if(clpcf[i]) match.fun(catFUNc[i])(X[nnu], g, ...) else DTagg(as.name(catFUNc[i]), nnu, ...))
#           } else
#             res[[2]] <- if(tcatFUNc) lapply(catFUN, DTagg, nnu, ...) else lapply(catFUNc, function(x) DTagg(as.name(x), nnu, ...))
#           res <- .Internal(unlist(res, FALSE, FALSE))
#           names(res) <- c(FUNc, catFUNc)
#         } else { # Single Function Aggregation
#           clpf <- any(FUNc == clpfun) # good ?? ????
#           clpcf <- any(catFUNc == clpfun)
#
#           if(clpf || clpcf)  # some fast functions
#             g <- if(listby) GRP(by, seq_along(by), keep.groups = !drop.by, order = sort.row) else GRP(X, num, keep.groups = !drop.by, order = sort.row)
#           res <- list(NULL,NULL)
#           if(clpf) res[[1]] <- match.fun(FUNc)(X[nu], g, use.g.names = FALSE, ...) else { # without fast functions
#             if(!isDT) {
#               setDT(X)
#               isDT = TRUE
#             }
#             res[[1]] <- if(tFUNc) DTagg(FUN, nu, ...) else DTagg(as.name(FUNc), nu, ...)
#           }
#           if(clpcf) res[[2]] <- match.fun(catFUNc)(X[nnu], g, use.g.names = FALSE, ...) else { # without fast functions
#             if(!isDT) setDT(X)
#             res[[2]] <- if(tcatFUNc) DTagg(catFUN, nnu, ...) else DTagg(as.name(catFUNc), nnu, ...)
#           }
#           res <- .Internal(unlist(res, FALSE, FALSE))
#           names(res) <- c(FUNc, catFUNc)
#         }
#       } else { # Only numeric or categorical variables
#         if (nnul) {
#           FUN <- catFUN
#           nu <- nnu
#         }
#         FUNc <- FUNchar(FUN)
#         tFUNc <- FUNc[[2]]
#         FUNc <- FUNc[[1]]
#         if(length(FUNc) > 1) { # Multi-function aggregation
#           clpf <- FUNc %in% clpfun
#           if(any(clpf)) { # with some fast functions
#             g <- if(listby) GRP(by, seq_along(by), keep.groups = !drop.by, order = sort.row) else GRP(X, num, keep.groups = !drop.by, order = sort.row)
#             res <- if(all(clpf)) lapply(FUNc, function(x) match.fun(x)(X[nu], g, use.g.names = FALSE, ...)) else {
#               setDT(X)
#               if(tFUNc) lapply(seq_along(FUNc), function(i) if(clpf[i]) match.fun(FUNc[i])(X[nu], g, use.g.names = FALSE, ...) else DTagg(FUN[[i]], nu, ...)[,-1:lby]) else
#                 lapply(seq_along(FUNc), function(i) if(clpf[i]) match.fun(FUNc[i])(X[nu], g, use.g.names = FALSE, ...) else DTagg(as.name(FUNc[i]), nu, ...)[,-1:lby])
#             }
#             if(reshape.long) {
#               res[clpf] <- lapply(res[clpf],function(x)c(g[["groups"]],x))
#               res <- data.table::rbindlist(res, fill = FALSE, idcol = ifelse(nul,"FUN","catFUN"))
#               class(res) <- clx
#             } else if (!list.out) {
#               res[!clpf] <- lapply(res[!clpf],function(x)x[,-1:lby])
#               res <- do.call(cbind.data.frame, c(g[["groups"]], res)) # faster ?? using c() ??
#               class(res) <- clx
#             }
#           } else { # without fast functions
#             setDT(X)
#             res <- if(tFUNc) lapply(FUN, DTagg, nu, ...) else lapply(FUNc, function(x) DTagg(as.name(x), nu, ...))
#             names(res) <- FUNc
#             if(reshape.long) {
#               res <- data.table::rbindlist(res, fill = FALSE, idcol = ifelse(nul,"FUN","catFUN"))
#               class(res) <- clx
#             } else if (!list.out) {
#               res[-1] <- lapply(res[-1], function(x)x[,-1:lby])
#               res <- do.call(cbind.data.frame, res) # faster ?? using c() ??
#               class(res) <- clx
#             }
#           }
#         } else { # Simple Aggregation
#           if(any(FUNc == clpfun)) { # with a fast function
#             g <- if(listby) GRP(by, seq_along(by), keep.groups = !drop.by, order = sort.row) else GRP(X, num, keep.groups = !drop.by, order = sort.row)
#             res <- match.fun(FUNc)(X[nu], g, use.g.names = FALSE, ...)
#           } else { # without fast functions
#             setDT(X)
#             res <- if(tFUNc) DTagg(FUN, nu, ...) else DTagg(as.name(FUNc), nu, ...)
#           }
#         }
#       }
#     } else { # Custom Aggregation -> still implement !!
#
#     }
#   } else { # non-grouped aggregation !!!
#
#     if(is.null(custom)) { # Default mode !!
#
#       # Characterizing the variables
#       cols <- seq_along(X)
#       nu <- which(vapply(X, is.numeric, FUN.VALUE = logical(length(X))))
#       nnu <- setdiff(cols, nu)
#
#       if(length(nu) && length(nnu)) {
#         res <- list(NULL, NULL)
#         if(is.function(FUN)) res[[1]] <- lapply(X[nu], FUN, ...) else if(is.character(FUN))
#       }
#       # in general: just make funlist and use lapply, not DT!!
#
#       if(length(nu) && length(nnu)) { # Both numeric and categorical variables !!
#
#       }
#
#     } else { # Custom Mode
#
#     }
#   }
# }
#
#
# # Original collap:
# # If Xcols and 2-sided formula give error message and say: please choose!!
# # I'd also remove the factors argument, and use factor.vars and as.numeric.factor, or simply aggregate as.numeric.factor(X)
# collap.data.frame <- function(X, by = NULL, FUN = mean, catFUN = Mode, Xcols = NULL,
#                               custom = NULL, custom.names = TRUE, sort.col = TRUE, # renamed sort to sort.col
#                               sort.row = TRUE, # added argument -> if true use keyby with data.table, else by
#                               reshape.long = FALSE,
#                               list.out = FALSE, # renamed as.list to list.out -> more intuitive
#                               drop.cat = FALSE, drop.by = FALSE, show.statistic = TRUE, # renamed dropcat and dropby to drop.cat and drop.by -> more consistent with the other arguments
#                               parallel = FALSE, ...) # data.table is on now by default,if you think it is stable, we can do only data.table
# {
#
#   # All Arguments writen out with default settings (for testing)
#   # X; by = NULL; FUN = mean; catFUN = Mode; factors = "as.categorical"
#   # custom = NULL; custom.names = TRUE; collapse = TRUE; sort.col = TRUE
#   # sort.row = TRUE; ndigits = NULL; simplify = FALSE; reshape.long = FALSE
#   # na.rm = TRUE; replace.nan = TRUE; list.out = FALSE; drop.cat = FALSE
#   # drop.by = FALSE; show.statistic = TRUE; data.table = TRUE; parallel = FALSE
#
#   # Identifying the inputs: X
#   cl = class(X)
#   isDT = any(cl == "data.table")
#
#   # Identifying the inputs: FUN, catFUN and custom
#   customl = !is.null(custom)
#   if(customl) {
#     namFUN = names(FUN)
#     FUN = lapply(namFUN, match.fun)
#     lFUN = length(FUN)
#     iFUN = seq_len(lFUN)
#   } else {
#     nlf = !is.list(FUN)
#     if(is.character(FUN)) {
#       FUN = unlist(strsplit(FUN,",",fixed = TRUE), use.names = FALSE)
#       FUN = sapply(FUN, match.fun, descend = FALSE)
#     } else if (nlf) {
#       FUN = setNames(list(FUN), deparse(substitute(FUN)))
#     }
#     namFUN = names(FUN)
#     lFUN = length(FUN)
#     iFUN = seq_len(lFUN)
#
#     nlcf = !is.list(catFUN)
#     if(is.character(catFUN)) {
#       catFUN = unlist(strsplit(catFUN,",",fixed = TRUE), use.names = FALSE)
#       catFUN = sapply(catFUN, match.fun, descend = FALSE)
#     } else if (nlcf) {
#       catFUN = setNames(list(catFUN), deparse(substitute(catFUN)))
#     }
#     namcatFUN = names(catFUN)
#     lcatFUN = length(catFUN)
#     icatFUN = seq_len(lcatFUN)
#   }
#
#   # Identifying the inputs: by
#   byl = !is.null(by)
#   if(byl) {
#       if(is.list(by)) {
#         if(any(lengths(by) != nrow(X))) stop("arguments must have same length")
#           naml = names(by)
#         if(is.null(naml))
#           names(by) = paste0("Group.",seq_along(by)) else {
#           indl <- which(!nzchar(naml))
#           names(by)[indl] <- paste0("Group.", indl)
#         }
#         X = cbind(by,X) # needed ?? use GRP !!!
#         by = names(by)
#         num = match(by, names(X))
#       } else if(is.call(by)) {
#         if(length(by)>2) X = get_all_vars(by,X)
#         by = attr(terms(by),"term.labels") # faster way using all.vars
#         num = match(by, names(X))
#       } else if(is.character(by)) {
#       by = unlist(strsplit(by,",",fixed = TRUE), use.names = FALSE)
#       num = match(by, names(X))
#      } else if (is.numeric(by)) {
#       num = by
#       by = names(X)[num]
#      } else stop("'by' must be either a formula, a numeric or character vector signifying the columns over which to aggregate, or a list of columns")
#     lby = length(by)
#     iby = seq_len(lby)
#   }
#
#   # Characterizing the variables
#   cols = setdiff(seq_along(X), num)
#   nu = setdiff(which(vapply(X, is.numeric, TRUE)), num)
#   nnu = setdiff(cols,nu)
#
#   # Preprocessing
#   Xby = X[num] # use GRP !!
#
#   # Functions to perform the aggregation
#     # nlfs = all(nlf,ifelse(drop.cat,TRUE,nlcf)) # not needed anymore !!
#     agg <- function(df, by, FUN, nam, ...)
#       if(any(nam == clpfun)) FUN(df, by, ...) else BY.data.frame(df, by, FUN, ...)
#
#   # Exchanging arguments in special cases
#   if (!drop.cat) {
#     if (length(nu)==0 && length(nnu)>0) {
#       nu = nnu; drop.cat = TRUE; FUN = catFUN; lFUN = lcatFUN; iFUN = icatFUN; namFUN = namcatFUN
#     } else if (lFUN<=1 && lcatFUN>1) {
#       a = nu; nu = nnu; nnu = a; b = FUN; FUN = catFUN; catFUN = b
#       a = lFUN; lFUN = lcatFUN; lcatFUN = a; b = namFUN; namFUN = namcatFUN; namcatFUN = b
#       a = iFUN; iFUN = icatFUN; icatFUN = a
#     } else if (lcatFUN>1) {
#       reshape.long = FALSE; list.out = FALSE
#     }
#   }
#
#   # Computing output
#   if(!customl) {
#     if(parallel && lFUN>1) { # Numeric Variables
#       no_cores <- parallel::detectCores() - 1
#       cl <- parallel::makeCluster(no_cores)
#       # parallel::clusterExport(cl,ifelse(data.table,"data.table","aggregate.data.frame"))
#       res = parallel::parLapply(cl,iFUN,function(i) agg(df=X[nu],by=Xby,FUN[[i]],namFUN[i], ...))
#       stopCluster(cl)
#     } else if (lFUN>1) {
#       res = lapply(iFUN,function(i) agg(df=X[nu],by=Xby,FUN[[i]],namFUN[i], ...))
#     } else {
#       res = agg(df=X[nu],by=Xby,FUN[[1]],namFUN, ...)
#     }
#     if (!drop.cat && length(nnu)>0) { # Categorical Variables
#       if(parallel && lcatFUN>1) {
#         no_cores <- parallel::detectCores() - 1
#         cl <- parallel::makeCluster(no_cores)
#         # parallel::clusterExport(cl,ifelse(data.table,"data.table","aggregate.data.frame"))
#         Catres = parallel::parLapply(cl,icatFUN,function(i) agg(df=X[nnu],by=Xby,catFUN[[i]],namcatFUN[i], ...))
#         stopCluster(cl)
#       } else if (lcatFUN>1) {
#         Catres = lapply(icatFUN,function(i) agg(df=X[nnu],by=Xby,catFUN[[i]],namcatFUN[i], ...))
#       } else {
#         Catres = agg(df=X[nnu],by=Xby,catFUN[[1]],namcatFUN, ...)
#       }
#     }
#   } else { # Custom Mode
#     unlistFUN = unlist(FUN, use.names = FALSE)
#     ivf = !is.list(unlistFUN) # If FUN is a vector that means it contains colindices or colnames
#     if (ivf) {
#       if (is.null(namFUN)) stop("Names of list must be valid function names")
#       if (is.character(unlistFUN)) { # FUN contains column names as character strings
#         indl = lapply(FUN,function(x)match(unlist(strsplit(x,",",fixed = TRUE), use.names = FALSE), names(X)))
#         unf = namFUN
#         FUN = sapply(namFUN, match.fun)
#       } else { # FUN contains column indices
#         indl = FUN
#         unf = namFUN
#         FUN = sapply(namFUN, match.fun)
#       }
#     } else { # # If FUN was a character string or a list of functions
#       ind = sort(union(nu,nnu)); fFUN = FUN; FUN = unique(FUN)
#       unf = if (is.null(namFUN)) seq_along(FUN) else unique(namFUN)
#       indl = lapply(FUN,function(x)ind[sapply(fFUN,identical,x)])
#       if (length(ind) != lFUN) stop("Vector of custom functions needs to match columns of data")
#     }
#     if (!is.null(namFUN) && custom.names) { res = lapply(seq_along(unf),function(i){r = agg(df=X[indl[[i]]],by=Xby,FUN[[i]],unf[i])
#     names(r)[-iby] = paste0(names(r)[-iby],".",unf[i]); r})
#     } else res = lapply(seq_along(unf),function(i) agg(df=X[indl[[i]]],by=Xby,FUN[[i]],unf[i]))
#     if (list.out != "FUN") {
#       res = data.frame(res[[1]][iby],cbind.data.frame(lapply(res,function(x)x[-iby]), stringsAsFactors = FALSE), stringsAsFactors = FALSE)
#       if (!collapse) res = merge.data.frame(Xby,res)[order(ordr), , drop = FALSE]
#       rownames(res) = NULL
#       if (sort.col) res = res[order(c(num,unlist(indl, use.names = FALSE)))]
#     } else { names(res) = unf
#     if (drop.by) res = lapply(res, function(x)x[-iby])
#     }
#     lFUN = 1
#   }
#
#   # Preparing Output
#   if (is.null(custom)) {
#     if (reshape.long || lFUN == 1 || list.out == "FUN") { # Anything but multiple functions columns in parallel
#       if (lFUN>1) { # if multiple functions
#         if (!drop.cat && length(nnu)>0) { # if categorical variables
#           if (list.out != "FUN" && show.statistic) { # show a statistic
#             res = lapply(iFUN,function(x){
#               data.frame(res[[x]][iby],Statistic=namFUN[x],res[[x]][-iby],Catres[-iby], stringsAsFactors = FALSE)})
#             if (!collapse) res = lapply(res,function(x)merge.data.frame(Xby, x)[order(ordr), , drop = FALSE])
#             if (sort.col) { ord = order(c(num,-Inf,nu,nnu))
#             res = lapply(res,function(x)x[ord])}
#           } else { # no statistic
#             res = lapply(iFUN,function(x){ # only if show.statistic argument
#               data.frame(res[[x]],Catres[-iby], stringsAsFactors = FALSE)})
#             if (!collapse) res = lapply(res,function(x)merge.data.frame(Xby, x)[order(ordr), , drop = FALSE])
#             if (sort.col) { ord = order(c(num,nu,nnu))
#             res = lapply(res,function(x)x[ord])}
#           }
#         } else {
#           if (list.out != "FUN" && show.statistic) {
#             res = lapply(iFUN, function(x){ # No categorical but with statistic
#               data.frame(res[[x]][iby],Statistic=namFUN[x],res[[x]][-iby], stringsAsFactors = FALSE)})
#             if (!collapse) res = lapply(res,function(x)merge.data.frame(Xby, x)[order(ordr), , drop = FALSE])
#             if (sort.col) { ord = order(c(num,-Inf,nu))
#             res = lapply(res,function(x)x[ord])}
#           } else { # No categorical no statistic
#             if (!collapse) res = lapply(res,function(x)merge.data.frame(Xby, x)[order(ordr), , drop = FALSE])
#             if (sort.col) { ord = order(c(num,nu))
#             res = lapply(res,function(x)x[ord])}
#           }
#         }
#         names(res) = namFUN
#         if (list.out != "FUN") { # Combine unless FUN
#           res = if (show.statistic) { do.call(rbind.data.frame, c(res, make.row.names = FALSE, stringsAsFactors = FALSE))
#           } else do.call(rbind.data.frame, c(res, stringsAsFactors = FALSE))
#         }
#       } else if (!drop.cat && length(nnu)>0) { # Only one function but with categorical variables
#         res = data.frame(res,Catres[-iby], stringsAsFactors = FALSE)
#         if (!collapse) res = merge.data.frame(Xby,res)[order(ordr), , drop = FALSE]
#         if (sort.col) res = res[order(c(num,nu,nnu))]
#       } else { # if one function and no categorical variables
#         if (!collapse) res = merge.data.frame(Xby,res)[order(ordr), , drop = FALSE]
#         if (sort.col) res = res[order(c(num,nu))]
#       }
#     } else { # If multiple functions all in a row
#       res = lapply(iFUN,function(x) {
#         names(res[[x]])[-iby] = paste0(names(res[[x]])[-iby],".",namFUN[x]); res[[x]]})
#       res[-1] = lapply(iFUN[-1], function(x) res[[x]][-iby])
#       if (!drop.cat && length(nnu)>0) { # If categorical variables
#         if (lcatFUN>1) { # Multiple categorical functions
#           Catres = lapply(icatFUN,function(x) {
#             names(Catres[[x]])[-iby] = paste0(names(Catres[[x]])[-iby],".",namcatFUN[x]); Catres[[x]]})
#           Catres[-1] = lapply(icatFUN[-1], function(x) Catres[[x]][-iby])
#           Catres = do.call(cbind,Catres)
#           nnu = rep(nnu,lcatFUN)
#         }
#         Catres = Catres[-iby]
#         res = do.call(cbind,res)
#         res = cbind.data.frame(res,Catres)
#         if (!collapse) res = merge.data.frame(Xby,res)[order(ordr), , drop = FALSE]
#         if (sort.col) res = res[order(c(num,rep(nu,lFUN),nnu))]
#       } else { # No categorical variables
#         res = do.call(cbind.data.frame,res)
#         if (!collapse) res = merge.data.frame(Xby,res)[order(ordr), , drop = FALSE]
#         if (sort.col) res = res[order(c(num,rep(nu,lFUN)))]
#       }
#     }
#   }
#   if (list.out==FALSE) {
#     if (lby == 1 && by == "nullID") res = res[-1]
#     if (drop.by) res = res[-match(by,names(res))]
#     if (isDT) setDT(res)
#     types = vapply(res, class, FUN.VALUE = character(1))
#     if (!is.null(ndigits)) {
#       numt = which(types %in% c("numeric","integer","double","atomic"))
#       res[numt] = round(res[numt],ndigits)
#     }
#     if (simplify) { # Probably an extended set of condition is faster than calling vapply again, to be comprehensive we also need to keep track of the Xby types to choose whether to simplify or not.
#       #if ((any(cl=="matrix") || (any(dim(res)==1) && (drop.cat || length(nnu)<=0))) && (lFUN==1 || !reshape.long))
#       if (length(unique(types))==1) res = drop(as.matrix(res))
#     }
#   } else {
#     if (lby == 1 && by == "nullID") lapply(res,function(x)x[-1])
#     if (list.out == "by" && by != "nullID") {
#       num = match(by,names(res))
#       res = split.data.frame(res[-num],res[num], drop = TRUE, lex.order = TRUE)
#     }
#     if (isDT) lapply(res,setDT)
#     inndigits = !is.null(ndigits)
#     res = lapply(res, function(x) {
#       types = vapply(x, class, FUN.VALUE = character(1))
#       numt = which(types %in% c("numeric","integer","double","atomic"))
#       if (inndigits) x[numt] = round(x[numt],ndigits)
#       if (simplify && length(unique(types))==1) x = drop(as.matrix(x))
#       x
#     }) # again, checking with conditions will be faster, therefore need to check Xby types above.
#     #if ((any(cl=="matrix") || (any(dim(res)==1) && (drop.cat || length(nnu)<=0))) && (lFUN==1 || !reshape.long))
#     #res = lapply(res,function(x)drop(as.matrix(x)))
#   }
#   return(res)
# }
