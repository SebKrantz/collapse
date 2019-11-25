# library(Rcpp)
# sourceCpp('src/fgrowth.cpp', rebuild = TRUE) # Todo: Thoroughly test these functions !!!
# sourceCpp('src/fgrowtha.cpp', rebuild = TRUE)
# sourceCpp('src/fgrowthl.cpp', rebuild = TRUE)
# source("R/GRP.R")
# source("R/small_helper.R")
# source("R/quick_conversion.R")

# For principle innovations of this code see flag.R and flag.cpp # stubs instead of give.names !!

# ,logdiff,stubs  # TRY   x, n = 1, diff = 1, ..., t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE
fgrowth <- function(x, n = 1, diff = 1, ...) { # , g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE
  UseMethod("fgrowth", x)
}
fgrowth.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fgrowth.default", dotstostr(...)))
  if(is.null(g))
    return(fgrowthCpp(x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),logdiff,stubs)) else if(is.atomic(g)) {
      if(is.factor(g)) nl <- fnlevels(g) else {
        g <- qG(g)
        nl <- attr(g, "N.groups")
      }
      return(fgrowthCpp(x,n,diff,fill,nl,g,NULL,G_t(t),logdiff,stubs))
    } else {
      if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
      return(fgrowthCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),logdiff,stubs))
    }
}
fgrowth.pseries <- function(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fgrowth.pseries", dotstostr(...)))
  index <- attr(x, "index")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  if(is.matrix(x))
  fgrowthmCpp(x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs) else
  fgrowthCpp(x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs)
}
fgrowth.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fgrowth.matrix", dotstostr(...)))
  if(is.null(g))
    return(fgrowthmCpp(x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),logdiff,stubs)) else if(is.atomic(g)) {
      if(is.factor(g)) nl <- fnlevels(g) else {
        g <- qG(g)
        nl <- attr(g, "N.groups")
      }
      fgrowthmCpp(x,n,diff,fill,nl,g,NULL,G_t(t),logdiff,stubs)
    } else {
      if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
      fgrowthmCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),logdiff,stubs)
    }
}
fgrowth.grouped_df <- function(x, n = 1, diff = 1, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, keep.ids = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fgrowth.grouped_df", dotstostr(...)))
  g <- GRP.grouped_df(x)
  tsym <- deparse(substitute(t))
  nam <- names(x)
  gn <- which(nam %in% g[[5L]])
  if(!(tsym == "NULL" || is.na(tn <- match(tsym, nam)))) {
    if(any(gn == tn)) stop("timevar coincides with grouping variables!")
    t <- x[[tn]]
    gn <- c(gn, tn)
  }
  if(length(gn)) {
    if(!keep.ids)
      return(fgrowthlCpp(x[-gn],n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),logdiff,stubs)) else {
        ax <- attributes(x)
        class(x) <- NULL # Works for multiple lags !!
        res <- c(x[gn],fgrowthlCpp(x[-gn],n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),logdiff,stubs))
        ax[["names"]] <- names(res)
        return(setAttributes(res, ax))
      }
  } else return(fgrowthlCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),logdiff,stubs))
}
fgrowth.data.frame <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fgrowth.data.frame", dotstostr(...)))
  if(is.null(g))
    return(fgrowthlCpp(x,n,diff,fill, t = G_t(t,FALSE), names = stubs)) else if(is.atomic(g)) {
      if(is.factor(g)) nl <- fnlevels(g) else {
        g <- qG(g)
        nl <- attr(g, "N.groups")
      }
      fgrowthlCpp(x,n,diff,fill,nl,g,NULL,G_t(t),logdiff,stubs)
    } else {
      if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
      fgrowthlCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),logdiff,stubs)
    }
}
fgrowth.pdata.frame <- function(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
  if(!missing(...)) stop(sprintf("Unknown argument %s passed to fgrowth.pdata.frame", dotstostr(...)))
  index <- attr(x, "index")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  fgrowthlCpp(x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs)
}

# x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...
G <- function(x, n = 1, diff = 1, ...) {
  UseMethod("G", x)
}
G.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
  fgrowth.default(x, n, diff, g, t, fill, logdiff, stubs, ...)
}
G.pseries <- function(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
  fgrowth.pseries(x, n, diff, fill, logdiff, stubs, ...)
}
G.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
  fgrowth.matrix(x, n, diff, g, t, fill, logdiff, stubs, ...)
}
G.grouped_df <- fgrowth.grouped_df
G.data.frame <- function(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric,
                         fill = NA, logdiff = FALSE, stubs = TRUE, keep.ids = TRUE, ...) {

  if(!missing(...)) stop(sprintf("Unknown argument %s passed to G.data.frame", dotstostr(...)))
  if(is.call(by) || is.call(t)) {
    ax <- attributes(x)
    class(x) <- NULL
    nam <- ax[["names"]]

    if(is.call(by)) {
      if(length(by) == 3L) {
        cols <- anyNAerror(match(all.vars(by[[2L]]), nam), "Unknown variables passed to by!")
        gn <- anyNAerror(match(all.vars(by[[3L]]), nam), "Unknown variables passed to by!")
      } else {
        gn <- anyNAerror(match(all.vars(by), nam), "Unknown variables passed to by!")
        cols <- if(is.null(cols)) seq_along(x)[-gn] else cols2int(cols, x, nam)
      }
      by <- if(length(gn) == 1L) at2GRP(x[[gn]]) else GRP(x, gn, return.groups = FALSE)
      if(!keep.ids) gn <- NULL
    } else {
      gn <- NULL
      if(!is.null(cols)) cols <- cols2int(cols, x, nam)
      if(!is.GRP(by)) by <- if(is.null(by)) list(0L, 0L, NULL) else if(is.atomic(by)) # Necessary for if by is passed externally !!
        at2GRP(by) else GRP(by, return.groups = FALSE)
    }

    if(is.call(t)) {
      t <- all.vars(t)
      tn <- anyNAerror(match(t, nam), "Unknown variables passed to t!")
      t <- x[[tn]]
      cols <- if(is.null(cols)) seq_along(x)[-tn] else cols[cols != tn]
      if(keep.ids) gn <- c(gn, tn)
    }

    res <- if(length(gn))
      c(x[gn], fgrowthlCpp(x[cols],n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),logdiff,stubs)) else
        fgrowthlCpp(x[cols],n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),logdiff,stubs)
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(!is.null(cols)) {
    ax <- attributes(x)
    x <- if(is.function(cols)) unclass(x)[vapply(x, cols, TRUE)] else unclass(x)[cols]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }

  if(is.null(by))
    return(fgrowthlCpp(x,n,diff,fill,0L,0L,NULL,G_t(t,FALSE),logdiff,stubs)) else if(is.atomic(by)) {
      if(is.factor(by)) nl <- fnlevels(by) else {
        by <- qG(by)
        nl <- attr(by, "N.groups")
      }
      fgrowthlCpp(x,n,diff,fill,nl,by,NULL,G_t(t),logdiff,stubs)
    } else {
      if(!is.GRP(by)) by <- GRP(by, return.groups = FALSE)
      fgrowthlCpp(x,n,diff,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),logdiff,stubs)
    }
}
G.pdata.frame <- function(x, n = 1, diff = 1, cols = is.numeric, fill = NA, logdiff = FALSE, stubs = TRUE, keep.ids = TRUE, ...) {

  if(!missing(...)) stop(sprintf("Unknown argument %s passed to G.pdata.frame", dotstostr(...)))
  ax <- attributes(x)
  nam <- ax[["names"]]
  index <- ax[["index"]]

  if(keep.ids) {
    gn <- which(nam %in% names(index))
    if(length(gn) && is.null(cols)) cols <- seq_along(x)[-gn]
  } else gn <- NULL

  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])

  if(!is.null(cols)) cols <- cols2int(cols, x, nam)

  if(length(gn) && !is.null(cols)) {
    class(x) <- NULL # Works for multiple lags !!
    res <- c(x[gn], fgrowthlCpp(x[cols],n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs))
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(!length(gn)) # could speed up ??
    return(fgrowthlCpp(x[cols],n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs)) else
      return(fgrowthlCpp(x,n,diff,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs))
}


# Previous versions were deleted -> see previous verstions of L.data.frame and L.pdata.frame in flag.R !!


# Old versions: several functions for t, and no Xcols...
# is_factor_G <- function(x) {
#   if(is.factor(x) || is.integer(x)) x else qG(x)
# }
# is_GRP_G <- function(l) {
#   if(all(class(l) == "GRP")) l[[2L]] else GRP(l, return.groups = FALSE)[[2L]]
# }
# fgrowth <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
#   UseMethod("fgrowth", x)
# }
# fgrowth.default <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
#   if(is.null(g)) {
#     if(is.null(t)) fgrowthCpp(x,n,diff,fill,0L,0L,NULL,NULL,logdiff,stubs) else if(is.atomic(t))
#       fgrowthCpp(x,n,diff,fill,0L,0L,NULL,is_factor_G(t),logdiff,stubs) else
#         fgrowthCpp(x,n,diff,fill,0L,0L,NULL,is_GRP_G(t),logdiff,stubs)
#   } else if(is.atomic(g)) {
#     if(!is.factor(g)) g = qF(g)
#     if(is.null(t)) {
#       message("Panel-growth-rate computed without timevar: Assuming ordered data")
#       fgrowthCpp(x,n,diff,fill,nlevels(g),g,NULL,NULL,logdiff,stubs)
#     } else if(is.atomic(t)) fgrowthCpp(x,n,diff,fill,nlevels(g),g,NULL,is_factor_G(t),logdiff,stubs) else
#       fgrowthCpp(x,n,diff,fill,nlevels(g),g,NULL,is_GRP_G(t),logdiff,stubs)
#   } else {
#     if(!all(class(g) == "GRP")) g = GRP(g, return.groups = FALSE)
#     if(is.null(t)) {
#       message("Panel-growth-rate computed without timevar: Assuming ordered data")
#       fgrowthCpp(x,n,diff,fill,g[[1L]],g[[2L]],NULL,NULL,logdiff,stubs)
#     } else if(is.atomic(t)) fgrowthCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],is_factor_G(t),logdiff,stubs) else
#       fgrowthCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],is_GRP_G(t),logdiff,stubs)
#   }
# }
# fgrowth.pseries <- function(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
#   index <- attr(x, "index")
#   fgrowthCpp(x,n,diff,fill,nlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs)
# }
# fgrowth.matrix <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
#   if(is.null(g)) {
#     if(is.null(t)) fgrowthmCpp(x,n,diff,fill,0L,0L,NULL,NULL,logdiff,stubs) else if(is.atomic(t))
#       fgrowthmCpp(x,n,diff,fill,0L,0L,NULL,is_factor_G(t),logdiff,stubs) else
#         fgrowthmCpp(x,n,diff,fill,0L,0L,NULL,is_GRP_G(t),logdiff,stubs)
#   } else if(is.atomic(g)) {
#     if(!is.factor(g)) g = qF(g)
#     if(is.null(t)) {
#       message("Panel-growth-rate computed without timevar: Assuming ordered data")
#       fgrowthmCpp(x,n,diff,fill,nlevels(g),g,NULL,NULL,logdiff,stubs)
#     } else if(is.atomic(t)) fgrowthmCpp(x,n,diff,fill,nlevels(g),g,NULL,is_factor_G(t),logdiff,stubs) else
#       fgrowthmCpp(x,n,diff,fill,nlevels(g),g,NULL,is_GRP_G(t),logdiff,stubs)
#   } else {
#     if(!all(class(g) == "GRP")) g = GRP(g, return.groups = FALSE)
#     if(is.null(t)) {
#       message("Panel-growth-rate computed without timevar: Assuming ordered data")
#       fgrowthmCpp(x,n,diff,fill,g[[1L]],g[[2L]],NULL,NULL,logdiff,stubs)
#     } else if(is.atomic(t)) fgrowthmCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],is_factor_G(t),logdiff,stubs) else
#       fgrowthmCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],is_GRP_G(t),logdiff,stubs)
#   }
# }
# fgrowth.data.frame <- function(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
#   if(is.null(g)) {
#     if(is.null(t)) fgrowthlCpp(x,n,diff,fill,0L,0L,NULL,NULL,logdiff,stubs) else if(is.atomic(t))
#       fgrowthlCpp(x,n,diff,fill,0L,0L,NULL,is_factor_G(t),logdiff,stubs) else
#         fgrowthlCpp(x,n,diff,fill,0L,0L,NULL,is_GRP_G(t),logdiff,stubs)
#   } else if(is.atomic(g)) {
#     if(!is.factor(g)) g = qF(g)
#     if(is.null(t)) {
#       message("Panel-growth-rate computed without timevar: Assuming ordered data")
#       fgrowthlCpp(x,n,diff,fill,nlevels(g),g,NULL,NULL,logdiff,stubs)
#     } else if(is.atomic(t)) fgrowthlCpp(x,n,diff,fill,nlevels(g),g,NULL,is_factor_G(t),logdiff,stubs) else
#       fgrowthlCpp(x,n,diff,fill,nlevels(g),g,NULL,is_GRP_G(t),logdiff,stubs)
#   } else {
#     if(!all(class(g) == "GRP")) g = GRP(g, return.groups = FALSE)
#     if(is.null(t)) {
#       message("Panel-growth-rate computed without timevar: Assuming ordered data")
#       fgrowthlCpp(x,n,diff,fill,g[[1L]],g[[2L]],NULL,NULL,logdiff,stubs)
#     } else if(is.atomic(t)) fgrowthlCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],is_factor_G(t),logdiff,stubs) else
#       fgrowthlCpp(x,n,diff,fill,g[[1L]],g[[2L]],g[[3L]],is_GRP_G(t),logdiff,stubs)
#   }
# }
# fgrowth.pdata.frame <- function(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, stubs = TRUE, ...) {
#   index <- attr(x, "index")
#   fgrowthlCpp(x,n,diff,fill,nlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],logdiff,stubs)
# }

