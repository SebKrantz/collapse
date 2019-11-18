# library(Rcpp)
# sourceCpp('src/flag.cpp', rebuild = TRUE) # Todo: Thoroughly check !!
# sourceCpp('src/flaga.cpp', rebuild = TRUE)
# sourceCpp('src/flagl.cpp', rebuild = TRUE)
# source("R/GRP.R")
# source("R/small_helper.R")
# source("R/quick_conversion.R")

# instead of more options for L.matrix and L.data.frame, could do without method dispatch??
# stubs instead of give.names
  # g = NULL, t = NULL, fill = NA, stubs = TRUE
flag <- function(x, n = 1, ...) {
  UseMethod("flag", x)
}
flag.default <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(is.null(g))
    return(flagleadCpp(x,n,fill,0L,0L,NULL,G_t(t,FALSE),stubs)) else if(is.atomic(g)) {
    if(is.factor(g)) nl <- fnlevels(g) else {
      g <- qG(g)
      nl <- attr(g, "N.groups")
    }
    return(flagleadCpp(x,n,fill,nl,g,NULL,G_t(t),stubs))
  } else {
    if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
    return(flagleadCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs))
  }
}
flag.pseries <- function(x, n = 1, fill = NA, stubs = TRUE, ...) {
  index <- attr(x, "index")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  if(is.matrix(x))
  flagleadmCpp(x,n,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs) else
  flagleadCpp(x,n,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs)
}
flag.matrix <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(is.null(g))
    return(flagleadmCpp(x,n,fill,0L,0L,NULL,G_t(t,FALSE),stubs)) else if(is.atomic(g)) {
    if(is.factor(g)) nl <- fnlevels(g) else {
      g <- qG(g)
      nl <- attr(g, "N.groups")
    }
    flagleadmCpp(x,n,fill,nl,g,NULL,G_t(t),stubs)
  } else {
    if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
    flagleadmCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs)
  }
}
flag.grouped_df <- function(x, n = 1, t = NULL, fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {
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
      return(flagleadlCpp(x[-gn],n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs)) else {
        ax <- attributes(x)
        class(x) <- NULL # Works for multiple lags !!
        res <- c(x[gn],flagleadlCpp(x[-gn],n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs))
        ax[["names"]] <- names(res)
        return(setAttributes(res, ax))
      }
  } else return(flagleadlCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs))
}
flag.data.frame <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  if(is.null(g))
    return(flagleadlCpp(x,n,fill,0L,0L,NULL,G_t(t,FALSE),stubs)) else if(is.atomic(g)) {
      if(is.factor(g)) nl <- fnlevels(g) else {
        g <- qG(g)
        nl <- attr(g, "N.groups")
      }
      flagleadlCpp(x,n,fill,nl,g,NULL,G_t(t),stubs)
    } else {
      if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
      flagleadlCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs)
    }
}
flag.pdata.frame <- function(x, n = 1, fill = NA, stubs = TRUE, ...) {
  index <- attr(x, "index")
  if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
  flagleadlCpp(x,n,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs)
}

# use xt instead of by ??? #   x, n = 1, g = NULL, by = NULL, t = NULL, cols = is.numeric, fill = NA, stubs = TRUE, keep.ids = TRUE, ...
L <- function(x, n = 1, ...) {
  UseMethod("L", x)
}
L.default <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  flag.default(x, n, g, t, fill, stubs, ...)
}
L.pseries <- function(x, n = 1, fill = NA, stubs = TRUE, ...) {
  flag.pseries(x, n, fill, stubs, ...)
}
L.matrix <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  flag.matrix(x, n, g, t, fill, stubs, ...)
}
L.grouped_df <- flag.grouped_df
L.data.frame <- function(x, n = 1, by = NULL, t = NULL, cols = is.numeric,
                         fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {

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
    c(x[gn], flagleadlCpp(x[cols],n,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),stubs)) else
    flagleadlCpp(x[cols],n,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),stubs)
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(!is.null(cols)) {
    ax <- attributes(x)
    x <- if(is.function(cols)) unclass(x)[vapply(x, cols, TRUE)] else unclass(x)[cols]
    ax[["names"]] <- names(x)
    setattributes(x, ax)
  }

  if(is.null(by))
    return(flagleadlCpp(x,n,fill,0L,0L,NULL,G_t(t,FALSE),stubs)) else if(is.atomic(by)) {
    if(is.factor(by)) nl <- fnlevels(by) else {
      by <- qG(by)
      nl <- attr(by, "N.groups")
    }
    flagleadlCpp(x,n,fill,nl,by,NULL,G_t(t),stubs)
  } else {
    if(!is.GRP(by)) by <- GRP(by, return.groups = FALSE)
    flagleadlCpp(x,n,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),stubs)
  }
}
L.pdata.frame <- function(x, n = 1, cols = is.numeric, fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {

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
    res <- c(x[gn], flagleadlCpp(x[cols],n,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs))
    ax[["names"]] <- names(res)
    return(setAttributes(res, ax))
  } else if(!length(gn)) # could speed up ??
    return(flagleadlCpp(x[cols],n,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs)) else
    return(flagleadlCpp(x,n,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],stubs))
}
 # , g = NULL, t = NULL, fill = NA, stubs = TRUE
flead <- function(x, n = 1, ...) {
  UseMethod("flead", x)
}
flead.default <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  flag.default(x, -n, g, t, fill, stubs, ...)
}
flead.pseries <- function(x, n = 1, fill = NA, stubs = TRUE, ...) {
  flag.pseries(x, -n, fill, stubs, ...)
}
flead.matrix <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  flag.matrix(x, -n, g, t, fill, stubs, ...)
}
flead.grouped_df <- function(x, n = 1, t = NULL, fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {
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
      return(flagleadlCpp(x[-gn],-n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs)) else {
        ax <- attributes(x)
        class(x) <- NULL # Works for multiple lags !!
        res <- c(x[gn],flagleadlCpp(x[-gn],-n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs))
        ax[["names"]] <- names(res)
        return(setAttributes(res, ax))
      }
  } else return(flagleadlCpp(x,-n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),stubs))
}
flead.data.frame <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  flag.data.frame(x, -n, g, t, fill, stubs, ...)
}
flead.pdata.frame <- function(x, n = 1, fill = NA, stubs = TRUE, ...) {
  flag.pdata.frame(x, -n, fill, stubs, ...)
}
 # x, n = 1, g = NULL, by = NULL, t = NULL, cols = is.numeric, fill = NA, stubs = TRUE, keep.ids = TRUE, keep.ids = TRUE, ...
F <- function(x, n = 1, ...) {
  UseMethod("F", x)
}
F.default <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  flag.default(x, -n, g, t, fill, stubs, ...)
}
F.pseries <- function(x, n = 1, fill = NA, stubs = TRUE, ...) {
  flag.pseries(x, -n, fill, stubs, ...)
}
F.matrix <- function(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, ...) {
  flag.matrix(x, -n, g, t, fill, stubs, ...)
}
F.grouped_df <- flead.grouped_df
F.data.frame <- function(x, n = 1, by = NULL, t = NULL, cols = is.numeric,
                         fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {
  L.data.frame(x, -n, by, t, cols, fill, stubs, keep.ids, ...)
}
F.pdata.frame <- function(x, n = 1, cols = is.numeric, fill = NA, stubs = TRUE, keep.ids = TRUE, ...) {
  L.pdata.frame(x, -n, cols, fill, stubs, keep.ids, ...)
}



# Previous versions: still allowing column names and indices:
# L.grouped_df <- function(x, n = 1, t = NULL, fill = NA, give.names = TRUE, drop.ids = FALSE, ...) {
#   g <- GRP.grouped_df(x)
#   tsym <- deparse(substitute(t))
#   nam <- names(x)
#   gn <- which(nam %in% g[[5L]])
#   if(!(tsym == "NULL" || is.na(tn <- match(tsym, nam)))) {
#     if(any(gn == tn)) stop("timevar coincides with grouping variables!")
#     t <- x[[tn]]
#     gn <- c(gn, tn)
#   }
#   if(length(gn)) {
#     if(drop.ids)
#       return(flagleadlCpp(x[-gn],n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),give.names)) else {
#         ax <- attributes(x)
#         class(x) <- NULL # Works for multiple lags !!
#         res <- c(x[gn],flagleadlCpp(x[-gn],n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),give.names))
#         ax[["names"]] <- names(res)
#         return(`attributes<-`(res, ax))
#       }
#   } else return(flagleadlCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),give.names))
# }
# L.data.frame <- function(x, n = 1, by = NULL, t = NULL, cols = is.numeric,
#                          fill = NA, give.names = TRUE, drop.ids = FALSE, ...) {
#   if(!is.null(t)) {
#     if(is.call(t)) t <- all.vars(t)
#     t1l <- length(t) == 1L
#     if(t1l) {
#       tv <- t
#       t <- x[[t]]
#       # if(drop.t) x[[tv]] <- NULL # not needed, drop.ids is good (drops all identifiers) !!
#     }
#   } else t1l <- FALSE
#   if(is.null(by)) {
#     if(is.null(cols)) return(flagleadlCpp(x,n,fill,0L,0L,NULL,G_t(t,FALSE),give.names)) else {
#       if(is.function(cols)) cols <- vapply(x, cols, TRUE)
#       return(flagleadlCpp(x[cols],n,fill,0L,0L,NULL,G_t(t,FALSE),give.names))
#     }
#   } else if (is.atomic(by) || is.call(by)) {
#     if(is.call(by) || length(by) != nrow(x)) {
#       nam <- names(x)
#       if(is.call(by) && length(by) == 3) {
#         v <- match(all.vars(by[[2L]]), nam)
#         by <- match(all.vars(by[[3L]]), nam)
#       } else {
#         if(is.call(by)) by <- match(all.vars(by), nam) else if(is.character(by)) by <- match(by, nam)
#         v <- if(is.null(cols)) seq_along(x)[-by] else if(is.function(cols))
#           setdiff(which(vapply(x, cols, TRUE)), by) else if(is.character(cols))
#             match(cols, nam) else cols
#       }
#       g <- GRP(x, by, return.groups = FALSE)
#       if(t1l) { # If time-variable supplied !!
#         if(is.character(tv)) tv <- match(tv, nam)
#         v <- setdiff(v, tv)
#         by <- c(by, tv)
#       }
#       if(drop.ids)
#         return(flagleadlCpp(x[v],n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),give.names)) else {
#           ax <- attributes(x)
#           class(x) <- NULL # Works for multiple lags !!
#           res <- c(x[by], flagleadlCpp(x[v],n,fill,g[[1L]],g[[2L]],g[[3L]],G_t(t),give.names))
#           ax[["names"]] <- names(res)
#           return(`attributes<-`(res, ax))
#         }
#     } else if(is.factor(by)) return(flagleadlCpp(x,n,fill,fnlevels(by),by,NULL,G_t(t),give.names)) else {
#       by <- qG(by, ordered = FALSE)
#       return(flagleadlCpp(x,n,fill,attr(by,"N.groups"),by,NULL,G_t(t),give.names))
#     }
#   } else {
#     if(!is.GRP(by)) by <- GRP(by, return.groups = FALSE)
#     if(is.null(cols)) return(flagleadlCpp(x,n,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),give.names)) else {
#       if(is.function(cols)) cols <- vapply(x, cols, TRUE)
#       return(flagleadlCpp(x[cols],n,fill,by[[1L]],by[[2L]],by[[3L]],G_t(t),give.names))
#     }
#   }
# }
# L.pdata.frame <- function(x, n = 1, cols = is.numeric, fill = NA, give.names = TRUE, drop.ids = FALSE, ...) {
#   ax <- attributes(x)
#   nam <- ax[["names"]]
#   index <- ax[["index"]]
#   gn <- match(names(index), nam)
#   gn2 = gn <- gn[!is.na(gn)]
#   if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
#   if(!is.null(cols)) {
#     if(is.function(cols)) cols <- seq_along(x)[!vapply(x, cols, TRUE)] else if(is.character(cols))
#       cols <- seq_along(x)[-match(cols, nam)] else if(is.logical(cols))
#         cols <- seq_along(x)[!cols] else if(is.numeric(cols))
#           cols <- seq_along(x)[-cols] else stop("cols needs to be a function, column names, indices or a logical vector")
#         if(drop.ids) gn <- c(gn, cols) else if(length(gn)) gn2 <- c(gn2, cols) else gn <- cols
#   }
#   if(length(gn)) {
#     if(drop.ids)
#       return(flagleadlCpp(x[-gn],n,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],give.names)) else {
#         class(x) <- NULL # Works for multiple lags !!
#         res <- c(x[gn], flagleadlCpp(x[-gn2],n,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],give.names))
#         ax[["names"]] <- names(res)
#         return(`attributes<-`(res, ax))
#       }
#   } else return(flagleadlCpp(x[-gn],n,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],give.names))
# }








# Older versions: several functions for t, and no Xcols...
# is_factor_G <- function(x) {
#   if(is.null(x) || is.factor(x) || is.integer(x)) x else qG(x)
# }
# is_GRP_G <- function(l) {
#   if(all(class(l) == "GRP")) l[[2L]] else GRP(l, return.groups = FALSE)[[2L]]
# }

# flag <- function(x, n = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...) {
#   UseMethod("flag", x)
# }
# flag.default <- function(x, n = 1, g = NULL, t = NULL, fill = NULL, give.names = TRUE, ...) {
#   if(is.null(g)) {
#     if(is.null(t)) flagleadCpp(x,n,fill, names = give.names) else if(is.atomic(t))
#       flagleadCpp(x,n,fill, t = is_factor_G(t), names = give.names) else
#         flagleadCpp(x,n,fill, t = is_GRP_G(t), names = give.names)
#   } else if(is.atomic(g)) {
#     if(!is.factor(g)) g = qF(g)
#     if(is.null(t)) {
#       message("Panel-lag computed without timevar: Assuming ordered data")
#       flagleadCpp(x,n,fill,fnlevels(g),g,NULL,NULL,give.names)
#     } else if(is.atomic(t)) flagleadCpp(x,n,fill,fnlevels(g),g,NULL,is_factor_G(t),give.names) else
#       flagleadCpp(x,n,fill,fnlevels(g),g,NULL,is_GRP_G(t),give.names)
#   } else {
#     if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
#     if(is.null(t)) {
#       message("Panel-lag computed without timevar: Assuming ordered data")
#       flagleadCpp(x,n,fill,g[[1L]],g[[2L]],NULL,NULL,give.names)
#     } else if(is.atomic(t)) flagleadCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],is_factor_G(t),give.names) else
#       flagleadCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],is_GRP_G(t),give.names)
#   }
# }
# flag.pseries <- function(x, n = 1, fill = NA, give.names = TRUE, ...) {
#   index <- attr(x, "index")
#   if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
#   flagleadCpp(x,n,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],give.names)
# }
# flag.matrix <- function(x, n = 1, g = NULL, t = NULL, fill = NULL, give.names = TRUE, ...) {
#   if(is.null(g)) {
#     if(is.null(t)) flagleadmCpp(x,n,fill, names = give.names) else if(is.atomic(t))
#       flagleadmCpp(x,n,fill, t = is_factor_G(t), names = give.names) else
#         flagleadmCpp(x,n,fill, t = is_GRP_G(t), names = give.names)
#   } else if(is.atomic(g)) {
#     if(!is.factor(g)) g = qF(g)
#     if(is.null(t)) {
#       message("Panel-lag computed without timevar: Assuming ordered data")
#       flagleadmCpp(x,n,fill,fnlevels(g),g,NULL,NULL,give.names)
#     } else if(is.atomic(t)) flagleadmCpp(x,n,fill,fnlevels(g),g,NULL,is_factor_G(t),give.names) else
#       flagleadmCpp(x,n,fill,fnlevels(g),g,NULL,is_GRP_G(t),give.names)
#   } else {
#     if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
#     if(is.null(t)) {
#       message("Panel-lag computed without timevar: Assuming ordered data")
#       flagleadmCpp(x,n,fill,g[[1L]],g[[2L]],NULL,NULL,give.names)
#     } else if(is.atomic(t)) flagleadmCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],is_factor_G(t),give.names) else
#       flagleadmCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],is_GRP_G(t),give.names)
#   }
# }
# flag.grouped_df <- function(x, n = 1, t = NULL, fill = NULL, give.names = TRUE, drop.groups = FALSE, drop.t = TRUE, ...) {
#   g <- GRP.grouped_df(x)
#   tsym <- deparse(substitute(t))
#   nam <- names(x)
#   gn <- match(names(g[[4]]), nam)
#   gn2 = gn <- gn[!is.na(gn)]
#   if(!(tsym == "NULL" || is.na(tn <- match(tsym, nam)))) {
#     if(any(gn == tn)) stop("timevar coincides with grouping variables!")
#     t <- x[[tn]]
#     if(drop.t) if(!length(gn)) x[[tn]] <- NULL else if(drop.groups) gn <- c(gn,tn) else gn2 <- c(gn2,tn)
#   }
#   if(length(gn)) {
#     if(drop.groups)
#       return(flagleadlCpp(x[-gn],n,fill,g[[1L]],g[[2L]],g[[3L]],is_factor_G(t),give.names)) else {
#         ax <- attributes(x)
#         attributes(x) <- NULL
#         ax[["names"]] <- c(nam[gn], if(give.names) paste0("W.",nam[-gn2]) else nam[-gn2])
#         return(`attributes<-`(c(x[gn],flagleadlCpp(x[-gn2],n,fill,g[[1L]],g[[2L]],g[[3L]],is_factor_G(t),give.names)), ax))
#       }
#   } else return(flagleadlCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],is_factor_G(t),give.names))
# }
#
# flag.data.frame <- function(x, n = 1, g = NULL, t = NULL, fill = NULL, give.names = TRUE, ...) {
#   if(is.null(g)) {
#     if(is.null(t)) flagleadlCpp(x,n,fill, names = give.names) else if(is.atomic(t))
#       flagleadlCpp(x,n,fill, t = is_factor_G(t), names = give.names) else
#         flagleadlCpp(x,n,fill, t = is_GRP_G(t), names = give.names)
#   } else if(is.atomic(g)) {
#     if(!is.factor(g)) g = qF(g)
#     if(is.null(t)) {
#       message("Panel-lag computed without timevar: Assuming ordered data")
#       flagleadlCpp(x,n,fill,fnlevels(g),g,NULL,NULL,give.names)
#     } else if(is.atomic(t)) flagleadlCpp(x,n,fill,fnlevels(g),g,NULL,is_factor_G(t),give.names) else
#       flagleadlCpp(x,n,fill,fnlevels(g),g,NULL,is_GRP_G(t),give.names)
#   } else {
#     if(!is.GRP(g)) g = GRP(g, return.groups = FALSE)
#     if(is.null(t)) {
#       message("Panel-lag computed without timevar: Assuming ordered data")
#       flagleadlCpp(x,n,fill,g[[1L]],g[[2L]],NULL,NULL,give.names)
#     } else if(is.atomic(t)) flagleadlCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],is_factor_G(t),give.names) else
#       flagleadlCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],is_GRP_G(t),give.names)
#   }
# }
# flag.pdata.frame <- function(x, n = 1, fill = NA, give.names = TRUE, ...) {
#   index <- attr(x, "index")
#   if(length(index) > 2L) index <- c(interaction(index[-length(index)], drop = TRUE), index[length(index)])
#   flagleadlCpp(x,n,fill,fnlevels(index[[1L]]),index[[1L]],NULL,index[[2L]],give.names)
# }


# L with no dispatch!! gives minor speed gain, but Xcols feature is better !!, and calling method directls is also still faster !!
# L <- function(X, n = 1, xt = NULL, t = NULL, fill = NULL, ...) {
#  if(is.atomic(X)) {
#    if(is.matrix(X)) flag.matrix(X, n, xt, t, fill, ...) else flag.default(X, n, xt, t, fill)
#  } else flag.data.frame(X, n, xt, t, fill, ...)
# }



# This function still needs work: keep and drop.xt need to be implemented, and speed comparable to flag !!
# L.data.frame <- function(X, n = 1, xt = NULL, t = NULL, fill = NA, Xcols = is.numeric, keep = 1, drop.xt = FALSE, ...) {
#   isDT <- inherits(X, "data.table") # compare using unclass??
#   if(!is.null(xt)) if(!(is.factor(xt) || inherits(xt, "GRP"))) {
#       if(is.list(xt)) {
#         if(!all(lengths(xt) == nrow(X))) stop("lengths(xt) need to match nrow(X)")
#         xt <- GRP(xt, return.groups = FALSE)
#       } else if(is.atomic(xt)) {
#         if(length(xt) == nrow(X)) xt <- GRP(xt, return.groups = FALSE) else xt <- GRP(X, xt, return.groups = FALSE) # what about two-sided formula
#       } else if(is.call(xt)) xt <- GRP(X, xt, return.groups = FALSE) else stop("xt needs to be a list, vector or formula")
#   }
#   if(!is.null(t)) if(!(is.factor(t) || inherits(t, "GRP"))) {
#     if(is.list(t)) {
#       if(!all(lengths(t) == nrow(X))) stop("lengths(xt) need to match nrow(X)")
#       t <- GRP(t, return.groups = FALSE)
#     } else if(is.atomic(t)) {
#       if(length(t) == nrow(X)) t <- GRP(t, return.groups = FALSE) else t <- GRP(X, t, return.groups = FALSE) # what about two-sided formula
#     } else if(is.call(t)) t <- GRP(X, t, return.groups = FALSE) else stop("xt needs to be a list, vector or formula")
#   }
#   if(!is.null(Xcols)) { # fastest way ?? making another copy ??
#     if(is.function(Xcols))
#       X <- base::Filter(is.numeric, X) else if(is.atomic(Xcols)) {
#         if(isDT) X <- X[ , Xcols, with = FALSE] else X <- X[ , Xcols, drop = FALSE]
#       } else stop("Xcols needs to be a function, column names or column indices")
#   }
#   flag.data.frame(X, n, xt, t, fill, ...)
# }




# Previous Versions: Only numeric data, without internal attribute copy and name generation, g can only be factor or GRP !!

# L <- function(X, n = 1, xt = NULL, t = NULL, fill = NA, ...) {
#   UseMethod("L", X)
# }
# L.default <- function(X, n = 1, xt = NULL, t = NULL, fill = NA, ...) {
#   if(!is.null(xt)) if(!(is.factor(xt) || inherits(xt, "GRP"))) xt <- GRP(xt, return.groups = FALSE)
#   if(!is.null(t)) if(!(is.factor(t) || inherits(t, "GRP"))) t <- GRP(t, return.groups = FALSE)
#   flag.default(X, n, xt, t, fill, ...)
# }
# L.matrix <- function(X, n = 1, xt = NULL, t = NULL, fill = NA, ...) {
#   if(!is.null(xt)) if(!(is.factor(xt) || inherits(xt, "GRP"))) xt <- GRP(xt, return.groups = FALSE)
#   if(!is.null(t)) if(!(is.factor(t) || inherits(t, "GRP"))) t <- GRP(t, return.groups = FALSE)
#   flag.matrix(X, n, xt, t, fill, ...)
# }


# lagnames <- function(nam, n) {
#   pos <- n > 0L
#   neg <- n < 0L
#   res <- character(length(n))
#   res[pos] <- paste0("L",n[pos],".")
#   res[neg] <- paste0("F",abs(n[neg]),".")
#   # or res <- c(paste0("L",n[pos],"."), paste0("F",abs(n[neg]),".")) ??
#   paste0(rep(res, length(nam)), rep(nam, each = length(n)))
# }

# flag <- function(x, n = 1, g = NULL, t = NULL, fill = NA, ...) {
#   UseMethod("flag", x)
# }
# flag.default <- function(x, n = 1, g = NULL, t = NULL, fill = NA, ...) {
#   res <- if(is.null(g)) {
#     if(is.null(t)) flagleadCpp(x,n,fill) else { # best way to organize this code ??
#       if(is.factor(t)) flagleadCpp(x,n,fill,0L,0L,NULL,t) else if (inherits(t,"GRP"))
#         flagleadCpp(x,n,fill,0L,0L,NULL,t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#     }
#   } else { # best way to organize this code ??
#     if(is.factor(g)) {
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         flagleadCpp(x,n,fill,fnlevels(g),g)
#       } else {
#         if(is.factor(t)) flagleadCpp(x,n,fill,fnlevels(g),g,NULL,t) else if (inherits(t,"GRP"))
#           flagleadCpp(x,n,fill,fnlevels(g),g,NULL,t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     } else {
#       if(!inherits(g,"GRP")) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         flagleadCpp(x,n,fill,g[[1L]],g[[2L]])
#       } else {
#         if(is.factor(t)) flagleadCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t) else if (inherits(t,"GRP"))
#           flagleadCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     }
#   }
#   if(length(n) == 1) attributes(res) <- attributes(x) else dimnames(res) <- list(names(x), n) # lagn(n)) # lagn(n)) # # paste0(ifelse(n>=0,"L","F"),abs(n)) # fastest way ?? or ax["dim"] <- attr(res,"dim") # list(as.character(ax[["names"]]), paste0(ifelse(n>=0,"L","F"),abs(n)))
#   res
# }
# flag.matrix <- function(x, n = 1, g = NULL, t = NULL, fill = NA, ...) {
#   res <- if(is.null(g)) {
#     if(any(abs(n) > nrow(x))) stop("lag-length exceeds nrow(x)") # efficient ?? -> yes !!
#     if(is.null(t)) flagleadmCpp(x,n,fill) else { # best way to organize this code ??
#       if(is.factor(t)) flagleadmCpp(x,n,fill,0L,0L,NULL,t) else if (inherits(t,"GRP"))
#         flagleadmCpp(x,n,fill,0L,0L,NULL,t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#     }
#   } else { # best way to organize this code ??
#     if (is.factor(g)) {
#       if(any(abs(n) > nrow(x)/fnlevels(g))) stop("lag-length exceeds average group size")
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         flagleadmCpp(x,n,fill,fnlevels(g),g)
#       } else {
#         if(is.factor(t)) flagleadmCpp(x,n,fill,fnlevels(g),g,NULL,t) else if (inherits(t,"GRP"))
#           flagleadmCpp(x,n,fill,fnlevels(g),g,NULL,t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     } else {
#       if(!inherits(g,"GRP")) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(any(abs(n) > nrow(x)/g[[1L]])) stop("lag-length exceeds average group size")
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         flagleadmCpp(x,n,fill,g[[1L]],g[[2L]])
#       } else {
#         if(is.factor(t)) flagleadmCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t) else if (inherits(t,"GRP"))
#           flagleadmCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     }
#   }
#   if(length(n)>1) {
#     dx <- dimnames(x)
#     dx[[2L]] <- lagnames(dx[[2L]], n)
#     dimnames(res) <- dx
#   } else dimnames(res) <- dimnames(x)
#   res
# }
# flag.data.frame <- function(x, n = 1, g = NULL, t = NULL, fill = NA, ...) {
#   ax <- attributes(x)
#   res <- if(is.null(g)) {
#     if(is.null(t)) flagleadlCpp(x,n,fill) else { # best way to organize this code ??
#       if(is.factor(t)) flagleadlCpp(x,n,fill,0L,0L,NULL,t) else if (inherits(t,"GRP"))
#         flagleadlCpp(x,n,fill,0L,0L,NULL,t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#     }
#   } else { # best way to organize this code ??
#     if (is.factor(g)) {
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         flagleadlCpp(x,n,fill,fnlevels(g),g)
#       } else {
#         if(is.factor(t)) flagleadlCpp(x,n,fill,fnlevels(g),g,NULL,t) else if (inherits(t,"GRP"))
#           flagleadlCpp(x,n,fill,fnlevels(g),g,NULL,t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     } else {
#       if(!inherits(g,"GRP")) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         flagleadlCpp(x,n,fill,g[[1L]],g[[2L]])
#       } else {
#         if(is.factor(t)) flagleadlCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t) else if (inherits(t,"GRP"))
#           flagleadlCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t[[2L]]) else stop("t must be a a factor, or a GRP() object, see ?GRP")
#       }
#     }
#   }
#   if(length(n)>1) ax[["names"]] <- lagnames(ax[["names"]], n)
#   attributes(res) <- ax
#   res
# }



# lagn <- function(n) {
#   pos <- n > 0
#   res <- character(length(n))
#   res[pos] <- paste0(".L",n[pos])
#   res[!pos] <- paste0(".F",abs(n[!pos]))
#   res[n == 0] <- ".--"
#   res
# } # Not really needed, but aesthetically pleasing !!


# flag.list <- function(x, n = 1, g = NULL, t = NULL, fill = NA, ...) {
#   ax <- attributes(x)
#   res <- if(is.null(g)) flagleadlCpp(x,n,fill) else { # best way to organize this code ??
#     if (is.factor(g)) {
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         flagleadlCpp(x,n,fill,fnlevels(g),g)
#       } else {
#         if(is.factor(t)) flagleadlCpp(x,n,fill,fnlevels(g),g,.Internal(tabulate(g,fnlevels(g))),t) else if (inherits(t,"GRP"))
#         flagleadlCpp(x,n,fill,fnlevels(g),g,.Internal(tabulate(g,fnlevels(g))),t[[2L]]) else stop("g must be a a factor, or a GRP() object, see ?GRP")
#       }
#     } else {
#       if(!inherits(g,"GRP")) stop("g must be a a factor, or a GRP() object, see ?GRP")
#       if(is.null(t)) {
#         warning("Panel-lag computed without timevar: Assuming ordered data")
#         flagleadlCpp(x,n,fill,g[[1L]],g[[2L]])
#       } else {
#         if(is.factor(t)) flagleadlCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t) else if (inherits(t,"GRP"))
#         flagleadlCpp(x,n,fill,g[[1L]],g[[2L]],g[[3L]],t[[2L]]) else stop("g must be a a factor, or a GRP() object, see ?GRP")
#       }
#     }
#   }
#   attributes(res) <- ax
#   res
# } # use rapply2d !!!
