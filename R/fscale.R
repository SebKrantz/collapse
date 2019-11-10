library(Rcpp)
sourceCpp('R/C++/fscale.cpp', rebuild = TRUE)
sourceCpp('R/C++/fscalea.cpp', rebuild = TRUE)
sourceCpp('R/C++/fscalel.cpp', rebuild = TRUE)
source("R/collapse R/GRP.R")
source("R/collapse R/small_helper.R")
source("R/collapse R/quick_conversion.R")

# w.type = "frequency"
# Todo: center and scale arguments -> link to fmean and fsd with TRA ??

fscale <- function(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) {
  UseMethod("fscale", x)
}
fscale.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) { 
  if(is.null(g)) return(fscaleCpp(x,0L,0L,NULL,w,na.rm,stable.algo)) else if (is.atomic(g)) {
    if(is.factor(g)) return(fscaleCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo)) else {
      g <- qG(g, ordered = FALSE)
      return(fscaleCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo))
    }
  } else {
    if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
    return(fscaleCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo))
  }
}
fscale.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) { 
  g <- if(length(effect) == 1L) attr(x, "index")[[effect]] else interaction(attr(x, "index")[effect], drop = TRUE)
  fscaleCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo) 
}
fscale.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) { 
  if(is.null(g)) return(fscalemCpp(x,0L,0L,NULL,w,na.rm,stable.algo)) else if (is.atomic(g)) {
    if(is.factor(g)) return(fscalemCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo)) else {
      g <- qG(g, ordered = FALSE)
      return(fscalemCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo))
    }
  } else {
    if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
    return(fscalemCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo))
  }
}
fscale.grouped_df <- function(x, w = NULL, na.rm = TRUE, stable.algo = TRUE, keep.groups = TRUE, keep.w = TRUE, ...) { 
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- names(x)
  gn2 <- match(names(g[[4]]), nam) 
  gn2 <- gn2[!is.na(gn2)]
  gn <- if(keep.groups) gn2 else NULL
  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- x[[wn]] 
    if(any(gn2 == wn)) stop("Weights coincide with grouping variables!")
    gn2 <- c(gn2,wn)
    if(keep.w) gn <- c(gn,wn) 
  }
  if(length(gn2)) {
    if(!length(gn)) 
      return(fscalelCpp(x[-gn2],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo)) else {
        ax <- attributes(x)
        attributes(x) <- NULL 
        ax[["names"]] <- c(nam[gn], nam[-gn2])
        return(`attributes<-`(c(x[gn],fscalelCpp(x[-gn2],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo)), ax))
      }
  } else return(fscalelCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo))
}
fscale.data.frame <- function(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) { 
  if(is.null(g)) return(fscalelCpp(x,0L,0L,NULL,w,na.rm,stable.algo)) else if (is.atomic(g)) {
    if(is.factor(g)) return(fscalelCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo)) else {
      g <- qG(g, ordered = FALSE)
      return(fscalelCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo))
    }
  } else {
    if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
    return(fscalelCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo))
  }
}
fscale.pdata.frame <- function(x, effect = 1L, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) {
  g <- if(length(effect) == 1L) attr(x, "index")[[effect]] else interaction(attr(x, "index")[effect], drop = TRUE)
  fscaleCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo) 
}


# Todo: Could still implement return = c("cols","all","add")

STD <- function(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) {
  UseMethod("STD", x)
}
STD.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) {
  fscale.default(x, g, w, na.rm, stable.algo)
}
STD.pseries <- function(x, effect = 1L, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) { 
  g <- if(length(effect) == 1L) attr(x, "index")[[effect]] else interaction(attr(x, "index")[effect], drop = TRUE)
  fscaleCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo) 
}
STD.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, give.names = TRUE, ...) {
  res <- fscale.matrix(x, g, w, na.rm, stable.algo)
  if(give.names) dimnames(res)[[2]] <- paste0("STD.",dimnames(res)[[2]])
  return(res)
}
STD.grouped_df <- function(x, w = NULL, na.rm = TRUE, stable.algo = TRUE, give.names = TRUE, keep.groups = TRUE, keep.w = TRUE, ...) { 
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- names(x)
  gn2 <- match(names(g[[4]]), nam) 
  gn2 <- gn2[!is.na(gn2)]
  gn <- if(keep.groups) gn2 else NULL
  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- x[[wn]] 
    if(any(gn2 == wn)) stop("Weights coincide with grouping variables!")
    gn2 <- c(gn2,wn)
    if(keep.w) gn <- c(gn,wn) 
  }
  if(length(gn2)) {
    if(!length(gn)) 
      return(give_nam(fscalelCpp(x[-gn2],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo), give.names, "STD.")) else {
        ax <- attributes(x)
        attributes(x) <- NULL 
        ax[["names"]] <- c(nam[gn], if(give.names) paste0("STD.",nam[-gn2]) else nam[-gn2])
        return(`attributes<-`(c(x[gn],fscalelCpp(x[-gn2],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo)), ax))
      }
  } else return(give_nam(fscalelCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo), give.names, "STD."))
}

# updated (best) version !!
STD.pdata.frame <- function(x, effect = 1L, w = NULL, cols = is.numeric, 
                            na.rm = TRUE, stable.algo = TRUE, give.names = TRUE, 
                            keep.ids = TRUE, keep.w = TRUE, ...) { 
  
  ax <- attributes(x)
  class(x) <- NULL
  nam <- ax[["names"]]
  g <- if(length(effect) == 1L) ax[["index"]][[effect]] else 
       interaction(ax[["index"]][effect], drop = TRUE)

  if(keep.ids) {
    gn <- match(names(ax[["index"]]), nam)
    gn <- gn[!is.na(gn)]
    if(length(gn) && is.null(cols)) cols <- seq_along(x)[-gn] 
  } else gn <- NULL
  
  if(!is.null(cols)) cols <- cols2int(cols, x, nam)

  if(is.call(w)) {
    w <- all.vars(w)
    wn <- match(w, nam)
    w <- x[[wn]]
    cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
    if(keep.w) gn <- c(gn, wn)
  } 
  
  if(length(gn) && !is.null(cols)) {
    ax[["names"]] <- c(nam[gn], if(give.names) paste0("STD.", nam[cols]) else nam[cols])
    return(`attributes<-`(c(x[gn], fscalelCpp(x[cols],fnlevels(g),g,NULL,w,na.rm,stable.algo)), ax))
  } else if(!length(gn)) {
    ax[["names"]] <- if(give.names) paste0("STD.", nam[cols]) else nam[cols]
    return(`attributes<-`(fscalelCpp(x[cols],fnlevels(g),g,NULL,w,na.rm,stable.algo), ax)) 
  } else {
    if(give.names) {
      ax[["names"]] <- paste0("STD.", nam)
      return(`attributes<-`(fscalelCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo), ax))
    } else
    return(fscalelCpp(`oldClass<-`(x, ax[["class"]]),fnlevels(g),g,NULL,w,na.rm,stable.algo))
  }
}
# updated, fast and data.table proof version !!! 
STD.data.frame <- function(x, by = NULL, w = NULL, cols = is.numeric, na.rm = TRUE, 
                           stable.algo = TRUE, give.names = TRUE, 
                           keep.by = TRUE, keep.w = TRUE, ...) {
  
  # fastest solution??
  if(is.call(by) || is.call(w)) {
    ax <- attributes(x)
    class(x) <- NULL
    nam <- ax[["names"]]
    
    if(is.call(by)) {
      if(length(by) == 3L) {
        cols <- match(all.vars(by[[2L]]), nam)
        gn <- match(all.vars(by[[3L]]), nam)
      } else {
        gn <- match(all.vars(by), nam)
        cols <- if(is.null(cols)) seq_along(x)[-gn] else cols2int(cols, x, nam)
      }
      by <- if(length(gn) == 1L) at2GRP(x[[gn]]) else GRP(x, gn, return.groups = FALSE)
      if(!keep.by) gn <- NULL
    } else {
      gn <- NULL
      if(!is.null(cols)) cols <- cols2int(cols, x, nam)
      if(!is.GRP(by)) by <- if(is.null(by)) list(0L, 0L, NULL) else if(is.atomic(by)) # Necessary for if by is passed externally !!
                            at2GRP(by) else GRP(by, return.groups = FALSE)
    }
    
    if(is.call(w)) {
      w <- all.vars(w)
      wn <- match(w, nam)
      w <- x[[wn]]
      cols <- if(is.null(cols)) seq_along(x)[-wn] else cols[cols != wn]
      if(keep.w) gn <- c(gn, wn)
    } 
    
    if(length(gn)) {
      ax[["names"]] <- c(nam[gn], if(give.names) paste0("STD.",nam[cols]) else nam[cols])
      return(`attributes<-`(c(x[gn], fscalelCpp(x[cols],by[[1]],by[[2]],by[[3]],w,na.rm,stable.algo)), ax))
    } else {
      ax[["names"]] <- if(give.names) paste0("STD.",nam[cols]) else nam[cols]
      return(`attributes<-`(fscalelCpp(x[cols],by[[1]],by[[2]],by[[3]],w,na.rm,stable.algo), ax))
    }
  } else if(!is.null(cols)) {
    ax <- attributes(x)
    x <- if(is.function(cols)) unclass(x)[vapply(x, cols, TRUE)] else unclass(x)[cols]
    ax[["names"]] <- names(x)
    attributes(x) <- ax
  } 
  if(give.names) names(x) <- paste0("STD.", names(x))
  
  if(is.null(by)) return(fscalelCpp(x,0L,0L,NULL,w,na.rm,stable.algo)) else if (is.atomic(by)) {
    if(is.factor(by)) return(fscalelCpp(x,fnlevels(by),by,NULL,w,na.rm,stable.algo)) else {
      by <- qG(by, ordered = FALSE)
      return(fscalelCpp(x,attr(by,"N.groups"),by,NULL,w,na.rm,stable.algo))
    }
  } else {
    if(!is.GRP(by)) by <- GRP(by, return.groups = FALSE)
    return(fscalelCpp(x,by[[1]],by[[2]],by[[3]],w,na.rm,stable.algo))
  }
}





# Previous Verions: Also allowing for names and indices !! Before making things more efficient following qsu enhancements !!
# STD.pdata.frame <- function(X, effect = 1L, w = NULL, Xcols = is.numeric, na.rm = TRUE, stable.algo = TRUE, give.names = FALSE, drop.xt = FALSE, drop.w = TRUE, ...) { 
#   ax <- attributes(X)
#   nam <- ax[["names"]]
#   gn <- match(names(ax[["index"]]), nam)
#   gn2 = gn <- gn[!is.na(gn)]
#   g <- if(length(effect) == 1L) ax[["index"]][[effect]] else interaction(ax[["index"]][effect], drop = TRUE)
#   if(!is.null(Xcols)) {
#     if(is.function(Xcols)) Xcols <- seq_along(X)[!vapply(X, Xcols, TRUE)] else if(is.character(Xcols)) 
#       Xcols <- seq_along(X)[-match(Xcols, nam)] else if(is.logical(Xcols)) 
#         Xcols <- seq_along(X)[!Xcols] else if(is.numeric(Xcols)) 
#           Xcols <- seq_along(X)[-Xcols] else stop("Xcols needs to be a function, column names, indices or a logical vector")
#         if(drop.xt) gn <- c(gn, Xcols) else if(length(gn)) gn2 <- c(gn2, Xcols) else gn <- Xcols 
#   }
#   if(!is.null(w)) {
#     if(is.call(w)) w <- all.vars(w) 
#     if(length(w) == 1) {
#       wn <- match(w, nam)
#       if(any(gn == wn)) stop("Weights coincide with grouping variables!")
#       w <- X[[wn]]
#       if(drop.w) if(!length(gn)) X[[wn]] <- NULL else if(drop.xt) gn <- c(gn,wn) else gn2 <- c(gn2,wn)
#     }
#   }
#   if(length(gn)) {
#     if(drop.xt) 
#       return(give_nam(fscalelCpp(X[-gn],fnlevels(g),g,NULL,w,na.rm,stable.algo), give.names, "STD.")) else {
#         attributes(X) <- NULL 
#         ax[["names"]] <- c(nam[gn], if(give.names) paste0("STD.",nam[-gn2]) else nam[-gn2])
#         return(`attributes<-`(c(X[gn],fscalelCpp(X[-gn2],fnlevels(g),g,NULL,w,na.rm,stable.algo)), ax))
#       }
#   } else return(give_nam(fscalelCpp(X,fnlevels(g),g,NULL,w,na.rm,stable.algo), give.names, "STD."))
# }
# STD.data.frame <- function(X, by = NULL, w = NULL, Xcols = is.numeric, na.rm = TRUE, 
#                            stable.algo = TRUE, give.names = FALSE, drop.by = FALSE, drop.w = TRUE, ...) { 
#   if(!is.null(w)) {
#     if(is.call(w)) w <- all.vars(w) 
#     if(length(w) == 1) {
#       v <- w
#       w <- X[[w]]
#       if(drop.w) X[[v]] <- NULL
#     }
#   }
#   if(is.null(by)) {
#     if(is.null(Xcols)) return(give_nam(fscalelCpp(X,0L,0L,NULL,w,na.rm,stable.algo), give.names, "STD.")) else {
#       if(is.function(Xcols)) Xcols <- vapply(X, Xcols, TRUE) 
#       return(give_nam(fscalelCpp(X[Xcols],0L,0L,NULL,w,na.rm,stable.algo), give.names, "STD."))
#     }
#   } else if (is.atomic(by) || is.call(by)) {
#     if(is.call(by) || length(by) != nrow(X)) {
#       nam <- names(X)
#       if(is.call(by) && length(by) == 3) {
#         v <- match(all.vars(by[[2]]), nam)
#         by <- match(all.vars(by[[3]]), nam)
#       } else {
#         if(is.call(by)) by <- match(all.vars(by), nam) else if(is.character(by)) by <- match(by, nam) 
#         v <- if(is.null(Xcols)) seq_along(X)[-by] else if(is.function(Xcols)) 
#           setdiff(which(vapply(X, Xcols, TRUE)), by) else if(is.character(Xcols)) 
#             match(Xcols, nam) else Xcols
#       }
#       g <- GRP(X, by, return.groups = FALSE)
#       if(drop.by) 
#         return(give_nam(fscalelCpp(X[v],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo), give.names, "STD.")) else {
#           ax <- attributes(X)
#           ax[["names"]] <- c(nam[by], if(give.names) paste0("STD.", nam[v]) else nam[v]) 
#           attributes(X) <- NULL
#           return(`attributes<-`(c(X[by],fscalelCpp(X[v],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo)), ax))
#         }
#     } else if(is.factor(by)) return(give_nam(fscalelCpp(X,fnlevels(by),by,NULL,w,na.rm,stable.algo), give.names, "STD.")) else {
#       by <- qG(by, ordered = FALSE)
#       return(give_nam(fscalelCpp(X,attr(by,"N.groups"),by,NULL,w,na.rm,stable.algo), give.names, "STD."))
#     }
#   } else {
#     if(!is.GRP(by)) by <- GRP(by, return.groups = FALSE)
#     if(is.null(Xcols)) return(give_nam(fscalelCpp(X,by[[1]],by[[2]],by[[3]],w,na.rm,stable.algo), give.names, "STD.")) else {
#       if(is.function(Xcols)) Xcols <- vapply(X, Xcols, TRUE) 
#       return(give_nam(fscalelCpp(X[Xcols],by[[1]],by[[2]],by[[3]],w,na.rm,stable.algo), give.names, "STD."))
#     }
#   }
# }







# OLD: using tabulate for factor group sizes, and before Xcols and formula for data.frame
# fscale <- function(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) {
#   UseMethod("fscale", x)
# }
# fscale.default <- function(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) { 
#   if(is.null(g)) fscaleCpp(x,0L,0L,0L,w,na.rm,stable.algo) else {
#     if(is.atomic(g)) {
#       if(is.factor(g)) ng <- nlevels(g) else {
#         g <- qG(g)
#         ng <- attr(g, "N.groups")
#       }
#       gs <- if(!stable.algo && is.null(w) && !na.rm) .Internal(tabulate(g, ng)) else 0L        
#       fscaleCpp(x,ng,g,gs,w,na.rm,stable.algo)
#     } else {
#       if(is.GRP(g))  
#       fscaleCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo) else if(.Internal(islistfactor(g, FALSE))) {
#         g <- interaction(g)
#         ng <- nlevels(g)
#         gs <- if(!stable.algo && is.null(w) && !na.rm) .Internal(tabulate(g, ng)) else 0L        
#         fscaleCpp(x,ng,g,gs,w,na.rm,stable.algo)
#       } else {
#         g <- GRP(g, return.groups = FALSE)
#         fscaleCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo)
#       } 
#     }
#   }
# }
# fscale.pseries <- function(x, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) { # good ?? 
#   g <- attr(x, "index")[[1]]
#   nl <- nlevels(g)
#   gs <- if(!stable.algo && is.null(w) && !na.rm) .Internal(tabulate(g, nl)) else 0L
#   fscaleCpp(x,nl,g,gs,w,na.rm,stable.algo) 
# }
# fscale.matrix <- function(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) { 
#   if(is.null(g)) fscalemCpp(x,0L,0L,0L,w,na.rm,stable.algo) else {
#     if(is.atomic(g)) {
#       if(is.factor(g)) ng <- nlevels(g) else {
#         g <- qG(g)
#         ng <- attr(g, "N.groups")
#       }
#       gs <- if(!stable.algo && is.null(w) && !na.rm) .Internal(tabulate(g, ng)) else 0L        
#       fscalemCpp(x,ng,g,gs,w,na.rm,stable.algo)
#     } else {
#       if(is.GRP(g))  
#         fscalemCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo) else if(.Internal(islistfactor(g, FALSE))) {
#           g <- interaction(g)
#           ng <- nlevels(g)
#           gs <- if(!stable.algo && is.null(w) && !na.rm) .Internal(tabulate(g, ng)) else 0L        
#           fscalemCpp(x,ng,g,gs,w,na.rm,stable.algo)
#         } else {
#           g <- GRP(g, return.groups = FALSE)
#           fscalemCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo)
#         } 
#     }
#   }
# }
# fscale.data.frame <- function(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) { 
#   if(is.null(g)) fscalelCpp(x,0L,0L,0L,w,na.rm,stable.algo) else {
#     if(is.atomic(g)) {
#       if(is.factor(g)) ng <- nlevels(g) else {
#         g <- qG(g)
#         ng <- attr(g, "N.groups")
#       }
#       gs <- if(!stable.algo && is.null(w) && !na.rm) .Internal(tabulate(g, ng)) else 0L        
#       fscalelCpp(x,ng,g,gs,w,na.rm,stable.algo)
#     } else {
#       if(is.GRP(g))  
#         fscalelCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo) else if(.Internal(islistfactor(g, FALSE))) {
#           g <- interaction(g)
#           ng <- nlevels(g)
#           gs <- if(!stable.algo && is.null(w) && !na.rm) .Internal(tabulate(g, ng)) else 0L        
#           fscalelCpp(x,ng,g,gs,w,na.rm,stable.algo)
#         } else {
#           g <- GRP(g, return.groups = FALSE)
#           fscalelCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo)
#         } 
#     }
#   }
# }
# fscale.pdata.frame <- function(x, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...) { 
#   g <- attr(x, "index")[[1]]
#   nl <- nlevels(g)
#   gs <- if(!stable.algo && is.null(w) && !na.rm) .Internal(tabulate(g, nl)) else 0L
#   fscalelCpp(x,nl,g,gs,w,na.rm,stable.algo) 
# }
# fscale.grouped_df <- function(x, w = NULL, na.rm = TRUE, drop.groups = FALSE, stable.algo = TRUE, ...) {
#   g <- GRP.grouped_df(x)
#   gn <- match(names(g[[4]]), names(x))
#   gn <- gn[!is.na(gn)]
#   if(length(gn)) {
#     if(drop.groups) fscalelCpp(x[-gn],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo) else {
#       clx <- class(x)
#       class(x) <- NULL
#       x[-gn] <- fscalelCpp(x[-gn],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo) 
#       class(x) <- clx
#       x
#     }
#   } else fscalelCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo)
# }
