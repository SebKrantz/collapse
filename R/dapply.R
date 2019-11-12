library(Rcpp)
sourceCpp("C:/Users/Sebastian Krantz/Documents/R/mrtl_type_dispatch.cpp")
# todo, Make setdapply(), or Xcols option??. Should be consistent with the rest. all set-functions should have copy option, or: if assigned to something make copy, else no!!
# definitely make add option for MARGIN = 1?? -> nah, rather setdapply, with Xcols and Add option
# -> Do just like B and W, both dapply and setdapply, and each has Xcols and add = 0,1,2 option.
# Note: in setdapply, 0 just means replace columns?? -> yes!!

# This version does not always return data.frame !!, but currently needs properly defined matrices and data.fram input !!
dapply3 <- function(X, FUN, ..., MARGIN = 2, parallel = FALSE,
                    mc.cores = 1L, return = c("same","matrix","data.frame")) {
  ax <- attributes(X)
  arl <- is.array(X)
  rowwl <- MARGIN == 1
  retmatl <- switch(return[1L], same = arl, matrix = TRUE, data.frame = FALSE, stop("Unknown return option!"))
  aplyfun <- if(parallel) function(...) parallel::mclapply(..., mc.cores = mc.cores) else lapply
  if(arl) {
    dX <- dim(X)
    if(length(dX) > 2L) stop("dapply cannot handle higher-dimensional arrays")
    if(!retmatl) { # without checking is.null(ax), can only input matrices and data.frames !!
      dn <- dimnames(X) # faster than ax[["dimnames"]] !!
      ax <- c(list(names = dn[[2L]], row.names = dn[[1L]], class = "data.frame"),
              ax[!(names(ax) %in% c("dim","dimnames","class"))])
    }
    X <- if(rowwl) aplyfun(mrtl(X), FUN, ...) else aplyfun(mctl(X), FUN, ...)
  } else {
    attributes(X) <- NULL
    dX <- c(length(X[[1L]]), length(X))
    X <- if(rowwl) aplyfun(mrtl(do.call(cbind, X)), FUN, ...) else aplyfun(X, FUN, ...)
    if(retmatl) ax <- c(list(dim = dX, dimnames = list(ax[["row.names"]], ax[["names"]])),
                        ax[!(names(ax) %in% c("names","row.names","class"))])
  }
  lx1 <- length(X[[1L]])
  # diffl <- lx1 != switch(MARGIN, dX[2L], dX[1L]) # further input check ??
  # if(diffl) nx1 <- names(X[[1L]])

  if(retmatl) {
    if(rowwl) {
      if(lx1 != dX[2L] && !is.null(nx1 <- names(X[[1L]]))) ax[["dimnames"]][[2L]] <- nx1
      X <- matrix(unlist(X, use.names = FALSE), ncol = lx1, byrow = TRUE)
    } else {
      if(lx1 != dX[1L] && !is.null(nx1 <- names(X[[1L]]))) ax[["dimnames"]][[1L]] <- nx1
      X <- do.call(cbind, X)
    }
  } else {
    if(rowwl) {
      if(lx1 != dX[2L] && !is.null(nx1 <- names(X[[1L]]))) ax[["names"]] <- nx1
      X <- mctl(matrix(unlist(X, use.names = FALSE), ncol = lx1, byrow = TRUE))
    } else if(lx1 != dX[1L] && !is.null(nx1 <- names(X[[1L]]))) ax[["row.names"]] <- nx1
  }
  return(setAttributes(X, ax))
}


# adding parallelism: No speed loss !!
dapply <- function(X, FUN, ..., MARGIN = 2, parallel = FALSE, mc.cores = 1L) {
  ax <- attributes(X)
  aplyfun <- if(parallel) function(...) parallel::mclapply(..., mc.cores = mc.cores) else lapply
  if(is.array(X)) {
    dX <- dim(X) # ax[["dim"]] ?? -> nope, slower !!
    if(length(dX) > 2L) stop("dapply cannot handle higher-dimensional arrays")
    lXo <- dX[1L]
    if(MARGIN == 1) {
      X <- aplyfun(mrtl(X), FUN, ...)
      lx1 <- length(X[[1L]])
      nx1 <- names(X[[1L]])
      X <- mctl(matrix(unlist(X, use.names = FALSE), ncol = lx1, byrow = TRUE)) # data.table transposelist?? -> slower!!
      if(lx1 != dX[2L]) {
        ax[["names"]] <- if(!is.null(nx1)) nx1 else if(lx1 == 1L)
          deparse(substitute(FUN)) else paste(deparse(substitute(FUN)), seq_len(lx1), sep = ".")
      } else ax[["names"]] <- ax[["dimnames"]][[2L]]
    } else {
      X <- aplyfun(mctl(X), FUN, ...)
      ax[["names"]] <- ax[["dimnames"]][[2L]]
    }
    ax[["row.names"]] <- ax[["dimnames"]][[1L]]
    ax[c("dim","dimnames")] <- NULL
  } else {
    attributes(X) <- NULL # faster !!
    lXo <- length(X[[1L]])
    if(MARGIN == 1L) {
      lX <- length(X)
      X <- aplyfun(mrtl(do.call(cbind, X)), FUN, ...) # X <- mrtl(matrix(unlist(X, use.names = FALSE), ncol = lX)) # slower!
      lx1 <- length(X[[1L]])
      nx1 <- names(X[[1L]])
      X <- mctl(matrix(unlist(X, use.names = FALSE), ncol = lx1, byrow = TRUE)) # mctl(do.call(rbind, X)) # But this is lower (GGDC) !!
      if(lx1 != lX) ax[["names"]] <- if(!is.null(nx1)) nx1 else if(lx1 == 1L)
        deparse(substitute(FUN)) else paste(deparse(substitute(FUN)), seq_len(lx1), sep = ".")
    } else X <- aplyfun(X, FUN, ...)
  }
  lXn <- length(X[[1L]])
  if(is.null(ax)) {
    class(X) <- "data.frame"
    attr(X, "row.names") <- .set_row_names(lXn)
  } else {
    if(is.null(ax[["row.names"]]) || lXo != lXn) ax[["row.names"]] <- .set_row_names(lXn)
    if(!any(ax[["class"]] == "data.frame")) ax[["class"]] <- c(ax[["class"]], "data.frame") # right order??
    setattributes(X, ax) # a lot faster for large data !! # attributes(X) <- ax is a tiny bit faster for small data !!
  }
  return(X)
}


# # Idea of implementing Xcols, but not really necessary I think!! data.table can do that!! this package should not substitute for data.table
# # also B() can do this, so you don't need dapply to do it! keep this function simple, it is already perfect the way it is, otherwise the programming gets really messy!!
# # also consider that B() is based on data.table, so it will be faster for those kinds of adding column operations!!
# indXcols <- function(X, xc, il) { # Faster way??
#   if (il) {
#     if (is.function(xc)) which(unlist(lapply(X,xc), use.names = FALSE)) else if (is.numeric(xc))
#       xc else if (is.character(xc)) match(xc,names(X)) else # could do pmatch
#         stop("Xcols needs to be a function, or a vector of colnames or indices")
#   } else {
#     if (is.function(xc)) stop("Cannot subset array using a function") else if (is.numeric(xc))
#       xc else if (is.character(xc)) match(xc,colnames(X)) else # could do pmatch
#         stop("Xcols needs to be a function, or a vector of colnames or indices")
#   }
# }
# dapply <- function(X, FUN, ..., MARGIN = 2, Xcols = NULL, return = 0L) {
#   # for efficiency, do first if by return option!!. perhaps even define dapply as below, and then do the different options. But tbh you could also just screw Xcols and return for this function, you can simple replace stuff,a nd for all else use data.table
#   ar = is.array(X)
#   if (ar) {
#     dX = dim(X)
#     if (length(dX)>2) stop("dapply cannot handle higher-dimensional arrays")
#     if (length(Xcols)) X = X[, indXcols(X,Xcols,FALSE), drop = FALSE] # no adding columns to arrays!! ? or convert to DF..-> do that, same treatment as X
#     ax = attributes(X) # good here??
#     lXo = dX[1]
#     if (MARGIN == 1) {
#       X = lapply(mrtl(X), FUN, ...)
#       lx1 = length(X[[1]])
#       nx1 = names(X[[1]])
#       X = mctl(matrix(unlist(X, use.names = FALSE), ncol = lx1, byrow = TRUE))
#       if (lx1 != dX[2]) {
#         ax[["names"]] = if (!is.null(nx1)) nx1 else if (lx1==1)
#           deparse(substitute(FUN)) else
#             paste(deparse(substitute(FUN)), seq_len(lx1), sep = ".")
#       } else ax[["names"]] = ax[["dimnames"]][[2]]
#     } else {
#       X = lapply(mctl(X), FUN, ...)
#       ax[["names"]] = ax[["dimnames"]][[2]]
#     }
#     ax[["row.names"]] = ax[["dimnames"]][[1]]
#     ax[c("dim","dimnames")] = NULL
#   } else {
#     if (length(Xcols)) if (return == 0L) X = X[indXcols(X,Xcols,FALSE)] else {
#       Y = X; ind = indXcols(X,Xcols,FALSE); X = X[ind] # faster? wrap into function?? do without copy?? i.e. just take apart and reorder later?? See what is faster with large data
#     }
#     ax = attributes(X) # good here??
#     lXo = length(X[[1]])
#     if (MARGIN == 1) {
#       lX = length(X)
#       X = mrtl(matrix(unlist(X, use.names = FALSE), ncol = lX))
#       X = lapply(X, FUN, ...)
#       lx1 = length(X[[1]])
#       nx1 = names(X[[1]])
#       X = mctl(matrix(unlist(X, use.names = FALSE), ncol = lx1, byrow = TRUE))
#       if (lx1 != lX) ax[["names"]] = if (!is.null(nx1)) nx1 else if (lx1==1)
#         deparse(substitute(FUN)) else
#           paste(deparse(substitute(FUN)), seq_len(lx1), sep = ".")
#     } else
#       X = lapply(X, FUN, ...)
#   }
#   lXn = length(X[[1]])
#   if (is.null(ax)) {
#     class(X) <- "data.frame"
#     attr(X, "row.names") <- .set_row_names(lXn)
#   } else {
#     if (is.null(ax[["row.names"]]) || lXo != lXn) ax[["row.names"]] <- .set_row_names(lXn)
#     if (!any(ax[["class"]] == "data.frame")) ax[["class"]] <- c(ax[["class"]],"data.frame") # right order??
#     attributes(X) = ax
#   }
#   if (length(Xcols) && !ar) {
#     if (return == 0L) X  else if (return == 1L) {Y[ind] = X; Y} else if (return == 2L)
#       cbind(Y,setNames(X,paste(deparse(substitute(FUN)),names(X), sep ="."))) # faster? wrap into function?? do without copy?? i.e. just take apart and reorder later?? See what is faster with large data
#   } else X
# }
#


# Old version:
# DO new rownames??
# dapply <- function(X, FUN, ..., MARGIN = 2) { # , possibly add margins argument?? -> nah..
#   # This if would disable seq_along() calls..
#   #if (!is.list(X)) stop("dapply only takes data.frames or lists that can be coerced to data.frame")
#   ax = attributes(X)
#   if (is.array(X)) {
#     #  lXo = nrow(X) # right?? what about HD Arrays??
#     # if (MARGIN != 2) {
#     #   dX = dim(X)
#     #   X = drop(t(apply(X, MARGIN, FUN, ...))) # See View(apply), take central part!!
#     #   dX2 = dim(X)
#     #   lXo = NROW(X) # right?? what about HD Arrays??
#     #   if (is.null(dX2) || !all(dX2==dX)) { # better solution??
#     #     ax = attributes(X) # necesary?? This will not coerce to data.table again...
#     #     n = names(ax) == "names"
#     #     if (any(n)) names(ax)[n] = "row.names"
#     #   }
#     #   #X = as.data.frame(X) # as.data.frame.array faster way??, separate out??
#     #   X = if (is.null(dX2)) list(as.vector(X)) else lapply(1:ncol(X), function(i) as.vector(X[,i]))
#     dX = dim(X)
#     if (length(dX)>2) stop("dapply cannot handle higher-dimensional arrays")
#     lXo = dX[1]
#     if (MARGIN == 1) {
#       X = lapply(mrtl(X), FUN, ...)
#       lx1 = length(X[[1]])
#       nx1 = names(X[[1]])
#       X = mctl(matrix(unlist(X, use.names = FALSE), ncol = lx1, byrow = TRUE))
#       if (lx1 != dX[2]) {
#         ax[["names"]] = if (!is.null(nx1)) nx1 else if (lx1==1) deparse(substitute(FUN)) else paste(deparse(substitute(FUN)), seq_len(lx1), sep = ".")
#       } else ax[["names"]] = ax[["dimnames"]][[2]]
#     } else {
#       X = lapply(mctl(X), FUN, ...) #lapply(1:ncol(X), function(i) FUN(as.vector(X[,i]), ...))
#       ax[["names"]] = ax[["dimnames"]][[2]]
#     }
#     ax[["row.names"]] = ax[["dimnames"]][[1]]
#     ax[c("dim","dimnames")] = NULL
#   } else {
#     lXo = length(X[[1]])
#     if (MARGIN == 1) {
#       lX = length(X)
#       X = mrtl(matrix(unlist(X, use.names = FALSE), ncol = lX))
#       X = lapply(X, FUN, ...)
#       lx1 = length(X[[1]])
#       nx1 = names(X[[1]])
#       X = mctl(matrix(unlist(X, use.names = FALSE), ncol = lx1, byrow = TRUE))
#       if (lx1 != lX) ax[["names"]] = if (!is.null(nx1)) nx1 else if (lx1==1) deparse(substitute(FUN)) else paste(deparse(substitute(FUN)), seq_len(lx1), sep = ".") # ax = ax[c("class","row.names")]
#     } else {
#       X = lapply(X, FUN, ...)
#     }
#   }
#   lXn = length(X[[1]])
#   if (is.null(ax)) {
#     class(X) <- "data.frame"
#     attr(X, "row.names") <- .set_row_names(lXn)
#   } else {
#     if (is.null(ax[["row.names"]]) || lXo != lXn) ax[["row.names"]] <- .set_row_names(lXn)
#     # if (is.null(ax[["class"]])) ax[["class"]] <- "data.frame"
#     if (!any(ax[["class"]] == "data.frame")) ax[["class"]] <- c(ax[["class"]],"data.frame") # right order??
#     attributes(X) = ax
#   }
#   #if (simplify) simplify2array(X, higher = FALSE) else X
#   X
# }
# }
# ldply(FEVD, function(x)x) This is what you want!!
# see functions is.atomic, is.list, is.recursive, prehaps they can make checking conditions easier!!, perhaps better than checking dim??
# also if check = FALSE, it should still work with nested lists of the same dimensions!! -> Works!!
# ID is still factor under certain corcumstances!!
