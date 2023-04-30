

# TODO: keep argument? -> not needed, can use fselect beforehand...
fcount_core <- function(x, g, w = NULL, name = "N", add = FALSE) {
  # TODO: don't need integer group sizes if this is the case....
  if(length(w)) g$group.sizes <- .Call(C_fwtabulate, g$group.id, w, g$N.groups, FALSE) # na.rm in g is not needed (FALSE)
  # if(is.atomic(x)) { # what about factors and sort argument?? and dropping levels??
  #   if(add) {
  #     res <- list(x, .Call(C_subsetVector, g$group.sizes, g$group.id, FALSE))
  #     names(res) <- c(g$group.vars, name[1L])
  #   } else {
  #     res <- g$groups
  #     res[[name[1L]]] <- g$group.sizes
  #   }
  #   attr(res, "row.names") <- .set_row_names(.Call(C_fnrow, res))
  #   oldClass(res) <- "data.frame"
  #   return(res)
  # }
  if(add) {
    gs <- .Call(C_subsetVector, g$group.sizes, g$group.id, FALSE)
    # return(`add_vars<-`(x, "end", `names<-`(list(gs), name[1L])))
    if(add == 2L) {
      x <- # if(inherits(x, "grouped_df")) fgroup_vars(x) else # Better keep groups, does no harm... can use fungroup()
        .Call(C_subsetCols, x, ckmatch(g$group.vars, attr(x, "names")), TRUE)
    }
    res <- c(x, `names<-`(list(gs), name[1L]))
    return(condalc(copyMostAttributes(res, x), inherits(x, "data.table")))
  }
  res <- g$groups
  if(!is.object(res) && is.object(x)) { # inherits(x, c("grouped_df", "indexed_frame"))
    res[[name[1L]]] <- g$group.sizes
    return(condCopyAttrib(res, x))
  }
  condalc(copyMostAttributes(c(res, `names<-`(list(g$group.sizes), name[1L])), res), inherits(x, "data.table"))
}

fcount <- function(x, ..., w = NULL, name = "N", add = FALSE, sort = FALSE, decreasing = FALSE) {
  if(is.list(x)) w <- eval(substitute(w), x, parent.frame())
  else x <- qDF(x)
  if(is.character(add)) add <- switch(add, gv =, group_vars = 2L, stop("add must be TRUE, FALSE or group_vars (gv)")) # add = "g", "groups" or "group_vars"
  # Note: this code duplication with GRP() is needed for GRP() to capture x (using substitute) if x is atomic.
  # if(is.atomic(x)) `names<-`(list(x), l1orlst(as.character(substitute(x)))) else
  g <- if(missing(...)) GRP(x, sort = sort, decreasing = decreasing, return.groups = !add, return.order = FALSE, call = FALSE) else
    GRP.default(fselect(x, ...), sort = sort, decreasing = decreasing, return.groups = !add, return.order = FALSE, call = FALSE)
  fcount_core(x, g, w, name, add)
}

fcountv <- function(x, cols = NULL, w = NULL, name = "N", add = FALSE, sort = FALSE, ...) {
  # Safe enough ? or only allow character ? what about collapv() ?, extra option ?
  # if(length(w) == 1L && is.list(x) && length(unclass(x)) > 1L && (is.character(w) || is.integer(w) || (is.numeric(w) && w %% 1 < 1e-6)))
  if(is.atomic(x)) x <- qDF(x)
  if(length(w) == 1L && is.character(w)) {
    w <- .subset2(x, w) # Problem: if w is wrong character: NULL
    if(is.null(w)) stop("Unknown column: ", w)
  }
  if(is.character(add)) add <- switch(add, gv =, group_vars = 2L, stop("add must be TRUE, FALSE or group_vars (gv)")) # add = "g", "groups" or "group_vars"
  g <- if(is.null(cols)) GRP(x, sort = sort, return.groups = !add, return.order = FALSE, call = FALSE, ...) else
    GRP.default(colsubset(x, cols), sort = sort, return.groups = !add, return.order = FALSE, call = FALSE, ...)
  fcount_core(x, g, w, name, add)
}
