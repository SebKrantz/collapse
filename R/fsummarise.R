fFUN_add_groups <- function(x) {
  x$g <- quote(.g_)  # Faster than [["g"]]
  x$use.g.names <- FALSE
  x
}

othFUN_compute <- function(x) {
  if(length(x) == 2L) # No additional function arguments
    return(substitute(.copyMostAttributes_(unlist(lapply(gsplit(a, .g_), b), FALSE, FALSE), a),
                      list(a = x[[2L]], b = x[[1L]])))
  # With more arguments, things become more complex..
  lapply_call <- as.call(c(list(quote(lapply), substitute(gsplit(a, .g_), list(a = x[[2L]]))), as.list(x[-2L])))
  substitute(.copyMostAttributes_(unlist(a, FALSE, FALSE), b),
             list(a = lapply_call, b = x[[2L]]))
}


fsummarise <- function(.data, ..., keep.group_vars = TRUE) {
  if(!is.list(.data)) stop(".data needs to be a list of equal length columns or a data.frame")
  e <- substitute(list(...))
  cld <- oldClass(.data)
  if(any(cld == "grouped_df")) {
    g <- GRP.grouped_df(.data, call = FALSE)
    ax <- attributes(fungroup(.data))
    for(i in seq_along(e)[-1L]) { # This is good and very fast
      ei <- e[[i]]
      e[[i]] <- if(any(startsWith(as.character(ei[[1L]]), .FAST_STAT_FUN_POLD))) # could pass collapse::flast.default etc..
                  fFUN_add_groups(ei) else othFUN_compute(ei)
    }
    res <- eval(e, c(list(.g_ = g, .copyMostAttributes_ = copyMostAttributes), .data), parent.frame())
    if(keep.group_vars) res <- c(g[["groups"]], res)
    ax[["names"]] <- names(res)
    ax[["row.names"]] <- .set_row_names(g[[1L]])
    return(condalcSA(res, ax, any(cld == "data.table")))
  }
  ax <- attributes(.data)
  res <- eval(e, .data, parent.frame())
  ax[["names"]] <- names(res)
  ax[["row.names"]] <- 1L
  return(condalcSA(res, ax, any(cld == "data.table")))
}

smr <- fsummarise

# sumr # -> yes, but above is more consistent with other shortcuts

# Some speed improvement before gsplit:
#
# othFUN_compute <- function(x) {
#   if(length(x) == 2L) # No additional function arguments
#     return(substitute(copyMostAttrib(unlist(lapply(split(a, .g_f), b), FALSE, FALSE), a),
#                       list(a = x[[2L]], b = x[[1L]])))
#   # With more arguments, things become more complex..
#   lapply_call <- as.call(c(list(quote(lapply), substitute(split(a, .g_f), list(a = x[[2L]]))), as.list(x[-2L])))
#   substitute(copyMostAttrib(unlist(a, FALSE, FALSE), b),
#              list(a = lapply_call, b = x[[2L]]))
# }
#
# fsummarise <- function(.data, ..., keep.group_vars = TRUE) {
#   if(!is.list(.data)) stop(".data needs to be a list of equal length columns or a data.frame")
#   e <- substitute(list(...))
#   cld <- oldClass(.data)
#   dm <- c(list(.g_ = g), .data)
#   ofl <- TRUE
#   if(any(cld == "grouped_df")) {
#     g <- GRP.grouped_df(.data, call = FALSE)
#     ax <- attributes(fungroup(.data))
#     # FUNs <- vapply(e[-1L], function(x) as.character(x[[1L]]), character(1L), USE.NAMES = FALSE)
#     # if(any(FUNs %!in% .FAST_STAT_FUN_POLD)) ...
#     for(i in seq_along(e)[-1L]) { # This is good and very fast
#       ei <- e[[i]]
#       if(any(startsWith(as.character(ei[[1L]]), .FAST_STAT_FUN_POLD))) { # could pass collapse::flast.default etc..
#         e[[i]] <- fFUN_add_groups(ei)
#       } else {
#         if(ofl) {
#           dm$.g_f <- as_factor_GRP(g)
#           ofl <- FALSE
#         }
#         e[[i]] <- othFUN_compute(ei)
#       }
#     }
#     res <- eval(e, dm, parent.frame())
#     if(keep.group_vars) res <- c(g[["groups"]], res)
#     ax[["names"]] <- names(res)
#     ax[["row.names"]] <- .set_row_names(g[[1L]])
#     return(condalcSA(res, ax, any(cld == "data.table")))
#   }
#   ax <- attributes(.data)
#   res <- eval(e, .data, parent.frame())
#   ax[["names"]] <- names(res)
#   ax[["row.names"]] <- 1L
#   return(condalcSA(res, ax, any(cld == "data.table")))
# }



