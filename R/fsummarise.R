fFUN_add_groups <- function(x) {
  x$g <- quote(.g_)  # Faster than [["g"]]
  x$use.g.names <- FALSE
  x
}

othFUN_compute <- function(x) {
  if(length(x) == 2L) # No additional function arguments
    return(substitute(copyMostAttrib(unlist(lapply(split(a, as_factor_GRP(.g_)), b), FALSE, FALSE), a),
                      list(a = x[[2L]], b = x[[1L]])))
  # With more arguments, things become more complex..
  lapply_call <- as.call(c(list(quote(lapply), substitute(split(a, as_factor_GRP(.g_)), list(a = x[[2L]]))), as.list(x[-2L])))
  substitute(copyMostAttrib(unlist(a, FALSE, FALSE), b),
             list(a = lapply_call, b = x[[2L]]))
}

fsummarise <- function(.data, ..., keep.group_vars = TRUE) {
  if(!is.list(.data)) stop(".data needs to be a list of equal length columns or a data.frame")
  e <- substitute(list(...))
  cld <- oldClass(.data)
  if(any(cld == "grouped_df")) {
    g <- GRP.grouped_df(.data, call = FALSE)
    ax <- attributes(fungroup(.data))
      # FUNs <- vapply(e[-1L], function(x) as.character(x[[1L]]), character(1L), USE.NAMES = FALSE)
      # if(any(FUNs %!in% .FAST_STAT_FUN)) ...
    for(i in seq_along(e)[-1L]) { # This is good and very fast
      ei <- e[[i]]
      e[[i]] <- if(any(as.character(ei[[1L]]) %in% .FAST_STAT_FUN)) # need %in% here because could pass collapse::fmean etc..
                  fFUN_add_groups(ei) else othFUN_compute(ei)
    }
    res <- eval(e, c(list(.g_ = g), .data), parent.frame())
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
