# Old, simple version:
# fFUN_add_groups <- function(x) {
#   x$g <- quote(.g_)  # Faster than [["g"]]
#   x$use.g.names <- FALSE
#   x
# }

fFUN_smr_add_groups <- function(z) {
  if(!is.call(z)) return(z)
  cz <- l1orlst(as.character(z[[1L]])) # needed if collapse::fmean etc..
  if(any(cz == .FAST_FUN_MOPS)) {
    z$g <- quote(.g_)
    if(any(cz == .FAST_STAT_FUN_POLD)) z$use.g.names <- FALSE
  } # This works for nested calls (nothing more required, but need to put at the end..)
  if(length(z) > 2L || is.call(z[[2L]])) return(as.call(lapply(z, fFUN_smr_add_groups)))
  z
}
# Works: fFUN_smr_add_groups(quote(mean(fmax(min(fmode(mpg))))/fmean(mpg) + e + f + 1 + fsd(hp) + sum(bla) / 20))
# Also: quote(sum(x) + fmean(x) + e - 1 / fmedian(z))
# Also: quote(sum(z)/2+4+e+g+h+(p/sum(u))+(q-y))
# Also: quote(b-c/i(u))
# Also: quote(i(u)-b/p(z-u/log(a)))
# Also: q/p

# Old version...
# othFUN_compute <- function(x) {
#   if(length(x) == 2L) # No additional function arguments
#     return(substitute(.copyMostAttributes_(unlist(lapply(gsplit(a, .g_), b), FALSE, FALSE), a),
#                       list(a = x[[2L]], b = x[[1L]])))
#   # With more arguments, things become more complex..
#   lapply_call <- as.call(c(list(quote(lapply), substitute(gsplit(a, .g_), list(a = x[[2L]]))), as.list(x[-2L])))
#   substitute(.copyMostAttributes_(unlist(a, FALSE, FALSE), b),
#              list(a = lapply_call, b = x[[2L]]))
# }

# Note: Need unclass here because of t_list() in do_across(), which only works if also the interior of the list is a list!
smr_funi_simple <- function(i, data, .data_, funs, aplvec, ce, ...) {
  # return(list(i = i, data = data, .data_ = .data_, funs = funs, aplvec = aplvec, ce = ce))
  .FUN_ <- funs[[i]]
  nami <- names(funs)[i]
  if(aplvec[i]) {
    value <- if(missing(...)) lapply(unattrib(.data_), .FUN_) else
      do.call(lapply, c(list(unattrib(.data_), .FUN_), eval(substitute(list(...)), data, ce)), envir = ce)
    names(value) <- names(.data_)
  } else if(any(nami == .FAST_STAT_FUN_POLD)) {
    if(missing(...)) return(unclass(.FUN_(.data_, drop = FALSE)))
    fcal <- as.call(c(list(as.name(nami), quote(.data_)), as.list(substitute(list(...))[-1L])))
    fcal$drop <- FALSE
    return(unclass(eval(fcal, c(list(.data_ = .data_), data), ce)))
  } else {
    value <- if(missing(...)) .FUN_(.data_) else
      do.call(.FUN_, c(list(.data_), eval(substitute(list(...)), data, ce)), envir = ce)
    oldClass(value) <- NULL
  }
  return(value)
  # Check is already done at the end...
  # if(all_eq(vlengths(value, FALSE))) stop("All computations must result in data values of equal length")
}

smr_funi_grouped <- function(i, data, .data_, funs, aplvec, ce, ...) {
  g <- data[[".g_"]]
  .FUN_ <- funs[[i]]
  nami <- names(funs)[i]
  if(aplvec[i]) {
    value <- if(missing(...)) lapply(unattrib(.data_), copysplaplfun, g, .FUN_) else
      dots_apply_grouped(.data_, g, .FUN_, eval(substitute(list(...)), data, ce))
    names(value) <- names(.data_)
  } else if(any(nami == .FAST_STAT_FUN_POLD)) {
    if(missing(...)) return(unclass(.FUN_(.data_, g = g, use.g.names = FALSE)))
    fcal <- as.call(c(list(as.name(nami), quote(.data_), g = quote(.g_)), as.list(substitute(list(...))[-1L])))
    fcal$use.g.names <- FALSE
    return(unclass(eval(fcal, c(list(.data_ = .data_), data), ce)))
  } else {
    value <- dots_apply_grouped_bulk(.data_, g, .FUN_, if(missing(...)) NULL else eval(substitute(list(...)), data, ce))
    value <- .Call(C_rbindlist, unclass(value), FALSE, FALSE, NULL)
    oldClass(value) <- NULL
  }
  return(value) # Again checks are done below
}




fsummarise <- function(.data, ..., keep.group_vars = TRUE) {
  if(!is.list(.data)) stop(".data needs to be a list of equal length columns or a data.frame")
  e <- substitute(list(...))
  nam <- names(e)
  nullnam <- is.null(nam)
  pe <- parent.frame()
  cld <- oldClass(.data) # This needs to be called cld, because across fetches it from here !!

  if(any(cld == "grouped_df")) {
    oldClass(.data) <- NULL
    g <- GRP.grouped_df(.data, call = FALSE)
    attr(.data, "groups") <- NULL
    ax <- attributes(.data)
    ax[["class"]] <- fsetdiff(cld, c("GRP_df", "grouped_df"))
    .data[[".g_"]] <- g
    .data[[".gsplit_"]] <- gsplit
    res <- vector("list", length(e))
    for(i in 2:length(e)) { # This is good and very fast
      ei <- e[[i]]
      if(nullnam || nam[i] == "") { # Across
        if(ei[[1L]] != quote(across) && ei[[1L]] != quote(acr)) stop("expressions need to be named or start with across(), or its shorthand acr().")
        ei[[1L]] <- quote(do_across)
        ei$.eval_funi <- quote(smr_funi_grouped)
        # return(eval(ei, list(do_across = do_across, smr_funi_grouped = smr_funi_grouped), pe))
        res[[i]] <- eval(ei, list(do_across = do_across, smr_funi_grouped = smr_funi_grouped), pe)
      } else { # Tagged vector expressions
        eif <- all_funs(ei)
        res[[i]] <- list(if(any(eif %in% .FAST_STAT_FUN_POLD))  # startsWith(eif, .FAST_STAT_FUN_POLD) Note: startsWith does not reliably capture expressions e.g. e <- quote(list(b = fmean(log(mpg)) + max(qsec))) does not work !!
                         eval(fFUN_smr_add_groups(ei), .data, pe) else
                         do_grouped_expr(ei, length(eif), .data, g, pe))
      }
    }
    names(res) <- nam
    res[[1L]] <- if(keep.group_vars) g$groups else NULL
    res <- unlist(res, FALSE, use.names = TRUE)
    # replicating groups if more rows per computation...
    if(!all_eq(lr <- vlengths(res, FALSE))) {
      if(!keep.group_vars) stop("all computations need to result in vectors of equal length")
      gi <- seq_along(g$group.vars)
      ef <- lr[length(gi)+1L] / g[[1L]]
      if(!all_eq(lr[-gi]) || ef %% 1 > 0) stop("all computations need to result in vectors of equal length")
      res[gi] <- .Call(C_subsetDT, res[gi], rep(seq_len(g[[1L]]), each = ef), gi, FALSE) # Using C_subsetvector is not really faster... (1-2 microseconds gain)
    }
  } else {
    # Without groups...
    ax <- attributes(.data)
    oldClass(.data) <- NULL # Not strictrly needed but just to make sure execution is efficient in across etc..
    if(nullnam || bsum(!nzchar(nam)) > 1L) { # Likely Across statement...
      for(i in 2:length(e)) {
        ei <- e[[i]]
        if(nullnam || nam[i] == "") {
          if(ei[[1L]] != quote(across) && ei[[1L]] != quote(acr)) stop("expressions need to be named or start with across(), or its shorthand acr().")
          ei[[1L]] <- quote(.do_across)
          ei$.eval_funi <- quote(.smr_funi_simple)
          e[[i]] <- ei
        } else e[[i]] <- as.call(list(quote(list), ei))
      }
      # return(eval(e, c(.data, list(.do_across = do_across, .smr_funi_simple = smr_funi_simple)), pe))
      res <- unlist(eval(e, c(.data, list(.do_across = do_across, .smr_funi_simple = smr_funi_simple)), pe), FALSE, use.names = TRUE)
    } else res <- eval(e, .data, pe)
    # return(res)
    if(!all_eq(vlengths(res, FALSE))) stop("all computations need to result in vectors of equal length")
  }
  ax[["names"]] <- names(res)
  ax[["row.names"]] <- .set_row_names(length(res[[1L]]))
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



