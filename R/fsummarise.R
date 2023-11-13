# Old, simple version:
# fFUN_add_groups <- function(x) {
#   x$g <- quote(.g_)  # Faster than [["g"]]
#   x$use.g.names <- FALSE
#   x
# }

fFUN_smr_add_groups <- function(z) {
  if(!is.call(z)) return(z)
  cz <- as.character(z[[1L]])
  if(length(cz) > 1L) cz <- if(any(cz == "collapse")) cz[length(cz)] else "" # needed if collapse::fmean etc..
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




fsummarise <- function(.data, ..., keep.group_vars = TRUE, .cols = NULL) {
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
    .data[c(".g_", ".gsplit_")] <- list(g, gsplit)
    res <- vector("list", length(e))
    for(i in 2:length(e)) { # This is good and very fast
      ei <- e[[i]]
      if(nullnam || nam[i] == "") { # Across
        if(ei[[1L]] == quote(across) || ei[[1L]] == quote(acr)) {
          ei[[1L]] <- quote(do_across)
          ei$.eval_funi <- quote(smr_funi_grouped)
          # return(eval(ei, list(do_across = do_across, smr_funi_grouped = smr_funi_grouped), pe))
          res[[i]] <- eval(ei, list(do_across = do_across, smr_funi_grouped = smr_funi_grouped), pe)
        } else res[[i]] <- do_grouped_expr_list(ei, .data, g, pe, .cols, ax)
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
      # if(!keep.group_vars) stop("all computations need to result in vectors of equal length")
      # gi <- seq_along(g$group.vars)
      # ef <- lr[length(gi)+1L] / g[[1L]]
      rnglr <- .range(lr)
      ef <- rnglr / g[[1L]]
      if(ef[1L] < 1) stop("An expression did not return a value for some groups. Please ensure that a value is returned for each group")
      ef <- ef[2L]
      # if(!all_eq(lr[-gi]) || ef %% 1 > 0) stop("all computations need to result in vectors of equal length")
      gi <- whichv(lr, rnglr[2L], invert = TRUE)
      if(ef != as.integer(ef) || !all_eq(lr[gi])) stop("all computations need to result in vectors of length 1 or the maximum length of any expression")
      res[gi] <- .Call(C_subsetDT, res, rep(seq_len(g[[1L]]), each = ef), gi, FALSE) # Using C_subsetvector is not really faster... (1-2 microseconds gain)
    }
  } else {
    # Without groups...
    ax <- attributes(.data)
    oldClass(.data) <- NULL # Not strictrly needed but just to make sure execution is efficient in across etc..
    if(nullnam || bsum(!nzchar(nam)) > 1L) { # Likely Across statement...
      for(i in 2:length(e)) {
        ei <- e[[i]]
        if(nullnam || nam[i] == "") {
          if(ei[[1L]] == quote(across) || ei[[1L]] == quote(acr)) { # stop("expressions need to be named or start with across(), or its shorthand acr().")
            ei[[1L]] <- quote(.do_across)
            ei$.eval_funi <- quote(.smr_funi_simple)
          }
          e[[i]] <- ei
        } else e[[i]] <- as.call(list(quote(list), ei))
      }
      # return(eval(e, c(.data, list(.do_across = do_across, .smr_funi_simple = smr_funi_simple)), pe))
      res <- unlist(eval(e, c(.data, list(.do_across = do_across, .smr_funi_simple = smr_funi_simple)), pe), FALSE, use.names = TRUE)
    } else res <- eval(e, .data, pe)
    # return(res)
    if(!all_eq(lr <- vlengths(res, FALSE))) {
      maxlr <- bmax(lr)
      gi <- whichv(lr, maxlr, invert = TRUE)
      if(!allv(lr[gi], 1L)) stop("all computations need to result in vectors of length 1 or the maximum length of any expression")
      res[gi] <- .Call(C_subsetDT, res, rep.int(1L, maxlr), gi, FALSE)
    }
  }
  ax[c("names", "row.names")] <- list(names(res), .set_row_names(.Call(C_fnrow, res)))
  return(condalcSA(res, ax, any(cld == "data.table")))
}

fsummarize <- fsummarise

smr <- fsummarise
