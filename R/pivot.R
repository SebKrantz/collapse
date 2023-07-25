

proc_names_longer <- function(x) {
  if(is.null(x)) return(list("variable", "value"))
  if(is.list(x)) { # is.character(x) : list is not necessary but clearer (also regarding multiple casts etc.) !!!
    if(is.null(names(x))) {
      if(length(x) != 2L) stop("If how = 'longer', 'names' needs to be a length-2 list or a named list. You specified a list length: ", length(x))
      return(x)
    }
    if(length(x) > 2L) stop("If how = 'longer', 'names' needs to be a length-2 list or a named list length-1 or -2. You specified a list length: ", length(x))
    res <- list(variable = "variable", value = "value")
    ind <- ckmatch(names(x), names(res), e = "Unknown keywords (must be 'variable' and/or 'value'):")
    res[ind] <- x
    return(res)
  }
  stop("If how = 'longer', 'names' needs to be a (named) list. You supplied a vector of type: ", typeof(x))
}

proc_names_recast <- function(x, data) {
  if(is.null(x)) {
    nam_col <- data[["variable"]]
    if(is.null(nam_col)) stop("Need to provide 'names'. The default name 'variable' was not found in the data.")
    return(list(nam_col, "variable"))
  }
  if(is.list(x)) {
    if(is.null(names(x))) {
      if(length(x) != 2L) stop("If how = 'recast', 'names' needs to be a length-2 list or a named list. You specified a list length: ", length(x))
    } else {
      if(length(x) > 2L) stop("If how = 'recast', 'names' needs to be a length-2 list or a named list length-1 or -2. You specified a list length: ", length(x))
      res <- list(from = "variable", to = "variable")
      ind <- ckmatch(names(x), names(res), e = "Unknown keywords (must be 'from' and/or 'to'):")
      res[ind] <- x
      x <- res
    }
    ind <- cols2int(x[[1L]], data, names(data))
    nam_col <- if(length(ind) == 1L) data[[ind]] else finteraction(data[ind], sort = FALSE, sep = "_")
    return(list(nam_col, x[[2L]]))
  }
  stop("If how = 'recast', 'names' needs to be a (named) list. You supplied a vector of type: ", typeof(x))
}


# Crbindlist <- function(x) .Call(C_rbindlist, x, FALSE, FALSE, NULL)
# Faster than do.call(c, unattrib(data[values])):
# c_to_vec <- function(l) .Call(C_rbindlist, lapply(unattrib(l), list), FALSE, FALSE, NULL)[[1L]]
# Same thing (also same speed), a bit less cumbersome...
# c_to_vec2 <- function(l) .Call(C_pivot_long, l, NULL, FALSE)

# Special case: no ids supplied
melt_all <- function(vd, names, factor, na.rm, labels, check.dups) {
  if(check.dups && fnrow(vd) > 1L) warning("duplicates detected: you have supplied no ids and the data has ", fnrow(vd), " rows. Consider supplying ids so that that records in the long format data frame are identified.")
  if(length(labels)) labs <- vlabels(vd, use.names = FALSE)
  # 6 cases: label or not, factor or not (either id or label)
  if(length(labels) || factor[1L]) { # if labels: generate id to expand vectors: faster than rep...
    nam <- names(vd)
    attributes(vd) <- NULL
  }
  if(na.rm) vd <- lapply(vd, na_rm) # Note: beforehand is faster, I tested it...
  res <- .Call(C_pivot_long, vd, NULL, TRUE) # rbindlist gives factor value: .Call(C_rbindlist, lapply(unattrib(vd), list), FALSE, FALSE, "id")
  names(res) <- names
  if(length(labels)) {
    if(factor[2L]) {
      label_col <- res[[1L]]
      attr(label_col, "levels") <- labs
      oldClass(label_col) <- c("factor", "na.included")
    } else label_col <- Csv(labs, res[[1L]])
    label_col <- list(label_col)
    names(label_col) <- if(is.character(labels)) labels else "label"
    res <- c(res[1L], label_col, res[2L])
  }
  if(factor[1L]) {
    attr(res[[1L]], "levels") <- nam
    oldClass(res[[1L]]) <- c("factor", "na.included")
  } else if(length(labels)) res[[1L]] <- Csv(nam, res[[1L]])
  res
}

# TODO: multiple pivots
# TODO: set new labels!!
# TODO: Think about: values could be list input, names only atomic. that would make more sense...
# Or: allow for both options... needs to be consisntent with "labels" though...
# TODO: Multithreading, SIMD
# TODO: option sort = TRUE | FALSE
# TODO: Option transpose = c(names = FALSE, columns = FALSE)
# TODO: leave args empty...

# Check labels and attributes..
pivot <- function(data,
                  ids = NULL,
                  values = NULL,
                  names = NULL, # list is better
                  labels = NULL, # todo: should also be list, e.g. for recast problem!!
                  how = c("longer", "wider", "recast"), # Better to only have one?, because the other arguments use multiple??
                  na.rm = FALSE,
                  factor = c("names", "labels"),
                  check.dups = FALSE,
                  fill = NULL, # Fill is for pivot_wider
                  drop = TRUE, # Same as with dcast()
                  sort = c(ids = FALSE, names = FALSE),
                  transpose = c(cols = FALSE, names = FALSE)) {

    if(!is.list(data)) stop("pivot only supports data.frame-like objects")
    ad <- attributes(data)
    oldClass(data) <- NULL
    nam <- names(data)
    if(length(ids)) ids <- cols2int(ids, data, nam)
    if(length(values)) values <- cols2int(values, data, nam)
    factor <- c("names", "labels") %in% factor # TODO: needed here, or in switch??


    res <- switch(how[1L],

      l = , longer = { # TODO: multiple output columns
        names <- proc_names_longer(names)
        if(is.null(ids)) melt_all(if(is.null(values)) data else data[values],
                                  names, factor, na.rm, labels, check.dups)
        else {
          if(is.null(values)) values <- seq_along(data)[-ids]
          vd <- data[values]
          if(length(labels) || factor[1L]) attributes(vd) <- NULL
          if(check.dups && force(ng <- fnunique(data[ids])) < fnrow(data))
             warning("duplicated id values detected: there are ", ng, " unique id-combinations, but the data has ", fnrow(data),
                     " rows. This means you have on average ", round(fnrow(data)/ng, 1), " duplicates per id-combination. ",
                     "Consider adding additional ids or aggregating your data (e.g. using collap()) before applying pivot().")
          if(na.rm) {
            cc <- lapply(vd, whichNA, invert = TRUE) # TODO: could do this all internally using a single vector
            # cc_vec <- c_to_vec(cc)
            # id_cols <- .Call(C_subsetDT, data, cc_vec, ids, FALSE)
            id_cols <- lapply(data[ids], function(x) .Call(C_pivot_long, alloc(x, length(cc), FALSE), cc, FALSE))
            value_cols <- .Call(C_pivot_long, vd, cc, TRUE)
            # value_col <- .Call(C_pivot_long, vd, cc, FALSE) # Csv(c_to_vec(data[values]), cc_vec)
            # variable_col <- rep(if(factor[1L]) seq_along(values) else nam[values], vlengths(cc))
          } else {
            id_cols <- .Call(C_rbindlist, alloc(data[ids], length(values)), FALSE, FALSE, NULL) # .Call(C_subsetDT, data, rep.int(seq_len(n), length(values)), ids, FALSE)
            # This is faster than .Call(C_pivot_long, vd, NULL) because rep() is slow...
            value_cols <- .Call(C_pivot_long, vd, NULL, TRUE) # .Call(C_rbindlist, lapply(vd, list), FALSE, FALSE, "id")
            # value_col <- .Call(C_pivot_long, vd, NULL)   # c_to_vec(data[values])
            # variable_col <- rep(if(factor[1L]) seq_along(values) else nam[values], each = fnrow(data))
          }
          names(value_cols) <- names # TODO: multiple pivots this does not work...
          if(length(labels)) { # TODO: what is list etc?? Need to accommodate this case??
            labs <- vlabels(vd, use.names = FALSE)
            if(factor[2L]) {
              label_col <- value_cols[[1L]]
              attr(label_col, "levels") <- labs
              oldClass(label_col) <- c("factor", "na.included")
            } else label_col <- Csv(labs, value_cols[[1L]])
            label_col <- list(label_col)
            names(label_col) <- if(is.character(labels)) labels else "label"
            value_cols <- c(value_cols[1L], label_col, value_cols[2L])
          }
          if(factor[1L]) {
            attr(value_cols[[1L]], "levels") <- nam[values]
            oldClass(value_cols[[1L]]) <- c("factor", "na.included")
          } else if(length(labels)) value_cols[[1L]] <- Csv(nam[values], value_cols[[1L]])
          c(id_cols, value_cols)
        }
      },

      # In general: names specifies where variable names are coming from. If multiple then interact them using "_"
      # Same for labels. drop specifies that factor levels should be dropped if a single factor column is passed to names
      w = , wider = {
        # (1) Preprocessing Arguments
        if(is.null(names)) {
          names <- whichv(nam, "variable")
          if(!length(names)) stop("Need to provide 'names' if how = 'wider'. The default name 'variable' was not found in the data.")
        } else names <- cols2int(names, data, nam)
        if(length(labels)) labels <- cols2int(labels, data, nam)
        if(is.null(values)) {
          if(is.null(ids)) {
            values <- whichv(nam, "value")
            if(!length(values)) stop("Need to provide values if how = 'wider' and is.null(ids). The default name 'value' was not found in the data.")
          } else values <- seq_along(data)[-c(ids, names, labels)]
        }
        if(is.null(ids)) ids <- seq_along(data)[-c(names, labels, values)]
        # (2) Missing Value Removal
        if(na.rm) { # TODO: better way?
          data <- na_omit(data[c(ids, names, values, labels)], cols = values)
          ids <- seq_along(ids)
          names <- seq_along(names) + length(ids)
          values <- seq_along(values) + length(ids) + length(names)
          if(length(labels)) labels <- seq_along(labels) + length(ids) + length(names) + length(values)
        }
        # (3) Compute ID Columns
        if(sort["ids"]) {
          g <- GRP.default(data[ids], sort = TRUE, return.order = FALSE, call = FALSE)
          id_cols <- g[[4L]]
          g <- g[[2L]]
        } else { # Could also use GRP(), but this avoids computing a potentially large and redundant group sizes vector
          g <- group(data[ids], starts = TRUE)
          id_cols <- .Call(C_subsetDT, data, attr(g, "starts"), ids, FALSE)
        }
        # (4) Compute Names and Labels Columns
        names_g <- GRP(if(length(names) == 1L && is.null(labels)) data[[names]] else data[names],
                       sort = sort["names"], group.sizes = check.dups, drop = drop, call = FALSE)
        names <- GRPnames(names_g, sep = "_")
        if(length(labels)) {
          if(check.dups && any(vary <- varying(data[labels], names_g)))  # See if there are duplicate labels
            stop("The following 'labels' columns vary by 'names': ", paste(names(vary)[vary], collapse = ", "))
          labels <- if(length(labels) == 1L) Csv(data[[labels]], names_g$group.starts) else
            do.call(paste, c(.Call(C_subsetDT, data, names_g$group.starts, labels, FALSE), list(sep = " - ")))
        }
        g_v <- names_g[[2L]]
        attr(g_v, "N.groups") <- names_g[[1L]]
        # (5) Optional duplicates check
        if(check.dups) {
          # Old way of doing it:
          # if(force(ng <- fnunique(list(g, g_v))) < fnrow(data))
          #   warning("duplicates detected: there are ", ng, " unique combinations of id- and name-columns, but the data has ", fnrow(data),
          #           " rows. This means you have on average ", round(fnrow(data)/ng, 1), " duplicates per id-name-combination. If how = 'wider', pivot() will take the last of those duplicates in first-appearance-order. Consider aggregating your data e.g. using collap() before applying pivot().")
          # This should be faster... but check !!
          ndg <- fndistinct.default(g, names_g, use.g.names = FALSE, na.rm = FALSE)
          if(!identical(ngd, names_g[[3L]])) {
            ng <- bsum(ndg)
            warning("duplicates detected: there are ", ng, " unique combinations of id- and name-columns, but the data has ", fnrow(data),
                    " rows. This means you have on average ", round(fnrow(data)/ng, 1), " duplicates per id-name-combination. If how = 'wider', pivot() will take the last of those duplicates in first-appearance-order. Consider aggregating your data e.g. using collap() before applying pivot().")
          }
        }
        # (6) Compute Reshaped Values
        if(length(values) > 1L) { # Multiple columns, as in dcast... TODO: check pivot_wider
          namv <- names(data)[values]
          attributes(data) <- NULL
          value_cols <- lapply(data[values], function(x) .Call(C_pivot_wide, g, g_v, x))
          if(length(labels)) value_cols <- lapply(value_cols, `vlabels<-`, value = labels)
          value_cols <- unlist(if(transpose["cols"]) t_list2(value_cols) else value_cols, FALSE, FALSE)
          namv_res <- if(transpose["names"]) t(outer(names, namv, paste, sep = "_")) else outer(namv, names, paste, sep = "_")
          names(value_cols) <- if(transpose["cols"]) namv_res else t(namv_res)
        } else {
          value_cols <- .Call(C_pivot_wide, g, g_v, data[[values]])
          names(value_cols) <- names
          if(length(labels)) vlabels(value_cols) <- labels
        }
        c(id_cols, value_cols)
      },

      # TODO: na.rm option implementation!!
      r = , recast = {
        # The optimization applied here is to avoid materialization of the "long" id-columns
        # There are two ways to do it, first the long value cast and then wide cast, or many wide casts and row-biding.
        # The complication is that the long cast requires construction of an id-column, which probably can only be efficiently
        # done by creating yet another C-function. Thus I try the wide option first.
        # -> initial benchmarks show that this is also definitely faster than recast from long frame...
        # but presumably because grouping is much faster. If an id is constructed we don't need to group a long frame though...
        g <- group(data[ids], starts = TRUE)
        id_cols <- .Call(C_subsetDT, data, attr(g, "starts"), ids, FALSE)
        if(!(is.list(names) && length(names) == 2L)) stop("'names' must be a list of length 2. See Documentation.")
        names1 <- names[[1L]] # TODO: could be integer etc??
        names_col <- if(length(names1) == 1L) data[[names1]] else finteraction(data[names1], sort = FALSE, sep = "_")
        if(is.factor(names_col)) {
          g_v <- names_col
          if(!inherits(g_v, "na.included")) g_v <- addNA2(g_v)
          if(drop && length(names1) == 1L) g_v <- fdroplevels.factor(g_v)
          lev <- attr(g_v, "levels")
          attr(g_v, "N.groups") <- length(lev)
        } else g_v <- group(names_col, starts = TRUE)
        if(check.dups && force(ng <- fnunique(list(g, g_v))) < fnrow(data))
          warning("duplicates detected: there are ", ng, " unique combinations of id- and name-columns, but the data has ", fnrow(data),
                  " rows. This means you have on average ", round(fnrow(data)/ng, 1), " duplicates per id-name-combination. If how = 'recast', pivot() will take the last of those duplicates in first-appearance-order. Consider aggregating your data e.g. using collap() before applying pivot().")
        if(factor[1L]) attributes(data) <- NULL
        value_cols <- lapply(data[values], function(x) .Call(C_pivot_wide, g, g_v, x))
        names(value_cols[[1L]]) <- if(is.factor(names_col)) lev else names_col[attr(g_v, "starts")]
        if(length(value_cols) > 1L) vlabels(value_cols[[1L]]) <- NULL
        id_cols <- .Call(C_rbindlist, alloc(id_cols, length(value_cols)), FALSE, FALSE, NULL)
        value_cols <- .Call(C_rbindlist, value_cols, FALSE, FALSE, names[[2L]]) # Final column is "variable" name
        if(factor[1L]) {
          attr(value_cols[[1L]], "levels") <- nam[values]
          oldClass(value_cols[[1L]]) <- c("factor", "na.included")
        }
        c(id_cols, value_cols)
      },

      stop("Unknown pivoting method: ", how[1L])
    ) # end of switch

    if(is.null(ad)) return(res) # Redundant ??
    if(any(ad$class == "data.frame")) ad$row.names <- .set_row_names(fnrow(res))
    ad$names <- names(res)
    .Call(C_setattributes, res, ad)
    if(any(ad$class == "data.table")) return(alc(res))
    return(res)
}
