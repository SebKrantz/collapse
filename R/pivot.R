

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
    ind <- whichv(names(data), "variable")
    if(!length(ind)) stop("Need to provide 'names'. The default name 'variable' was not found in the data.")
    return(list(ind, "variable"))
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
    # nam_col <- if(length(ind) == 1L) data[[ind]] else finteraction(data[ind], sort = sort, sep = "_")
    return(list(ind, x[[2L]]))
  }
  stop("If how = 'recast', 'names' needs to be a (named) list. You supplied a vector of type: ", typeof(x))
}

proc_labels_recast <- function(x, data) {
  if(is.list(x)) {
    if(is.null(names(x))) {
      if(length(x) != 2L && length(x) != 3L) stop("If how = 'recast', 'labels' needs to be a length-2 list or a named list. You specified a list length: ", length(x))
    } else {
      if(length(x) > 3L) stop("If how = 'recast', 'labels' needs to be a length-2 list or a named list length-1 or -2. You specified a list length: ", length(x))
      res <- list(from = NULL, to = NULL, new = NULL)
      ind <- ckmatch(names(x), names(res), e = "Unknown keywords (must be 'from', 'to' or 'new'):")
      res[ind] <- x
      x <- res
    }
    ind <- if(length(x[[1L]])) cols2int(x[[1L]], data, names(data)) else NULL
    return(list(ind, x[[2L]], x[[3L]]))
  }
  stop("If how = 'recast', 'labels' needs to be a (named) list. You supplied a vector of type: ", typeof(x))
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
    if(is.list(labels)) stop("Since no ids are specified, please just use setLabels() or relabel() following pivot to assign new variable labels")
    if(factor[2L]) {
      label_col <- res[[1L]]
      attr(label_col, "levels") <- labs
      oldClass(label_col) <- "factor" # c("factor", "na.included")
    } else label_col <- Csv(labs, res[[1L]])
    label_col <- list(label_col)
    names(label_col) <- if(is.character(labels)) labels else "label"
    res <- c(res[1L], label_col, res[2L])
  }
  if(factor[1L]) {
    attr(res[[1L]], "levels") <- nam
    oldClass(res[[1L]]) <- "factor" # c("factor", "na.included")
  } else if(length(labels)) res[[1L]] <- Csv(nam, res[[1L]])
  res
}

# Retain labels in wider reshaping
add_labels <- function(l, labs) {
  ll <- .Call(C_vlabels, l, "label", FALSE)
  if(!allNA(ll)) labs <- paste(ll, labs, sep = " - ")
  .Call(C_setvlabels, l, "label", labs, NULL)
}

apply_external_FUN <- function(data, g, FUN, args, name) {
  FUN <- match.fun(FUN)
  if(is.null(args)) {
    if(any(name == .FAST_STAT_FUN)) return(FUN(data, g = g, TRA = "fill"))
    return(TRA(data, BY(data, g, FUN, use.g.names = FALSE, reorder = FALSE), "fill", g))
  }
  if(any(name == .FAST_STAT_FUN)) return(do.call(FUN, c(list(x = data, g = g, TRA = "fill"), args)))
  TRA(data, do.call(BY, c(list(x = data, g = g, FUN = FUN, use.g.names = FALSE, reorder = FALSE), args)), "fill", g)
}

# TODO: Think about: values could be list input, names only atomic. that would make more sense...
# Or: allow for both options... needs to be consistent with "labels" though...

# Transposition Example:
# pivot(BWA, names = list(from = c("Variable", "Year"), to = "Sectorcode"), how = "r")

# data = BWA
# ids = NULL
# names = list(from = c("Variable", "Year"), to = "Sectorcode")
# labels = NULL
# values = NULL
# how = "r"
# na.rm = FALSE
# factor = c("names", "labels")
# check.dups = FALSE
# fill = NULL
# drop = TRUE
# sort = FALSE
# nthreads = 1L
# transpose = FALSE



# Check labels and attributes..
pivot <- function(data,
                  ids = NULL,
                  values = NULL,
                  names = NULL,   # list is better
                  labels = NULL,
                  how = "longer", # Better to only have one?, because the other arguments use multiple??
                  na.rm = FALSE,
                  factor = c("names", "labels"),
                  check.dups = FALSE,
                  FUN = "last",
                  FUN.args = NULL,
                  nthreads = .op[["nthreads"]],
                  fill = NULL, # Fill is for pivot_wider
                  drop = TRUE, # Same as with dcast()
                  sort = FALSE, # c("ids", "names")
                  transpose = FALSE) # c(columns = FALSE, names = FALSE))
{

  if(!is.list(data)) stop("pivot only supports data.frame-like objects")
  ad <- attributes(data)
  oldClass(data) <- NULL
  nam <- names(data)
  if(length(ids)) ids <- cols2int(ids, data, nam)
  if(length(values)) values <- cols2int(values, data, nam)
  factor <- c("names", "labels") %in% factor
  how <- switch(how, l = , longer = 1L, w = , wider = 2L, r = , recast = 3L,
                stop("Unknown pivoting method: ", how))

  if(how == 1L) { # TODO: multiple output columns
        names <- proc_names_longer(names)
        if(is.null(ids) && is.null(values)) res <- melt_all(if(is.null(values)) data else data[values],
                                  names, factor, na.rm, labels, check.dups)
        else {
          if(is.null(values)) values <- seq_along(data)[-ids]
          else if(is.null(ids)) ids <- seq_along(data)[-values]
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
          if(length(values) > 1L) vlabels(value_cols) <- NULL # Could solve at C-level with additional argument...
          names(value_cols) <- names # TODO: multiple pivots this does not work...
          if(length(labels)) {
            labs <- vlabels(vd, use.names = FALSE)
            if(factor[2L]) {
              label_col <- value_cols[[1L]]
              attr(label_col, "levels") <- labs
              oldClass(label_col) <- "factor" # c("factor", "na.included")
            } else label_col <- Csv(labs, value_cols[[1L]])
            label_col <- list(label_col)
            if(is.list(labels)) { # Setting new labels...
              if(is.null(names(labels))) {
                new_labels <- labels[[2L]]
                label <- labels[[1L]]
              } else {
                new_labels <- labels[["new"]]
                label <- labels[["name"]]
                if(is.null(label)) label <- "label"
              }
              if(!is.character(label)) stop("label column name supplied in a list needs to be character typed, you passed an object of type: ", typeof(labels))
              if(!is.character(new_labels)) stop("new labels need to be specified as a character vector, you passed an object of type: ", typeof(new_labels))
              names(label_col) <- label
              value_cols <- c(value_cols[1L], label_col, value_cols[2L])
              if(is.null(names(new_labels))) {
                if(length(new_labels) != length(value_cols)) stop("Number of new labels supplied must match number of new columns in long format frame. There are ", length(value_cols), " new columns in the molten frame, and you supplied ", length(new_labels), " new labels")
                vlabels(value_cols) <- new_labels
              } else vlabels(value_cols)[names(new_labels)] <- new_labels
            } else {
              names(label_col) <- if(is.character(labels)) labels else "label"
              value_cols <- c(value_cols[1L], label_col, value_cols[2L])
            }
          }
          if(factor[1L]) {
            attr(value_cols[[1L]], "levels") <- nam[values]
            oldClass(value_cols[[1L]]) <- "factor" # c("factor", "na.included")
          } else if(length(labels)) value_cols[[1L]] <- duplAttributes(Csv(nam[values], value_cols[[1L]]), value_cols[[1L]])
          res <- c(id_cols, value_cols)
        }

    } else {

     sort <- if(is.logical(sort)) rep(sort, length.out = 2L) else c("ids", "names") %in% sort
     transpose <- if(is.logical(transpose)) rep(transpose, length.out = 2L) else c("columns", "names") %in% transpose

     if (how == 2L) { # Wide Pivot

      # Note: No Complete Pivoting (no ids and values) supported! This does not make a lot of sense!

      # In general: names specifies where variable names are coming from. If multiple then interact them using "_"
      # Same for labels. drop specifies that factor levels should be dropped if a single factor column is passed to names
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
          data <- data[c(ids, names, values, labels)]
          ids <- seq_along(ids)
          names <- seq_along(names) + length(ids)
          values <- seq_along(values) + length(ids) + length(names)
          if(length(labels)) labels <- seq_along(labels) + length(ids) + length(names) + length(values)
          data <- na_omit(data, cols = values, prop = 1)
        }
        # (3) Compute ID Columns
        if(sort[1L]) {
          g <- GRP.default(data[ids], sort = TRUE, return.order = FALSE, call = FALSE)
          id_cols <- g[[4L]]
          ng <- g[[1L]]
          g <- g[[2L]]
          attr(g, "N.groups") <- ng
        } else { # Could also use GRP(), but this avoids computing a potentially large and redundant group sizes vector
          g <- group(data[ids], starts = TRUE)
          id_cols <- .Call(C_subsetDT, data, attr(g, "starts"), ids, FALSE)
        }
        # (4) Compute Names and Labels Columns
        names_g <- GRP(if(length(names) == 1L && is.null(labels)) data[[names]] else data[names],
                       sort = sort[2L], group.sizes = check.dups, drop = drop, call = FALSE)
        names <- GRPnames(names_g, sep = "_")
        if(length(labels)) {
          if(check.dups && any(vary <- varying(data[labels], names_g)))  # See if there are duplicate labels
            stop("The following 'labels' columns vary by 'names': ", paste(names(vary)[vary], collapse = ", "))
          labels <- if(length(labels) == 1L) tochar(Csv(data[[labels]], names_g$group.starts)) else
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

          # With 10 million obs, 1 million id groups (g), and 100 names groups, this is 2x faster than the fnunique() option + could multithread
          ndg <- fndistinct.default(g, names_g, use.g.names = FALSE, na.rm = FALSE, nthreads = nthreads)
          attributes(ndg) <- NULL
          if(!identical(ndg, names_g[[3L]])) {
            ng <- fsumC(ndg, narm = FALSE)
            warning("duplicates detected: there are ", ng, " unique combinations of id- and name-columns, but the data has ", fnrow(data),
                    " rows. This means you have on average ", round(fnrow(data)/ng, 1), " duplicates per id-name-combination. If how = 'wider', pivot() will take the last of those duplicates in first-appearance-order. Consider aggregating your data e.g. using collap() before applying pivot().")
          }
        }
        # (6) Compute Reshaped Values
        if(length(values) > 1L) { # Multiple columns, as in dcast... TODO: check pivot_wider
          namv <- names(data)[values]
          attributes(data) <- NULL
          if(!is.character(FUN)) {
            data[values] <- apply_external_FUN(data[values], group(list(g, g_v)), FUN, FUN.args, l1orlst(as.character(substitute(FUN))))
            FUN <- "last"
          }
          value_cols <- lapply(data[values], function(x) .Call(C_pivot_wide, g, g_v, x, fill, nthreads, FUN, na.rm))
          if(length(labels)) value_cols <- lapply(value_cols, add_labels, labels)
          value_cols <- unlist(if(transpose[1L]) t_list2(value_cols) else value_cols, FALSE, FALSE)
          namv_res <- if(transpose[2L]) t(outer(names, namv, paste, sep = "_")) else outer(namv, names, paste, sep = "_")
          names(value_cols) <- if(transpose[1L]) namv_res else t(namv_res)
        } else {
          if(!is.character(FUN)) {
            data[[values]] <- apply_external_FUN(data[[values]], group(list(g, g_v)), FUN, FUN.args, l1orlst(as.character(substitute(FUN))))
            FUN <- "last"
          }
          value_cols <- .Call(C_pivot_wide, g, g_v, data[[values]], fill, nthreads, FUN, na.rm)
          names(value_cols) <- names
          if(length(labels)) vlabels(value_cols) <- labels
        }
        res <- c(id_cols, value_cols)

      } else { # Recast Pivot

        # The optimization applied here is to avoid materialization of the "long" id-columns
        # There are two ways to do it, first the long value cast and then wide cast, or many wide casts and row-biding.
        # The complication is that the long cast requires construction of an id-column, which probably can only be efficiently
        # done by creating yet another C-function. Thus I try the wide option first.
        # -> initial benchmarks show that this is also definitely faster than recast from long frame...
        # but presumably because grouping is much faster. If an id is constructed we don't need to group a long frame though...

        # TODO: multiple recast?? -> I think in such cases it would be justifyable to call pivot() 2 times,
        # the syntax with recast could become very complicated

        # (1) Preprocessing Arguments
        names <- proc_names_recast(names, data) # List of 2 elements...
        names1 <- names[[1L]]
        if(length(labels)) {
          labels <- proc_labels_recast(labels, data)
          labels1 <- labels[[1L]]
        } else labels1 <- NULL

        if(is.null(values)) values <- seq_along(data)[-c(ids, names1, labels1)]
        else if(is.null(ids)) ids <- seq_along(data)[-c(names1, labels1, values)]
        # (2) Compute ID Columns
        if(length(ids)) {
          if(sort[1L]) {
            g <- GRP.default(data[ids], sort = TRUE, return.order = FALSE, call = FALSE)
            id_cols <- g[[4L]]
            ng <- g[[1L]]
            g <- g[[2L]]
            attr(g, "N.groups") <- ng
          } else { # Could also use GRP(), but this avoids computing a potentially large and redundant group sizes vector
            g <- group(data[ids], starts = TRUE)
            id_cols <- .Call(C_subsetDT, data, attr(g, "starts"), ids, FALSE)
          }
        } else {
          g <- alloc(1L, fnrow(data)) # TODO: Better create a C-level exemption?? but this is inefficient anyway (row-binding single rows...)
          attr(g, "N.groups") <- 1L
          id_cols <- NULL
        }
        # (3) Compute Names and Labels Columns
        names_g <- GRP(if(length(names1) == 1L && is.null(labels1)) data[[names1]] else data[names1],
                       sort = sort[2L], group.sizes = check.dups, drop = drop, call = FALSE)
        if(length(labels1)) {
          if(check.dups && any(vary <- varying(data[labels1], names_g)))  # See if there are duplicate labels
            stop("The following 'labels' columns vary by 'names': ", paste(names(vary)[vary], collapse = ", "))
          labels1 <- if(length(labels1) == 1L) tochar(Csv(data[[labels1]], names_g$group.starts)) else
            do.call(paste, c(.Call(C_subsetDT, data, names_g$group.starts, labels1, FALSE), list(sep = " - ")))
        }
        g_v <- names_g[[2L]]
        attr(g_v, "N.groups") <- names_g[[1L]]
        names1 <- GRPnames(names_g, sep = "_")

        # (4) Optional duplicates check...
        if(check.dups) {
          ndg <- fndistinct.default(g, names_g, use.g.names = FALSE, na.rm = FALSE, nthreads = nthreads)
          attributes(ndg) <- NULL
          if(!identical(ndg, names_g[[3L]])) {
            ng <- fsumC(ndg, narm = FALSE)
            warning("duplicates detected: there are ", ng, " unique combinations of id- and name-columns, but the data has ", fnrow(data),
                    " rows. This means you have on average ", round(fnrow(data)/ng, 1), " duplicates per id-name-combination. If how = 'recast', pivot() will take the last of those duplicates in first-appearance-order. Consider aggregating your data e.g. using collap() before applying pivot().")
          }
        }

        # (5) Compute Reshaped Values
        save_labels <- !is.null(labels[[2L]])
        vd <- data[values]
        if(save_labels || factor[1L]) {
          namv <- names(vd)
          attributes(vd) <- NULL
        }
        if(!is.character(FUN)) {
          vd <- apply_external_FUN(vd, group(list(g, g_v)), FUN, FUN.args, l1orlst(as.character(substitute(FUN))))
          FUN <- "last"
        }
        value_cols <- lapply(vd, function(x) .Call(C_pivot_wide, g, g_v, x, fill, nthreads, FUN, na.rm))
        if(length(id_cols)) id_cols <- .Call(C_rbindlist, alloc(id_cols, length(value_cols)), FALSE, FALSE, NULL)
        value_cols <- .Call(C_rbindlist, value_cols, FALSE, FALSE, names[[2L]]) # Final column is "variable" name

        names(value_cols) <- c(names[[2L]], names1)
        if(length(labels1)) vlabels(value_cols)[-1L] <- labels1
        else if(length(vd) > 1L) vlabels(value_cols) <- NULL

        # (6) Missing Value Removal
        if(na.rm) { # TODO: better way???
          cc <- whichv(missing_cases(value_cols, prop = 1), FALSE)
          if(length(cc) != fnrow(value_cols)) {
            value_cols <- .Call(C_subsetDT, value_cols, cc, seq_along(value_cols), FALSE)
            id_cols <- .Call(C_subsetDT, id_cols, cc, seq_along(id_cols), FALSE)
          }
        }

        # (7) Properly deal with variable names and labels
        if(save_labels) {
          if(!is.character(labels[[2L]])) stop("label column name supplied in a list needs to be character typed, you passed an object of type: ", typeof(labels[[2L]]))
          labs <- vlabels(vd, use.names = FALSE)
          if(factor[2L]) {
            label_col <- value_cols[[1L]]
            attr(label_col, "levels") <- labs
            oldClass(label_col) <- "factor" # c("factor", "na.included")
          } else label_col <- Csv(labs, value_cols[[1L]])
          label_col <- list(label_col)
          names(label_col) <- labels[[2L]]
          value_cols <- c(value_cols[1L], label_col, value_cols[-1L])
        }

        if(factor[1L]) {
          attr(value_cols[[1L]], "levels") <- namv
          oldClass(value_cols[[1L]]) <- "factor" # c("factor", "na.included")
        } else if(save_labels) value_cols[[1L]] <- Csv(namv, value_cols[[1L]])

        if(length(new_labels <- labels[[3L]])) {
          if(is.null(names(new_labels))) {
            if(length(new_labels) == length(value_cols)) vlabels(value_cols) <- new_labels
            else if(length(new_labels) == 1L+save_labels) vlabels(value_cols)[seq_len(1L+save_labels)] <- new_labels
            else stop("Number of new labels supplied must match either number of new ids (names/label-columns) or total number of new columns in recasted frame. There are ", length(value_cols), " new columns in the frame, of which ", 1L+save_labels, " are ids, and you supplied ", length(new_labels), " new labels. Alternatively, please provide a named vector matching labels to columns.")
          } else vlabels(value_cols)[names(new_labels)] <- new_labels
        }

        res <- if(length(id_cols)) c(id_cols, value_cols) else value_cols
      }
    }

  if(is.null(ad)) return(res) # Redundant ??
  if(any(ad$class == "data.frame")) ad$row.names <- .set_row_names(fnrow(res))
  ad$names <- names(res)
  .Call(C_setattributes, res, ad)
  if(any(ad$class == "data.table")) return(alc(res))
  return(res)
}
