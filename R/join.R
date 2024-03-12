################################
# Implementation of Table Joins
################################

sort_merge_join <- function(x_sorted, table, count = FALSE) {
  ot <- radixorderv(table, decreasing = FALSE, na.last = TRUE)
  .Call(C_sort_merge_join, x_sorted, table, ot, count)
}

multi_match <- function(m, g) .Call(C_multi_match, m, g)

# Modeled after Pandas/Polars:
# https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.join.html
# https://pola-rs.github.io/polars/py-polars/html/reference/dataframe/api/polars.DataFrame.join.html
join <- function(x, y,
                 on = NULL, # union(names(x), names(y)),
                 how = "left",
                 suffix = NULL, # c("_x", "_y")
                 validate = "m:m",  # NULL,
                 multiple = FALSE,
                 sort = FALSE,
                 keep.col.order = TRUE,
                 drop.dup.cols = FALSE,
                 verbose = .op[["verbose"]],
                 column = NULL,
                 attr = NULL, ...) { # method = c("hash", "radix") -> implicit to sort...

  # Initial checks
  if(!is.list(x)) stop("x must be a list")
  if(!is.list(y)) stop("y must be a list")

  # Get names and attributes
  ax <- attributes(x)
  x_name <- as.character(substitute(x))
  if(length(x_name) != 1L || x_name == ".") x_name <- "x" # Piped use
  y_name <- as.character(substitute(y))
  if(length(y_name) != 1L || y_name == ".") y_name <- "y" # Piped use
  oldClass(x) <- NULL
  oldClass(y) <- NULL
  xnam <- names(x)
  ynam <- names(y)
  how <- switch(how, l = "left", r = "right", i = "inner", f = "full", s = "semi", a = "anti", how)

  # Get join columns
  if(is.null(on)) {
    xon <- on <- xnam[xnam %in% ynam]
    if(length(on) == 0L) stop("No matching column names between x and y, please specify columns to join 'on'.")
    if(anyDuplicated.default(on) > 0L) stop("Duplicated join columns: ", paste(on[fduplicated(on)], collapse = ", "), ". Please supply 'on' columns and ensure that each data frame has unique column names.")
    ixon <- match(on, xnam)
    iyon <- match(on, ynam)
  } else {
    if(!is.character(on)) stop("need to provide character 'on'")
    xon <- names(on)
    if(is.null(xon)) xon <- on
    else if(any(miss <- !nzchar(xon))) xon[miss] <- on[miss]
    ixon <- ckmatch(xon, xnam, "Unknown x columns:")
    iyon <- ckmatch(on, ynam, "Unknown y columns:")
  }

  # Matching step
  rjoin <- switch(how, right = TRUE, FALSE)
  count <- verbose || validate != "m:m"

  if(sort) {
    if(rjoin) {
      y <- roworderv(y, cols = iyon, decreasing = FALSE, na.last = TRUE)
      m <- sort_merge_join(y[iyon], x[ixon], count = count)
    } else {
      x <- roworderv(x, cols = ixon, decreasing = FALSE, na.last = TRUE)
      m <- sort_merge_join(x[ixon], y[iyon], count = count)
      if(how == "left" && length(ax[["row.names"]])) ax[["row.names"]] <- attr(x, "row.names")
    }
  } else {
    m <- if(rjoin) fmatch(y[iyon], x[ixon], nomatch = NA_integer_, count = count, ...) else
                   fmatch(x[ixon], y[iyon], nomatch = NA_integer_, count = count, ...)
  }

  # TODO: validate full join...
  switch(validate,
    "m:m" = TRUE,
    "1:1" = {
      c1 <- attr(m, "N.distinct") != length(m) - attr(m, "N.nomatch")
      c2 <- attr(m, "N.groups") != attr(m, "N.distinct") && any_duplicated(if(rjoin) x[ixon] else y[iyon])
      if(rjoin) {
        tmp <- c2
        c2 <- c1
        c1 <- tmp
      }
      if(c1 || c2) stop("Join is not 1:1: ", x_name, " (x) is ", if(c1) "not " else "", "unique on the join columns; ", y_name, " (y) is ", if(c2) "not " else "", "unique on the join columns")
    },
    "1:m" = {
      cond <- if(rjoin) attr(m, "N.groups") != attr(m, "N.distinct") && any_duplicated(x[ixon]) else
              attr(m, "N.distinct") != length(m) - attr(m, "N.nomatch")
      if(cond) stop("Join is not 1:m: ", x_name, " (x) is not unique on the join columns")
    },
    "m:1" = {
      cond <- if(rjoin) attr(m, "N.distinct") != length(m) - attr(m, "N.nomatch") else
              attr(m, "N.groups") != attr(m, "N.distinct") && any_duplicated(y[iyon])
      if(cond) stop("Join is not m:1: ", y_name, " (y) is not unique on the join columns")
    },
    stop("validate must be one of '1:1', '1:m', 'm:1' or 'm:m'")
  )

  if(multiple) {
    g <- group(if(rjoin) x[ixon] else y[iyon], group.sizes = TRUE)
    m <- multi_match(m, g)
    if(is.list(m)) {
      if(rjoin) y <- .Call(C_subsetDT, y, m[[1L]], seq_along(y), FALSE)
      else x <- .Call(C_subsetDT, x, m[[1L]], seq_along(x), FALSE)
      m <- m[[2L]]
      if(how == "left" && length(ax[["row.names"]])) ax[["row.names"]] <- .set_row_names(length(m))
    }
  }

  if(verbose) {
    nx <- length(m) - attr(m, "N.nomatch")
    ny <- attr(m, "N.distinct")
    Ny <- attr(m, "N.groups")
    if(verbose == 2L) {
      cin_x <- paste0(xon, ":", vclasses(x[ixon], FALSE))
      cin_y <- paste0(on, ":", vclasses(y[iyon], FALSE))
    } else {
      cin_x <- xon
      cin_y <- on
    }
    xstat <- paste0(nx, "/", length(m), " (", signif(nx/length(m)*100, 3), "%)")
    ystat <- paste0(ny, "/", Ny, " (", signif(ny/Ny*100, 3), "%)")
    if(rjoin) {
      tmp <- ystat
      ystat <- xstat
      xstat <- tmp
    }
    cat(how, " join: ",
        x_name, "[", paste(cin_x, collapse = ", "), "] ",
        xstat, " <", validate ,"> ",
        y_name, "[", paste(cin_y, collapse = ", "), "] ",
        ystat, "\n", sep = "")
  }

  # Check for duplicate columns and suffix as needed
  if(any(nm <- match(ynam[-iyon], xnam, nomatch = 0L))) {
    nnm <- nm != 0L
    nam <- xnam[nm[nnm]]
    if(is.character(drop.dup.cols) || drop.dup.cols) {
      switch(drop.dup.cols,
        y = {
          rmyi <- logical(length(ynam))
          rmyi[-iyon][nnm] <- TRUE
          y[rmyi] <- NULL
          ynam <- names(y)
          tmp <- rmyi
          tmp[iyon] <- TRUE
          iyon <- which(tmp[!rmyi])
          if(verbose) cat("duplicate columns: ", paste(nam, collapse = ", "), " => dropped from y\n", sep = "")
        },
        x = {
          x[nm[nnm]] <- NULL
          tmp <- logical(length(xnam))
          xnam <- names(x)
          tmp[ixon] <- TRUE
          ixon <- which(tmp[-nm[nnm]])
          if(verbose) cat("duplicate columns: ", paste(nam, collapse = ", "), " => dropped from x\n", sep = "")
        },
        stop("drop.dup.cols needs to be 'y', 'x', or TRUE")
      )
    } else {
      if(length(suffix) <= 1L) { # Only appends y with name
        if(is.null(suffix)) suffix <- paste0("_", y_name)
        names(y)[-iyon][nnm] <- paste0(nam, suffix)
      } else {
        names(x)[nm[nnm]] <- paste0(nam, suffix[[1L]]) # if(suffix[[1L]] != "") ??
        names(y)[-iyon][nnm] <- paste0(nam, suffix[[2L]])
      }
      if(verbose) cat("duplicate columns: ", paste(nam, collapse = ", "), " => renamed using suffix ",
          if(length(suffix) == 1L) paste0("'", suffix, "' for y") else paste0("'", suffix[[1L]], "' for x and '", suffix[[2L]], "' for y"), "\n", sep = "")
    }
  }

  # Core: do the joins
  res <- switch(how,
    left = {
      y_res <- .Call(C_subsetDT, y, m, seq_along(y)[-iyon], if(count) attr(m, "N.nomatch") else TRUE)
      c(x, y_res)
    },
    inner = {
      anyna <- if(count) attr(m, "N.nomatch") > 0L else anyNA(m)
      if(anyna) {
        x_ind <- whichNA(m, invert = TRUE)
        x <- .Call(C_subsetDT, x, x_ind, seq_along(x), FALSE)
        m <- na_rm(m)
        # rn <- ax[["row.names"]] # TODO: Works inside switch??
        # if(length(rn)) ax[["row.names"]] <- if(is.numeric(rn) || is.null(rn) || rn[1L] == "1")
        #             .set_row_names(length(x_ind)) else Csv(rn, x_ind)
      }
      y_res <- .Call(C_subsetDT, y, m, seq_along(y)[-iyon], FALSE)
      c(x, y_res)
    },
    full = {
      cond <- !count || attr(m, "N.distinct") != attr(m, "N.groups")
      if(cond) {
        um <- if(!count || length(m)-attr(m, "N.distinct")-attr(m, "N.nomatch") != 0L)
          .Call(C_funique, m) else m # This gets the rows of table matched
        if(!count || attr(m, "N.nomatch")) um <- na_rm(um)
        if(count) tsize <- attr(m, "N.groups")
        else {
          tsize <- fnrow(y)
          cond <- length(um) != tsize
        }
      }
      if(cond) { # TODO: special case ? 1 distinct value etc.??
        tind <- seq_len(tsize)[-um] # TODO: Table may not be unique.
        res_nrow <- length(m) + length(tind)
        x_res <- .Call(C_subsetDT, x, seq_len(res_nrow), seq_along(x)[-ixon], TRUE)  # Need check here because oversize indices !!
        y_res <- .Call(C_subsetDT, y, vec(list(m, tind)), seq_along(y)[-iyon], TRUE) # Need check here because oversize indices !!
        on_res <- .Call(C_rbindlist, list(x[ixon], .Call(C_subsetDT, y, tind, iyon, FALSE)), FALSE, FALSE, NULL)
        # if(length(ax[["row.names"]])) ax[["row.names"]] <- .set_row_names(res_nrow)
        if(keep.col.order) {
          if(length(x_res)) add_vars(x_res, pos = ixon) <- on_res
          else x_res <- on_res
          c(x_res, y_res)
        } else {
          keep.col.order <- 2L # has global effects !!
          c(on_res, x_res, y_res)
        }
      } else { # If all elements of table are matched, this is simply a left join
        how <- "left"
        y_res <- .Call(C_subsetDT, y, m, seq_along(y)[-iyon], if(count) attr(m, "N.nomatch") else TRUE) # anyNA(um) ??
        c(x, y_res)
      }
    },
    right = {
      x_res <- .Call(C_subsetDT, x, m, seq_along(x)[-ixon], if(count) attr(m, "N.nomatch") else TRUE)
      # if(length(ax[["row.names"]])) ax[["row.names"]] <- .set_row_names(length(m))
      y_on <- y[iyon]
      names(y_on) <- xon
      if(keep.col.order) {
        if(length(x_res)) add_vars(x_res, pos = ixon) <- y_on
        else x_res <- y_on
        c(x_res, y[-iyon])
      } else {
        keep.col.order <- 2L # has global effects !!
        c(y_on, x_res, y[-iyon])
      }
    },
    semi = { # = return rows in x that have matching values in y
      anyna <- if(count) attr(m, "N.nomatch") > 0L else anyNA(m)
      if(anyna) {
        x_ind <- whichNA(m, invert = TRUE)
        # rn <- ax[["row.names"]] # TODO: Works inside switch??
        # if(length(rn)) ax[["row.names"]] <- if(is.numeric(rn) || is.null(rn) || rn[1L] == "1")
        #            .set_row_names(x_ind) else Csv(rn, x_ind)
        .Call(C_subsetDT, x, x_ind, seq_along(x), FALSE)
      } else x
    },
    # = return rows in x that have no matching values in y
    anti = .Call(C_subsetDT, x, whichNA(m), seq_along(x), FALSE),
    stop("Unknown join method: ", how)
  )

  # Join column and reordering
  if(length(column)) {
    if(is.list(column)) {
      lev <- column[[2L]]
      column <- column[[1L]]
      x_name <- lev[[1L]]
      y_name <- lev[[2L]]
      matched <- lev[[3L]]
    } else matched <- "matched"
    # TODO: better?
    # matched <- paste0(y_name, "_", y_name)
    mc <- switch(how,
                 left = structure(is.na(m) + 1L, levels = c(matched, x_name), class = c("factor", "na.included")),
                 right = structure(is.na(m) + 1L, levels = c(matched, y_name), class = c("factor", "na.included")),
                 full = structure(vec(list(is.na(m) + 1L, alloc(3L, fnrow(res)-length(m)))), levels = c(matched, x_name, y_name), class = c("factor", "na.included")),
                 inner =, semi = structure(alloc(1L, fnrow(res)), levels = matched, class = c("factor", "na.included")),
                 anti = structure(alloc(1L, fnrow(res)), levels = x_name, class = c("factor", "na.included")))
    attr(mc, "on.cols") <- `names<-`(list(xon, `names<-`(on, NULL)), c(x_name, y_name))
    mc_name <- if(is.character(column)) column else ".join"
    if(keep.col.order == 1L) res[[mc_name]] <- mc
    else res <- c(res[ixon], `names<-`(list(mc), mc_name), res[-ixon])
  } else if(!keep.col.order) res <- c(res[ixon], res[-ixon])

  # Final steps
  if(length(attr)) ax[[if(is.character(attr)) attr else "join.match"]] <- list(call = match.call(),
                                                                               on.cols = list(x = xon, y = `names<-`(on, NULL)),
                                                                               match = m) # TODO: sort merge join also report o?
  if(sort && how == "full") res <- roworderv(res, cols = xon)
  if(how != "left" && length(ax[["row.names"]])) ax[["row.names"]] <- .set_row_names(fnrow(res))
  ax[["names"]] <- names(res)
  .Call(C_setattributes, res, ax)
  if(any(ax$class == "data.table")) return(alc(res))
  return(res)
}
