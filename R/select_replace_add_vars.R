

fselect <- function(x, ..., return = "data") { # This also takes names and indices ....
  ax <- attributes(x)
  oldClass(x) <- NULL # attributes ?
  nam <- names(x)
  nl <- `names<-`(as.vector(seq_along(x), "list"), nam)
  vars <- eval(substitute(c(...)), nl, parent.frame())
  # if(!is.integer(vars)) stop(paste0("Unknown columns: ", .c(...))) # if(!is.integer(vars) || max(vars) > length(nam)) # nah, a bit redundant..
  nam_vars <- names(vars)
  if(is.character(vars)) vars <- ckmatch(vars, nam)
  if(length(nam_vars)) { # Allow renaming during selection
    nonmiss <- nzchar(nam_vars)
    nam[vars[nonmiss]] <- nam_vars[nonmiss]
  }
  # if(!is.numeric(vars)) stop("... needs to be column names, or character / integer / logical vectors")
  switch(return,
         data = setAttributes(x[vars], `[[<-`(ax, "names", nam[vars])), # Also Improvements in code below ?
         names = nam[vars],
         indices = vars,
         named_indices = `names<-`(vars, nam[vars]),
         logical = `[<-`(logical(length(x)), vars, TRUE),
         named_logical = `names<-`(`[<-`(logical(length(x)), vars, TRUE), nam))
}

# or slt sel, selt, sct -> shortcut ?
slt <- fselect # good, consistent

# fselect(GGDC10S, Country, AGR:SUM)
# fselect(GGDC10S, Variable == "VA" & Year > 1990, Country, Year, AGR:SUM) -> why no error ?? first argument is just ignored ... ??



"fselect<-" <- function(x, ..., value) {
  clx <- oldClass(x)
  oldClass(x) <- NULL
  nl <- `names<-`(as.vector(seq_along(x), "list"), names(x))
  vars <- eval(substitute(c(...)), nl, parent.frame())
  if(is.character(vars)) vars <- ckmatch(vars, names(x))
  # if(!is.numeric(vars)) stop("... needs to be column names, or character / integer / logical vectors")
  # if(!is.integer(vars)) stop(paste0("Unknown columns: ", .c(...)))
  if(is.null(value)) {
    if(!length(vars)) return(`oldClass<-`(x, clx))
    ax <- attributes(x)
    ax[["names"]] <- ax[["names"]][-vars]
    ax[["class"]] <- clx
    return(setAttributes(.subset(x, -vars), ax))
  }
  if(is.list(value)) {
    oldClass(value) <- NULL
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(value) != length(vars)) stop("NCOL(value) must match selected variables")
    x[vars] <- value
    if(length(nam <- names(value))) names(x)[vars] <- nam
  } else {
    if(NROW(unclass(value)) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(vars) != 1L) stop("NCOL(value) must match selected variables")
    x[[vars]] <- value
  }
  return(`oldClass<-`(x, clx))
}

"slt<-" <- `fselect<-`


# STD(fselect(GGDC10S, Country, Variable, Year, AGR:SUM))

# Idea: also do this for replacement functions, replacing characters renames, replacong number reorders, replacing 3 does renaming and reordering?

get_vars_FUN <- function(x, FUN, return = "data")
  switch(return, data = colsubsetFUN(x, FUN),
         names = attr(x, "names")[vapply(`attributes<-`(x, NULL), FUN, TRUE)],
         indices = which(vapply(`attributes<-`(x, NULL), FUN, TRUE)),
         named_indices = which(`names<-`(vapply(`attributes<-`(x, NULL), FUN, TRUE), attr(x, "names"))),
         logical = vapply(`attributes<-`(x, NULL), FUN, TRUE),
         named_logical = `names<-`(vapply(`attributes<-`(x, NULL), FUN, TRUE), attr(x, "names")),
         stop("Unknown return option!"))

"get_vars_FUN<-" <- function(x, FUN, value) {
  ind <- which(vapply(`attributes<-`(x, NULL), FUN, TRUE))
  if(is.null(value)) {
    if(!length(ind)) return(x)
    ax <- attributes(x)
    ax[["names"]] <- ax[["names"]][-ind] # fastest ? -> yes !
    return(setAttributes(.subset(x, -ind), ax))
  }
  clx <- oldClass(x)
  oldClass(x) <- NULL
  if(is.list(value)) {
    oldClass(value) <- NULL # fastest ?? if(is.object(value)) oldClass(value) <- NULL ??
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(value) != length(ind)) stop("NCOL(value) must match selected variables") # length(num_vars(x))
    x[ind] <- value
    if(length(nam <- names(value))) names(x)[ind] <- nam #  == length(ind)
  } else {
    if(NROW(unclass(value)) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(ind) != 1L) stop("NCOL(value) must match selected variables") # length(num_vars(x))
    x[[ind]] <- value
  }
  return(`oldClass<-`(x, clx))
}

num_vars <- function(x, return = "data") get_vars_FUN(x, is.numeric, return)
nv <- num_vars

"num_vars<-" <- function(x, value) `get_vars_FUN<-`(x, is.numeric, value)
"nv<-" <- `num_vars<-`

char_vars <- function(x, return = "data") get_vars_FUN(x, is.character, return)
"char_vars<-" <- function(x, value) `get_vars_FUN<-`(x, is.character, value)

fact_vars <- function(x, return = "data") get_vars_FUN(x, is.factor, return)
"fact_vars<-" <- function(x, value) `get_vars_FUN<-`(x, is.factor, value)

logi_vars <- function(x, return = "data") get_vars_FUN(x, is.logical, return)
"logi_vars<-" <- function(x, value) `get_vars_FUN<-`(x, is.logical, value)

Date_vars <- function(x, return = "data") get_vars_FUN(x, is.Date, return)
"Date_vars<-" <- function(x, value) `get_vars_FUN<-`(x, is.Date, value)

cat_vars <- function(x, return = "data") {
  switch(return, data = fcolsubset(x, !vapply(`attributes<-`(x, NULL), is.numeric, TRUE)),
         names = attr(x, "names")[!vapply(`attributes<-`(x, NULL), is.numeric, TRUE)],
         indices = which(!vapply(`attributes<-`(x, NULL), is.numeric, TRUE)),
         named_indices = which(`names<-`(!vapply(`attributes<-`(x, NULL), is.numeric, TRUE), attr(x, "names"))),
         logical = !vapply(`attributes<-`(x, NULL), is.numeric, TRUE),
         named_logical = `names<-`(!vapply(`attributes<-`(x, NULL), is.numeric, TRUE), attr(x, "names")),
         stop("Unknown return option!"))

}
"cat_vars<-" <- function(x, value) {
  ind <- which(!vapply(`attributes<-`(x, NULL), is.numeric, TRUE))
  if(is.null(value)) {
    if(!length(ind)) return(x)
    ax <- attributes(x)
    ax[["names"]] <- ax[["names"]][-ind]
    return(setAttributes(.subset(x, -ind), ax))
  }
  clx <- oldClass(x)
  oldClass(x) <- NULL
  if(is.list(value)) {
    oldClass(value) <- NULL
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(value) != length(ind)) stop("NCOL(value) must match length(cat_vars(x))")
    x[ind] <- value
    if(length(nam <- names(value))) names(x)[ind] <- nam
  } else {
    if(NROW(unclass(value)) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(ind) != 1L) stop("NCOL(value) must match length(cat_vars(x))")
    x[[ind]] <- value
  }
  return(`oldClass<-`(x, clx))
}


get_vars <- function(x, vars, return = "data", regex = FALSE, ...) {
 if(!regex && !missing(...)) unused_arg_action(match.call(), ...)
 if(is.function(vars)) return(get_vars_FUN(x, vars, return))
 switch(return,
        data = colsubset(x, if(regex) rgrep(vars, attr(x, "names"), ...) else vars),
        names = if(regex) rgrep(vars, attr(x, "names"), value = TRUE, ...) else attr(x, "names")[vars], # error?
        indices = if(regex) rgrep(vars, attr(x, "names"), ...) else if(is.character(vars)) ckmatch(vars, attr(x, "names")) else
          stop("For indices, vars must be a function, character names or a regular expression"),
        named_indices = if(is.character(vars)) {
                           nam <- attr(x, "names")
                           ind <- if(regex) rgrep(vars, nam, ...) else ckmatch(vars, nam)
                           `names<-`(ind, nam[ind])
                        } else stop("For named indices, vars must be a function, character names or a regular expression"),
        logical = if(regex) rgrepl(vars, attr(x, "names"), ...) else if(is.character(vars)) attr(x, "names") %in% vars else # rgrepl ?
          stop("For logical, vars must be a function, character names or a regular expression"),
        named_logical = if(is.character(vars)) {
          nam <- attr(x, "names")
          `names<-`(if(regex) rgrepl(vars, nam, ...) else nam %in% vars, nam)
        } else stop("For named logical, vars must be a function, character names or a regular expression"),
        stop("Unknown return option!"))
}

gv <- function(x, vars, return = "data", ...) {
  if(!missing(...)) {
    warning("Please use the new shortcut 'gvr' for regex column selection.")
    return(get_vars(x, vars, return, ...))
  }
  if(is.function(vars)) return(get_vars_FUN(x, vars, return))
  switch(return,
         data = colsubset(x, vars),
         names = attr(x, "names")[vars], # error?
         indices = if(is.character(vars)) ckmatch(vars, attr(x, "names")) else stop("For indices, vars must be a function or character names"),
         named_indices = if(is.character(vars)) {
           nam <- attr(x, "names")
           ind <- ckmatch(vars, nam)
           `names<-`(ind, nam[ind])
         } else stop("For named indices, vars must be a function, character names or a regular expression"),
         logical = if(is.character(vars)) attr(x, "names") %in% vars else stop("For logical, vars must be a function, character names or a regular expression"),
         named_logical = if(is.character(vars)) {
           nam <- attr(x, "names")
           `names<-`(nam %in% vars, nam)
         } else stop("For named logical, vars must be a function, character names or a regular expression"),
         stop("Unknown return option!"))
}

gvr <- function(x, vars, return = "data", ...) {
    switch(return,
           data = fcolsubset(x, rgrep(vars, attr(x, "names"), ...)),
           names = rgrep(vars, attr(x, "names"), value = TRUE, ...),
           indices = rgrep(vars, attr(x, "names"), ...),
           named_indices = {
             nam <- attr(x, "names")
             ind <- rgrep(vars, nam, ...)
             `names<-`(ind, nam[ind])
           },
           logical = rgrepl(vars, attr(x, "names"), ...),
           named_logical = {
             nam <- attr(x, "names")
             `names<-`(rgrepl(vars, nam, ...), nam)
           },
           stop("Unknown return option!"))
  }



"get_vars<-" <- function(x, vars, regex = FALSE, ..., value) {
  clx <- oldClass(x)
  oldClass(x) <- NULL
  if(is.numeric(vars)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(max(abs(vars)) > length(x)) stop("Column index out of range abs(1:length(x))")
    if(any(vars < 0)) vars <- seq_along(x)[vars]
  } else if(is.logical(vars)) {
    if(!missing(...)) unused_arg_action(match.call(), ...)
    if(length(vars) != length(x)) stop("Logical subsetting vector must match length(x)")
    vars <- which(vars)
  } else {
    if(!regex && !missing(...)) unused_arg_action(match.call(), ...)
    vars <- if(is.function(vars)) which(vapply(unattrib(x), vars, TRUE)) else if(regex) # if(is.character(vars) && regex): redundant...
            rgrep(vars, names(x), ...) else ckmatch(vars, names(x))
  }
  if(is.null(value)) {
    if(!length(vars)) return(`oldClass<-`(x, clx))
    ax <- attributes(x)
    ax[["names"]] <- names(x)[-vars]
    ax[["class"]] <- clx
    return(setAttributes(x[-vars], ax))
  }
  if(is.list(value)) {
    oldClass(value) <- NULL
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(value) != length(vars)) stop("NCOL(value) must match length(vars)")
    x[vars] <- value
    if(length(nam <- names(value))) names(x)[vars] <- nam #  == length(vars)
  } else {
    if(NROW(unclass(value)) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(vars) != 1L) stop("NCOL(value) must match length(vars)")
    x[[vars]] <- value
  }
  return(`oldClass<-`(x, clx))
}

"gv<-" <- function(x, vars, ..., value) {
  if(!missing(...)) {
    warning("Please use the new shortcut 'gvr<-' for regex column replacement.")
    return(`get_vars<-`(x, vars, ..., value = value))
  }
  clx <- oldClass(x)
  oldClass(x) <- NULL
  if(is.numeric(vars)) {
    if(max(abs(vars)) > length(x)) stop("Column index out of range abs(1:length(x))")
    if(any(vars < 0)) vars <- seq_along(x)[vars]
  } else if(is.logical(vars)) {
    if(length(vars) != length(x)) stop("Logical subsetting vector must match length(x)")
    vars <- which(vars)
  } else vars <- if(is.function(vars)) which(vapply(unattrib(x), vars, TRUE)) else ckmatch(vars, names(x))

  if(is.null(value)) {
    if(!length(vars)) return(`oldClass<-`(x, clx))
    ax <- attributes(x)
    ax[["names"]] <- names(x)[-vars]
    ax[["class"]] <- clx
    return(setAttributes(x[-vars], ax))
  }
  if(is.list(value)) {
    oldClass(value) <- NULL
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(value) != length(vars)) stop("NCOL(value) must match length(vars)")
    x[vars] <- value
    if(length(nam <- names(value))) names(x)[vars] <- nam #  == length(vars)
  } else {
    if(NROW(unclass(value)) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(vars) != 1L) stop("NCOL(value) must match length(vars)")
    x[[vars]] <- value
  }
  return(`oldClass<-`(x, clx))
}
"gvr<-" <- function(x, vars, ..., value) {
  clx <- oldClass(x)
  oldClass(x) <- NULL
  vars <- rgrep(vars, names(x), ...)

  if(is.null(value)) {
    if(!length(vars)) return(`oldClass<-`(x, clx))
    ax <- attributes(x)
    ax[["names"]] <- names(x)[-vars]
    ax[["class"]] <- clx
    return(setAttributes(x[-vars], ax))
  }
  if(is.list(value)) {
    oldClass(value) <- NULL
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(value) != length(vars)) stop("NCOL(value) must match length(vars)")
    x[vars] <- value
    if(length(nam <- names(value))) names(x)[vars] <- nam #  == length(vars)
  } else {
    if(NROW(unclass(value)) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(vars) != 1L) stop("NCOL(value) must match length(vars)")
    x[[vars]] <- value
  }
  return(`oldClass<-`(x, clx))
}

# Make faster ?
"add_vars<-" <- function(x, pos = "end", value) {
  ax <- attributes(x)
  attributes(x) <- NULL
  lx <- length(x)
  if(is.list(value)) {
    oldClass(value) <- NULL # fastest ?
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    # res <- c(x, value)  # FASTER than commented out below
    if(is.character(pos)) {
      if(pos == "end") {
        ax[["names"]] <- if(length(nam <- names(value)))  c(ax[["names"]], nam) else
          c(ax[["names"]], paste0("V", seq(lx+1L, lx+length(value))))
        return(setAttributes(c(x, value), ax))
      } else if(pos != "front") stop("pos needs to be 'end', 'front' or a suitable numeric / integer vector of positions!")
        ax[["names"]] <- if(length(nam <- names(value)))  c(nam, ax[["names"]]) else
          c(paste0("V", seq_along(value)), ax[["names"]])
        return(setAttributes(c(value, x), ax))
    }
    lv <- length(value)
    tl <- lv+lx
    if(!is.numeric(pos) || length(pos) != lv || max(pos) > tl) stop("pos needs to be 'end', 'front' or a suitable numeric / integer vector of positions!")
    o <- forder.int(c(seq_len(tl)[-pos], pos))
    ax[["names"]] <- if(length(nam <- names(value)))  c(ax[["names"]], nam)[o] else
        c(ax[["names"]], paste0("V", pos))[o] # FASTER THIS WAY? -> It seems so...
    return(setAttributes(c(x, value)[o], ax)) # fastest ?? use setcolorder ? (probably not )
    # ind <- seq(lx+1L, lx+length(value))
    # x[ind] <- value  # FASTER than simply using x[names(value)] <- value ? -> Yes !
    # ax[["names"]] <- if(length(nam <- names(value)))  c(ax[["names"]], nam) else
    #   c(ax[["names"]], paste0("V", ind))
  } else {
    if(NROW(value) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    # res <- c(x, list(value)) # FASTER than below ? -> Nope
    # ax[["names"]] <- c(ax[["names"]], paste0("V", lx+1L))
    nam <- l1orlst(as.character(substitute(value)))
    if(is.character(pos)) {
      if(pos == "end") {
        x[[lx+1L]] <- value
        ax[["names"]] <- c(ax[["names"]], nam) # paste0("V", lx+1L)
        return(setAttributes(x, ax))
      } else if(pos != "front") stop("pos needs to be 'end', 'front' or a suitable numeric / integer vector of positions!")
      ax[["names"]] <- c(nam, ax[["names"]])
      return(setAttributes(c(list(value), x), ax))
    }
    if(!is.numeric(pos) || length(pos) > 1L || pos > lx+1L) stop("pos needs to be 'end', 'front' or a suitable numeric / integer vector of positions!")
    o <- forder.int(c(1:lx, pos-1L))
    ax[["names"]] <- c(ax[["names"]], nam)[o]
    return(setAttributes(c(x, list(value))[o], ax))
  }
}
"av<-" <- `add_vars<-`

add_vars <- function(x, ..., pos = "end") {
  if(...length() == 1L) return(`add_vars<-`(x, pos, ...))
  l <- c(...)
  if(!all(fnrow2(x) == lengths(l, FALSE))) stop("if multiple arguments are passed to '...', each needs to be a data.frame/list with column-lengths matching nrow(x)")
  return(`add_vars<-`(x, pos, l)) # very minimal ! Doesn't work for vectors etc !
}
av <- add_vars



# Previous Versions:

# fselect(GGDC10S, Country, AGR:SUM) -> doesn't work !!
# "fselect<-" <- function(x, ..., value) { # This also takes names and indices ....
#   e <- substitute(list(...))
#   ax <- attributes(x)
#   if(length(e) == 2L) { # only one ... could be a sequence with :
#     oldClass(x) <- NULL
#     nam <- names(x)
#     nl <- `names<-`(as.vector(seq_along(x), "list"), nam)
#     vars <- eval(e[[2L]], nl, parent.frame())
#     ax[["names"]] <- nam[vars]
#     return(setAttributes(x[vars], ax))
#   } else {
#     res <- eval(e, x, parent.frame()) # gives errors for unknown !! great !!
#     ax[["names"]] <- all.vars(e)
#     return(setAttributes(res, ax))
#   }
# }
#

# "num_vars<-" <- function(x, value) {
#   ax <- attributes(x)
#   attributes(x) <- NULL
#   ind <- which(vapply(x, is.numeric, TRUE, USE.NAMES = FALSE))
#   if(is.list(value)) {
#     oldClass(value) <- NULL
#     if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
#     if(length(value) != length(ind)) stop("NCOL(value) must match length(char_vars(x))")
#     x[ind] <- value
#     if(length(nam <- names(value))) ax[["names"]][ind] <- nam
#   } else if(is.null(value)) {
#     if(!length(ind)) return(setAttributes(x, ax))
#     ax[["names"]] <- ax[["names"]][-ind]
#     return(setAttributes(x[-ind], ax))
#   } else {
#     if(NROW(value) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
#     if(length(ind) != 1L) stop("NCOL(value) must match length(char_vars(x))")
#     x[[ind]] <- value
#   }
#   return(setAttributes(x, ax))
# }

# "num_vars<-" <- function(x, value) {
#     ax <- attributes(x)
#     ilv <- is.list(value)
#     d <- dim(value)
#     if(NROW2(value, d) != nrow(x)) stop("NROW(value) must match nrow(x)")
#     attributes(x) <- NULL # vapply without attributes is faster !!!. if Fail, data unchanged !!
#     ind <- which(vapply(x, is.numeric, TRUE, USE.NAMES = FALSE))
#     if(NCOL2(d, ilv) != length(ind)) stop("NCOL(value) must match length(num_vars(x))")
#     if(ilv) x[ind] <- value else x[[ind]] <- value #   x[ind] <- if(ilv) value else list(value)
#     if(ilv && length(nam <- attr(value, "names")) == length(ind)) # does stop execution of if statement if ilv is false !!!!!
#        ax[["names"]][ind] <- nam
#     return(setAttributes(x, ax))
# }


# "get_vars<-" <- function(x, vars, regex = FALSE, ..., value) {
#   ax <- attributes(x)
#   attributes(x) <- NULL
#   if(is.numeric(vars)) {
#     if(!missing(...)) stop("Unknown argument ", dotstostr(...))
#     if(max(abs(vars)) > length(x)) stop("Column index out of range abs(1:length(x))")
#     if(any(vars < 0)) vars <- seq_along(x)[vars]
#   } else if(is.logical(vars)) {
#     if(!missing(...)) stop("Unknown argument ", dotstostr(...))
#     if(length(vars) != length(x)) stop("Logical subsetting vector must match length(x)")
#     vars <- which(vars)
#   } else {
#     if(!regex && !missing(...)) stop("Unknown argument ", dotstostr(...))
#     vars <- if(is.function(vars)) which(vapply(x, vars, TRUE, USE.NAMES = FALSE)) else if(regex) # if(is.character(vars) && regex): redundant...
#       rgrep(vars, ax[["names"]], ...) else ckmatch(vars, ax[["names"]])
#   }
#   if(is.list(value)) {
#     oldClass(value) <- NULL # fastest ??
#     if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
#     if(length(value) != length(vars)) stop("NCOL(value) must match length(vars)")
#     x[vars] <- value
#     if(length(nam <- names(value))) ax[["names"]][vars] <- nam #  == length(vars)
#   } else if(is.null(value)) {
#     if(!length(vars)) return(setAttributes(x, ax))
#     ax[["names"]] <- ax[["names"]][-vars] # fastest ??? -> Yes !! This is slower: x[vars] <- NULL
#     return(setAttributes(x[-vars], ax))
#   } else {
#     if(NROW(unclass(value)) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
#     if(length(vars) != 1L) stop("NCOL(value) must match length(vars)")
#     x[[vars]] <- value
#   }
#   return(setAttributes(x, ax))
# }

# "get_vars<-" <- function(x, vars, regex = FALSE, ..., value) {
#   ax <- attributes(x)
#   ilv <- is.list(value)
#   d <- dim(value)
#   # Slightly faster but a bit less secure..
#   # attributes(x) <- NULL
#   # if(NROW2(value, d) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
#   if(NROW2(value, d) != nrow(x)) stop("NROW(value) must match nrow(x)")
#   attributes(x) <- NULL
#   if(is.numeric(vars)) {
#     if(!missing(...)) stop("Unknown argument ", dotstostr(...))
#     if(max(abs(vars)) > length(x)) stop("Column index out of range abs(1:length(x))")
#     if(any(vars < 0)) vars <- seq_along(x)[vars]
#   } else if(is.logical(vars)) {
#     if(!missing(...)) stop("Unknown argument ", dotstostr(...))
#     if(length(vars) != length(x)) stop("Logical subsetting vector must match length(x)")
#     vars <- which(vars)
#   } else {
#     if(!regex && !missing(...)) stop("Unknown argument ", dotstostr(...))
#     vars <- if(is.function(vars)) which(vapply(x, vars, TRUE, USE.NAMES = FALSE)) else if(is.character(vars) && regex)
#            rgrep(vars, ax[["names"]], ...) else ckmatch(vars, ax[["names"]])
#   }
#   if(NCOL2(d, ilv) != length(vars)) stop("NCOL(value) must match length(vars)")
#   if(ilv) x[vars] <- value else x[[vars]] <- value
#   if(ilv && length(nam <- attr(value, "names")) == length(vars)) ax[["names"]][vars] <- nam
#   return(setAttributes(x, ax))
# }

# "add_vars<-" <- function(x, value) {
#   ax <- attributes(x)
#   ilv <- is.list(value)
#   d <- dim(value)
#   if(NROW2(value, d) != nrow(x)) stop("NROW(value) must match nrow(x)")
#   attributes(x) <- NULL # vapply without attributes is faster !!!. if Fail, data unchanged !!
#   lx <- length(x)
#   if(ilv) {
#     ind <- seq(lx+1L, lx+length(value))
#     x[ind] <- value  # FATER than simply using x[names(value)} <- value ?
#     ax[["names"]] <- if(length(nam <- attr(value, "names")))  c(ax[["names"]], nam) else
#                      c(ax[["names"]], paste0("V", ind))
#   } else {
#     x[[lx+1L]] <- value
#     ax[["names"]] <- c(ax[["names"]], paste0("V", lx+1L))
#   }
#   return(setAttributes(x, ax))
# }


# Exercises:
# repl <- function(x)x
# `repl<-` <- function(x, value) {
#   x <- value
#   x
# }
# repl(x)[2] <- 4 # Works!!
# http://adv-r.had.co.nz/Functions.html#special-calls

# This works because the expression names(x)[2] <- "two" is evaluated as if you had written:

#`*tmp*` <- names(x)
#`*tmp*`[2] <- "two"
#names(x) <- `*tmp*`
