# Possibly option return.names onstead of data or return = 0,1,2? for data,names,indices But names(num.vars(...)) is a.so good, compare speed!!
# also do this for replacement functions!!!, replacing characters renames, replacong number reorders, replacing 3 does renaming and reordering!!
# use vapply ??

# Todo: Could make methods for dplyr grouped_df and pdata.frame which keep identifier columns..
# -> NaH, for dplyr there is group_keys

# Use _ instead of . (because of classes) or no gap at all ??
# _ is good !! looks like dplyr and good practice !! (maybe there will be a class 'vars' one day )

# Note: using setAttributes gives dramatic speed-up on large data!!
# Also:: match.arg takes 10 microseconds !!!, the solution you have now is better and also allows numeric indices !!
# # could do error if index is out of range or logical vector too long !!

num_vars <- function(x, return = c("data","names","indices","named_indices")) {
  switch(return[1L], data = colsubsetFUN(x, is.numeric), # switch(match.arg(return)
         names = names(which(vapply(x, is.numeric, TRUE))),
         indices = which(vapply(x, is.numeric, TRUE, USE.NAMES = FALSE)),
         named_indices = which(vapply(x, is.numeric, TRUE)),
         stop("Unknown return option!"))
}
nv <- num_vars
"num_vars<-" <- function(x, value) {
  ax <- attributes(x)
  attributes(x) <- NULL
  ind <- which(vapply(x, is.numeric, TRUE, USE.NAMES = FALSE))
  if(is.list(value)) {
    class(value) <- NULL # fastest ??
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(value) != length(ind)) stop("NCOL(value) must match length(num_vars(x))")
    x[ind] <- value
    if(length(nam <- names(value))) ax[["names"]][ind] <- nam #  == length(ind)
  } else if(is.null(value)) {
    ax[["names"]] <- ax[["names"]][-ind] # fastest ??? -> yes !!!
    return(setAttributes(x[-ind], ax))
  } else {
    if(NROW(value) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(ind) != 1L) stop("NCOL(value) must match length(num_vars(x))")
    x[[ind]] <- value
  }
  return(setAttributes(x, ax))
}

"nv<-" <- `num_vars<-`

# Previous Version:
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

cat_vars <- function(x, return = c("data","names","indices","named_indices")) {
  switch(return[1L], data = colsubset(x, !vapply(x, is.numeric, TRUE)), # a lot faster than using is.categorical
         names = names(which(!vapply(x, is.numeric, TRUE))),
         indices = which(!vapply(x, is.numeric, TRUE, USE.NAMES = FALSE)),
         named_indices = which(!vapply(x, is.numeric, TRUE)),
         stop("Unknown return option!"))
}
"cat_vars<-" <- function(x, value) {
  ax <- attributes(x)
  attributes(x) <- NULL
  ind <- which(!vapply(x, is.numeric, TRUE, USE.NAMES = FALSE))
  if(is.list(value)) {
    class(value) <- NULL
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(value) != length(ind)) stop("NCOL(value) must match length(cat_vars(x))")
    x[ind] <- value
    if(length(nam <- names(value))) ax[["names"]][ind] <- nam
  } else if(is.null(value)) {
    ax[["names"]] <- ax[["names"]][-ind]
    return(setAttributes(x[-ind], ax))
  } else {
    if(NROW(value) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(ind) != 1L) stop("NCOL(value) must match length(cat_vars(x))")
    x[[ind]] <- value
  }
  return(setAttributes(x, ax))
}

char_vars <- function(x, return = c("data","names","indices","named_indices")) {
  switch(return[1L], data = colsubsetFUN(x, is.character),
         names = names(which(vapply(x, is.character, TRUE))),
         indices = which(vapply(x, is.character, TRUE, USE.NAMES = FALSE)),
         named_indices = which(vapply(x, is.character, TRUE)),
         stop("Unknown return option!"))
}
"char_vars<-" <- function(x, value) {
  ax <- attributes(x)
  attributes(x) <- NULL
  ind <- which(vapply(x, is.character, TRUE, USE.NAMES = FALSE))
  if(is.list(value)) {
    class(value) <- NULL
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(value) != length(ind)) stop("NCOL(value) must match length(char_vars(x))")
    x[ind] <- value
    if(length(nam <- names(value))) ax[["names"]][ind] <- nam
  } else if(is.null(value)) {
    ax[["names"]] <- ax[["names"]][-ind]
    return(setAttributes(x[-ind], ax))
  } else {
    if(NROW(value) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(ind) != 1L) stop("NCOL(value) must match length(char_vars(x))")
    x[[ind]] <- value
  }
  return(setAttributes(x, ax))
}

fact_vars <- function(x, return = c("data","names","indices","named_indices")) {
  switch(return[1L], data = colsubsetFUN(x, is.factor),
         names = names(which(vapply(x, is.factor, TRUE))),
         indices = which(vapply(x, is.factor, TRUE, USE.NAMES = FALSE)),
         named_indices = which(vapply(x, is.factor, TRUE)),
         stop("Unknown return option!"))
}
"fact_vars<-" <- function(x, value) {
  ax <- attributes(x)
  attributes(x) <- NULL
  ind <- which(vapply(x, is.factor, TRUE, USE.NAMES = FALSE))
  if(is.list(value)) {
    class(value) <- NULL
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(value) != length(ind)) stop("NCOL(value) must match length(fact_vars(x))")
    x[ind] <- value
    if(length(nam <- names(value))) ax[["names"]][ind] <- nam
  } else if(is.null(value)) {
    ax[["names"]] <- ax[["names"]][-ind]
    return(setAttributes(x[-ind], ax))
  } else {
    if(NROW(value) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(ind) != 1L) stop("NCOL(value) must match length(fact_vars(x))")
    x[[ind]] <- value
  }
  return(setAttributes(x, ax))
}

logi_vars <- function(x, return = c("data","names","indices","named_indices")) {
  switch(return[1L], data = colsubsetFUN(x, is.logical),
         names = names(which(vapply(x, is.logical, TRUE))),
         indices = which(vapply(x, is.logical, TRUE, USE.NAMES = FALSE)),
         named_indices = which(vapply(x, is.logical, TRUE)),
         stop("Unknown return option!"))
}
"logi_vars<-" <- function(x, value) {
  ax <- attributes(x)
  attributes(x) <- NULL
  ind <- which(vapply(x, is.logical, TRUE, USE.NAMES = FALSE))
  if(is.list(value)) {
    class(value) <- NULL
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(value) != length(ind)) stop("NCOL(value) must match length(logi_vars(x))")
    x[ind] <- value
    if(length(nam <- names(value))) ax[["names"]][ind] <- nam
  } else if(is.null(value)) {
    ax[["names"]] <- ax[["names"]][-ind]
    return(setAttributes(x[-ind], ax))
  } else {
    if(NROW(value) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(ind) != 1L) stop("NCOL(value) must match length(logi_vars(x))")
    x[[ind]] <- value
  }
  return(setAttributes(x, ax))
}

Date_vars <- function(x, return = c("data","names","indices","named_indices")) {
  switch(return[1L], data = colsubsetFUN(x, is.Date),
         names = names(which(vapply(x, is.Date, TRUE))),
         indices = which(vapply(x, is.Date, TRUE, USE.NAMES = FALSE)),
         named_indices = which(vapply(x, is.Date, TRUE)),
         stop("Unknown return option!"))
}
"Date_vars<-" <- function(x, value) {
  ax <- attributes(x)
  attributes(x) <- NULL
  ind <- which(vapply(x, is.Date, TRUE, USE.NAMES = FALSE))
  if(is.list(value)) {
    class(value) <- NULL
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(value) != length(ind)) stop("NCOL(value) must match length(Date_vars(x))")
    x[ind] <- value
    if(length(nam <- names(value))) ax[["names"]][ind] <- nam
  } else if(is.null(value)) {
    ax[["names"]] <- ax[["names"]][-ind]
    return(setAttributes(x[-ind], ax))
  } else {
    if(NROW(value) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(ind) != 1L) stop("NCOL(value) must match length(Date_vars(x))")
    x[[ind]] <- value
  }
  return(setAttributes(x, ax))
}


get_vars <- function(x, vars, return = c("data","names","indices","named_indices"), regex = FALSE, ...) { # vars is a function or regex call. perhaps also do like this for lsubset and has.elem??, because if ... is in the second place, you cannot put other things
 if(!regex && !missing(...)) stop("Unknown argument ", dotstostr(...))
 switch(return[1L],
        data = colsubset(x, if(is.character(vars) && regex) rgrep(vars, attr(x, "names"), ...) else vars),
        names = if(is.function(vars)) names(which(vapply(x, vars, TRUE))) else if(is.character(vars) && regex)
                                     attr(x, "names")[rgrep(vars, attr(x, "names"), ...)] else attr(x, "names")[vars],
        indices = if(is.function(vars)) which(vapply(x, vars, TRUE, USE.NAMES = FALSE)) else if(is.character(vars) && regex)
                   rgrep(vars, attr(x, "names"), ...) else anyNAerror(match(vars, attr(x, "names")), "Unknown column names!"), # put error !!
        named_indices = if(is.function(vars)) which(vapply(x, vars, TRUE)) else if(is.character(vars) && regex) {
                       nam <- attr(x, "names")
                       ind <- rgrep(vars, nam, ...)
                      `names<-`(ind, nam[ind])} else stop("vars must be a function or a regular expression"),
        stop("Unknown return option!"))
}
gv <- get_vars


"get_vars<-" <- function(x, vars, regex = FALSE, ..., value) {
  ax <- attributes(x)
  attributes(x) <- NULL
  if(is.numeric(vars)) {
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    if(max(abs(vars)) > length(x)) stop("Column index out of range abs(1:length(x))")
    if(any(vars < 0)) vars <- seq_along(x)[vars]
  } else if(is.logical(vars)) {
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    if(length(vars) != length(x)) stop("Logical subsetting vector must match length(x)")
    vars <- which(vars)
  } else {
    if(!regex && !missing(...)) stop("Unknown argument ", dotstostr(...))
    vars <- if(is.function(vars)) which(vapply(x, vars, TRUE, USE.NAMES = FALSE)) else if(is.character(vars) && regex)
      rgrep(vars, ax[["names"]], ...) else anyNAerror(match(vars, ax[["names"]]), "Unknown column names!")
  }
  if(is.list(value)) {
    class(value) <- NULL # fastest ??
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(value) != length(vars)) stop("NCOL(value) must match length(vars)")
    x[vars] <- value
    if(length(nam <- names(value))) ax[["names"]][vars] <- nam #  == length(vars)
  } else if(is.null(value)) {
    ax[["names"]] <- ax[["names"]][-vars] # fastest ??? -> Yes !! This is slower: x[vars] <- NULL
    return(setAttributes(x[-vars], ax))
  } else {
    if(NROW(value) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    if(length(vars) != 1L) stop("NCOL(value) must match length(vars)")
    x[[vars]] <- value
  }
  return(setAttributes(x, ax))
}

"gv<-" <- `get_vars<-`

# Previous Version:
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
#            rgrep(vars, ax[["names"]], ...) else anyNAerror(match(vars, ax[["names"]]), "Unknown column names!")
#   }
#   if(NCOL2(d, ilv) != length(vars)) stop("NCOL(value) must match length(vars)")
#   if(ilv) x[vars] <- value else x[[vars]] <- value
#   if(ilv && length(nam <- attr(value, "names")) == length(vars)) ax[["names"]][vars] <- nam
#   return(setAttributes(x, ax))
# }


"add_vars<-" <- function(x, value) {
  ax <- attributes(x)
  attributes(x) <- NULL
  lx <- length(x)
  if(is.list(value)) {
    class(value) <- NULL # fastest ??
    if(length(value[[1L]]) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    # res <- c(x, value)  # FASTER than commented out below !!!
    ax[["names"]] <- if(length(nam <- names(value)))  c(ax[["names"]], nam) else
      c(ax[["names"]], paste0("V", seq(lx+1L, lx+length(value))))
    return(setAttributes(c(x, value), ax))
    # ind <- seq(lx+1L, lx+length(value))
    # x[ind] <- value  # FASTER than simply using x[names(value)] <- value ??????????? -> Yes !!!
    # ax[["names"]] <- if(length(nam <- names(value)))  c(ax[["names"]], nam) else
    #   c(ax[["names"]], paste0("V", ind))
  } else {
    if(NROW(value) != length(x[[1L]])) stop("NROW(value) must match nrow(x)")
    # res <- c(x, list(value)) # FASTER than below ??? -> Nope !!!!
    # ax[["names"]] <- c(ax[["names"]], paste0("V", lx+1L))
    x[[lx+1L]] <- value
    ax[["names"]] <- c(ax[["names"]], paste0("V", lx+1L))
    return(setAttributes(x, ax))
  }
  # return(setAttributes(res, ax))
  # return(setAttributes(x, ax))
}

"av<-" <- `add_vars<-`

# Previous Version:
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







# Note: This does not yet work properly with data.tables !!! (perhaps can also use unclass and reclass???)

# "get.vars<-" # still do!! # already exists in some packages, perhaps select.vars. nah, thats too much dpler??



# Neat Examples:
#unlist2d(rapply2d(IRF,dim))
#list.elem(IRF) %>% rapply2d(setRownames) %>% unlist2d(c("type","shock"),TRUE,"time") %>% head
#unlist2d(rapply2d(list.elem(IRF),setRownames),c("type","shock"),TRUE,"time")


# Exercizes:
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
