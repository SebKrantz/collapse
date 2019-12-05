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
  switch(return[1L], data = colsubset(x, is.numeric), # switch(match.arg(return)
         names = names(which(vapply(x, is.numeric, TRUE))),
         indices = which(vapply(x, is.numeric, TRUE, USE.NAMES = FALSE)),
         named_indices = which(vapply(x, is.numeric, TRUE)),
         stop("Unknown return option!"))
}
"num_vars<-" <- function(x, value) {
    ax <- attributes(x)
    ilv <- is.list(value)
    d <- dim(value)
    if(NROW2(value, d) != nrow(x)) stop("NROW(value) must match nrow(x)")
    attributes(x) <- NULL # vapply without attributes is faster !!!. if Fail, data unchanged !!
    ind <- which(vapply(x, is.numeric, TRUE, USE.NAMES = FALSE))
    if(NCOL2(d, ilv) != length(ind)) stop("NCOL(value) must match length(num.vars(x))")
    if(ilv) x[ind] <- value else x[[ind]] <- value #   x[ind] <- if(ilv) value else list(value)
    if(ilv && length(nam <- names(value)) == length(ind)) # does stop execution of if statement if ilv is false !!!!!
       ax[["names"]][ind] <- nam
    return(setAttributes(x, ax))
}

cat_vars <- function(x, return = c("data","names","indices","named_indices")) {
  switch(return[1L], data = colsubset(x, !vapply(x, is.numeric, TRUE)), # a lot faster than using is.categorical
         names = names(which(!vapply(x, is.numeric, TRUE))),
         indices = which(!vapply(x, is.numeric, TRUE, USE.NAMES = FALSE)),
         named_indices = which(!vapply(x, is.numeric, TRUE)),
         stop("Unknown return option!"))
}
"cat_vars<-" <- function(x, value) {
  ax <- attributes(x)
  ilv <- is.list(value)
  d <- dim(value)
  if(NROW2(value, d) != nrow(x)) stop("NROW(value) must match nrow(x)")
  attributes(x) <- NULL
  ind <- which(!vapply(x, is.numeric, TRUE, USE.NAMES = FALSE))
  if(NCOL2(d, ilv) != length(ind)) stop("NCOL(value) must match length(cat.vars(x))")
  if(ilv) x[ind] <- value else x[[ind]] <- value
  if(ilv && length(nam <- names(value)) == length(ind)) ax[["names"]][ind] <- nam
  return(setAttributes(x, ax))
}

char_vars <- function(x, return = c("data","names","indices","named_indices")) {
  switch(return[1L], data = colsubset(x, is.character),
         names = names(which(vapply(x, is.character, TRUE))),
         indices = which(vapply(x, is.character, TRUE, USE.NAMES = FALSE)),
         named_indices = which(vapply(x, is.character, TRUE)),
         stop("Unknown return option!"))
}
"char_vars<-" <- function(x, value) {
  ax <- attributes(x)
  ilv <- is.list(value)
  d <- dim(value)
  if(NROW2(value, d) != nrow(x)) stop("NROW(value) must match nrow(x)")
  attributes(x) <- NULL
  ind <- which(vapply(x, is.character, TRUE, USE.NAMES = FALSE))
  if(NCOL2(d, ilv) != length(ind)) stop("NCOL(value) must match length(char.vars(x))")
  if(ilv) x[ind] <- value else x[[ind]] <- value
  if(ilv && length(nam <- names(value)) == length(ind)) ax[["names"]][ind] <- nam
  return(setAttributes(x, ax))
}

fact_vars <- function(x, return = c("data","names","indices","named_indices")) {
  switch(return[1L], data = colsubset(x, is.factor),
         names = names(which(vapply(x, is.factor, TRUE))),
         indices = which(vapply(x, is.factor, TRUE, USE.NAMES = FALSE)),
         named_indices = which(vapply(x, is.factor, TRUE)),
         stop("Unknown return option!"))
}
"fact_vars<-" <- function(x, value) {
  ax <- attributes(x)
  ilv <- is.list(value)
  d <- dim(value)
  if(NROW2(value, d) != nrow(x)) stop("NROW(value) must match nrow(x)")
  attributes(x) <- NULL
  ind <- which(vapply(x, is.factor, TRUE, USE.NAMES = FALSE))
  if(NCOL2(d, ilv) != length(ind)) stop("NCOL(value) must match length(fact.vars(x))")
  if(ilv) x[ind] <- value else x[[ind]] <- value
  if(ilv && length(nam <- names(value)) == length(ind)) ax[["names"]][ind] <- nam
  return(setAttributes(x, ax))
}

logi_vars <- function(x, return = c("data","names","indices","named_indices")) {
  switch(return[1L], data = colsubset(x, is.logical),
         names = names(which(vapply(x, is.logical, TRUE))),
         indices = which(vapply(x, is.logical, TRUE, USE.NAMES = FALSE)),
         named_indices = which(vapply(x, is.logical, TRUE)),
         stop("Unknown return option!"))
}
"logi_vars<-" <- function(x, value) {
  ax <- attributes(x)
  ilv <- is.list(value)
  d <- dim(value)
  if(NROW2(value, d) != nrow(x)) stop("NROW(value) must match nrow(x)")
  attributes(x) <- NULL
  ind <- which(vapply(x, is.logical, TRUE, USE.NAMES = FALSE))
  if(NCOL2(d, ilv) != length(ind)) stop("NCOL(value) must match length(logi.vars(x))")
  if(ilv) x[ind] <- value else x[[ind]] <- value
  if(ilv && length(nam <- names(value)) == length(ind)) ax[["names"]][ind] <- nam
  return(setAttributes(x, ax))
}

Date_vars <- function(x, return = c("data","names","indices","named_indices")) {
  switch(return[1L], data = colsubset(x, is.Date),
         names = names(which(vapply(x, is.Date, TRUE))),
         indices = which(vapply(x, is.Date, TRUE, USE.NAMES = FALSE)),
         named_indices = which(vapply(x, is.Date, TRUE)),
         stop("Unknown return option!"))
}
"Date_vars<-" <- function(x, value) {
  ax <- attributes(x)
  ilv <- is.list(value)
  d <- dim(value)
  if(NROW2(value, d) != nrow(x)) stop("NROW(value) must match nrow(x)")
  attributes(x) <- NULL
  ind <- which(vapply(x, is.Date, TRUE, USE.NAMES = FALSE))
  if(NCOL2(d, ilv) != length(ind)) stop("NCOL(value) must match length(Date.vars(x))")
  if(ilv) x[ind] <- value else x[[ind]] <- value
  if(ilv && length(nam <- names(value)) == length(ind)) ax[["names"]][ind] <- nam
  return(setAttributes(x, ax))
}


get_vars <- function(x, vars, return = c("data","names","indices","named_indices"), regex = FALSE, ...) { # vars is a function or regex call. perhaps also do like this for lsubset and has.elem??, because if ... is in the second place, you cannot put other things
 if(!regex && !missing(...)) stop("Unknown argument ", dotstostr(...))
 switch(return[1L],
        data = colsubset(x, if(is.character(vars) && regex) rgrep(vars, names(x), ...) else vars),
        names = if(is.function(vars)) names(which(vapply(x, vars, TRUE))) else if(is.character(vars) && regex)
                                     names(x)[rgrep(vars, names(x), ...)] else names(x)[vars],
        indices = if(is.function(vars)) which(vapply(x, vars, TRUE, USE.NAMES = FALSE)) else if(is.character(vars) && regex)
                   rgrep(vars, names(x), ...) else anyNAerror(match(vars, names(x)), "Unknown column names!"), # put error !!
        named_indices = if(is.function(vars)) which(vapply(x, vars, TRUE)) else if(is.character(vars) && regex) {
                       nam <- names(x)
                       ind <- rgrep(vars, nam, ...)
                      `names<-`(ind, nam[ind])} else stop("vars must be a function or a regular expression"),
        stop("Unknown return option!"))
}
"get_vars<-" <- function(x, vars, regex = FALSE, ..., value) {
  ax <- attributes(x)
  ilv <- is.list(value)
  d <- dim(value)
  if(NROW2(value, d) != nrow(x)) stop("NROW(value) must match nrow(x)")
  attributes(x) <- NULL
  if(is.numeric(vars)) {
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    if(max(abs(vars)) > length(x)) stop("Column index out of range abs(1:length(x))")
  } else if(is.logical(vars)) {
    if(!missing(...)) stop("Unknown argument ", dotstostr(...))
    if(length(vars) != length(x)) stop("Logical subsetting vector must match length(x)")
  } else {
    if(!regex && !missing(...)) stop("Unknown argument ", dotstostr(...))
    vars <- if(is.function(vars)) which(vapply(x, vars, TRUE, USE.NAMES = FALSE)) else if(is.character(vars) && regex)
           rgrep(vars, ax[["names"]], ...) else anyNAerror(match(vars, ax[["names"]]), "Unknown column names!")
  }
  if(NCOL2(d, ilv) != length(vars)) stop("NCOL(value) must match length(num.vars(x))")
  if(ilv) x[vars] <- value else x[[vars]] <- value
  if(ilv && length(nam <- names(value)) == length(vars)) ax[["names"]][vars] <- nam
  return(setAttributes(x, ax))
}


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
