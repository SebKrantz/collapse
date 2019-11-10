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

get_vars <- function(x, FoR, return = c("data","names","indices","named_indices"), regex = FALSE, ...) { # FoR is a function or regex call. perhaps also do like this for lsubset and has.elem??, because if ... is in the second place, you cannot put other things
     switch(return[1L], 
            data = colsubset(x, if(is.character(FoR) && regex) rgrep(FoR, names(x), ...) else FoR), 
            names = if(is.function(FoR)) names(which(vapply(x, FoR, TRUE))) else if(is.character(FoR) && regex) 
                                         names(x)[rgrep(FoR, names(x), ...)] else names(x)[FoR],
            indices = if(is.function(FoR)) which(vapply(x, FoR, TRUE, USE.NAMES = FALSE)) else if(is.character(FoR) && regex)
                       rgrep(FoR, names(x), ...) else anyNAerror(match(FoR, names(x)), "Unknown column names!"), # put error !!
            named_indices = if(is.function(FoR)) which(vapply(x, FoR, TRUE)) else if(is.character(FoR) && regex) {
                           nam <- names(x)
                           ind <- rgrep(FoR, nam, ...)
                          `names<-`(ind, nam[ind])} else stop("FoR must be a function or a regular expression"),
            stop("Unknown return option!"))
}
"get_vars<-" <- function(x, FoR, regex = FALSE, ..., value) {
  ax <- attributes(x)
  ilv <- is.list(value)
  d <- dim(value)
  if(NROW2(value, d) != nrow(x)) stop("NROW(value) must match nrow(x)")
  attributes(x) <- NULL
  # if(!(is.numeric(FoR) || is.logical(FoR))) # old version without checks !!!
  #   FoR <- if(is.function(FoR)) which(vapply(x, FoR, TRUE, USE.NAMES = FALSE)) else if(is.character(FoR) && regex)
  #          rgrep(FoR, ax[["names"]], ...) else match(FoR, ax[["names"]])
  if(is.numeric(FoR)) { 
    if(max(abs(FoR)) > length(x)) stop("Column index out of range abs(1:length(x))") 
  } else if(is.logical(FoR)) { 
    if(length(FoR) != length(x)) stop("Logical subsetting vector must match length(x)") 
  } else {
    FoR <- if(is.function(FoR)) which(vapply(x, FoR, TRUE, USE.NAMES = FALSE)) else if(is.character(FoR) && regex)
           rgrep(FoR, ax[["names"]], ...) else anyNAerror(match(FoR, ax[["names"]]), "Unknown column names!")
  }
  if(NCOL2(d, ilv) != length(FoR)) stop("NCOL(value) must match length(num.vars(x))")
  if(ilv) x[FoR] <- value else x[[FoR]] <- value
  if(ilv && length(nam <- names(value)) == length(FoR)) ax[["names"]][FoR] <- nam 
  return(setAttributes(x, ax))
}


# Note: This does not yet work properly with data.tables !!! (perhaps can also use unclass and reclass???)

# "get.vars<-" # still do!! # already exists in some packages, perhaps select.vars. nah, thats too much dpler??


# Previous Versions:
# num.vars <- function(X, return = 0L) if (is.list(X)) {if (return == 0L) base::Filter(is.numeric,X) else if (return == 1L) names(which(unlist(lapply(X,is.numeric)))) else if (return == 2L) which(unlist(lapply(X,is.numeric), use.names = FALSE)) else which(unlist(lapply(X,is.numeric))) } else if (is.numeric(X)) X # X[unlist(lapply(X,is.numeric), use.names = FALSE, recursive = FALSE)]
# "num.vars<-" <- function(X, value) {
#   if (is.list(X)) { 
#     if (inherits(X, "data.table")) # make a function myfilter that does this for data.table??
#       X[,which(unlist(lapply(X,is.numeric), use.names = FALSE, recursive = FALSE))] = value else
#         X[unlist(lapply(X,is.numeric), use.names = FALSE, recursive = FALSE)] = value
#   } else if (is.numeric(X)) X = value
#   X
# }
# cat.vars <- function(X, return = 0L) {
#   is.categorical <- function(x) !is.numeric(x)
#   if (is.list(X)) {if (return == 0L) base::Filter(is.categorical,X) else if (return == 1L) names(which(unlist(lapply(X,is.categorical)))) else if (return == 2L)  which(unlist(lapply(X,is.categorical), use.names = FALSE)) else which(unlist(lapply(X,is.categorical))) } else if (is.categorical(X)) X 
# }
# "cat.vars<-" <- function(X, value) {
#   if (is.list(X)) { 
#     if (inherits(X, "data.table")) 
#       X[,which(!unlist(lapply(X,is.numeric), use.names = FALSE, recursive = FALSE))] = value else
#         X[!unlist(lapply(X,is.numeric), use.names = FALSE, recursive = FALSE)] = value
#   } else if (!is.numeric(X)) X = value
#   X
# }
# char.vars <- function(X, return = 0L) if (is.list(X)) {if (return == 0L) base::Filter(is.character,X) else if (return == 1L) names(which(unlist(lapply(X,is.character)))) else if (return == 2L) which(unlist(lapply(X,is.character), use.names = FALSE)) else which(unlist(lapply(X,is.character))) } else if (is.character(X)) X
# "char.vars<-" <- function(X, value) {
#   if (is.list(X)) { 
#     if (inherits(X, "data.table")) 
#       X[,which(unlist(lapply(X,is.character), use.names = FALSE, recursive = FALSE))] = value else
#         X[unlist(lapply(X,is.character), use.names = FALSE, recursive = FALSE)] = value
#   } else if (is.character(X)) X = value
#   X
# }
# fact.vars <- function(X, return = 0L) if (is.list(X)) {if (return == 0L) base::Filter(is.factor,X) else if (return == 1L) names(which(unlist(lapply(X,is.factor)))) else if (return == 2L) which(unlist(lapply(X,is.factor), use.names = FALSE)) else which(unlist(lapply(X,is.factor)))} else if (is.factor(X)) X
# "fact.vars<-" <- function(X, value) {
#   if (is.list(X)) { 
#     if (inherits(X, "data.table")) 
#       X[,which(unlist(lapply(X,is.factor), use.names = FALSE, recursive = FALSE))] = value else
#         X[unlist(lapply(X,is.factor), use.names = FALSE, recursive = FALSE)] = value
#   } else if (is.factor(X)) X = value
#   X
# }
# logi.vars <- function(X, return = 0L) if (is.list(X)) {if (return == 0L) base::Filter(is.logical,X) else if (return == 1L) names(which(unlist(lapply(X,is.logical)))) else if (return == 2L)  which(unlist(lapply(X,is.logical), use.names = FALSE)) else which(unlist(lapply(X,is.logical))) } else if (is.logical(X)) X
# "logi.vars<-" <- function(X, value) {
#   if (is.list(X)) { 
#     if (inherits(X, "data.table")) 
#       X[,which(unlist(lapply(X,is.logical), use.names = FALSE, recursive = FALSE))] = value else
#         X[unlist(lapply(X,is.logical), use.names = FALSE, recursive = FALSE)] = value
#   } else if (is.logical(X)) X = value
#   X
# }
# get.vars <- function(X, FoR, return = 0L, regex = TRUE, ...) { # FoR is a function or regex call. perhaps also do like this for lsubset and has.elem??, because if ... is in the second place, you cannot put other things
#   if (is.function(FoR)) { # what about ... to functions?? Filter doesn't allow
#     if (is.list(X)) {if (return == 0L) base::Filter(FoR,X) else if (return == 1L) names(which(unlist(lapply(X,FoR)))) else if (return == 2L)  which(unlist(lapply(X,FoR), use.names = FALSE)) else which(unlist(lapply(X,FoR))) } else if (FoR(X)) X
#   } else if (is.character(FoR)) {
#     if (regex) { # could do pmatch or charmatch, but less options!!
#       rgrep <- function(exp, nam, ...) if (length(exp)>1) sort(unique(unlist(lapply(exp,grep,nam, ...), use.names = FALSE))) else grep(exp,nam,...)
#       if (is.list(X)) { 
#         if (return == 0L) X[rgrep(FoR,names(X),...)] else if (return == 1L) names(X)[rgrep(FoR,names(X),...)] else if (return == 2L)  rgrep(FoR,names(X),...) else {ind <- rgrep(FoR,names(X),...); setNames(ind,names(X)[ind]) } 
#       } else if (is.matrix(X)) {
#         if (return == 0L) X[,rgrep(FoR,colnames(X),...)] else if (return == 1L) colnames(X)[rgrep(FoR,colnames(X),...)] else if (return == 2L) rgrep(FoR,colnames(X),...) else {ind <- rgrep(FoR,colnames(X),...); setNames(ind,colnames(X)[ind]) }    
#       }
#     } else {
#       if (is.list(X)) { 
#         if (return == 0L) X[FoR] else if (return == 1L) FOR else if (return == 2L) match(FoR,names(X)) else setNames(match(FoR,names(X)),FoR)  
#       } else if (is.matrix(X)) {
#         if (return == 0L) X[,FoR, drop = FALSE] else if (return == 1L) FoR else if (return == 2L) match(FoR,colnames(X)) else setNames(match(FoR,colnames(X)),FoR)    
#       }
#     }
#   } else if (is.numeric(FoR)) {
#     if (is.list(X)) { 
#       if (return == 0L) X[FoR] else if (return == 1L) names(X)[FoR] else if (return == 2L)  FoR else setNames(FoR,names(X)[FoR])  
#     } else if (is.matrix(X)) {
#       if (return == 0L) X[,FoR, drop = FALSE] else if (return == 1L) colnames(X)[FoR] else if (return == 2L) FoR else setNames(FoR,colnames(X)[FoR])    
#     }
#   } else stop("FoR needs to be a function, character vector, numeric vector or regular expression (character)")
# }


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