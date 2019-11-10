# or getGRP: more speed using tabulate ?? or code your own??
GRP <- function(X, ...) UseMethod("GRP", X)

# To do: use data.table:::CsubsetDT and CsubsetVector !!
GRP.default <- function(X, by = NULL, sort = TRUE, order = 1L, 
                na.last = FALSE, return.groups = TRUE, return.order = FALSE, ...) { # , gs = TRUE # o
  if(is.list(X) && is.null(by)) {
    by <- seq_along(X)
    namby <- names(X)
  } else if(is.call(by)) 
    namby <- by <- all.vars(by) else if(is.null(by)) 
    namby <- deparse(substitute(X)) else namby <- names(X)[by] # is.call is fastest!!!  attr(terms.formula(by), "term.labels")
  o = data.table:::forderv(X, by, TRUE, sort, order, na.last) # faster calling c??
  f = attr(o, "starts") # , exact = TRUE -> slightly slower 
  if (length(o)) { # if ordered, returns 0 !!! perhaps rewrite Cforder function??
    len = .Call(data.table:::Cuniqlengths, f, length(o)) # len = data.table:::uniqlengths(f, length(o))
    # grpuo = rep.int(seq_along(len),len)[data.table:::forderv(o)] # faster way?? -> try to speed this up.. 
    grpuo = .Call(data.table:::Cfrank, o, f, len, "dense") #-1L
    res <- if (return.groups) {
      # if (gs) 
      if (return.order) { 
                   list(N.groups = length(f), group.id = grpuo, group.sizes = len, # speed up this part??
                        groups = if(inherits(X, "data.table")) X[o[f], by, with = FALSE] else if (is.atomic(X)) X[o[f]] else if(inherits(X, "data.frame"))
                        X[o[f], by, drop = FALSE] else qDF(X)[o[f], by, drop = FALSE], ordered = c(GRP.sort = sort, initially.ordered = FALSE), group.vars = namby, order = o, call = match.call()) 
        } else
                        list(N.groups = length(f), group.id = grpuo, group.sizes = len, # speed up this part??
                             groups = if (inherits(X, "data.table")) X[o[f], by, with = FALSE] else if (is.atomic(X)) X[o[f]] else if(inherits(X, "data.frame"))
                             X[o[f], by, drop = FALSE] else qDF(X)[o[f], by, drop = FALSE], ordered = c(GRP.sort = sort, initially.ordered = FALSE), group.vars = namby, call = match.call())
      # else 
                       # list(ng = length(f), g = grpuo,  # speed up this part??
                       #      groups = if (inherits(X, "data.table")) X[o[f], by, with = FALSE] else if (is.atomic(X)) X[o[f]] else
                       #        X[o[f], by, drop = FALSE], ordered = order)
    } else {
      # if (gs) 
      if (return.order) list(N.groups = length(f), group.id = grpuo, group.sizes = len, 
                             ordered = c(GRP.sort = sort, initially.ordered = FALSE), group.vars = namby, order = o, call = match.call()) else # else list(ng = length(f), g = grpuo, ordered = order)
                        list(N.groups = length(f), group.id = grpuo, group.sizes = len, 
                             ordered = c(GRP.sort = sort, initially.ordered = FALSE), group.vars = namby, call = match.call())
    }
  } else { 
    len = .Call(data.table:::Cuniqlengths, f, NROW(X)) # or cumsubtract !!!
    grpuo = rep.int(seq_along(len),len)
    # o = f # for unique. right???
    res <- if (return.groups) {
     # if (gs) 
      if (return.order) {
        # o2 = .Call(data.table:::Cfrank, seq_len(NROW(X)), f, len, "sequence") # good ???? This does not really give order, but unique row id's !!!
        # attributes(o2) = attributes(o) # good ????
        list(N.groups = length(f), group.id = grpuo, group.sizes = len,  # speed up this part??
                             groups = if (inherits(X, "data.table")) X[f, by, with = FALSE] else if (is.atomic(X)) X[f] else if(inherits(X, "data.frame"))
                             X[f, by, drop = FALSE] else qDF(X)[f, by, drop = FALSE], ordered = c(GRP.sort = sort, initially.ordered = TRUE), group.vars = namby, order = o, call = match.call()) 
        } else 
                        list(N.groups = length(f), group.id = grpuo, group.sizes = len,  # speed up this part??
                             groups = if (inherits(X, "data.table")) X[f, by, with = FALSE] else if (is.atomic(X)) X[f] else if(inherits(X, "data.frame"))
                             X[f, by, drop = FALSE] else qDF(X)[f, by, drop = FALSE], ordered = c(GRP.sort = sort, initially.ordered = TRUE), group.vars = namby, call = match.call())
      
      # else 
                       # list(ng = length(f), g = grpuo,  # speed up this part??
                       #      groups = if (inherits(X, "data.table")) X[f, by, with = FALSE] else if (is.atomic(X)) X[f] else
                       #        X[f, by, drop = FALSE], ordered = order)
    } else {
      #if (gs)
      if (return.order) {
        attributes(o) <- NULL
        list(N.groups = length(f), group.id = grpuo, group.sizes = len, 
                             ordered = c(GRP.sort = sort, initially.ordered = TRUE), group.vars = namby, order = o, call = match.call()) 
        } else # else list(ng = length(f), g = grpuo, ordered = order)
                               list(N.groups = length(f), group.id = grpuo, group.sizes = len, 
                                    ordered = c(GRP.sort = sort, initially.ordered = TRUE), group.vars = namby, call = match.call())
    }
  }
  class(res) <- "GRP"
  res
  # or only group: frankv(mtcars, cols = c("cyl","vs","am"), ties.method = "dense")
}

is.GRP <- function(x) inherits(x, "GRP")

group.names.GRP <- function(g) { # fastest !!!
  groups <- g[["groups"]]
  if(is.atomic(groups)) groups else do.call(paste,c(groups, list(sep = "."))) 
}

print.GRP <- function(g, n = 6) {
  # fs <- function(x) {
  #   x[1] = paste0("   ",x[1])
  #   x
  # }
  ord <- g[["ordered"]]
  cat(paste("collapse grouping object of length",length(g[[2]]),"with",g[[1]],ifelse(any(ord),"ordered","unordered"),"groups"), fill = TRUE)
  #if(!is.null(g[["gs"]])) {
  cat("\nCall: ", paste0(deparse(g[["call"]]),", ",ifelse(ord[2],"ordered","unordered")), "\n\n", sep = "")
    cat("Distribution of group sizes: ", fill = TRUE) 
    print.summaryDefault(summary.default(g[["group.sizes"]]))
  #}
  if(!is.null(g[["groups"]])) {
    ug <- g[["groups"]]
    cat("\nGroups with sizes: ", fill = TRUE) 
    if(is.atomic(ug)) {
      if(length(ug)>2*n) {
        ind <- seq.int(g[[1]]-n+1,g[[1]])
        print.default(setNames(g[["group.sizes"]][1:n],ug[1:n])) # cat(, sep = ", ")
        cat("  ---", fill = TRUE) 
        print.default(setNames(g[["group.sizes"]][ind],ug[ind]))
      } else {
        print.default(setNames(g[["group.sizes"]],ug))
      } 
    } else {
      if(length(ug[[1]])>2*n) {
        ind <- seq.int(g[[1]]-n+1,g[[1]])
        print.default(setNames(g[["group.sizes"]][1:n], do.call(paste, c(lapply(ug,function(x)x[1:n]), list(sep = ".")))))
        cat("  ---", fill = TRUE) 
        print.default(setNames(g[["group.sizes"]][ind], do.call(paste, c(lapply(ug,function(x)x[ind]), list(sep = ".")))))
      } else print.default(setNames(g[["group.sizes"]], do.call(paste, c(ug, list(sep = ".")))))
    }
  }
}

plot.GRP <- function(g, breaks = "auto", type = "s", horizontal = FALSE) {
  settings <- par(c("mfrow","mar","mgp"))
  par(mfrow = if(horizontal) 1:2 else 2:1, mar = c(3.9,4.1,2.1,1), mgp = c(2.5,1,0))
  if(breaks == "auto") {
    ugs <- length(unique.default(g[[3]]))
    breaks <- if(ugs > 80) 80 else ugs
  }
  #if(breaks == "auto") {
    plot(seq_len(g[[1]]), g[[3]], type = type, xlab = "Group id", ylab = "Group Size", main = paste0("Sizes of ",g[[1]]," ",ifelse(any(g[["ordered"]]),"Ordered","Unordered")," Groups"), frame.plot = FALSE)
    if(breaks == 1) plot(g[[3]][1], g[[1]], type = "h", ylab = "Frequency", xlab = "Group Size", main = "Histogram of Group Sizes", frame.plot = FALSE) else 
    hist(g[[3]], breaks, xlab = "Group Size", main = "Histogram of Group Sizes")
  #} else {
    #if(length(breaks) != 2) stop("'breaks' must supply a vector of list of the breaks settings for the two histograms")
   # plot(seq_len(g[[1]]), g[[3]], type = type, xlab = "Group-id", main = paste0("Histogram of ",g[[1]]," ",ifelse(g[["ordered"]][1],"Ordered","Unordered")," Groups"))
  #  hist(g[[3]], breaks[[1]], xlab = "Group Size", main = paste0("Histogram of Group-Sizes"))
  #}
  par(settings)
}

as.factor.GRP <- function(g) { 
  # if(is.factor(g)) return(g)
  # if(class(g2) != "GRP") stop("g must be a 'GRP' object")
  f <- g[[2]]
  gr <- g[["groups"]] # faster !!!
  if(is.null(gr)) # .set_row_names(g[[1]]) # as.character ?? faster ?? // formatC ?? 
    attr(f, "levels") <- as.character(seq_len(g[[1]])) else {
    attr(f, "levels") <- if(is.atomic(gr)) as.character(gr) else # or formatC ??
      do.call(paste, c(gr, list(sep = ".")))
  }
  class(f) <- if(any(g[["ordered"]])) c("ordered","factor") else "factor"
  f
}

GRP.qG <- function(X) { # need to class !!
  ng <- attr(X, "N.groups")
  gs <- .Internal(tabulate(X, ng))
  ordered <- if(is.ordered(X)) c(TRUE,TRUE) else c(FALSE,FALSE)
  attributes(X) <- NULL
  `class<-`(list(N.groups = nl, group.id = X, group.sizes = gs, 
                 ordered = ordered, call = match.call()), "GRP")
}

GRP.factor <- function(X) { 
  lev <- attr(X, "levels")
  nl <- length(lev)
  gs <- .Internal(tabulate(X, nl))
  ordered <- if(is.ordered(X)) c(TRUE,TRUE) else c(FALSE,FALSE)
  attributes(X) <- NULL
  `class<-`(list(N.groups = nl, group.id = X, 
                 group.sizes = gs, groups = lev, 
                 ordered = ordered, call = match.call()), "GRP")
}

GRP.pseries <- function(X) { 
  g <- attr(X, "index")[[1]]
  lev <- attr(g, "levels")
  nl <- length(lev)
  gs <- .Internal(tabulate(g, nl))
  ordered <- if(is.ordered(g)) c(TRUE,TRUE) else c(FALSE,FALSE)
  attributes(g) <- NULL
  `class<-`(list(N.groups = nl, group.id = g, 
                 group.sizes = gs, groups = lev, 
                 ordered = ordered, call = match.call()), "GRP")
}

GRP.pdata.frame <- function(X) { 
  g <- attr(X, "index")[[1]]
  lev <- attr(g, "levels")
  nl <- length(lev)
  gs <- .Internal(tabulate(g, nl))
  ordered <- if(is.ordered(g)) c(TRUE,TRUE) else c(FALSE,FALSE)
  attributes(g) <- NULL
  `class<-`(list(N.groups = nl, group.id = g, 
                 group.sizes = gs, groups = lev, 
                 ordered = ordered, call = match.call()), "GRP")
}

GRP.grouped_df <- function(X) { 
  g <- unclass(attr(X, "groups"))
  lg <- length(g)
  gr <- g[[lg]]
  ng <- length(gr)
  gs <- lengths(gr) # faster sorting still ?? qsort ?? or data.table forder -> nope, slower than order !!
  `class<-`(list(N.groups = ng, group.id = rep(seq_len(ng), gs)[.Internal(radixsort(TRUE, FALSE, FALSE, TRUE, .Internal(unlist(gr, FALSE, FALSE))))], #[order(unlist(gr, use.names = FALSE, recursive = FALSE))], 
                 group.sizes = gs, groups = g[-lg], ordered = c(TRUE,TRUE), call = match.call()), "GRP")
}




# More parsimonious but slower !!!!
# GRP2 <- function(X, by = if(is.atomic(X)) NULL else seq_along(X), sort = TRUE, order = 1L, 
#                  na.last = FALSE, return.groups = TRUE, return.order = FALSE, ...) { 
#   if(is.call(by)) by = all.vars(by) 
#   o = data.table:::forderv(X, by, TRUE, sort, order, na.last) 
#   f = attr(o, "starts") 
#   mg = attr(o, "maxgrpn")
#   if (length(o)) { 
#     len = .Call(data.table:::Cuniqlengths, f, length(o)) 
#     grpuo = .Call(data.table:::Cfrank, o, f, len, "dense") 
#     ordered = c(GRP.sort = sort, initially.ordered = FALSE)
#     if (return.groups) {
#       groups = if (inherits(X, "data.table")) X[o[f], by, with = FALSE] else if (is.atomic(X)) X[o[f]] else X[o[f], by, drop = FALSE]
#       if (return.order) { 
#         attributes(o) <- NULL
#         res <- list(length(f), grpuo, len, f, mg, groups, ordered, o, match.call()) 
#         attr(res,"names") = c("N.groups","group.id","group.sizes","starts","maxgrpn","groups","ordered","order","call")
#       } else {
#         res <- list(length(f), grpuo, len, f, mg, groups, ordered, match.call()) 
#         attr(res,"names") = c("N.groups","group.id","group.sizes","starts","maxgrpn","groups","ordered","call")
#       }
#     } else {
#       if (return.order) {
#         attributes(o) <- NULL
#         res <- list(length(f), grpuo, len, f, mg, ordered, o, match.call()) 
#         attr(res,"names") = c("N.groups","group.id","group.sizes","starts","maxgrpn","ordered","order","call")
#       } else {
#         res <- list(length(f), grpuo, len, f, mg, ordered, match.call()) 
#         attr(res,"names") = c("N.groups","group.id","group.sizes","starts","maxgrpn","ordered","call")
#       }
#     }
#   } else { 
#     len = .Call(data.table:::Cuniqlengths, f, NROW(X)) # or cumsubtract !!
#     grpuo = rep.int(seq_along(len),len)
#     ordered = c(GRP.sort = sort, initially.ordered = TRUE)
#     if (return.groups) {
#       groups = if (inherits(X, "data.table")) X[f, by, with = FALSE] else if (is.atomic(X)) X[f] else X[f, by, drop = FALSE]
#       if (return.order) {
#         attributes(o) <- NULL
#         res <- list(length(f), grpuo, len, f, mg, groups, ordered, o, call = match.call()) 
#         attr(res,"names") = c("N.groups","group.id","group.sizes","starts","maxgrpn","groups","ordered","order","call")
#       } else 
#         res <- list(length(f), grpuo, len, f, mg, groups, ordered, call = match.call())       
#       attr(res,"names") = c("N.groups","group.id","group.sizes","starts","maxgrpn","groups","ordered","call")
#     } else {
#       if (return.order) {
#         attributes(o) <- NULL
#         res <- list(length(f), grpuo, len, f, mg, ordered, o, call = match.call()) 
#         attr(res,"names") = c("N.groups","group.id","group.sizes","starts","maxgrpn","ordered","order","call")
#       } else {
#         res <- list(length(f), grpuo, len, f, mg, ordered, call = match.call()) 
#         attr(res,"names") = c("N.groups","group.id","group.sizes","starts","maxgrpn","ordered","call")
#       }
#     }
#   }
#   class(res) <- "GRP"
#   res
# }

# with starts and maxgrpn
# # or getGRP: more speed using tabulate ?? or code your own??
# GRP <- function(X, by = if(is.atomic(X)) NULL else seq_along(X), sort = TRUE, order = 1L, 
#                 na.last = FALSE, return.groups = TRUE, return.order = FALSE, ...) { # , gs = TRUE # o
#   if(is.call(by)) by = all.vars(by) # is.call is fastest!!!  attr(terms.formula(by), "term.labels")
#   o = data.table:::forderv(X, by, TRUE, sort, order, na.last) # faster calling c??
#   f = attr(o, "starts") # , exact = TRUE -> slightly slower 
#   if (length(o)) { # if ordered, returns 0 !!! perhaps rewrite Cforder function??
#     len = .Call(data.table:::Cuniqlengths, f, length(o)) # len = data.table:::uniqlengths(f, length(o))
#     # grpuo = rep.int(seq_along(len),len)[data.table:::forderv(o)] # faster way?? -> try to speed this up.. 
#     grpuo = .Call(data.table:::Cfrank, o, f, len, "dense") #-1L
#     res <- if (return.groups) {
#       # if (gs) 
#       if (return.order) { 
#         attributes(o) <- NULL
#         list(N.groups = length(f), group.id = grpuo, group.sizes = len, starts = f, maxgrpn = attr(o, "maxgrpn"), # speed up this part??
#              groups = if (inherits(X, "data.table")) X[o[f], by, with = FALSE] else if (is.atomic(X)) X[o[f]] else
#                X[o[f], by, drop = FALSE], ordered = c(GRP.sort = sort, initially.ordered = FALSE), order = o, call = match.call()) 
#       } else
#         list(N.groups = length(f), group.id = grpuo, group.sizes = len, starts = f, maxgrpn = attr(o, "maxgrpn"), # speed up this part??
#              groups = if (inherits(X, "data.table")) X[o[f], by, with = FALSE] else if (is.atomic(X)) X[o[f]] else
#                X[o[f], by, drop = FALSE], ordered = c(GRP.sort = sort, initially.ordered = FALSE), call = match.call())
#       # else 
#       # list(ng = length(f), g = grpuo,  # speed up this part??
#       #      groups = if (inherits(X, "data.table")) X[o[f], by, with = FALSE] else if (is.atomic(X)) X[o[f]] else
#       #        X[o[f], by, drop = FALSE], ordered = order)
#     } else {
#       # if (gs) 
#       if (return.order) list(N.groups = length(f), group.id = grpuo, group.sizes = len, 
#                              ordered = c(GRP.sort = sort, initially.ordered = FALSE), order = o, call = match.call()) else # else list(ng = length(f), g = grpuo, ordered = order)
#                                list(N.groups = length(f), group.id = grpuo, group.sizes = len, 
#                                     ordered = c(GRP.sort = sort, initially.ordered = FALSE), call = match.call())
#     }
#   } else { 
#     len = .Call(data.table:::Cuniqlengths, f, NROW(X)) # or cumsubtract !!!
#     grpuo = rep.int(seq_along(len),len)
#     # o = f # for unique. right???
#     res <- if (return.groups) {
#       # if (gs) 
#       if (return.order) {
#         # o2 = .Call(data.table:::Cfrank, seq_len(NROW(X)), f, len, "sequence") # good ???? This does not really give order, but unique row id's !!!
#         # attributes(o2) = attributes(o) # good ????
#         attributes(o) <- NULL
#         list(N.groups = length(f), group.id = grpuo, group.sizes = len, starts = f, maxgrpn = attr(o, "maxgrpn"), # speed up this part??
#              groups = if (inherits(X, "data.table")) X[f, by, with = FALSE] else if (is.atomic(X)) X[f] else
#                X[f, by, drop = FALSE], ordered = c(GRP.sort = sort, initially.ordered = TRUE), order = o2, call = match.call()) 
#       } else 
#         list(N.groups = length(f), group.id = grpuo, group.sizes = len, starts = f, maxgrpn = attr(o, "maxgrpn"), # speed up this part??
#              groups = if (inherits(X, "data.table")) X[f, by, with = FALSE] else if (is.atomic(X)) X[f] else
#                X[f, by, drop = FALSE], ordered = c(GRP.sort = sort, initially.ordered = TRUE), call = match.call())
#       
#       # else 
#       # list(ng = length(f), g = grpuo,  # speed up this part??
#       #      groups = if (inherits(X, "data.table")) X[f, by, with = FALSE] else if (is.atomic(X)) X[f] else
#       #        X[f, by, drop = FALSE], ordered = order)
#     } else {
#       #if (gs)
#       if (return.order) {
#         attributes(o) <- NULL
#         list(N.groups = length(f), group.id = grpuo, group.sizes = len, starts = f, maxgrpn = attr(o, "maxgrpn"),
#              ordered = c(GRP.sort = sort, initially.ordered = TRUE), order = o, call = match.call()) 
#       } else # else list(ng = length(f), g = grpuo, ordered = order)
#         list(N.groups = length(f), group.id = grpuo, group.sizes = len, starts = f, maxgrpn = attr(o, "maxgrpn"),
#              ordered = c(GRP.sort = sort, initially.ordered = TRUE), call = match.call())
#     }
#   }
#   class(res) <- "GRP"
#   res
#   # or only group: frankv(mtcars, cols = c("cyl","vs","am"), ties.method = "dense")
# }
