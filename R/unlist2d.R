# library(Rcpp)
# sourceCpp("C:/Users/Sebastian Krantz/Documents/R/C++/mrtl_type_dispatch.cpp")
# Crbindlist <- data.table:::Crbindlist
# Csetcolorder <- data.table:::Csetcolorder

unlist2d <- function(l, idcols = ".id", row.names = FALSE, recursive = TRUE, id.factor = FALSE, DT = FALSE) {
  if (!is.list(l)) return(l) #stop("l is not a list")
  make.ids <- idcols[1L] != FALSE
  if(make.ids) id.names <- if(idcols[1L] == TRUE) ".id" else idcols[1L]
  keep.row.names <- row.names[1L] != FALSE
  if(keep.row.names) row.names <- if(row.names[1L] == TRUE) "row.names" else row.names[1L]

  DFDTl <- function(l) {
    attr(l, "row.names") <- .set_row_names(length(l[[1L]]))
    class(l) <- if(DT) c("data.table","data.frame") else "data.frame"
    l
  }
  idf <- function(x) if(inherits(x, "data.frame")) 2L else if (is.null(x)) 1L else 3L*is.atomic(x) # faster way ??
  addrn <- function(x) if(any(names(x) == row.names)) x else c(`names<-`(list(attr(x, "row.names")),row.names),x) # faster way??
  attol <- function(x) {
    if (is.array(x)) {
      d <- dim(x)
      if (length(d) > 2L) { # breaking down HDA
        dn <- dimnames(x)
        dim(x) <- c(d[1L], prod(d[-1L]))
        if (!is.null(dn)) {
          for (i in 2L:length(d)) if (is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
          dimnames(x) <- list(dn[[1L]], interaction(expand.grid(dn[-1L]))) # Good??
        }
      }
      dn <- dimnames(x)
      if(keep.row.names && !is.null(dn[[1L]]))
        x <- `names<-`(c(list(dn[[1L]]), mctl(x)), c(row.names, dn[[2L]])) else
        x <- `names<-`(mctl(x), dn[[2L]])
    } else x <- as.vector(x, "list")
    if (is.null(names(x))) names(x) <- paste0("V", seq_along(x))     # it seems this is not yet working for all (i.e. model objects..), also perhaps not start at V1, sepending on what other columsn there are.. i.e start at the right position??
    return(x)
  }
  ul2d <- function(y) {
    if(inherits(y, "data.frame") || is.atomic(y)) return(y)
    ident <- vapply(y, idf, 1L)
    if(is.list(y) && all(ident > 0L)) {
      at <- ident == 3L
      if(any(at)) y[at] <- lapply(y[at], attol)
      if(keep.row.names && any(df <- ident == 2L)) y[df] <- lapply(y[df], addrn) # better cbind for data.table???? or x[["row.names"]] =.. and the sort later??
        if(make.ids) {
          if(id.factor) {
            y <- y[ident != 1L] # better way?? y[ident!=1L] = NULL??
            nam <- names(y)
            names(y) <- NULL
            y <- DFDTl(.Call(C_rbindlist, y, TRUE, TRUE, id.names))
            # attributes(y[[1L]]) <- list(levels = nam, class = c("ordered", "factor"))
            setattributes(y[[1L]], pairlist(levels = nam, class = c("ordered", "factor"))) # a lot faster !!
            return(y)
          } else return(DFDTl(.Call(C_rbindlist, y[ident != 1L], TRUE, TRUE, id.names)))
        } else return(DFDTl(.Call(C_rbindlist, y[ident != 1L], TRUE, TRUE, NULL)))
    } else lapply(y, ul2d)
  }

  l <- ul2d(l)
  if(recursive) {
    while(!inherits(l, "data.frame")) l <- ul2d(l)
    if(make.ids) {
      nams <- names(l)
      ids <- which(nams == id.names)
      nid <- length(ids)
      if(nid > 1L) { # Still make sure row.names comes after ids, even if only one id!!
        nids <- seq_len(nid)
        names(l)[ids] <- if(length(idcols) == nid) idcols else paste(id.names, nids, sep = ".")
        if(keep.row.names) { # New!! It seems it lost a bot of speed through this part!!
          rn <- which(nams == row.names) # New!!
          if(!all(ids == nids) || rn != nid + 1L) .Call(C_setcolorder, l, c(ids, rn, seq_along(l)[-c(ids,rn)]))  # l <- l[c(ids,rn,seq_along(l)[-c(ids,rn)])] # New!! efficient? could replace only rownames if one of the conditions holds
        } else if (!all(ids == nids)) .Call(C_setcolorder, l, c(ids, seq_along(l)[-ids])) # l <- l[c(ids,seq_along(l)[-ids])] # Old!! before row.names!!
      }
    } else if (keep.row.names) { # New!!
      rn <- which(names(l) == row.names) # New!!
      if(rn != 1L) .Call(C_setcolorder, l, c(rn,seq_along(l)[-rn]))  # l <- l[c(rn,seq_along(l)[-rn])] # New!!
    }
  }
  attr(l, ".internal.selfref") <- NULL
  # setattrib(l, ".internal.selfref", NULL) # what if recursive = FALSE -> doesn't work.. but takes more time!!
  return(l)
}
# STill do a lot of checking !!
# perhaps process unnames vectors as columns, and names vectors as rows?? -> nah, it's rbinding

# Examples:

# # example:
# l = lapply(split(iris[1:4], iris[5]), function(x)list(N = fNobs(x), mean = fmean(x), sd = fsd(x)))
#
# # neat example:
# l = list(a = mtcars[1:8],b = list(c = mtcars[4:11], d = list(e = mtcars[2:10])))
# unlist2d(rapply2d(l,colMeans), recursive = FALSE)
# unlist2d(rapply2d(l,colMeans)) # now error again!!!
# unlist2d(l)
# unlist2d(rapply2d(l,dapply,sd))
#
# unlist2d(split(iris,iris["Species"]))
#
# nl = lapply(split(mtcars,mtcars[2]),function(x)split(x,x["vs"]))
# str(nl)
# View(unlist2d(nl))
# View(unlist2d(nl, idcols = FALSE))
# View(unlist2d(nl, idcols = c(".cyl",".vs")))
# str(unlist2d(nl, recursive = FALSE)) # why is .id a character string?? -> names!!
#
# # Neat example:
# # list.elem(IRF) %>% rapply2d(colSums) %>% unlist2d(c("type","shock")) %>% filter(type == "irf") %>% num.vars %>% dapply(function(x)sum(abs(x)),MARGIN = 1)
#
# unlist2d(qsu(mtcars,~cyl,~vs, data.out = TRUE)) # not dim, but is.data.frame
#
# unlist2d(rapply2d(reg.elem(SV),dim)) # error!! lots of nested NULL's .. doesn't perform
#

# Version befor: Above I added some efficiency !!
# unlist2dOld <- function(X, idcols = ".id", row.names = FALSE, recursive = TRUE, id.factor = FALSE, DT = FALSE) {
#   if (!is.list(X)) return(X) #stop("X is not a list")
#   make.ids <- idcols[1L] != FALSE
#   id.names <- if (idcols[1L] == TRUE) ".id" else idcols[1L] # New!! Made faster by making id.names size 1, and referring to idcols later on. But is it foolproof??
#   row.names <- if (row.names[1L] == TRUE) "row.names" else row.names[1L]
#   keep.row.names <- row.names != FALSE
#   idf <- function(x) if(is.data.frame(x)) 2L else if (is.null(x)) 1L else 3L*is.atomic(x) # faster way ??
#   addrn <- function(x) if(any(names(x) == row.names)) x else c(`names<-`(list(row.names.data.frame(x)),row.names),x) # faster way??
#   ul2d <- function(y) {
#     if (is.data.frame(y) || is.atomic(y)) return(y)
#     ident <- unlist(lapply(y,idf), use.names = FALSE)
#     if (is.list(y) && all(ident>0L)) {
#       ar <- ident == 3L
#       if (any(ar)) { # Faster!!
#         y[ar] <- lapply(y[ar], function(x){
#           if (is.array(x)) { # fast??
#             d <- dim(x)
#             if (length(d) > 2L) {
#               dn <- dimnames(x)
#               dim(x) <- c(d[1L], prod(d[-1L]))
#               if (!is.null(dn)) {
#                 if (length(dn[[1L]]))
#                   rownames(x) <- dn[[1L]]
#                 for (i in 2L:length(d)) if (is.null(dn[[i]]))
#                   dn[[i]] <- seq_len(d[i])
#                 colnames(x) <- interaction(expand.grid(dn[-1L])) # Good??
#               }
#             }
#             cnx <- colnames(x)
#             if (keep.row.names && !is.null(rnx <- rownames(x))) { # added this!! assignment in if is a tiny bit slower, but here it makes sense
#               x <- c(list(rnx),mctl(x)) # faster way??
#               names(x) <- c(row.names,cnx)
#             } else {
#               x <- mctl(x)
#               names(x) <- cnx
#             }
#           } else x <- as.vector(x,"list")
#           # it seems this is not yet working for all (i.e. model objects..), also perhaps not start at V1, sepending on what other columsn there are.. i.e start at the right position??
#           if (is.null(names(x))) names(x) <- paste0("V",seq_along(x))
#           x
#         })
#       }
#       if (keep.row.names && any(df <- ident == 2L)) y[df] <- lapply(y[df],addrn) # better cbind for data.table???? or x[["row.names"]] =.. and the sort later??
#       #if (make.ids) dplyr::bind_rows(y, .id = "id.id.id") else dplyr::bind_rows(y) # the former stillneeds to be changed # data.table is faster for larger tasks!!
#       if (DT) {
#         if (make.ids) {
#           if (id.factor) {
#             y = y[ident!=1L] # better way?? y[ident!=1L] = NULL??
#             nam = names(y)
#             names(y) = NULL
#             y = data.table::rbindlist(y, fill = recursive, idcol = idcols[1L])
#             data.table::setattr(y[[1]], "levels", nam)
#             data.table::setattr(y[[1]], "class", c("ordered", "factor"))
#             y
#           } else data.table::rbindlist(y[ident!=1L], fill = recursive, idcol = idcols[1L])
#         } else data.table::rbindlist(y[ident!=1L], fill = recursive)
#       } else {
#         if (make.ids) {
#           if (id.factor) {
#             y = y[ident!=1L] # better way?? y[ident!=1L] = NULL
#             nam = names(y)
#             names(y) = NULL
#             y = data.table::setattr(data.table::rbindlist(y, fill = recursive, idcol = idcols[1L]), "class", "data.frame")
#             data.table::setattr(y[[1]], "levels", nam)
#             data.table::setattr(y[[1]], "class", c("ordered", "factor"))
#             y
#           } else data.table::setattr(data.table::rbindlist(y[ident!=1L], fill = recursive, idcol = idcols[1L]), "class", "data.frame")
#         } else data.table::setattr(data.table::rbindlist(y[ident!=1L], fill = recursive), "class", "data.frame")
#       }
#     } else lapply(y,ul2d)
#   }
#   X <- ul2d(X)
#   if (recursive) {
#     while(!is.data.frame(X)) X <- ul2d(X)
#     if (make.ids) {
#       nams <- names(X)
#       ids <- which(nams == id.names)
#       if (length(ids)>1L) { # Still make sure row.names comes after ids, even if only one id!!
#         nids <- seq_along(ids)
#         if (length(idcols) == length(ids))
#           names(X)[ids] <- idcols else
#             names(X)[ids] <- paste0(id.names,".",nids)
#           if (keep.row.names) { # New!! It seems it lost a bot of speed through this part!!
#             rn <- which(nams == row.names) # New!!
#             if (!all(ids==nids) || rn!=max(nids)+1L) data.table::setcolorder(X,c(ids,rn,seq_along(X)[-c(ids,rn)])) # X <- X[c(ids,rn,seq_along(X)[-c(ids,rn)])] # New!! efficient? could replace only rownames if one of the conditions holds
#           } else if (!all(ids==nids)) data.table::setcolorder(X,c(ids,seq_along(X)[-ids])) # X <- X[c(ids,seq_along(X)[-ids])] # Old!! before row.names!!
#       }
#     } else if (keep.row.names) { # New!!
#       rn <- which(names(X) == row.names) # New!!
#       if (rn!=1L) data.table::setcolorder(X,c(rn,seq_along(X)[-rn])) # X <- X[c(rn,seq_along(X)[-rn])] # New!!
#     }
#   }
#   data.table::setattr(X, ".internal.selfref", NULL) # what if recursive = FALSE -> doesn't work.. but takes more time!!
#   X
# }

# Prior version before sorting row.names and efficiency removement in idcols not subsetting id.names
# unlist2d <- function(X, idcols = ".id", recursive = TRUE, row.names = FALSE) {
#   if (!is.list(X)) return(X) #stop("X is not a list")
#   make.ids <- idcols[1L] != FALSE
#   id.names <- if (idcols[1L] == TRUE) ".id" else idcols
#   row.names <- if (row.names[1L] == TRUE) "row.names" else row.names[1L]
#   keep.row.names <- row.names != FALSE
#   idf <- function(x) if(is.data.frame(x)) 2L else if (is.null(x)) 1L else 3L*is.atomic(x) # faster way ??
#   addrn <- function(x) if(any(names(x) == row.names)) x else c(`names<-`(list(row.names.data.frame(x)),row.names),x) # faster way??
#   ul2d <- function(y) {
#     if (is.data.frame(y) || is.atomic(y)) return(y)
#     ident <- unlist(lapply(y,idf), use.names = FALSE)
#     if (is.list(y) && all(ident>0L)) {
#       ar <- ident == 3L
#       if (any(ar)) { # Faster!!
#         y[ar] <- lapply(y[ar], function(x){
#           if (is.array(x)) { # fast??
#             d <- dim(x)
#             if (length(d) > 2L) {
#               dn <- dimnames(x)
#               dim(x) <- c(d[1L], prod(d[-1L]))
#               if (!is.null(dn)) {
#                 if (length(dn[[1L]]))
#                   rownames(x) <- dn[[1L]]
#                 for (i in 2L:length(d)) if (is.null(dn[[i]]))
#                   dn[[i]] <- seq_len(d[i])
#                 colnames(x) <- interaction(expand.grid(dn[-1L])) # Good??
#               }
#             }
#             cnx <- colnames(x)
#             if (keep.row.names && !is.null(rnx <- rownames(x))) { # added this!! assignment in if is a tiny bit slower, but here it makes sense
#               x <- c(list(rnx),mctl(x)) # faster way??
#               names(x) <- c(row.names,cnx)
#             } else {
#               x <- mctl(x)
#               names(x) <- cnx
#             }
#           } else x <- as.vector(x,"list")
#           if (is.null(names(x))) names(x) <- paste0("V",seq_along(x))
#           x
#         })
#       }
#       if (keep.row.names && any(df <- ident == 2L)) y[df] <- lapply(y[df],addrn) # better cbind for data.table???? or x[["row.names"]] =.. and the sort later??
#       #if (make.ids) dplyr::bind_rows(y, .id = "id.id.id") else dplyr::bind_rows(y) # the former stillneeds to be changed # data.table is faster for larger tasks!!
#       res <- if (make.ids) data.table::rbindlist(y[ident!=1L], fill = recursive, idcol = idcols[1L]) else data.table::rbindlist(y[ident!=1L], fill = recursive)
#       class(res) <- "data.frame" # Faster using setDF???
#       res
#     } else lapply(y,ul2d)
#   }
#   X <- ul2d(X)
#   if (recursive) {
#     while(!is.data.frame(X)) X <- ul2d(X)
#     if (make.ids) {
#       nams <- names(X)
#       ids <- which(nams == id.names[1L])
#       if (length(ids)>1L) {
#         nids <- seq_along(ids)
#         if (length(id.names) == length(ids))
#           names(X)[ids] <- id.names else
#             names(X)[ids] <- paste0(id.names[1L],".",nids)
#           if (!all(ids==nids)) X <- X[c(ids,seq_along(X)[-ids])]
#       }
#     }
#   }
#   X
# }


# Prior Version without row.names:
# unlist2d <- function(X, idcols = ".id", recursive = TRUE) {
#   if (!is.list(X)) return(X) #stop("X is not a list")
#   make.ids = idcols[1L] != FALSE
#   id.names = if (idcols[1L] == TRUE) ".id" else idcols
#   idf <- function(x) if(is.data.frame(x)) 2L else if (is.null(x)) 1L else 3L*is.atomic(x) # faster way ??
#   ul2d <- function(y) {
#     if (is.data.frame(y) || is.atomic(y)) return(y)
#     ident = unlist(lapply(y,idf), use.names = FALSE)
#     if (is.list(y) && all(ident>0L)) {
#       ar = ident == 3L
#       if (any(ar)) { # Faster!!
#         y[ar] = lapply(y[ar], function(x){
#           if (is.array(x)) { # fast??
#             d = dim(x)
#             if (length(d) > 2L) {
#               dn <- dimnames(x)
#               dim(x) <- c(d[1L], prod(d[-1L]))
#               if (!is.null(dn)) {
#                 if (length(dn[[1L]]))
#                   rownames(x) <- dn[[1L]]
#                 for (i in 2L:length(d)) if (is.null(dn[[i]]))
#                   dn[[i]] <- seq_len(d[i])
#                 colnames(x) <- interaction(expand.grid(dn[-1L])) # Good??
#               }
#             }
#             cnx = colnames(x)
#             x = mctl(x)
#             names(x) = cnx
#           } else x = as.vector(x,"list")
#           if (is.null(names(x))) names(x) = paste0("V",seq_along(x))
#           x
#         })
#       } #if (make.ids) dplyr::bind_rows(y, .id = "id.id.id") else dplyr::bind_rows(y) # the former stillneeds to be changed # data.table is faster for larger tasks!!
#       res = if (make.ids) data.table::rbindlist(y[ident!=1L], fill = recursive, idcol = idcols[1L]) else data.table::rbindlist(y[ident!=1L], fill = recursive)
#       class(res) = "data.frame"
#       res
#     } else lapply(y,ul2d)
#   }
#   X = ul2d(X)
#   if (recursive) {
#     while(!is.data.frame(X)) X = ul2d(X)
#     if (make.ids) {
#       nams = names(X)
#       ids = which(nams == id.names[1L])
#       if (length(ids)>1L) {
#         nids = seq_along(ids)
#         if (length(id.names) == length(ids))
#           names(X)[ids] = id.names else
#             names(X)[ids] = paste0(id.names[1L],".",nids)
#           if (!all(ids==nids)) X = X[c(ids,seq_along(X)[-ids])]
#       }
#     }
#   }
#   X
# }



# Old Version:
# unlist2dold <- function(X, make.IDs = TRUE, ID.names = "LID", make.row.names = FALSE, recursive = TRUE, check = TRUE) {
#   # Still make LID a fully optional choice i.e. a variable may be named LID!!
#   if (!is.list(X)) stop("x is not a list")
#   ns <- function(x) if (is.array(x)) colnames(x) else if (is.null(names(x))) seq_along(x) else names(x) # right?? -> yup, good!!
#   hasdimorvec <- function(x) is.data.frame(x) || is.atomic(x) #!is.null(dim(x))
#   nrowor1 <- function(x) if (!is.null(dim(x))) nrow(x) else 1
#   mn <- function(x) { # make names
#     rows = unlist(lapply(x,nrowor1), use.names = FALSE)
#     nam = if (is.null(names(x))) {
#       seq_along(x)
#     } else { nx = names(x)
#     nas = is.na(nx) | nx == ""
#     nx[nas] = which(nas)
#     nx
#     }
#     rep(nam, rows)
#   }
#   ul2d <- function(y) { df = unlist(lapply(y,is.data.frame), use.names = FALSE)
#   ar = unlist(lapply(y,is.atomic), use.names = FALSE)
#   if (is.list(y) && !is.data.frame(y) && all(df+ar>0)) { #is.null(dim(y))
#     if (!all(df)) { # Only if this is not true, do this code, else use bind_rows
#       if (check) {
#         nams = lapply(y,ns) # do.call doesn't work here
#         maxid = max(unlist(lapply(nams,function(x)sum(x=="id.id.id")), use.names = FALSE))
#         nams = unlist(nams, use.names = FALSE)
#         if (maxid!=0) nams = c(rep("id.id.id",maxid),unique.default(setdiff(nams,"id.id.id"))) else nams = unique.default(nams)
#         #y[df] = lapply(y[df], function(j){}) doing in two parts for ar and df could be faster
#         y = lapply(y, function(j){
#           namj = ns(j) # rbind seems to be capable of ordering, otherwise uncomment all below!!
#           v = nams[!nams %in% namj] #setdiff(nams,ns(j)) # problem here: set only takes unique # can speed up? i.e. use the above??
#           namj = c(v,namj)
#           nj = which(namj!="id.id.id")
#           nms = match(nams,namj)[nj] # correct?? prove!!
#           # add condition if names not the same??
#           if (is.array(j)) { # necessary here?? -> yes
#             if (length(v)) {
#               vmat = matrix(NA, ncol = length(v), nrow = nrow(j))
#               colnames(vmat) = v
#               j = cbind(vmat,j)
#             }
#             j[, nj] = j[, nms] # drop? # make this faster? conditional??
#             colnames(j)[nj] = nams[nj]
#           } else {
#             if (length(v)) j[v] = NA
#             j[nj] = j[nms] # make this faster? conditional??
#             names(j)[nj] = nams[nj]
#           }
#           j
#         })
#       } else if (recursive && make.IDs) { # possibility to make faster i.e. keeping track of list somehow??
#         nids = unlist(lapply(y,function(j)sum(ns(j)=="id.id.id")), use.names = FALSE)
#         if (!all(nids==nids[1])) {
#           maxnids = max(nids)
#           #wm = which(nids!=maxnids)
#           y = lapply(seq_along(y),function(j){ # faster way??
#             dLID = maxnids - nids[j]
#             while (dLID>0) { # not needed, the whole function is recursive -> FALSE, it is needed!!
#               y[[j]] = cbind(id.id.id = NA, y[[j]]) # faster way??
#               dLID = dLID-1
#             }
#             y[[j]]
#           })
#           #y[wm] = lapply(y[wm],function(j)cbind(id.id.id = NA,j)) # faster way?
#         }
#       }
#       #cbind(LID = mn(y),do.call(plyr::rbind.fill,y))
#       if (make.IDs) {  # the !is.null(names(y)) condition below is because otherwise it will produce a character matrix if y has names
#         anylist = any(df) #unlist(lapply(y,is.list), use.names = FALSE))
#         if (anylist || !is.null(names(y))) { # faster calling method directly??
#           if (anylist) { cbind(id.id.id = mn(y),do.call(rbind.data.frame,c(y, make.row.names = make.row.names, stringsAsFactors = FALSE))) # Stringasfactors is still not executing properly!!
#           } else cbind.data.frame(id.id.id = mn(y), do.call(rbind,y), stringsAsFactors = FALSE) # Set global option for stringasfactors
#         } else cbind(id.id.id = mn(y),do.call(rbind,y))
#       } else {
#         if (any(df)) { #unlist(lapply(y,is.list), use.names = FALSE)
#           do.call(rbind.data.frame,c(y, make.row.names = make.row.names, stringsAsFactors = FALSE))
#         } else do.call(rbind,y)
#       }
#     } else {
#       #if (make.IDs) dplyr::bind_rows(y, .id = "id.id.id") else dplyr::bind_rows(y) # the former stillneeds to be changed
#       # data.table is faster for larger tasks!!
#       if (make.IDs) data.table::setDF(data.table::rbindlist(y, fill = check || recursive, idcol = "id.id.id")) else data.table::setDF(data.table::rbindlist(y, fill = check || recursive))
#     }
#   } else if (hasdimorvec(y)) {
#     y
#   } else {
#     lapply(y,ul2d)
#   }
#   }
#   X = ul2d(X)
#   if (recursive) {
#     while(!hasdimorvec(X)) X = ul2d(X) #is.null(dim(x))
#     if (make.IDs) {
#       nams = ns(X)
#       ids = which(nams=="id.id.id")
#       if (length(ids)>1) {
#         if (!all(ids==seq_along(ids))) { # This is only because of dplyr
#           names(X)[ids] = nams[ids] = paste0("id.id.id",ids)
#           X = X[c(nams[ids],nams[-ids])]  #x = cbind(x[,ids, drop = FALSE],x[,-ids, drop = FALSE]) # I altered the check syntax so this is no longer necessary. Hopefully it is faster!! -> a bit, yes!!
#           ids = seq_along(ids)
#         }
#         if (!all(ID.names==ID.names[1])) {
#           if (length(ID.names) == length(ids)) {
#             colnames(X)[ids] = ID.names # seq_along(ids)
#           } else colnames(X)[ids] = paste0(ID.names[1],".",ids) #seq_along(ids)
#         } else colnames(X)[ids] = paste0(ID.names[1],".",ids) #seq_along(ids)
#       } else colnames(X)[ids] = ID.names[1]
#     }
#   } else if (make.IDs) { # && ID.names[1]!="LID"
#     X = rapply2d(X, function(j){ # speed up?
#       if (is.array(j)) colnames(j)[colnames(j)=="id.id.id"] = ID.names[1] else names(j)[names(j)=="id.id.id"] = ID.names[1]
#       j
#     })
#   }
#   X
# }
# neat example:
#l = list(a = mtcars[1:8],b = list(c = mtcars[4:11], d = list(e = mtcars[2:10])))
#unlist2d(rapply2d(l,colMeans), recursive = FALSE)
#unlist2d(rapply2d(l,colMeans)) # now error again!!!
#unlist2d(l)
#unlist2d(rapply2d(l,dapply,sd))
# error:
#unlist2d(qsu(mtcars,~cyl,~vs, data.out = TRUE)) # not dim, but is.data.frame
