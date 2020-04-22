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
  addrn <- function(x) if(any(attr(x, "names") == row.names)) x else c(`names<-`(list(attr(x, "row.names")),row.names),x) # faster way??
  attol <- function(x) {
    # class(x) <- NULL # tables are also arrays, although only 1D, not because of the class but because they have a dimension attribute.
    if (length(d <- dim(x)) > 1L) { # is.array(x) # length could also be 0...
      # d <- dim(x)
      if (length(d) > 2L) { # breaking down HDA
        dn <- dimnames(x)
        dim(x) <- c(d[1L], prod(d[-1L]))
        if (!is.null(dn)) {
          for (i in 2L:length(d)) if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
          dimnames(x) <- list(dn[[1L]], interact_names(dn[-1L])) # Good??
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
    class(y) <- NULL # perhaps unclassing y would put more safety ?? -> yes !!
    ident <- vapply(y, idf, 1L) # `attributes<-`(y, NULL) # possibly you can still get a few microseconds in the apply commands, but beware, this removes names in output!!
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
      nams <- attr(l, "names")
      ids <- which(nams == id.names)
      nid <- length(ids)
      if(nid > 1L) { # Still make sure row.names comes after ids, even if only one id!!
        nids <- seq_len(nid)
        attr(l, "names")[ids] <- if(length(idcols) == nid) idcols else paste(id.names, nids, sep = ".")
        if(keep.row.names) { # New!! It seems it lost a bot of speed through this part!!
          rn <- which(nams == row.names) # New!!
          if(!all(ids == nids) || rn != nid + 1L) .Call(C_setcolorder, l, c(ids, rn, seq_along(l)[-c(ids,rn)]))  # l <- l[c(ids,rn,seq_along(l)[-c(ids,rn)])] # New!! efficient? could replace only rownames if one of the conditions holds
        } else if (!all(ids == nids)) .Call(C_setcolorder, l, c(ids, seq_along(l)[-ids])) # l <- l[c(ids,seq_along(l)[-ids])] # Old!! before row.names!!
      }
    } else if (keep.row.names) { # New!!
      rn <- which(attr(l, "names") == row.names) # New!!
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
