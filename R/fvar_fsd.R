library(Rcpp)
sourceCpp('R/C++/fvarsd.cpp')
sourceCpp('R/C++/fvarsda.cpp')
sourceCpp('R/C++/fvarsdl.cpp')
sourceCpp('R/C++/TRA.cpp')
sourceCpp('R/C++/TRAl.cpp')
sourceCpp('R/C++/TRAa.cpp')

# w.type = "frequency"
# Note: for principal innovations of this code see fsum.R !!

fsd <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) {
  UseMethod("fsd", x)
}
fsd.default <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fvarsdCpp(x,0L,0L,NULL,w,na.rm,stable.algo)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fvarsdCpp(x,length(lev),g,NULL,w,na.rm,stable.algo), lev))
      } else {
        if(is.factor(g)) return(fvarsdCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo)) else {
          g <- qG(g)
          return(fvarsdCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fvarsdCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo), group.names.GRP(g))) else 
        return(fvarsdCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fvarsdCpp(x,0L,0L,NULL,w,na.rm,stable.algo),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fvarsdCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fvarsdCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fvarsdCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fsd.matrix <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fvarsdmCpp(x,0L,0L,NULL,w,na.rm,stable.algo,TRUE,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fvarsdmCpp(x,length(lev),g,NULL,w,na.rm,stable.algo), list(lev, dimnames(x)[[2]])))
      } else {
        if(is.factor(g)) return(fvarsdmCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo)) else {
          g <- qG(g)
          return(fvarsdmCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fvarsdmCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo), list(group.names.GRP(g), dimnames(x)[[2]]))) else
        return(fvarsdmCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fvarsdmCpp(x,0L,0L,NULL,w,na.rm,stable.algo),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fvarsdmCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fvarsdmCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fvarsdmCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fsd.data.frame <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fvarsdlCpp(x,0L,0L,NULL,w,na.rm,stable.algo,TRUE,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fvarsdlCpp(x,length(lev),g,NULL,w,na.rm,stable.algo), lev))
      } else {
        if(is.factor(g)) return(fvarsdlCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo)) else {
          g <- qG(g)
          return(fvarsdlCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group.names.GRP(g))) 
        return(setRow.names(fvarsdlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo), if(is.double(groups)) paste0(groups) else groups)) else 
          return(fvarsdlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo)) 
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fvarsdlCpp(x,0L,0L,NULL,w,na.rm,stable.algo),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fvarsdlCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fvarsdlCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fvarsdlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fsd.grouped_df <- function(x, w = NULL, TRA = FALSE, na.rm = TRUE, drop.groups = FALSE, drop.w = TRUE, stable.algo = TRUE, ...) { 
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- names(x)
  gn <- match(names(g[[4]]), nam) 
  gn2 = gn <- gn[!is.na(gn)]
  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- x[[wn]] 
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    if(drop.w) if(drop.groups) gn <- c(gn,wn) else gn2 <- c(gn2,wn)
  }
  if(length(gn)) {
    if(drop.groups) { 
      if(TRA == FALSE) return(fvarsdlCpp(x[-gn],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo)) else {
        x <- x[-gn] 
        return(TRAlCpp(x,fvarsdlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo),g[[2]],TRAtoInt(TRA)))
      }
    } else {
      ax <- attributes(x)
      attributes(x) <- NULL 
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn2])
      if(TRA == FALSE) {
        ax[["row.names"]] <- .set_row_names(g[[1]])
        return(`attributes<-`(c(g[[4]],fvarsdlCpp(x[-gn2],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo)), ax))
      } else 
        return(`attributes<-`(c(x[gn],TRAlCpp(x[-gn2],fvarsdlCpp(x[-gn2],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo),g[[2]],TRAtoInt(TRA))), ax))
    }
  } else {
    if(TRA == FALSE)
      return(fvarsdlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo)) else 
        return(TRAlCpp(x,fvarsdlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo),g[[2]],TRAtoInt(TRA)))
  }
}


fvar <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) {
  UseMethod("fvar", x)
}
fvar.default <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fvarsdCpp(x,0L,0L,NULL,w,na.rm,stable.algo,FALSE)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fvarsdCpp(x,length(lev),g,NULL,w,na.rm,stable.algo,FALSE), lev))
      } else {
        if(is.factor(g)) return(fvarsdCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE)) else {
          g <- qG(g)
          return(fvarsdCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fvarsdCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE), group.names.GRP(g))) else 
        return(fvarsdCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fvarsdCpp(x,0L,0L,NULL,w,na.rm,stable.algo,FALSE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fvarsdCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fvarsdCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fvarsdCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fvar.matrix <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fvarsdmCpp(x,0L,0L,NULL,w,na.rm,stable.algo,FALSE,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fvarsdmCpp(x,length(lev),g,NULL,w,na.rm,stable.algo,FALSE), list(lev, dimnames(x)[[2]])))
      } else {
        if(is.factor(g)) return(fvarsdmCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE)) else {
          g <- qG(g)
          return(fvarsdmCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fvarsdmCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE), list(group.names.GRP(g), dimnames(x)[[2]]))) else
        return(fvarsdmCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fvarsdmCpp(x,0L,0L,NULL,w,na.rm,stable.algo,FALSE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fvarsdmCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fvarsdmCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fvarsdmCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fvar.data.frame <- function(x, g = NULL, w = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, stable.algo = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fvarsdlCpp(x,0L,0L,NULL,w,na.rm,stable.algo,FALSE,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fvarsdlCpp(x,length(lev),g,NULL,w,na.rm,stable.algo,FALSE), lev))
      } else {
        if(is.factor(g)) return(fvarsdlCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE)) else {
          g <- qG(g)
          return(fvarsdlCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group.names.GRP(g))) 
        return(setRow.names(fvarsdlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE), if(is.double(groups)) paste0(groups) else groups)) else 
          return(fvarsdlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE)) 
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fvarsdlCpp(x,0L,0L,NULL,w,na.rm,stable.algo,FALSE),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fvarsdlCpp(x,fnlevels(g),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fvarsdlCpp(x,attr(g,"N.groups"),g,NULL,w,na.rm,stable.algo,FALSE),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fvarsdlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fvar.grouped_df <- function(x, w = NULL, TRA = FALSE, na.rm = TRUE, drop.groups = FALSE, drop.w = TRUE, stable.algo = TRUE, ...) { 
  g <- GRP.grouped_df(x)
  wsym <- deparse(substitute(w))
  nam <- names(x)
  gn <- match(names(g[[4]]), nam) 
  gn2 = gn <- gn[!is.na(gn)]
  if(!(wsym == "NULL" || is.na(wn <- match(wsym, nam)))) {
    w <- x[[wn]] 
    if(any(gn == wn)) stop("Weights coincide with grouping variables!")
    if(drop.w) if(drop.groups) gn <- c(gn,wn) else gn2 <- c(gn2,wn)
  }
  if(length(gn)) {
    if(drop.groups) { 
      if(TRA == FALSE) return(fvarsdlCpp(x[-gn],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE)) else {
        x <- x[-gn] 
        return(TRAlCpp(x,fvarsdlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE),g[[2]],TRAtoInt(TRA)))
      }
    } else {
      ax <- attributes(x)
      attributes(x) <- NULL 
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn2])
      if(TRA == FALSE) {
        ax[["row.names"]] <- .set_row_names(g[[1]])
        return(`attributes<-`(c(g[[4]],fvarsdlCpp(x[-gn2],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE)), ax))
      } else 
        return(`attributes<-`(c(x[gn],TRAlCpp(x[-gn2],fvarsdlCpp(x[-gn2],g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE),g[[2]],TRAtoInt(TRA))), ax))
    }
  } else {
    if(TRA == FALSE)
      return(fvarsdlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE)) else 
        return(TRAlCpp(x,fvarsdlCpp(x,g[[1]],g[[2]],g[[3]],w,na.rm,stable.algo,FALSE),g[[2]],TRAtoInt(TRA)))
  }
}
