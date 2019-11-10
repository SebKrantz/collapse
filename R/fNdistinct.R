library(Rcpp)
sourceCpp('R/C++/fNdistinct.cpp')
sourceCpp('R/C++/fNdistincta.cpp')
sourceCpp('R/C++/TRA.cpp')
sourceCpp('R/C++/TRAl.cpp')
sourceCpp('R/C++/TRAa.cpp')

# Better data.frame funtcion ??
# For foundational changes to this code see fsum.R !!
# also: matrix method needs memory equal to size of the object, while data.frame method does not need any memory !!

fNdistinct <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) {
  UseMethod("fNdistinct", x)
}
fNdistinct.default <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fNdistinctCpp(x, narm = na.rm)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`names<-`(fNdistinctCpp(x,length(lev),g,NULL,na.rm), lev))
      } else {
        if(is.factor(g)) return(fNdistinctCpp(x,fnlevels(g),g,NULL,na.rm)) else {
          g <- qG(g)
          return(fNdistinctCpp(x,attr(g,"N.groups"),g,NULL,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`names<-`(fNdistinctCpp(x,g[[1]],g[[2]],g[[3]],na.rm), group.names.GRP(g))) else 
        return(fNdistinctCpp(x,g[[1]],g[[2]],g[[3]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRACpp(x,fNdistinctCpp(x, narm = na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRACpp(x,fNdistinctCpp(x,fnlevels(g),g,NULL,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRACpp(x,fNdistinctCpp(x,attr(g,"N.groups"),g,NULL,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRACpp(x,fNdistinctCpp(x,g[[1]],g[[2]],g[[3]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fNdistinct.matrix <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, ...) { 
  if(TRA == FALSE) {
    if(is.null(g)) return(fNdistinctmCpp(x,0L,0L,NULL,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(`dimnames<-`(fNdistinctmCpp(x,length(lev),g,NULL,na.rm), list(lev, dimnames(x)[[2]])))
      } else {
        if(is.factor(g)) return(fNdistinctmCpp(x,fnlevels(g),g,NULL,na.rm)) else {
          g <- qG(g)
          return(fNdistinctmCpp(x,attr(g,"N.groups"),g,NULL,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names) return(`dimnames<-`(fNdistinctmCpp(x,g[[1]],g[[2]],g[[3]],na.rm), list(group.names.GRP(g), dimnames(x)[[2]]))) else
        return(fNdistinctmCpp(x,g[[1]],g[[2]],g[[3]],na.rm))
    }
  } else {
    if(is.null(g)) return(TRAmCpp(x,fNdistinctmCpp(x, narm = na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAmCpp(x,fNdistinctmCpp(x,fnlevels(g),g,NULL,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAmCpp(x,fNdistinctmCpp(x,attr(g,"N.groups"),g,NULL,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAmCpp(x,fNdistinctmCpp(x,g[[1]],g[[2]],g[[3]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fNdistinct.data.frame <- function(x, g = NULL, TRA = FALSE, na.rm = TRUE, drop = TRUE, use.g.names = TRUE, ...) {
  if(TRA == FALSE) {
    if(is.null(g)) return(fNdistinctlCpp(x,0L,0L,NULL,na.rm,drop)) else if (is.atomic(g)) {
      if(use.g.names && !inherits(x, "data.table")) {
        if(!is.factor(g)) g <- qF(g)
        lev <- attr(g, "levels")
        return(setRow.names(fNdistinctlCpp(x,length(lev),g,NULL,na.rm), lev))
      } else {
        if(is.factor(g)) return(fNdistinctlCpp(x,fnlevels(g),g,NULL,na.rm)) else {
          g <- qG(g)
          return(fNdistinctlCpp(x,attr(g,"N.groups"),g,NULL,na.rm))
        }
      }
    } else {
      if(!is.GRP(g)) g <- if(use.g.names) GRP(g) else GRP(g, return.groups = FALSE)
      if(use.g.names && !inherits(x, "data.table") && !is.null(groups <- group.names.GRP(g))) 
        return(setRow.names(fNdistinctlCpp(x,g[[1]],g[[2]],g[[3]],na.rm), if(is.double(groups)) paste0(groups) else groups)) else 
          return(fNdistinctlCpp(x,g[[1]],g[[2]],g[[3]],na.rm)) 
    }
  } else {
    if(is.null(g)) return(TRAlCpp(x,fNdistinctlCpp(x, narm = na.rm),0L,TRAtoInt(TRA))) else if (is.atomic(g)) {
      if(is.factor(g)) return(TRAlCpp(x,fNdistinctlCpp(x,fnlevels(g),g,NULL,na.rm),g,TRAtoInt(TRA))) else {
        g <- qG(g)
        return(TRAlCpp(x,fNdistinctlCpp(x,attr(g,"N.groups"),g,NULL,na.rm),g,TRAtoInt(TRA)))
      }
    } else {
      if(!is.GRP(g)) g <- GRP(g, return.groups = FALSE)
      return(TRAlCpp(x,fNdistinctlCpp(x,g[[1]],g[[2]],g[[3]],na.rm),g[[2]],TRAtoInt(TRA)))
    }
  }
}
fNdistinct.grouped_df <- function(x, TRA = FALSE, na.rm = TRUE, drop.groups = FALSE, ...) { 
  g <- GRP.grouped_df(x)
  gn <- match(names(g[[4]]), names(x))
  gn <- gn[!is.na(gn)]
  if(length(gn)) {
    if(drop.groups) {
      if(TRA == FALSE) return(fNdistinctlCpp(x[-gn],g[[1]],g[[2]],g[[3]],na.rm)) else {
        x <- x[-gn] 
        return(TRAlCpp(x,fNdistinctlCpp(x,g[[1]],g[[2]],g[[3]],na.rm),g[[2]],TRAtoInt(TRA)))
      }
    } else {
      ax <- attributes(x)
      attributes(x) <- NULL 
      ax[["names"]] <- c(ax[["names"]][gn], ax[["names"]][-gn])
      if(TRA == FALSE) {
        ax[["row.names"]] <- .set_row_names(g[[1]])
        return(`attributes<-`(c(g[[4]],fNdistinctlCpp(x[-gn],g[[1]],g[[2]],g[[3]],na.rm)), ax))
      } else 
        return(`attributes<-`(c(x[gn],TRAlCpp(x[-gn],fNdistinctlCpp(x[-gn],g[[1]],g[[2]],g[[3]],na.rm),g[[2]],TRAtoInt(TRA))), ax))
    }
  } else {
    if(TRA == FALSE)
      return(fNdistinctlCpp(x,g[[1]],g[[2]],g[[3]],na.rm)) else 
        return(TRAlCpp(x,fNdistinctlCpp(x,g[[1]],g[[2]],g[[3]],na.rm),g[[2]],TRAtoInt(TRA)))
  }
}
