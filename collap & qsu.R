# Preliminaries to Load 
Mode <- function(x, na.rm = FALSE) {
  ax = attributes(x)
  if (na.rm == TRUE) x = x[!is.na(x)]
  ux <- unique(x)
  y = ux[which.max(tabulate(match(x, ux)))]
  attributes(y) = ax
  return(y)
}
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}


# Collap
collap <- function(X, by = NULL, FUN = mean, catFUN = Mode, factors = "as.categorical", 
                   custom = NULL, custom.names = TRUE, collapse = TRUE, sort = TRUE, 
                   reshape.long = FALSE, na.rm = TRUE, replace.nan = TRUE, as.list = FALSE, 
                   dropcat = FALSE, dropby = FALSE, show.statistic = TRUE,  
                   data.table = FALSE, parallel = FALSE, ...) {
  
  # Identifying the inputs: X
  cl = class(X)
  if (!all(cl == "data.frame")) {
    iv = is.vector(X)  
    if (iv) namx = make.names(deparse(substitute(X)))
    X = as.data.frame(X, stringsAsFactors = FALSE) 
    if (iv) names(X) = namx 
    if (any(cl == "data.table")) {
      data.table = TRUE
      isDT = TRUE
    } else isDT = FALSE
  } else isDT = FALSE
  
  
  # Some functions
  nantona <- function(x) {x[is.nan(x)] = NA; x} 
  quickdf <- function(l) { 
    class(l) <- "data.frame"
    attr(l, "row.names") <- .set_row_names(length(l[[1]]))
    l
  }
  ident <- function(x) { 
    y <- as.factor(x)
    l <- length(levels(y))
    s <- as.character(seq_len(l)) 
    n <- nchar(s) 
    levels(y) <- paste0(strrep("0", n[l] - n), s)
    as.character(y)
  }
  
  # Identifying the inputs: FUN, catFUN and custom
  if (!is.null(custom)) FUN = custom
  nlf = !is.list(FUN)
  if (is.character(FUN)) { 
    FUN = unlist(strsplit(FUN,",",fixed = TRUE), use.names = FALSE)
    FUN = sapply(FUN, match.fun, descend = FALSE) 
  } else if (nlf) {
    nFUN = deparse(substitute(FUN)) 
    FUN = list(FUN)
    names(FUN) = nFUN
  }
  namFUN = names(FUN)
  lFUN = length(FUN)
  iFUN = seq_len(lFUN)
  nlcf = !is.list(catFUN)
  if (is.character(catFUN)) {
    catFUN = unlist(strsplit(catFUN,",",fixed = TRUE), use.names = FALSE)
    catFUN = sapply(catFUN, match.fun, descend = FALSE) 
  } else if (nlcf) { 
    nFUN = deparse(substitute(catFUN)) 
    catFUN = list(catFUN)
    names(catFUN) = nFUN
  }
  namcatFUN = names(catFUN)
  lcatFUN = length(catFUN)
  icatFUN = seq_len(lcatFUN)
  
  # Identifying the inputs: by
  if (!is.null(by)) {
    clby = class(by)
    if (is.list(by) || clby == "formula" || length(by) == NROW(X)) {
      if (is.list(by)) {
        if (any(lengths(by) != NROW(X))) stop("arguments must have same length")
        naml = names(by)
        if (is.null(naml)) { 
          names(by) = paste0("Group.",seq_along(by)) 
        } else {
          indl <- which(!nzchar(naml)) 
          names(by)[indl] <- paste0("Group.", indl)
        }
        X = cbind(by,X)
        by = names(by)
      } else if (!any(clby == c("matrix","formula"))) {
        namby = make.names(deparse(substitute(by)))
        X = data.frame(Groupby = by, X, stringsAsFactors = FALSE)
        by = names(X)[1] = namby 
      } else {
        if (length(by)>2) X = get_all_vars(by,X)
        by = attr(terms(by),"term.labels")      
      }
    }
    if (is.character(by)) { 
      by = unlist(strsplit(by,",",fixed = TRUE), use.names = FALSE)
      num = match(by,names(X))
    } else if (is.numeric(by)) { 
      num = by 
      by = names(X)[num]
    } else { 
      stop("'by' must be either a formula, a numeric or character vector signifying the columns over which to aggregate, or a vector, list or data.frame")
    }
    lby = length(by)
    iby = seq_len(lby)
  } else {
    by = "nullID"; num = 1; lby = 1; iby=1
    X = data.frame(nullID = 1, X, stringsAsFactors = FALSE)
  }
  
  # Characterizing the variables  
  cols = setdiff(seq_along(X),num)
  nu = setdiff(which(vapply(X, is.numeric, FUN.VALUE=logical(1))),num)
  if (!is.character(factors) || factors!="as.categorical") {
    fc = setdiff(which(vapply(X, is.factor, FUN.VALUE=logical(1))),num)
    ff = match.fun(factors)
    for (i in fc) X[[i]] = ff(X[[i]]) 
    nu = sort(c(nu,fc)) 
  }
  nnu = setdiff(cols,nu)
  
  # Preprocessing 
  cc = which(!complete.cases(X[num]))
  if (length(cc)) X = X[-cc, , drop = FALSE] 
  Xby = X[num]
  if (!collapse) ordr = do.call(order,Xby)   
  
  # Functions to perform the aggregation
  if (na.rm) {
    FUNwrap <- function(x,f, ...){x = x[!is.na(x)]; f(x, ...)}
  } else 
    FUNwrap <- function(x,f, ...)f(x, ...)
  if (data.table) { # data.table solution -> maximum speed
    nlfs = all(nlf,ifelse(dropcat,TRUE,nlcf))
    bstats = c("mean","median","Mode","sd","var","min","max","skewness","kurtosis")
    suppressPackageStartupMessages(library(data.table))
    agg <- function(df, by, FUN, nam, ...) { 
      if (na.rm && any(nam == bstats)) {
        eval(parse(text = paste0("df = setDT(df)[, lapply(.SD, ",nam,", na.rm = TRUE",ifelse(missing(...),"",paste0(", ",deparse(substitute(...)))),"), keyby = by]"))) 
      } else if (!na.rm && nlfs) { 
        eval(parse(text = paste0("df = setDT(df)[, lapply(.SD, ",nam,ifelse(missing(...),"",paste0(", ",deparse(substitute(...)))),"), keyby = by]"))) 
      } else {
        df = setDT(df)[, lapply(.SD, FUNwrap, f = FUN, ...), keyby = by]
      }
      setDF(df)   
      if (na.rm && replace.nan) df = quickdf(lapply(df,nantona))
      return(df)
    }
  } else { # base R solution (default)
    agg <- function(df, by, FUN, nam, ...) {
      grp <- lapply(by, ident) 
      names(grp) <- NULL
      grp <- do.call(paste, c(grp, list(sep = "."))) 
      y <- by[match(sort(unique(grp)), grp, 0L), , drop = FALSE] 
      nry <- NROW(y)
      z <- lapply(df, function(e) { 
        ans <- lapply(X = split(e, grp), FUN = FUNwrap, f = FUN, ...)
        unlist(ans, use.names = FALSE, recursive = FALSE) }) 
      if (na.rm && replace.nan) z = lapply(z, nantona) 
      len <- length(y)
      for (i in seq_along(z)) y[[len + i]] <- z[[i]] 
      names(y) <- c(names(by), names(df))
      row.names(y) <- NULL
      return(y)
    }  
  }
  
  # Exchanging arguments in special cases
  if (!dropcat) {
    if (length(nu)==0 && length(nnu)>0) {
      nu = nnu; dropcat = TRUE; FUN = catFUN; lFUN = lcatFUN; iFUN = icatFUN; namFUN = namcatFUN 
    } else if (lFUN<=1 && lcatFUN>1) {
      a = nu; nu = nnu; nnu = a; b = FUN; FUN = catFUN; catFUN = b
      a = lFUN; lFUN = lcatFUN; lcatFUN = a; b = namFUN; namFUN = namcatFUN; namcatFUN = b
      a = iFUN; iFUN = icatFUN; icatFUN = a
    } else if (lcatFUN>1) { 
      reshape.long = FALSE; as.list = FALSE
    }
  }
  
  # Computing output
  if (is.null(custom)) {
    if (parallel && lFUN>1) { # Numeric Variables
      suppressPackageStartupMessages(library(parallel))
      no_cores <- detectCores() - 1
      cl <- makeCluster(no_cores)
      clusterExport(cl,ifelse(data.table,"data.table","aggregate.data.frame"))
      res = parLapply(cl,iFUN,function(i) agg(df=X[nu],by=Xby,FUN[[i]],namFUN[i], ...)) 
      stopCluster(cl)
    } else if (lFUN>1) {
      res = lapply(iFUN,function(i) agg(df=X[nu],by=Xby,FUN[[i]],namFUN[i], ...))
    } else {
      res = agg(df=X[nu],by=Xby,FUN[[1]],namFUN, ...)
    }
    if (!dropcat && length(nnu)>0) { # Categorical Variables
      if (parallel && lcatFUN>1) {
        library(parallel) 
        no_cores <- detectCores() - 1
        cl <- makeCluster(no_cores)
        clusterExport(cl,ifelse(data.table,"data.table","aggregate.data.frame"))
        Catres = parLapply(cl,icatFUN,function(i) agg(df=X[nnu],by=Xby,catFUN[[i]],namcatFUN[i], ...)) 
        stopCluster(cl)
      } else if (lcatFUN>1) {
        Catres = lapply(icatFUN,function(i) agg(df=X[nnu],by=Xby,catFUN[[i]],namcatFUN[i], ...))
      } else {
        Catres = agg(df=X[nnu],by=Xby,catFUN[[1]],namcatFUN, ...)
      }
    }
  } else { # Custom Mode
    icf = is.character(unlist(FUN, use.names = FALSE)) 
    if (!nlf && !is.null(namFUN)) {
      indl = lapply(namFUN,function(x)match(unlist(strsplit(x,",",fixed = TRUE), use.names = FALSE), names(X))) 
      if (icf) {namFUN = unlist(FUN, use.names = FALSE); FUN = sapply(FUN, match.fun) }
      unf = namFUN 
    } else {
      ind = sort(union(nu,nnu)); fFUN = FUN; FUN = unique(FUN)
      unf = if (is.null(namFUN)) seq_along(FUN) else unique(namFUN) 
      indl = lapply(FUN,function(x)ind[sapply(fFUN,identical,x)])
      if (length(ind) != lFUN) stop("Vector of custom functions needs to match columns of data")
    }
    if ((nlf || icf) && custom.names) { res = lapply(seq_along(unf),function(i){r = agg(df=X[indl[[i]]],by=Xby,FUN[[i]],unf[i]) 
    names(r)[-iby] = paste0(names(r)[-iby],".",unf[i]); r})
    } else res = lapply(seq_along(unf),function(i) agg(df=X[indl[[i]]],by=Xby,FUN[[i]],unf[i]))
    if (as.list != "FUN") {
      res = data.frame(res[[1]][iby],cbind.data.frame(lapply(res,function(x)x[-iby]), stringsAsFactors = FALSE), stringsAsFactors = FALSE)
      if (!collapse) res = merge.data.frame(Xby,res)[order(ordr), , drop = FALSE] 
      rownames(res) = NULL
      if (sort) res = res[order(c(num,unlist(indl, use.names = FALSE)))]  
    } else names(res) = unf
    lFUN = 1
  }
  
  # Preparing Output
  if (is.null(custom)) {
    if (reshape.long || lFUN == 1 || as.list == "FUN") { # Anything but multiple functions columns in parallel
      if (lFUN>1) { # if multiple functions 
        if (!dropcat && length(nnu)>0) { # if categorical variables 
          if (as.list != "FUN" && show.statistic) { # show a statistic
            res = lapply(iFUN,function(x){
              data.frame(res[[x]][iby],Statistic=namFUN[x],res[[x]][-iby],Catres[-iby], stringsAsFactors = FALSE)})
            if (!collapse) res = lapply(res,function(x)merge.data.frame(Xby, x)[order(ordr), , drop = FALSE])
            if (sort) { ord = order(c(num,-Inf,nu,nnu))
            res = lapply(res,function(x)x[ord])}
          } else { # no statistic 
            res = lapply(iFUN,function(x){ # only if show.statistic argument
              data.frame(res[[x]],Catres[-iby], stringsAsFactors = FALSE)})
            if (!collapse) res = lapply(res,function(x)merge.data.frame(Xby, x)[order(ordr), , drop = FALSE])
            if (sort) { ord = order(c(num,nu,nnu))
            res = lapply(res,function(x)x[ord])}
          }
        } else { 
          if (as.list != "FUN" && show.statistic) { 
            res = lapply(iFUN, function(x){ # No categorical but with statistic
              data.frame(res[[x]][iby],Statistic=namFUN[x],res[[x]][-iby], stringsAsFactors = FALSE)})
            if (!collapse) res = lapply(res,function(x)merge.data.frame(Xby, x)[order(ordr), , drop = FALSE])
            if (sort) { ord = order(c(num,-Inf,nu))
            res = lapply(res,function(x)x[ord])}
          } else { # No categorical no statistic
            if (!collapse) res = lapply(res,function(x)merge.data.frame(Xby, x)[order(ordr), , drop = FALSE])
            if (sort) { ord = order(c(num,nu))
            res = lapply(res,function(x)x[ord])}
          }
        }
        names(res) = namFUN
        if (as.list != "FUN") { # Combine unless FUN
          res = if (show.statistic) { do.call(rbind.data.frame, c(res, make.row.names = FALSE, stringsAsFactors = FALSE))
          } else do.call(rbind.data.frame, c(res, stringsAsFactors = FALSE))  
        }
      } else if (!dropcat && length(nnu)>0) { # Only one function but with categorical variables
        res = data.frame(res,Catres[-iby], stringsAsFactors = FALSE)
        if (!collapse) res = merge.data.frame(Xby,res)[order(ordr), , drop = FALSE]
        if (sort) res = res[order(c(num,nu,nnu))]
      } else { # if one function and no categorical variables
        if (!collapse) res = merge.data.frame(Xby,res)[order(ordr), , drop = FALSE]
        if (sort) res = res[order(c(num,nu))]
      }
    } else { # If multiple functions all in a row
      res = lapply(iFUN,function(x) {
        names(res[[x]])[-iby] = paste0(names(res[[x]])[-iby],".",namFUN[x]); res[[x]]})
      res[-1] = lapply(iFUN[-1], function(x) res[[x]][-iby])
      if (!dropcat && length(nnu)>0) { # If categorical variables
        if (lcatFUN>1) { # Multiple categorical functions 
          Catres = lapply(icatFUN,function(x) {
            names(Catres[[x]])[-iby] = paste0(names(Catres[[x]])[-iby],".",namcatFUN[x]); Catres[[x]]})
          Catres[-1] = lapply(icatFUN[-1], function(x) Catres[[x]][-iby])
          Catres = do.call(cbind,Catres)
          nnu = rep(nnu,lcatFUN)
        }
        Catres = Catres[-iby]
        res = do.call(cbind,res)
        res = cbind.data.frame(res,Catres)
        if (!collapse) res = merge.data.frame(Xby,res)[order(ordr), , drop = FALSE]
        if (sort) res = res[order(c(num,rep(nu,lFUN),nnu))] 
      } else { # No categorical variables 
        res = do.call(cbind.data.frame,res)
        if (!collapse) res = merge.data.frame(Xby,res)[order(ordr), , drop = FALSE]
        if (sort) res = res[order(c(num,rep(nu,lFUN)))] 
      }
    }
  }
  if (as.list != "FUN") { 
    if (lby == 1 && by == "nullID") res = res[-1]
    if (as.list == "by" && by != "nullID") {
      num = match(by,names(res))
      res = split.data.frame(res[-num],res[num], drop = TRUE, lex.order = TRUE)
    } else if (dropby) {
      res = res[-match(by,names(res))]
    }
    if ((any(cl=="matrix") || (any(dim(res)==1) && (dropcat || length(nnu)<=0))) && (lFUN==1 || !reshape.long))   
      res = drop(as.matrix(res))
    if (isDT) setDT(res)
  } else {
    if (lby == 1 && by == "nullID") lapply(res,function(x)x[-1])
    if ((any(cl=="matrix") || (any(dim(res)==1) && (dropcat || length(nnu)<=0))) && (lFUN==1 || !reshape.long)) 
      res = lapply(res,function(x)drop(as.matrix(x)))
    if (isDT) lapply(res,setDT)
  }
  return(res)
}


# Qsu
qsu <- function(X, by = NULL, xt = NULL, FUN = NULL, Q = FALSE, Ext = FALSE, 
                trans = NULL, trans.by = FALSE, ndigits = 2, na.rm = TRUE, pretty = FALSE,
                labels = FALSE, factors = "as.categorical", combine.by = FALSE,
                combine.xt = TRUE, within.add.mean = TRUE, show.trans = TRUE,
                data.out = FALSE, data.out.drop = FALSE, xt.data.table = FALSE) {
  
  if (!all(class(X) == "data.frame")) X = as.data.frame(X, stringsAsFactors = FALSE) 
  
  # Find Labels
  if (labels != FALSE) {
    if (is.character(labels)) {
      labs = labels
    } else {
      labs = lapply(X,attr,"label")
      labs[lengths(labs) == 0] = NA
      labs = unlist(labs, use.names = FALSE, recursive = FALSE) 
    } 
  } else {
    labs = NULL
  }
  
  # Identifying the inputs: by
  if (!is.null(by)) {
    clby = class(by)
    if (is.list(by) || clby == "formula" || length(by) == NROW(X)) {
      if (is.list(by)) {
        if (any(lengths(by) != NROW(X))) stop("arguments must have same length")
        namby = names(by)
        if (is.null(namby)) { 
          names(by) = paste0("Groupby.",seq_along(by)) 
        } else {
          indl <- which(!nzchar(namby)) 
          names(by)[indl] <- paste0("Groupby.", indl)
        }
        X = cbind(by,X)
        by = names(by)
      } else if (!any(clby == c("matrix","formula"))) { 
        namby = make.names(deparse(substitute(by)))
        X = data.frame(Groupby.1 = by, X, stringsAsFactors = FALSE)
        by = names(X)[1] = namby 
      } else {
        if (length(by)>2 || length(xt)>2) 
          X = X[unique(unlist(strsplit(c(as.character(by)[-1],as.character(xt)[-1])," + ", fixed = TRUE), use.names = FALSE))] 
        by = attr(terms(by),"term.labels") 
      }
    }
    if (is.character(by)) { 
      by = unlist(strsplit(by,",",fixed = TRUE), use.names = FALSE)
      numby = match(by,names(X))
    } else if (is.numeric(by)) { 
      numby = by 
      by = names(X)[numby]
    } else { 
      stop("'by' must be either a formula, a numeric or character vector signifying the columns over which to summarize, or a vector, list or data.frame")    
    }
  } else {
    numby = 0
  }
  
  # Identifying the inputs: xt
  if (!is.null(xt)) {
    clxt = class(xt)
    if (is.list(xt) || clxt == "formula" || length(xt) == NROW(X)) {
      if (is.list(xt)) {
        if (any(lengths(xt) != NROW(X))) stop("arguments must have same length")
        namxt = names(xt)
        if (is.null(namxt)) { 
          names(xt) = paste0("Groupxt.",seq_along(xt)) 
        } else {
          indl <- which(!nzchar(namxt)) 
          names(xt)[indl] <- paste0("Groupxt.", indl)
        }
        X = cbind(X,xt)
        xt = names(xt)
      } else if (!any(clxt == c("matrix","formula"))) { 
        X = data.frame(X, Groupxt.1 = xt, stringsAsFactors = FALSE)
        xt = "Groupxt.1"
      } else {
        if (length(xt)>2 && is.null(by)) X = get_all_vars(xt,X)
        xt = attr(terms(xt),"term.labels") 
      }
    }
    if (is.character(xt)) { 
      xt = unlist(strsplit(xt,",",fixed = TRUE), use.names = FALSE)
      numxt = match(xt,names(X))
    } else if (is.numeric(xt)) { 
      numxt = xt 
      xt = names(X)[numxt]
    } else { 
      stop("'xt' must be either a formula, a numeric or character vector signifying the columns over which to summarize, or a vector, list or data.frame")    
    }
  } else {
    numxt = 0
  }
  
  # Some functions:
  quickdf <- function(l) { 
    class(l) <- "data.frame"
    attr(l, "row.names") <- .set_row_names(length(l[[1]]))
    l
  }
  
  # Character variables:
  cols = setdiff(seq_along(X),c(numby,numxt)) 
  nu = vapply(X, is.numeric, FUN.VALUE=logical(1))
  if (!is.character(factors) || factors!="as.categorical") {
    fc = vapply(X, is.factor, FUN.VALUE=logical(1))
    fcn = which(fc)
    ff = match.fun(factors)
    for (i in fcn) X[[i]] = ff(X[[i]]) 
    nu = nu+fc==1 }
  nu = intersect(which(nu),cols)
  nnu = setdiff(cols,nu)
  
  # Preprocessing DATA: 
  if (is.character(trans)) trans = match.fun(trans)
  transfun <- function(D) {
    if (na.rm) {
      quickdf(lapply(D,function(x){nna = which(!is.na(x)); x[nna] = trans(x[nna]); x}))
    } else 
      quickdf(lapply(D,trans))
  }
  if (!is.null(trans) && !trans.by) X[nu] = transfun(X[nu])
  
  # Preprocessing vector of functions:
  if (is.null(FUN)) {
    fun = list(N=length,D=function(x)length(unique(x)),Mean=mean,SD=sd,Min=min,Max=max)
    if (Q) fun = c(fun[-(5:6)],myQ=function(x)quantile(x, probs = seq(0,1,1/Q)))
    if (Ext) 
      if (Q) { fun = c(fun,Skew=function(x){n <- length(x); (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)},Kurt=function(x){n <- length(x); n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)})
      } else fun = c(fun[1:3],Median=median,fun[-(1:3)],Skew=function(x){n <- length(x); (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)},Kurt=function(x){n <- length(x); n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)})
  } else {
    if (is.character(FUN)) { fun = Reduce(c,strsplit(FUN,",",fixed = TRUE))
    fun = sapply(fun, match.fun, descend = FALSE) 
    } else if (!is.list(FUN)) {
      fun = list(FUN) 
      if (!is.null(xt) && is.null(names(fun))) names(fun) = ""
    } else {
      fun = FUN
      if (!is.null(xt) && is.null(names(fun))) names(fun) = ""
    }
  }
  
  # Functions to compute the summary:
  sumnum<-function(x){
    if (na.rm) x<-x[!is.na(x)]
    unlist(lapply(fun,function(j)j(x)), recursive = FALSE)
  } 
  sumchr<-function(x, l){
    if (na.rm) x<-x[!is.na(x)]
    c(length(x),length(unique(x)),rep(NA,l-2))
  } 
  su <- function(X, cust = !is.null(FUN)) {
    if (!data.out) {
      if (!cust) {
        resnu = t(sapply(X[nu],sumnum))
        if (length(nnu)) {
          sn1 = sumnum(1)
          resch = t(sapply(X[nnu],sumchr,l=length(sn1)))
          colnames(resch) = names(sn1)
        } else { resch = NULL }
        res = rbind.data.frame(resnu,resch)[order(c(nu,nnu)),]
        tryCatch({colnames(res) = gsub("myQ.","",colnames(res))},error = function(e)cat(""))
        if (nrow(res)==1) {if (is.null(labs) && is.null(xt) && combine.by==FALSE) res = drop(as.matrix(res)) else rownames(res) = ""} 
      } else {
        if (length(fun)==1 && (!is.null(xt) || !is.null(labs))) { 
          res = do.call(rbind,lapply(X[cols],sumnum))
        } else if (combine.by == TRUE || !is.null(xt)) { 
          res = as.data.frame(t(sapply(X[cols],sumnum)), stringsAsFactors = FALSE)
        } else {
          res = drop(t(sapply(X[cols],sumnum)))
          if (is.vector(res) && length(res)==1) names(res) = NULL
        }
      }
      tryCatch({res = round(res,ndigits)},error = function(e)cat(""))
      if (pretty==TRUE) {
        NAres = is.na(res)
        res = format(res, drop0trailing = TRUE, scientific = FALSE) 
        res[NAres] = "-"
      }
      return(res)
    } else {
      if (data.out.drop) X = drop(X[sort(c(nu,nnu))])
      return(X)
    }
  }
  xtsu <- function(X) {
    Xxt = X[numxt] 
    ord = do.call(order,Xxt) 
    Xxt = Xxt[ord, , drop = FALSE]
    Be = collap(X, numxt, data.table = xt.data.table)
    W1 = merge.data.frame(Xxt,Be, sort = FALSE)[order(ord), names(Be), drop = FALSE]
    names(Be) = names(X) 
    W = X
    if (within.add.mean==TRUE) { # The overall mean is added back in STATA
      W[nu] = quickdf(lapply(nu,function(i){W[[i]] - W1[[i]] + mean(W[[i]], na.rm = TRUE)}))    
    } else {
      W[nu] = W[nu] - W1[nu]
    }
    if (data.out==FALSE) {
      trans = c("overall","between","within")
      # Normal
      All = su(X=X)
      # Between    
      B = su(X=Be)
      # Within
      W = su(X=W)
      if (is.null(FUN)) { 
        if (pretty) { W[,1] = as.character(round(as.numeric(All[,1])/as.numeric(B[,1])))
        } else W[,1] = round(All[,1]/B[,1],2)
        if (combine.xt==FALSE) {
          names(W)[1] = "T"
        } else names(All)[1] = names(B)[1] = names(W)[1] = "N/T"
      }
      if (combine.xt==FALSE) {
        res = list(overall = drop(All), between = drop(B), within = drop(W))
      } else if (length(sumnum(1))==1) { 
        rnres = if (names(fun)=="") trans  else paste0(trans,".",names(fun))
        res = cbind.data.frame(All[,1],B[,1],W[,1]); colnames(res) = rnres 
        if (nrow(res)==1 || combine.by==FALSE) res = drop(as.matrix(res))
      } else {
        colnam = colnames(All)
        rownames(B) = paste0(rownames(B),".B")
        rownames(W) = paste0(rownames(W),".W")
        res = rbind.data.frame(All,B,W, stringsAsFactors = FALSE) 
        res = res[order(rep(seq_len(nrow(All)),3)),]
        if (show.trans) {
          res = cbind.data.frame(Trans = rep(trans,nrow(All)), res, stringsAsFactors = FALSE)
          colnames(res)[-1] = colnam
        } else colnames(res) = colnam      
      }
    } else {
      if (data.out.drop) { 
        sortc = sort(c(nu,nnu))
        res = list(overall = drop(X[sortc]), between = drop(Be[sortc]), within = drop(W[sortc]))
      } else 
        res = list(overall = X, between = Be, within = W)
    }
    return(res)
  }
  # Computing the Summary:
  if (is.null(by) && is.null(xt)) { # Normal Summary
    suppressWarnings({res = su(X=X)})
    if (!is.null(labs)) {res = as.data.frame(res, stringsAsFactors = FALSE); res$Label = labs[cols]}
  } else if (is.null(xt)) { # Summary by
    sumfun = ifelse(!is.null(trans) && trans.by, function(x){x[nu] = transfun(x[nu]); su(x)}, su)
    if (combine.by) {
      X = split.data.frame(X,X[numby], drop = TRUE, lex.order = TRUE) 
      suppressWarnings({res = lapply(X, sumfun)})
    } else {
      suppressWarnings({res = by(X, X[numby], sumfun, simplify = FALSE)})
    }
    if (!is.null(labs)) res = lapply(res,function(x){x = as.data.frame(x, stringsAsFactors = FALSE); x$Label = labs[cols]; x}) 
    if (combine.by) {
      res = do.call(rbind,res) 
      if (ncol(res)==1) res = drop(as.matrix(res)) 
    }
  } else if (is.null(by)) { # Panel Summary
    suppressWarnings({res = xtsu(X)})
    if (!is.null(labs)) res$Label = unlist(lapply(labs[cols],function(x)c(x,rep("",2))), use.names = FALSE, recursive = FALSE)
  } else { # Panel Summary by
    sumfun = ifelse(!is.null(trans) && trans.by, function(x){x[nu] = transfun(x[nu]); xtsu(x)}, xtsu)
    if (combine.by) {
      X = split.data.frame(X,X[numby], drop = TRUE, lex.order = TRUE) 
      suppressWarnings({res = lapply(X, sumfun)})
    } else {
      suppressWarnings({res = by(X, X[numby], sumfun, simplify = FALSE)})
    }
    if (!is.null(labs)) res = lapply(res,function(x){x$Label = unlist(lapply(labs[cols],function(z)c(z,rep("",2))), use.names = FALSE, recursive = FALSE); x}) 
    if (combine.by) res = do.call(rbind,res) 
  }
  return(res) 
}
