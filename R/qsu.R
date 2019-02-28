#' Quick summary for multi-type and multi-level data
#'
#' This function allows you to create quick summaries for data.
#'
#' @param X A vector, matrix, data.frame or data.table to summarize (anything that can be coerced to data.frame)
#' @param by Groups to summarize by, \bold{either contained in X} and indicated using a one-or two sided formula (two-sided if only certain columns in X are to be aggregated), column indices, a vector of column names, or a string of comma-separated column names, \bold{or externally supplied} in form of a vector, list of vectors or data.frame, with the number of elements/rows matching that of X.
#' @param xt Groups to compute statistics overall, between and within. The same flexibility as with the 'by' argument applies. If used together with 'by', a subgroup of 'by' should be used. If a two-sided formula is used together with 'by', it does not matter whether the LHS variables are specified in the 'by', 'xt' or in both arguments.
#' @param FUN Custom function(s) to apply to all columns in X apart from columns in the 'by' or 'xt' arguments. Functions must take a vector and return a vector of statistics. A single function can be supplied without quotes. Multiple functions can be supplied as a character vector, string of comma-separated function names, or as a named list of functions. Ad-hoc functions can be supplied. 'FUN' when it is used overrides the default set of statistics and the 'Q' and 'Ext' arguments.
#' @param Q Number of quantiles to compute.
#' @param Ext Request an Extended set of statistics including the \emph{median}, the \emph{skewness} and the \emph{kurtosis}
#' @param trans A transformation function applied to the numeric columns of the data (for example \emph{log}, \emph{scale}, \emph{diff} or growth rates)
#' @param trans.by If the 'by' option is used, 'trans' can be applied to groups separately (i.e. one could use it to obtain growth rates for multiple countries in a long country-time $\times$ variables dataset)
#' @param ndigits Number of digits to show. If set to NULL, all digits will be shown.
#' @param na.rm Internally removes missing values before applying any functions or transformations. It is not required for functions to have a 'na.rm' argument.
#' @param pretty Returns result as a character matrix where trailing zeros are eliminated and large numbers are written in standard (as opposed to scientific) notation.
#' @param labels Show variable labels next to statistics. If labels = TRUE, X must be a data.frame with variable labels stored as attributes [attr(X\$var1,"label")<-"label1"] etc. Alternatively, a character vector of labels of length ncol(X) can be passed to the labels argument.
#' @param factors Specifies the treatment of factor variables. Default is treatment as categorical variables. Alternatively factors can be coerced to numerical variables by spcifying "as.numeric", or the factor levels can be extracted and coerced to a numerical variable by specifying "as.numeric.fractor" (internally defined as: as.numeric.factor <- function(x) \{as.numeric(levels(x))[x]\})
#' @param combine.by If the 'by' option is used, combine.by = TRUE gives a compact output instead of a list.
#' @param combine.xt If the 'xt' option is used combine.xt = FALSE returns a list with overall, between group and within group statistics.
#' @param within.add.mean By default, within-group statistics are computed as $\bold{x}_{it}-\bar{\bold{x}}_i+\bar{\bar{\bold{x}}}$. If within.add.mean = FALSE, The within-transformed dataset is obtained as $\bold{x}_{it}-\bar{\bold{x}}_i$, which is a more classical within-transformation used i.e. for fixed-effects regression.
#' @param data.out Output transformed data used to compute the summary. If the 'xt' option is used, the output will be a named list of three datasets: An overall dataset (= the original dataset if trans = NULL), an aggregated dataset for the between-statistics, and a within-transformed dataset. All datasets come with the original column order, the aggregated dataset is sorted by the 'xt' identifiers, and the within-transformed dataset has the same row-order as the original dataset. In the aggregated dataset categorical variables were aggregated using the mode, while in the within-transformed dataset categorical variables are unaffected/untransformed.
#' @param data.out.drop Drop all identifiers supplied to 'by' or 'xt' before returning the dataset.
#' @param xt.data.table If the 'xt' option is used, \emph{qsu} internally utilizes \emph{collap} to aggregate the data and compute the within-transformed dataset. If xt.data.table = TRUE, \emph{collap} will internally use \emph{data.table}, yielding a much faster computation on large datasets.
#'
#' @export
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
