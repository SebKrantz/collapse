#' Aggregation of multi-type and multi-level data
#'
#' The function \code{collap} allows you to aggregate datasets with
#' multiple types (numeric, factor) in a single function call.
#'
#' @param X A vector, matrix, list, data.frame or data.table to aggregate (anything that can be coerced to data.frame)
#' @param by Columns to aggregate by, \bold{either contained in X} and indicated using a one-or two sided formula (two-sided if only certain columns in X are to be aggregated), column indices, a vector of column names, or a string of comma-separated column names, \bold{or externally supplied} in form of a vector, list of vectors or data.frame, with the number of elements/rows matching that of X. If 'by' is left empty, columns are fully aggregated.
#' @param FUN Function(s) to apply to numeric columns in X, defaults to the mean. A single function can be supplied without quotes. Multiple functions can be supplied as a character vector, string of comma-separated function names, or as a list of functions (preferably named). Ad-hoc functions can be supplied.
#' @param catFUN Function(s) to apply to categorical columns in X, defaults to the Mode. If all elements in a group defined by 'by' are distinct, the Mode defaults to the first element. Multiple functions can be supplied in the same manor as to 'FUN'.
#' @param factors Specifies treatment of factor variables. Default is treatment as categorical variables. Alternatively factors can be coerced to numerical variables by spcifying "as.numeric", or the factor levels can be extracted and coerced to a numerical variable by specifying "as.numeric.fractor" (internally defined as:
#' \code{as.numeric.factor <- function(x) \{as.numeric(levels(x))[x]\})}
#' @param custom Option to supply a custom vector or list of functions whose length must match the number of columns to be aggregated. Alternatively a named list can be provided with the names being the comma-separated names of the columns to be aggregated by different functions, i.e. list("var1,var2,var3" = mean, var4 = median, "var7,var8" = sd).
#' @param custom.names Interact the column names with the respective function names in 'custom'.
#' @param na.rm Removes missing values from all columns before applying any functions. This is done internally in \emph{collap}, thus it is not required for functions in 'FUN' or 'catFUN' to have a 'na.rm' argument.
#' @param replace.nan Replaces NaN values with NA values. NaN's are frequently generated if na.rm = TRUE, and aggregation takes place over an empty subset.
#' @param sort Sort restores the columns back to their original order after aggregation. If sort = FALSE, the dataset is returned with the 'by' columns in front, and the other columns following in the order of computation (first numeric columns and then categorical columns, or columns in the order they are passed to 'custom').
#' @param collapse If collapse = FALSE, the aggregated data will be matched with the original data in the 'by' argument and \emph{collap} will return a dataset that is aggregated but of the same dimensions and row-order as the original data, i.e. a between-transformed dataset.
#' @param reshape.long If multiple functions are supplied to either 'FUN' or 'catFUN', by default \emph{collap} returns a wider dataset. If reshape.long = TRUE, then a long form of the dataset is returned with an additional column 'Statistic' indicating the function used for aggregation.
#' @param show.statistic If multiple functions are called and reshape.long = TRUE, show.statistic = FALSE can be called to omit the 'Statistic' column and instead make appropriate row.names.
#' @param as.list Optionally the output can be requested as a list of vectors or data.frames. There are two options here: If as.list = "by", then a list will be returned whose elements are the aggregated output for each group in 'by'. If multiple functions are supplied to either 'FUN' or 'catFUN', calling as.list = "FUN" will return a list with the dataset aggregated by the different functions. as.list = "by" may come at some slight extra computational cost but as.list = "FUN" does not.
#' @param dropcat Drop all categorical variables apart from identifiers in 'by' (i.e. don't perform aggregation on them).
#' @param dropby Drop the columns in 'by' from the final output.
#' @param data.table By default \emph{collap} is built as a wrapper around \emph{aggregate.data.frame}. Calling this argument will internally use \emph{data.table} as workhorse function, yielding significant speed improvements for large datasets. Requires \code{data.table} package to be installed.
#' @param parallel If multiple functions are supplied to 'FUN' or 'catFUN', parallel = TRUE will automatically parallelize computation on $k-1$ of the available cores (using the \emph{parLapply} function from the \emph{parallel} package). The argument works together with data.table = TRUE to guarantee maximum performance on tasks involving large datsets and multiple functions.
#' @param ... Additional arguments supplied to 'FUN', 'catFUN' or to \emph{aggregate.data.frame} in the default mode.
#'
#' @export
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

  if (data.table && !requireNamespace("data.table", quietly=TRUE))
    stop("data.table package is not installed, install data.table or use `data.table=FALSE`")

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
    .SD = NULL
    agg <- function(df, by, FUN, nam, ...) {
      narmcalls = na.rm && any(nam == bstats)
      nonarmcalls = !na.rm && nlfs
      if (narmcalls || nonarmcalls) {
        if (narmcalls) {
          old = parse(text = paste0("df_old = data.table::setDT(df)[, lapply(.SD, ",nam,", na.rm = TRUE",ifelse(missing(...),"",paste0(", ",deparse(substitute(...)))),"), keyby = by]")) ## this line can be removed after fully moving from eval-parse to eval-subsitute
        } else {
          old = parse(text = paste0("df_old = data.table::setDT(df)[, lapply(.SD, ",nam,ifelse(missing(...),"",paste0(", ",deparse(substitute(...)))),"), keyby = by]"))
        }
        new = substitute(df <- data.table::setDT(df)[, .lapplyCall, keyby = by],
                         list(.nam=as.name(nam), .dots=ifelse(missing(...), substitute(), substitute(...)),
                              .lapplyCall = as.call(c(
                                list(as.name("lapply"), as.name(".SD"), as.name(nam)),
                                if (narmcalls) list(na.rm = TRUE) else list(), # here we handle narmcalls or nonarmcalls
                                if (missing(...)) list() else as.list(substitute(...)) # here we handle dots
                              ))))
        eval(old) ## this line can be removed after fully moving from eval-parse to eval-subsitute
        eval(new)
        if(!isTRUE(all.equal(df, df_old))) browser() ## this line can be removed after fully moving from eval-parse to eval-subsitute
        rm(df_old) ## this line can be removed after fully moving from eval-parse to eval-subsitute
      } else {
        df = setDT(df)[, lapply(.SD, FUNwrap, f = FUN, ...), keyby = by]
      }
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
