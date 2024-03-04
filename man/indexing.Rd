\name{indexing}
\alias{indexing}
\alias{findex_by}
\alias{iby}
\alias{findex}
\alias{ix}
\alias{unindex}
\alias{reindex}
\alias{is_irregular}
\alias{to_plm}
\alias{[.indexed_series}
\alias{[.indexed_frame}
\alias{$.indexed_frame}
\alias{[[.indexed_frame}
\alias{[.index_df}
\alias{print.index_df}
\alias{[<-.indexed_frame}
\alias{$<-.indexed_frame}
\alias{[[<-.indexed_frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Indexed Time Series and Panels
}
\description{
A fast and flexible indexed time series and panel data class that inherits from \emph{plm}'s 'pseries' and 'pdata.frame', but is more rigorous, natively handles irregularity, can be superimposed on any data.frame/list, matrix or vector, and supports ad-hoc computations inside data masking functions and model formulas.
}

\usage{
## Create an 'indexed_frame' containing 'indexed_series'
findex_by(.X, \dots, single = "auto", interact.ids = TRUE)
iby(.X, \dots, single = "auto", interact.ids = TRUE)  # Shorthand

## Retrieve the index ('index_df') from an 'indexed_frame' or 'indexed_series'
findex(x)
ix(x)     # Shorthand

## Remove index from 'indexed_frame' or 'indexed_series' (i.e. get .X back)
unindex(x)

## Reindex 'indexed_frame' or 'indexed_series' (or index vectors / matrices)
reindex(x, index = findex(x), single = "auto")

## Check if 'indexed_frame', 'indexed_series', index or time vector is irregular
is_irregular(x, any_id = TRUE)

## Convert 'indexed_frame'/'indexed_series' to normal 'pdata.frame'/'pseries'
to_plm(x, row.names = FALSE)

# Subsetting & replacement methods: [(<-) methods call NextMethod().
# Also methods for fsubset, funique and roworder(v), na_omit (internal).

\method{[}{indexed_series}(x, i, \dots, drop.index.levels = "id")

\method{[}{indexed_frame}(x, i, \dots, drop.index.levels = "id")

\method{[}{indexed_frame}(x, i, j) <- value

\method{$}{indexed_frame}(x, name)

\method{$}{indexed_frame}(x, name) <- value

\method{[[}{indexed_frame}(x, i, \dots)

\method{[[}{indexed_frame}(x, i) <- value

# Index subsetting and printing: optimized using ss()

\method{[}{index_df}(x, i, j, drop = FALSE, drop.index.levels = "id")

\method{print}{index_df}(x, topn = 5, ...)

}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.X}{a data frame or list-like object of equal-length columns.}
  \item{x}{an 'indexed_frame' or 'indexed_series'. \code{findex} also works with 'pseries' and 'pdata.frame's created with \emph{plm}. For \code{is_irregular} \code{x} can also be an index (inherits 'pindex') or a vector representing time. }
  \item{\dots}{for \code{findex_by}: variables identifying the individual (id) and/or time dimensions of the data. Passed either as unquoted comma-separated column names or (tagged) expressions involving columns, or as a vector of column names, indices, a logical vector or a selector function. The time variable must enter last. See Examples. Otherwise: further arguments passed to \code{\link[=NextMethod]{NextMethod()}}.}
  \item{single}{character. If only one indexing variable is supplied, this can be declared as \code{"id"} or \code{"time"} variable. \code{"auto"} chooses \code{"id"} if the variable has \code{\link{anyDuplicated}} values.  }
  \item{interact.ids}{logical. If \code{n > 2} indexing variables are passed, \code{TRUE} calls \code{\link{finteraction}} on the first \code{n-1} of them (\code{n}'th variable must be time). \code{FALSE} keeps all variables in the index. The latter slows down computations of lags / differences etc. because ad-hoc interactions need to be computed, but gives more flexibility for scaling / centering / summarising over different data dimensions.  }
  \item{index}{and index (inherits 'pindex'), or an atomic vector or list of factors matching the data dimensions. Atomic vectors or lists with 1 factor will must be declared, see \code{single}. Atomic vectors will additionally be grouped / turned into time-factors. See Details. }
  \item{drop.index.levels}{character. Subset methods also subset the index (= a data.frame of factors), and this argument regulates which factor levels should be dropped: either \code{"all"}, \code{"id"}, \code{"time"} or \code{"none"}. The default \code{"id"} only drops levels from id's. \code{"all"} or \code{"time"} should be used with caution because time-factors may contain levels for missing time periods (gaps in irregular sequences, or periods within a sequence removed through subsetting), and dropping those levels would create a variable that is ordinal but no longer represents time. The benefit of dropping levels is that it can speed-up subsequent computations by reducing the size of intermediate vectors created in C++.   }
  \item{any_id}{logical. For panel series: \code{FALSE} returns the irregularity check performed for each id, \code{TRUE} calls \code{\link{any}} on those checks.}
  \item{row.names}{logical. \code{TRUE} creates descriptive row-names (or names for pseries) as in \code{plm}. This can be expensive and is usually not required for \code{plm} models to work.}
  \item{topn}{integer. The number of first and last rows to print.}
  \item{i, j, name, drop, value}{Arguments passed to \code{\link{NextMethod}}, or as in the \link[=[.data.frame]{data.frame methods}. Note that for index subsetting to work, \code{i} needs to be integer or logical (or an expression evaluation to integer or logical if \code{x} is a \emph{data.table}).}
}
\details{
The first thing to note about these new 'indexed_frame', 'indexed_series' and 'index_df' classes is that they inherit \emph{plm}'s 'pdata.frame', 'pseries' and 'pindex' classes, respectively. They add, improve, and, in some cases, remove functionality offered by \emph{plm}, with the aim of striking an optimal balance of flexibility and performance. The inheritance means that all 'pseries' and 'pdata.frame' methods in \emph{collapse}, and also some methods in \emph{plm}, apply to them. Where compatibility or performance considerations allow for it, \emph{collapse} will continue to create methods for \emph{plm}'s classes instead of the new classes.

The use of these classes does not require much knowledge of \emph{plm}, but as a basic background: A 'pdata.frame' is a data.frame with an index attribute: a data.frame of 2 factors identifying the individual and time-dimension of the data. When pulling a variable out of the pdata.frame using a method like \code{$.pdata.frame} or \code{[[.pdata.frame} (defined in \emph{plm}), a 'pseries' is created by transferring the index attribute to the vector. Methods defined for functions like \code{\link{lag}} / \code{\link{flag}} etc. use the index for correct computations on this panel data, also inside \emph{plm}'s estimation commands.

\bold{Main Features and Enhancements}

The 'indexed_frame' and 'indexed_series' classes extend and enhance 'pdata.frame' and 'pseries' in a number of critical dimensions. Most notably they:

\itemize{
  \item Support \bold{both time series and panel data}, by allowing indexation of data with one, two or more variables.
  \item Are \bold{class-agnostic}: any data.frame/list (such as data.table, tibble, tsibble, sf etc.) can become an 'indexed_frame' and continue to function as usual for most use cases. Similarly, any vector or matrix (such as ts, mts, xts) can become an 'indexed_series'. This also allows for transient workflows e.g. \code{some_df |> findex_by(...) |> 'do something using collapse functions' |> unindex() |> 'continue working with some_df'}.
  \item Have a comprehensive and efficient set of \bold{methods for subsetting and manipulation}, including methods for \code{\link{fsubset}}, \code{\link{funique}}, \code{\link[=roworder]{roworder(v)}} (internal) and \code{\link{na_omit}} (internal, \code{\link{na.omit}} also works but is slower). It is also possible to group indexed data with \code{\link{fgroup_by}} for transformations e.g. using \code{\link{fmutate}}, but aggregation requires \code{unindex()}ing.
  \item \bold{Natively handle irregularity}: time objects (such as 'Date', 'POSIXct' etc.) are passed to \code{\link{timeid}}, which efficiently determines the temporal structure by finding the greatest common divisor (GCD), and creates a time-factor with levels corresponding to a complete time-sequence. The latter is also done with plain numeric vectors, which are assumed to represent unit time steps (GDC = 1) and coerced to integer (but can also be passed through \code{\link{timeid}} if non-unitary). Character time variables are converted to factor, which might also capture irregular gaps in panel series. Using this time-factor in the index, \emph{collapse}'s functions efficiently perform correct computations on irregular sequences and panels without the need to 'expand' the data / fill gaps. \code{is_irregular} can be used to check for irregularity in the entire sequence / panel or separately for each individual in panel data.
\item Support computations inside \bold{data-masking functions and formulas}, by virtue of "\bold{deep indexation}": Each variable inside an 'indexed_frame' is an 'indexed_series' which contains in its 'index_df' attribute an external pointer to the 'index_df' attribute of the frame. Functions operating on 'indexed_series' stored inside the frame (such as \code{with(data, flag(column))}) can fetch the index from this pointer. This allows worry-free application inside arbitrary data masking environments (\code{with}, \code{\%$\%}, \code{attach}, etc..) and estimation commands (\code{glm}, \code{feols}, \code{lmrob} etc..) without duplication of the index in memory. A limitation is that external pointers are only valid during the present R session, thus when saving an 'indexed_frame' and loading it again, you need to call \code{data = reindex(data)} before computing on it.
}

Indexed series also have simple \link[base]{Math} and \link[base]{Ops} methods, which apply the operation to the unindexed series and shallow copy the attributes of the original object to the result, unless the result it is a logical vector (from operations like \code{!}, \code{==} etc.). For \link[base]{Ops} methods, if the LHS object is an 'indexed_series' its attributes are taken, otherwise the attributes of the RHS object are taken.

\bold{Limits to plm Compatibility}

In contrast to 'pseries' and 'pdata.frame's, 'indexed_series' and 'indexed_frames' do not have descriptive "names" or "row.names" attributes attached to them, mainly for efficiency reasons.

Furthermore, the index is stored in an attribute named 'index_df' (same as the class name), not 'index' as in \emph{plm}, mainly to make these classes work with \emph{data.table}, \emph{tsibble} and \emph{xts}, which also utilize 'index' attributes. This for the most part poses no problem to plm compatibility because plm source code fetches the index using \code{attr(x, "index")}, and \code{\link{attr}} by default performs partial matching. %It however allows plm objects to be indexed again / doubly indexed with both 'index' and 'index_df' attributes, so care needs to be taken when working with \emph{plm}.

A much greater obstacle in working with \emph{plm} is that some internal \emph{plm} code is hinged on there being no \code{[.pseries} method, and the existence of \code{[.indexed_series} limits the use of these classes in most \emph{plm} estimation commands. Therefore the \code{to_plm} function is provided to efficiently coerce the classes to ordinary plm objects before estimation. See Examples.

Overall these classes don't really benefit \emph{plm}, especially given that collapse's plm methods also support native plm objects. However, they work very well inside other models and software, including \emph{stats} models, \emph{fixest} / \emph{lfe}, and a whole bunch of time series and ML models. See Examples.

\bold{Performance Considerations}

When indexing long time-series or panels with a single variable, setting \code{single = "id" or "time"} avoids a potentially expensive call to \code{\link{anyDuplicated}}. Note also that when panel-data are regular and sorted, omitting the time variable in the index can bring >= 2x performance improvements in operations like lagging and differencing (alternatively use \code{shift = "row"} argument to \code{\link{flag}}, \code{\link{fdiff}} etc.) .

When dealing with long Date or POSIXct time sequences, it may also be that the internal processing by \code{\link{timeid}} is slow simply because calling \code{\link{strftime}} on these sequences to create factor levels is slow. In this case you may choose to generate an index factor with integer levels by passing \code{timeid(t)} to \code{findex_by} or \code{reindex} (which by default generates a 'qG' object which is internally converted to factor using \code{as_factor_qG}. The lazy evaluation of expressions like \code{as.character(seq_len(nlev))} in modern R makes this extremely efficient).

With multiple id variables e.g. \code{findex_by(data, id1, id2, id3, time)}, the default call to \code{finteraction()} can be expensive because of pasting the levels together. In this case, users may gain performance by manually invoking \code{finteraction()} (or its shorthand \code{itn()}) with argument \code{factor = FALSE} e.g. \code{findex_by(data, ids = itn(id1, id2, id3, factor = FALSE), time)}. This will generate a factor with integer levels instead.



\bold{Print Method}

The print methods for 'indexed_frame' and 'indexed_series' first call \code{print(unindex(x), ...)}, followed by the index variables with the number of categories (index factor levels) in square brackets. If the time factor contains unused levels (= irregularity in the sequence), the square brackets indicate the number of used levels (periods), followed by the total number of levels (periods in the sequence) in parentheses.

}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{timeid}},
\link[=time-series-panel-series]{Time Series and Panel Series}, \link[=collapse-documentation]{Collapse Overview}}
\examples{
oldopts <- options(max.print = 70)
# Indexing panel data ----------------------------------------------------------

wldi <- findex_by(wlddev, iso3c, year)
wldi
wldi[1:100,1]                 # Works like a data frame
POP <- wldi$POP               # indexed_series
qsu(POP)                      # Summary statistics
G(POP)                        # Population growth
STD(G(POP, c(1, 10)))         # Within-standardized 1 and 10-year growth rates
psmat(POP)                    # Panel-Series Matrix
plot(psmat(log10(POP)))

POP[30:5000]                  # Subsetting indexed_series
Dlog(POP[30:5000])            # Log-difference of subset
psacf(identity(POP[30:5000])) # ACF of subset
L(Dlog(POP[30:5000], c(1, 10)), -1:1) # Multiple computations on subset

# Fast Statistical Functions don't have dedicated methods
# Thus for aggregation we need to unindex beforehand ...
fmean(unindex(POP))
wldi |> unindex() |>
  fgroup_by(iso3c) |> num_vars() |> fmean()

library(magrittr)
# ... or unindex after taking group identifiers from the index
fmean(unindex(fgrowth(POP)), ix(POP)$iso3c)
wldi |> num_vars() \%>\%
  fgroup_by(iso3c = ix(.)$iso3c) |>
  unindex() |> fmean()

# With matrix methods it is easier as most attributes are dropped upon aggregation.
G(POP, c(1, 10)) \%>\% fmean(ix(.)$iso3c)

# Example of index with multiple ids
GGDC10S |> findex_by(Variable, Country, Year) |> head() # default is interact.ids = TRUE
GGDCi <- GGDC10S |> findex_by(Variable, Country, Year, interact.ids = FALSE)
head(GGDCi)
findex(GGDCi)
# The benefit is increased flexibility for summary statistics and data transformation
qsu(GGDCi, effect = "Country")
STD(GGDCi$SUM, effect = "Variable")            # Standardizing by variable
STD(GGDCi$SUM, effect = c("Variable", "Year")) # ... by variable and year
# But time-based operations are a bit more expensive because of the necessary interactions
D(GGDCi$SUM)

# Panel-Data modelling ---------------------------------------------------------

# Linear model of 5-year annualized growth rates of GDP on Life Expactancy + 5y lag
lm(G(PCGDP, 5, p = 1/5) ~ L(G(LIFEEX, 5, p = 1/5), c(0, 5)), wldi) # p abbreviates "power"

# Same, adding time fixed effects via plm package: need to utilize to_plm function
plm::plm(G(PCGDP, 5, p = 1/5) ~ L(G(LIFEEX, 5, p = 1/5), c(0, 5)), to_plm(wldi), effect = "time")

# With country and time fixed effects via fixest
fixest::feols(G(PCGDP, 5, p=1/5) ~ L(G(LIFEEX, 5, p=1/5), c(0, 5)), wldi, fixef = .c(iso3c, year))
\dontrun{ % Not suggested packages
# Running a robust MM regression without fixed effects
robustbase::lmrob(G(PCGDP, 5, p = 1/5) ~ L(G(LIFEEX, 5, p = 1/5), c(0, 5)), wldi)

# Running a robust MM regression with country and time fixed effects
wldi |> fselect(PCGDP, LIFEEX) |>
  fgrowth(5, power = 1/5) |> ftransform(LIFEEX_L5 = L(LIFEEX, 5)) |>
  # drop abbreviates drop.index.levels (not strictly needed here but more consistent)
  na_omit(drop = "all") |> fhdwithin(na.rm = FALSE) |> # For TFE use fwithin(effect = "year")
  unindex() |> robustbase::lmrob(formula = PCGDP ~.)    # using lm() gives same result as fixest

# Using a random forest model without fixed effects
# ranger does not support these kinds of formulas, thus we need some preprocessing...
wldi |> fselect(PCGDP, LIFEEX) |>
  fgrowth(5, power = 1/5) |> ftransform(LIFEEX_L5 = L(LIFEEX, 5)) |>
  unindex() |> na_omit() |> ranger::ranger(formula = PCGDP ~.)
}

# Indexing other data frame based classes --------------------------------------

library(tibble)
wlditbl <- qTBL(wlddev) |> findex_by(iso3c, year)
wlditbl[,2] # Works like a tibble...
wlditbl[[2]]
wlditbl[1:1000, 10]
head(wlditbl)

library(data.table)
wldidt <- qDT(wlddev) |> findex_by(iso3c, year)
wldidt[1:1000]      # Works like a data.table...
wldidt[year > 2000]
wldidt[, .(sum_PCGDP = sum(PCGDP, na.rm = TRUE)), by = country] # Aggregation unindexes the result
wldidt[, lapply(.SD, sum, na.rm = TRUE), by = country, .SDcols = .c(PCGDP, LIFEEX)]
# This also works but is a bit inefficient since the index is subset and then dropped
# -> better unindex beforehand
wldidt[year > 2000, .(sum_PCGDP = sum(PCGDP, na.rm = TRUE)), by = country]
wldidt[, PCGDP_gr_5Y := G(PCGDP, 5, power = 1/5)]  # Can add Variables by reference
# Note that .SD is a data.table of indexed_series, not an indexed_frame, so this is WRONG!
wldidt[, .c(PCGDP_gr_5Y, LIFEEX_gr_5Y) := G(slt(.SD, PCGDP, LIFEEX), 5, power = 1/5)]
# This gives the correct outcome
wldidt[, .c(PCGDP_gr_5Y, LIFEEX_gr_5Y) := lapply(slt(.SD, PCGDP, LIFEEX), G, 5, power = 1/5)]
%# wldidt[, .c(PCGDP_growth_5Y, LIFEEX_growth_5Y) := G(slt(reindex(.SD, ix(wldidt)), PCGDP, LIFEEX), 5, power = 1/5)] # Works !!
\dontrun{
library(sf)
nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
nci <- findex_by(nc, SID74)
nci[1:10, "AREA"]
st_centroid(nci) # The geometry column is never indexed, thus sf computations work normally
st_coordinates(nci)
fmean(st_area(nci))

library(tsibble)
pedi <- findex_by(pedestrian, Sensor, Date_Time)
pedi[1:5, ]
findex(pedi) # Time factor with 17k levels from POSIXct
# Now here is a case where integer levels in the index can really speed things up
ix(iby(pedestrian, Sensor, timeid(Date_Time)))
library(microbenchmark)
microbenchmark(descriptive_levels = findex_by(pedestrian, Sensor, Date_Time),
               integer_levels = findex_by(pedestrian, Sensor, timeid(Date_Time)))
# Data has irregularity
is_irregular(pedi)
is_irregular(pedi, any_id = FALSE) # irregularity in all sequences
# Manipulation such as lagging with tsibble/dplyr requires expanding rows and grouping
# Collapse can just compute correct lag on indexed series or frames
library(dplyr)
microbenchmark(
  dplyr = fill_gaps(pedestrian) |> group_by_key() |> mutate(Lag_Count = lag(Count)),
  collapse = fmutate(pedi, Lag_Count = flag(Count)), times = 10)
}
# Indexing Atomic objects ---------------------------------------------------------

## ts
print(AirPassengers)
AirPassengers[-(20:30)]        # Ts class does not support irregularity, subsetting drops class
G(AirPassengers[-(20:30)], 12) # Annual Growth Rate: Wrong!
# Now indexing AirPassengers (identity() is a trick so that the index is named time(AirPassengers))
iAP <- reindex(AirPassengers, identity(time(AirPassengers)))
iAP
findex(iAP)    # See the index
iAP[-(20:30)]  # Subsetting
G(iAP[-(20:30)], 12)                # Annual Growth Rate: Correct!
L(G(iAP[-(20:30)], c(0,1,12)), 0:1) # Lagged level, period and annual growth rates...
\donttest{ % No code relying on suggested package
## xts
library(xts)
library(zoo) # Needed for as.yearmon() and index() functions
X <- wlddev |> fsubset(iso3c == "DEU", date, PCGDP:POP) \%>\% {
  xts(num_vars(.), order.by = as.yearmon(.$date))
  } |> ss(-(30:40)) \%>\% reindex(identity(index(.))) # Introducing a gap
# plot(G(unindex(X)))
diff(unindex(X))    # diff.xts gixes wrong result
fdiff(X)            # fdiff gives right result

# But xts range-based subsets do not work...
\dontrun{
X["1980/"]
}
# Thus a better way is not to index and perform ad-hoc omputations on the xts index
X <- unindex(X)
X["1980/"] \%>\% fdiff(t = index(.)) # xts index is internally processed by timeid()
}
## Of course you can also index plain vectors / matrices...
options(oldopts)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ manip }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
