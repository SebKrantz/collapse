\name{collapse-options}
\alias{collapse-options}

\title{\emph{collapse} Package Options}

\description{
\itemize{
\item \code{option("collapse_unused_arg_action")} regulates how generic functions (such as the \link[=A1-fast-statistical-functions]{Fast Statistical Functions}) in the package react when an unknown argument is passed to a method. The default action is \code{"warning"} which issues a warning. Other options are \code{"error"}, \code{"message"} or \code{"none"}, whereby the latter enables silent swallowing of such arguments.
\item \code{option("collapse_DT_alloccol")} sets how many empty columns \emph{collapse} data manipulation functions like \code{ftransform} allocate when taking a shallow copy of \emph{data.table}'s. The default is \code{100L}. Note that the \emph{data.table} default is \code{getOption("datatable.alloccol") = 1024L}.
}
}

\keyword{documentation}
