\name{collapse-options}
\alias{collapse-options}

\title{\emph{collapse} Package Options}

\description{
Currently \emph{collapse} only provides \code{option("collapse_unused_arg_action")}, which regulates how generic functions (such as the \link[=A1-fast-statistical-functions]{Fast Statistical Functions}) in the package react when an unknown argument is passed to a method. The default action is \code{"warning"} which issues a warning. Other options are \code{"error"}, \code{"message"} or \code{"none"}, whereby the latter enables silent swallowing of such arguments.
}

\keyword{documentation}
