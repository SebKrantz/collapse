\name{collapse-renamed}
\alias{collapse-renamed}
\alias{fNobs}
\alias{fNobs.default}
\alias{fNobs.matrix}
\alias{fNobs.data.frame}
\alias{fNobs.grouped_df}
\alias{fNdistinct}
\alias{fNdistinct.default}
\alias{fNdistinct.matrix}
\alias{fNdistinct.data.frame}
\alias{fNdistinct.grouped_df}
\alias{pwNobs}
\alias{fHDwithin}
\alias{fHDwithin.default}
\alias{fHDwithin.matrix}
\alias{fHDwithin.data.frame}
\alias{fHDwithin.pseries}
\alias{fHDwithin.pdata.frame}
\alias{fHDwithin.grouped_df}
\alias{fHDbetween}
\alias{fHDbetween.default}
\alias{fHDbetween.matrix}
\alias{fHDbetween.data.frame}
\alias{fHDbetween.pseries}
\alias{fHDbetween.pdata.frame}
\alias{fHDbetween.grouped_df}
\alias{as.factor_GRP}
\alias{as.factor_qG}
\alias{is.GRP}
\alias{is.qG}
\alias{is.unlistable}
\alias{is.categorical}
\alias{is.Date}
\alias{as.character_factor}
\alias{as.numeric_factor}
\alias{Date_vars}
\alias{Date_vars<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Renamed Functions
}
\description{
These functions were renamed moving from collapse 1.5.3 to 1.6.0 to make the namespace more consistent. With collapse 1.9.0, all of these functions are depreciated, and the core methods to \code{fNobs}, \code{fNdistinct}, \code{fHDbetween} and \code{fHDwithin} are no longer exported. I intend to remove all of these functions except for the S3 generics \code{fNobs}, \code{fNdistinct}, \code{fHDbetween} and \code{fHDwithin} from the package by end of 2023. %The S3 generics and the other functions will be depreciated in 2023 for the earliest. These all now give a message reminding you not to use them in fresh code.
}
\section{Renaming}{\preformatted{
fNobs -> fnobs
fNdistinct -> fndistinct
pwNobs -> pwnobs
fHDwithin -> fhdwithin
fHDbetween -> fhdbetween
as.factor_GRP -> as_factor_GRP
as.factor_qG -> as_factor_qG
is.GRP -> is_GRP
is.qG -> is_qG
is.unlistable -> is_unlistable
is.categorical -> is_categorical
is.Date -> is_date
as.numeric_factor -> as_numeric_factor
as.character_factor -> as_character_factor
Date_vars -> date_vars
`Date_vars<-` -> `date_vars<-`
}
}
