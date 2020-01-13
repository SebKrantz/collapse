## Resubmission
This is a resubmission. In this version I have:

* Renamed all.identical to all_identical

* Removed code suppressing pragma warnings in .c files and fixed the issues producing these warnings

* Replaced variable-sized arrays in C++ with vectors to eliminate compiler warnings. 


## License and Authorship
Since C-code from data.table and stats was copied/modified into src, data.table's
Mozilla Public License 2.0 was adopted and the authors listed as 'ctb'. 
Additionally, Simon Gaure, the author of lfe from which I import 
the demeanlist function (written in C), and Dirk Eddelbuettel, whose package
Rcpp I have thoroughly exploited in the creation of this one, were added as 'ctb'.
In addition, all authors have been mentioned on the collapse-package.Rd page.

(An alternative could be listing us all as 'aut'. My current understanding though 
 is that an 'aut' is someone with substantial stakes in the design, purpose or functionality 
 of a package, which is not the case here (none of these people know about or have actively 
 contributed to this package, I have not simply duplicated functionality of packages from which
 C-code was taken without substantial modifications, and the overall functionality of 
 collapse far outstrips that of the borrowed code).)

(The C-functions imported can be seen in src/ExportSymbols.cpp lines 5-15, 
 the rest is C++ of which I am the sole author).

## Test environments
* local Windows 8.1 install, R 3.6.1
* win-builder (devel and release)
* Ubuntu Linux 16.04 LTS, R-release, GCC (on Rhub)
* macOS 10.11 El Capitan, R-release (on Rhub)


## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

  * checking installed package size ... NOTE
    installed size is  6.7Mb
    sub-directories of 1Mb or more:
      libs   4.6Mb

This has to do with compiled files (.dll's). Pre-compilation, 
the size of all .R, .c, .cpp, .h, .man, .rda and .Rmd files 
together is about 2.4 Mb, of which 0.5 Mb is data (.rda).
