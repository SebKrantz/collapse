## Submission collapse 1.1.0
In this version I have:

* Fixed the 2 test failures on windows old release (r-oldrel-windows-ix86+x86_64) as demanded by Prof. Brian Ripley in a e-mail from 20th of March.

* Fixed some minor issues and added small improvements to the code. 

* Exported some small additional functions and reworked fbetween and fscale for more customized scaling and centering of data. 

## Test environments
* local Windows 8.1 install, R 3.6.1
* win-builder (devel and release)
* Ubuntu Linux 16.04 LTS, R-release, GCC (on Rhub)
* macOS 10.11 El Capitan, R-release (on Rhub)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

  * checking installed package size ... NOTE
    installed size is  7.9Mb
    sub-directories of 1Mb or more:
      doc    1.7Mb
      libs   4.7Mb

This has to do with compiled files. Data is 0.5 Mb. 
