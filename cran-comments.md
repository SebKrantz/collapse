## Resubmission
This is a resubmission (Nr. 4, original submission was on 13.01.2020). In this version I have:

* Added Martyn Plummer and 1999-2016 The R Core Team as cph in Authors@R as demanded by Swetlana Herbrandt in a response from 6th March 2020. In the previous submission R Core Team and contributors worldwide was just added as ctb.

* Added 2 additional vignettes. 

## Test environments
* local Windows 8.1 install, R 3.6.1
* win-builder (devel and release)
* Ubuntu Linux 16.04 LTS, R-release, GCC (on Rhub)
* macOS 10.11 El Capitan, R-release (on Rhub)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

  * checking installed package size ... NOTE
    installed size is  7.8Mb
    sub-directories of 1Mb or more:
      doc    1.7Mb
      libs   4.7Mb

This has to do with compiled files. Data is only 0.5 Mb. 
The note was already there in the previous submissions, 
the added vignettes are small.
