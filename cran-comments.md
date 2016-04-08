# mwaved v1.1.3 update candidate from v1.1.2 package. 

* Update to facilitate syntax changes in the `testthat` dependency. Otherwise errors are thrown during the `testthat` checks.

## Test environments
* locally Windows 10, R 3.2.4
* win-builder (devel 2016-03-12 r70325 and release)

## R CMD check results

* There were no ERRORs or WARNINGs. 
* There were 2 NOTEs:
    - Possibly mis-spelled words in DESCRIPTION:
        - Deconvolution (3:29)
        - Multichannel (3:8)
        - deconvolution (9:35)
    - Updated maintainer's email address.