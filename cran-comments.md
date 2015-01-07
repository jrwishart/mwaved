# mwaved v1.1.1 update candidate from v1.1.0 package. 

The updates are to address the new POSIX standard and address some other issues flagged in the R-devel version.

## Update description

* Removed ` += ` from the Makevars.in file
* Removed the `require(.)` commands from the source code and replaced with `requireNamespace( ., quietly = TRUE)`.
* Changed `NeedsCompilation:Yes` to `NeedsCompilation:yes` in the DESCRIPTION
* Misc: 
    + Changed code dropping suggested package `gridExtra` to prefer a leaner import a few functions from the `r-base` `grid` package.
    + Cleaned code for readability.
    + Converted the DESCRIPTION title to title case

## Test environments
* locally Ubuntu 14.04, R 3.1.2 and r-devel 2015-01-06 r67330
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

Possibly mis-spelled words in DESCRIPTION:
  Deconvolution (3:29)
  Multichannel (3:8)
  deconvolution (9:35)
