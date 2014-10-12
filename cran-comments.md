This is a candidate v1.1.0 update to the mwaved 1.0.1 package. 

Major changes

------------------------------------------------

* Changed interface slightly but adding adding a blur detection function, detectBlur, so the user is not required to specify this manually. Functions checks if input G matches structure of direct or box.car blur and returns character string 'direct' or 'box.car' respectively,  (returns 'smooth' otherwise). 
* Added resolution argument so that the user can specify the resolution selection method, warnings are thrown if the detectBlur function does not match the desired selection type.
* Added HTML vignette to document the package.
* Added unit tests via the testthat package.


There are a few minor patches fixed as well.

-------------------------------------------------

* Fixed plotting errors for smooth resolution method.
* Fixed ggplot labelling issue.
* Fixed bug in coarse coefficient computation in multiWaveD code.
* Made arguments in functions consistent.
