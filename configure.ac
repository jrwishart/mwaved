# The version set here will propagate to other files from here
AC_INIT([mwaved], 1.0.0)

# Checks for common programs using default macros

## Use Rscript to query Rcpp for compiler and linker flags
##
## look for Rscript, but use the one found via R_HOME to allow for multiple installations
AC_DEFUN(AC_PROG_RSCRIPT, [AC_CHECK_PROG(RSCRIPT,Rscript,yes)])
AC_PROG_RSCRIPT
if test x"${RSCRIPT}" == x"yes" ; then
     : ${R_HOME=`R RHOME`}
     if test -z "${R_HOME}"; then
         AC_MSG_ERROR([Could not determine R_HOME.])
     fi
     ## Rcpp compile flag providing header directory containing Rcpp.h
     RCPP_CXXFLAGS=`${R_HOME}/bin/Rscript -e 'Rcpp:::CxxFlags()'`
     ## link flag providing libary as well as path to library, and optionally rpath
     RCPP_LDFLAGS=`${R_HOME}/bin/Rscript -e 'Rcpp:::LdFlags()'`
else
    AC_MSG_ERROR([Rscript not found, is R installed?])
fi

# Check if the FFTW Library exists on the system
AC_MSG_NOTICE([Checking if FFTW3 library is available using pkg-config])
## --- FFTW 3 ----------------------------------------------------------------
PKG_CHECK_MODULES(FFTW, fftw3, [HAVE_FFTW=TRUE], [HAVE_FFTW=FALSE])

if test "$HAVE_FFTW" == TRUE
then
	AC_MSG_RESULT([FFTW3 found!])
	# See if all functions exist for correct compilation
	AC_MSG_NOTICE([Testing all required C level functions])
	# Insert all the compiler checks
	AC_CHECK_FUNCS([memset])
	# Math library checks
	AC_CHECK_HEADER(math.h)
	AC_CHECK_LIB([m], [cos])
	AC_CHECK_LIB([m], [sin])
	AC_CHECK_LIB([m], [sqrt])
	AC_CHECK_LIB([m], [floor])
	# FFTW library checks.
	AC_CHECK_LIB([fftw3], [fftw_malloc])
	AC_CHECK_LIB([fftw3], [fftw_free])
	AC_CHECK_LIB([fftw3], [fftw_execute])
	AC_CHECK_LIB([fftw3], [fftw_destroy_plan])
	AC_CHECK_LIB([fftw3], [fftw_plan_many_dft_r2c])
	AC_CHECK_LIB([fftw3], [fftw_plan_dft_c2r])
	AC_CHECK_LIB([fftw3], [fftw_plan_dft_r2c_1d])
else
  # Return message that package requires FFTW3 library
  AC_MSG_RESULT([FFTW library not found, please install fftw-3-3-3 or greater])
fi

# Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(FFTW_CFLAGS)
AC_SUBST(FFTW_LIBS)

AC_SUBST(CPPFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LIBS)

AC_SUBST(RCPP_CXXFLAGS)

AC_SUBST(RCPP_LDFLAGS)

AC_OUTPUT(src/Makevars)
