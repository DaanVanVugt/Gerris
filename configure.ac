AC_INIT(configure.ac)
AC_CONFIG_MACRO_DIR([m4])

dnl we to AC_DIVERT_PUSH/AC_DIVERT_POP these variable definitions so they
dnl are available for $ac_help expansion (don't we all *love* autoconf?)
#AC_DIVERT_PUSH(AC_DIVERSION_NOTICE)dnl
# Making releases:
#   GFS_MICRO_VERSION += 1;
#   GFS_INTERFACE_AGE += 1;
#   GFS_BINARY_AGE += 1;
# if any functions have been added, set GFS_INTERFACE_AGE to 0.
# if backwards compatibility has been broken,
# set GFS_BINARY_AGE and GFS_INTERFACE_AGE to 0.
#
# Remember to update rpm/gerris.spec when changing the version number.
GFS_MAJOR_VERSION=1
GFS_MINOR_VERSION=3
GFS_MICRO_VERSION=2
GFS_INTERFACE_AGE=0
GFS_BINARY_AGE=0
GFS_VERSION=$GFS_MAJOR_VERSION.$GFS_MINOR_VERSION.$GFS_MICRO_VERSION
GFS_COMPILATION_FLAGS=$CFLAGS
dnl
#AC_DIVERT_POP()dnl

AC_SUBST(GFS_MAJOR_VERSION)
AC_SUBST(GFS_MINOR_VERSION)
AC_SUBST(GFS_MICRO_VERSION)
AC_SUBST(GFS_VERSION)
AC_SUBST(GFS_COMPILATION_FLAGS)

# change this to the gerris web site root
WEBROOT=/var/www/site-gerris/gerris
AC_SUBST(WEBROOT)

# libtool versioning
LT_RELEASE=$GFS_MAJOR_VERSION.$GFS_MINOR_VERSION
LT_CURRENT=`expr $GFS_MICRO_VERSION - $GFS_INTERFACE_AGE`
LT_REVISION=$GFS_INTERFACE_AGE
LT_AGE=`expr $GFS_BINARY_AGE - $GFS_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

# For automake.
VERSION=$GFS_VERSION
PACKAGE=gerris

AM_INIT_AUTOMAKE($PACKAGE, $VERSION, no-define)

# Specify a configuration file
AM_CONFIG_HEADER(config.h)

AC_DEFINE_UNQUOTED(GFS_MAJOR_VERSION, $GFS_MAJOR_VERSION, [Major version])
AC_DEFINE_UNQUOTED(GFS_MINOR_VERSION, $GFS_MINOR_VERSION, [Minor version])
AC_DEFINE_UNQUOTED(GFS_MICRO_VERSION, $GFS_MICRO_VERSION, [Micro version])
AC_DEFINE_UNQUOTED(GFS_INTERFACE_AGE, $GFS_INTERFACE_AGE, [Interface age])
AC_DEFINE_UNQUOTED(GFS_BINARY_AGE, $GFS_BINARY_AGE,       [Binary age])
AC_DEFINE_UNQUOTED(GFS_VERSION, "$GFS_VERSION",           [Version])
AC_DEFINE_UNQUOTED(GFS_COMPILATION_FLAGS, "$GFS_COMPILATION_FLAGS", [Compilation flags])

dnl Initialize maintainer mode
AM_MAINTAINER_MODE

AC_ARG_WITH(mpicc,
              [  --with-mpicc=path       pick a specific path to mpicc],
              [with_mpicc=$withval],
              [with_mpicc=''])
AC_ARG_WITH(suffix,
              [  --with-suffix=.a        override the glib dynamic module suffix],
              [with_suffix=$withval],
              [with_suffix=default])

if test ! "x$with_suffix" = "xdefault"; then
   AC_DEFINE_UNQUOTED(MODULES_SUFFIX, "$with_suffix", [Dynamic module suffix])
fi

# check if MPI support is disabled
AC_ARG_ENABLE(mpi,
[  --disable-mpi           do not compile MPI support],
[ case "${enableval}" in
	yes) ;;
	*) with_mpicc=no ;;
  esac])

if test "x$with_mpicc" != "xno" ; then
  if test "x$with_mpicc" != "x" ; then
     if test -x "$with_mpicc" ; then
       use_mpicc=yes
       CC="$with_mpicc"
     fi
  else
     AC_CHECK_PROG(use_mpicc, mpicc, yes, no)
     if test "x$use_mpicc" = "xyes" ; then
       CC=mpicc
     fi
  fi
  if test "x$use_mpicc" = "xyes" ; then
     AC_DEFINE(HAVE_MPI, 1, [Define to 1 if you have MPI])
  else
     AC_MSG_WARN([MPI not found. No MPI support will be compiled in.])
  fi
else
  use_mpicc=no
fi

AC_PROG_CC

# -D_GNU_SOURCE is only necessary for old (< 2.10) glibc implementations of open_memstream()
CFLAGS="$CFLAGS -D_XOPEN_SOURCE=700 -D_BSD_SOURCE -D_GNU_SOURCE"
if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wmissing-prototypes -Wmissing-declarations -pipe -std=c99"
fi

dnl Initialize libtool
AC_LIBTOOL_WIN32_DLL
AM_PROG_LIBTOOL
AM_PROG_CC_C_O

dnl Enable verbose output from libtool (enable libtool's default)
AC_ARG_ENABLE(libtool-verbose,
	      [  --enable-libtool-verbose enable verbose libtool output (default disabled)],
	      [with_libtool_verbose=$enableval],
	      [with_libtool_verbose='no'])
if test "$with_libtool_verbose" = 'no'
then
  LIBTOOL="$LIBTOOL --silent"
fi

dnl Fix for MacOSX and cygwin libtool
case "$build" in
  *-apple-darwin*) 
  		   NO_UNDEFINED=""
		   CFLAGS="$CFLAGS -D_DARWIN_C_SOURCE"
  		   ;;
  *-pc-cygwin*)    NO_UNDEFINED="" ;;
  *)               NO_UNDEFINED="-no-undefined" ;;
esac
AC_SUBST(NO_UNDEFINED)

AC_PROG_AWK
AC_ISC_POSIX
AC_HEADER_STDC

# check that pointers can be stored in doubles
AC_MSG_CHECKING(whether pointers can be stored in doubles)
AC_TRY_RUN([
main () {
  if (sizeof (void *) > sizeof (double))
    return 1;
  return 0;
}
], can_store_pointers=yes, can_store_pointers=no,
[can_store_pointers=yes; echo $ac_n "cross compiling; assumed OK... $ac_c"])
AC_MSG_RESULT($can_store_pointers)
if test x$can_store_pointers = xno ; then
  AC_MSG_ERROR([
*** Pointers cannot be stored in doubles on this architecture.])
fi

PKG_CHECK_MODULES(GTS, [gts >= 0.7.4])

# check if we want to enable GTS casts checks
AC_ARG_ENABLE(gts-check,
[  --enable-gts-check      enable object type cast checks in GTS],
[ case "${enableval}" in
	*) gts_check_casts="--check" ;;
  esac])

if test -n "$gts_check_casts"; then
  GTS_CFLAGS="$GTS_CFLAGS -DGTS_CHECK_CASTS"
fi
AC_SUBST(GTS_CFLAGS)
AC_SUBST(GTS_LIBS)

# check whether GModules are supported
AC_MSG_CHECKING(whether modules are supported)
OLD_CFLAGS=$CFLAGS
CFLAGS=$GTS_CFLAGS
OLD_LIBS=$LIBS
LIBS=$GTS_LIBS
AC_TRY_RUN([
#include <glib.h>
#include <gmodule.h>

main () {
  if (!g_module_supported ())
    return 1;
  return 0;
}
], have_gmodule=yes, have_gmodule=no,
[have_gmodule=yes; echo $ac_n "cross compiling; assumed OK... $ac_c"])
AC_MSG_RESULT($have_gmodule)
CFLAGS=$OLD_CFLAGS
LIBS=$OLD_LIBS

# checking how to build modules
if test "$have_gmodule" = "yes"; then
   # try gcc first
   module_flags=""
   if test x$GCC = xyes; then
     gcc_flags="-fpic -x c"
     for flags in \
         "-shared" \
         "-dynamiclib -undefined dynamic_lookup" \
         "-bundle"
     do
	if test x"$module_flags" = x; then
	   if echo "void test(){}" | $CC $flags $gcc_flags - 2> /dev/null; then
	      module_flags="$flags"
	   fi
	fi
     done
   fi
   if test x"$module_flags" = x; then
       # use libtool if not using gcc
       sed -e "s|PREFIX|$prefix|g" -e "s|COMPILER|$CC|g" < build_function_libtool > build_function
   else
       sed -e "s|MODULE_FLAGS|$module_flags|" -e "s|COMPILER|$CC|g" < build_function_gcc > build_function
   fi
else
   # cannot use modules
   cat <<EOF > build_function
#!/bin/sh
echo "functions are not supported on this system"
exit 1
EOF
fi

# checks for pkg-config
if test "$have_gmodule" = "yes"; then
  AC_CHECK_PROG(have_pkg_config, pkg-config, yes, no)
  if test "$have_pkg_config" = "no"; then
    AC_MSG_WARN([pkg-config is required for functions in parameter files.])
  fi
fi
if test x"$have_pkg_config" = x"yes"; then
  AC_DEFINE(HAVE_PKG_CONFIG, 1, [Define to 1 if GModule are supported and if you have pkg-config.])
else
  AC_DEFINE(HAVE_PKG_CONFIG, 0)
fi
if test x"$have_pkg_config" = x""; then
   have_pkg_config=undefined
fi

# checks for m4
AC_CHECK_PROG(have_m4, m4, yes, no)
options_m4=""
if test "$have_m4" = "yes"; then
   AC_DEFINE(HAVE_M4, 1, [Define to 1 if m4 is available.])
   AC_MSG_CHECKING(check whether m4 accepts the -P option)
   if m4 -P < /dev/null 2> /dev/null; then
      options_m4="-P"
      accepts_P="yes"
   else
      accepts_P="no"
   fi
   AC_MSG_RESULT($accepts_P)
else
   AC_DEFINE(HAVE_M4, 0)
fi
AC_DEFINE_UNQUOTED(M4_OPTIONS, "$options_m4", [Options for m4])

AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)

# checks for libproj
AC_CHECK_LIB(proj, pj_fwd, proj="yes",
  AC_MSG_WARN([libproj not found. Map module will not be available.]), [-lm])
AC_CHECK_HEADERS(proj_api.h, proj="yes", proj="no")
AM_CONDITIONAL(HAS_LIBPROJ, test x$proj = xyes)

# checks for GNU Scientific Library
AM_PATH_GSL
AM_CONDITIONAL(HAS_GSL, test x$no_gsl != xyes)

# checks for netCDF
AC_CHECK_LIB(netcdf, nc_open, netcdf="true", netcdf="false", [-lm])
AC_CHECK_HEADERS(netcdf.h, netcdf="true", netcdf="false")
AM_CONDITIONAL(HAS_NETCDF, test x$netcdf = xtrue)

# Tide module requires both netCDF and GSL
if test x$netcdf = xfalse -o x$no_gsl = xyes; then
   AC_MSG_WARN([GSL and/or netCDF not found. Tide module will not be available.])
   tide=no
else
   tide=yes
fi
AM_CONDITIONAL(BUILD_TIDE, test x$netcdf = xtrue -a x$no_gsl = x)

# Stokes module requires a fortran compiler
AC_PROG_F77
if test x$F77 != x; then
   AC_F77_WRAPPERS
fi
AM_CONDITIONAL(BUILD_STOKES, test x$F77 != x)

# Enable the AGMG module (disabled by default due to non-free license)
AC_ARG_ENABLE(agmg,
	      [  --enable-agmg enable the AGMG module (beware of non-free license)],
	      [with_agmg=$enableval],
	      [with_agmg="no"])

# AGMG module requires the gfortran compiler and lapack/blas
if test x$with_agmg != xno; then
   AC_PROG_FC([gfortran f90 xlf90],90)
   AC_CHECK_LIB(lapack, sorg2l_, lapack="true", lapack="false", [-lblas -lm])
fi
AM_CONDITIONAL(BUILD_AGMG, test x$FC = xgfortran -a x$lapack = xtrue -a x$with_agmg != xno)

# Build wavewatch module if wavewatch is installed
AM_CONDITIONAL(BUILD_WAVEWATCH, test -r $HOME/.wwatch3.env)
if test -r $HOME/.wwatch3.env; then
   ww3_dir=`grep WWATCH3_DIR $HOME/.wwatch3.env | awk '{print $2}'`
   ww3_version=`grep "WWVER  = " $ww3_dir/ftn/*.ftn | \
                awk '{print substr($(NF-1),2)}' | \
                sed 's/\.//'`
   W3INIT="w3init$ww3_version"
   AC_SUBST(W3INIT)
   AC_DEFINE_UNQUOTED(WW3_VERSION, $ww3_version, [WaveWatch III version number])
fi

# Build ode module if ODE is installed
AC_CHECK_PROG(ode, ode-config, yes, no)
if test x$ode = xyes; then
  ODE_CFLAGS=`ode-config --cflags`
  ODE_LIBS="`ode-config --libs` -lstdc++"
  AC_SUBST(ODE_CFLAGS)
  AC_SUBST(ODE_LIBS)
else
  AC_MSG_WARN([ODE not found. ode module will not be available.])
fi
AM_CONDITIONAL(BUILD_ODE, test x$ode = xyes)

# Build hypre module if HYPRE is installed
hypre=no
if test "x$use_mpicc" = "xyes" ; then # HYPRE needs MPI
   AC_CHECK_LIB(HYPRE, HYPRE_IJMatrixCreate, hypre="yes", hypre="no", [-lm])
   if test x$hypre = xyes; then
      AC_CHECK_HEADERS(HYPRE.h, hypre="yes", hypre="no")
   fi
   if test x$hypre = xno; then
      AC_MSG_WARN([HYPRE not found. hypre Module won't be available.])
   fi
fi
AM_CONDITIONAL(BUILD_HYPRE, test x$hypre = xyes)

# Build lis module if lis is installed
AC_CHECK_LIB(lis, lis_initialize, lis="yes", lis="no", [-lgfortran -lm])
if test x$lis = xyes; then
   AC_CHECK_HEADERS(lis.h, lis="yes", lis="no")
fi
if test x$lis = xno; then
   AC_MSG_WARN([lis not found. lis Module won't be available.])
fi
AM_CONDITIONAL(BUILD_LIS, test x$lis = xyes)

# Build fft module if fttw3 is installed
PKG_CHECK_MODULES(FFTW3, [ fftw3 >= 3 ], [ fftw3="yes" ], [ fftw3="no" ])
AM_CONDITIONAL(HAS_FFTW3, test x$fftw3 = xyes)
AC_SUBST(FFTW3_CFLAGS)
AC_SUBST(FFTW3_LIBS)

# header file checks
AC_CHECK_HEADERS(fenv.h, AC_DEFINE(HAVE_FENV_H))
AC_CHECK_HEADERS(unistd.h, AC_DEFINE(HAVE_UNISTD_H))
AC_CHECK_HEADERS(getopt.h, AC_DEFINE(HAVE_GETOPT_H))

# functions checks
OLD_CFLAGS=$CFLAGS
CFLAGS=$GTS_CFLAGS
OLD_LIBS=$LIBS
LIBS=$GTS_LIBS
AC_CHECK_FUNCS(getopt_long g_mkdtemp open_memstream funopen)
CFLAGS=$OLD_CFLAGS
LIBS=$OLD_LIBS

# use alternative implementation of open_memstream() if the system
# (e.g. MacOS X 10.3, FreeBSD 6.0, NetBSD 3.0, OpenBSD 3.8.) does not
# support POSIX 2008 natively
AC_CHECK_FUNC(open_memstream, have_open_memstream=yes, have_open_memstream=no)
AM_CONDITIONAL(OPEN_MEMSTREAM, test x$have_open_memstream = xno)

# link flags
GFS2D_LIBS="\$(top_builddir)/src/libgfs2D.la $GTS_LIBS"
GFS3D_LIBS="\$(top_builddir)/src/libgfs3D.la $GTS_LIBS"
AC_SUBST(GFS2D_LIBS)
AC_SUBST(GFS3D_LIBS)

# check whether the code is under darcs
AM_CONDITIONAL(DARCS_CONTROLLED, test -d _darcs)

AC_CONFIG_FILES([
Makefile
src/Makefile
src/gerris2D.pc
src/gerris3D.pc
tools/Makefile
modules/Makefile
modules/RStarTree/Makefile
modules/kdt/Makefile
modules/fes2004/Makefile
modules/wavewatch/Makefile
modules/culvert/Makefile
test/Makefile
doc/Makefile
doc/tutorial/Makefile
doc/examples/Makefile
doc/manpages/Makefile
doc/examples/gfs2doc
doc/examples/gfs-highlight
doc/examples/crossref.sh
desktop/Makefile
debian/Makefile
])


AC_OUTPUT()

echo
echo "Configuration summary"
echo
echo "  System            = $build"
echo "  C   Compiler      = $CC"
echo "  C   Flags         = $CFLAGS"
echo "  MPI enabled       = $use_mpicc"
echo "  GModule support   = $have_gmodule"
echo "  pkg-config        = $have_pkg_config"
echo "  m4                = $have_m4"
if test x$have_m4 = xyes; then
   echo "  m4 options        = $options_m4"
fi
echo "  open_memstream    = $have_open_memstream"
echo
echo "Modules"
echo
echo "  map               = $proj"
echo "  tide              = $tide"
if test x$F77 != x; then
echo "  stokes            = yes"
else
echo "  stokes            = no"
fi
if test -r $HOME/.wwatch3.env; then
echo "  wavewatch         = yes"
echo "  wavewatch version = $ww3_version"
else
echo "  wavewatch         = no"
fi
echo "  hypre             = $hypre"
echo "  lis               = $lis"
if test x$FC = xgfortran -a x$lapack = xtrue -a x$with_agmg != xno; then
echo "  agmg              = yes"
else
echo "  agmg              = no"
fi
echo "  ode               = $ode"
if test x$no_gsl != xyes; then
echo "  bubbles           = yes"
else
echo "  bubbles           = no"
fi
if test x$fftw3 != xno; then
echo "  fft               = yes"
else
echo "  fft               = no"
fi
echo
