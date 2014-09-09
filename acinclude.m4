# Configure paths for GTS
# Stéphane Popinet  2001-10-4
#       adapted from
# Configure paths for GLIB
# Owen Taylor       97-11-3

dnl AM_PATH_GTS([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND ]]])
dnl Test for GTS, and define GTS_CFLAGS and GTS_LIBS
dnl
AC_DEFUN([AM_PATH_GTS],
[dnl 
dnl Get the cflags and libraries from the gts-config script
dnl
AC_ARG_WITH(gts-prefix,[  --with-gts-prefix=PFX   Prefix where GTS is installed (optional)],
            gts_config_prefix="$withval", gts_config_prefix="")
AC_ARG_WITH(gts-exec-prefix,[  --with-gts-exec-prefix=PFX Exec prefix where GTS is installed (optional)],
            gts_config_exec_prefix="$withval", gts_config_exec_prefix="")
AC_ARG_ENABLE(gtstest, [  --disable-gtstest       Do not try to compile and run a test GTS program],
		    , enable_gtstest=yes)

  if test x$gts_config_exec_prefix != x ; then
     gts_config_args="$gts_config_args --exec-prefix=$gts_config_exec_prefix"
     if test x${GTS_CONFIG+set} != xset ; then
        GTS_CONFIG=$gts_config_exec_prefix/bin/gts-config
     fi
  fi
  if test x$gts_config_prefix != x ; then
     gts_config_args="$gts_config_args --prefix=$gts_config_prefix"
     if test x${GTS_CONFIG+set} != xset ; then
        GTS_CONFIG=$gts_config_prefix/bin/gts-config
     fi
  fi

  for module in . $4
  do
      case "$module" in
         gmodule) 
             gts_config_args="$gts_config_args gmodule"
         ;;
         gthread) 
             gts_config_args="$gts_config_args gthread"
         ;;
      esac
  done

  AC_PATH_PROG(GTS_CONFIG, gts-config, no)
  min_gts_version=ifelse([$1], ,0.4.2,$1)
  AC_MSG_CHECKING(for GTS - version >= $min_gts_version)
  no_gts=""
  if test "$GTS_CONFIG" = "no" ; then
    no_gts=yes
  else
    GTS_CFLAGS=`$GTS_CONFIG $gts_config_args --cflags`
    GTS_LIBS=`$GTS_CONFIG $gts_config_args --libs`
    gts_config_major_version=`$GTS_CONFIG $gts_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    gts_config_minor_version=`$GTS_CONFIG $gts_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    gts_config_micro_version=`$GTS_CONFIG $gts_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`
    if test "x$enable_gtstest" = "xyes" ; then
      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
      CFLAGS="$CFLAGS $GTS_CFLAGS"
      LIBS="$GTS_LIBS $LIBS"
dnl
dnl Now check if the installed GTS is sufficiently new. (Also sanity
dnl checks the results of gts-config to some extent
dnl
      rm -f conf.gtstest
      AC_TRY_RUN([
#include <gts.h>
#include <stdio.h>
#include <stdlib.h>

int 
main ()
{
  int major, minor, micro;
  char *tmp_version;

  system ("touch conf.gtstest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = g_strdup("$min_gts_version");
  if (sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) != 3) {
     printf("%s, bad version string\n", "$min_gts_version");
     exit(1);
   }

  if ((gts_major_version != $gts_config_major_version) ||
      (gts_minor_version != $gts_config_minor_version) ||
      (gts_micro_version != $gts_config_micro_version))
    {
      printf("\n*** 'gts-config --version' returned %d.%d.%d, but GTS (%d.%d.%d)\n", 
             $gts_config_major_version, $gts_config_minor_version, $gts_config_micro_version,
             gts_major_version, gts_minor_version, gts_micro_version);
      printf ("*** was found! If gts-config was correct, then it is best\n");
      printf ("*** to remove the old version of GTS. You may also be able to fix the error\n");
      printf("*** by modifying your LD_LIBRARY_PATH enviroment variable, or by editing\n");
      printf("*** /etc/ld.so.conf. Make sure you have run ldconfig if that is\n");
      printf("*** required on your system.\n");
      printf("*** If gts-config was wrong, set the environment variable GTS_CONFIG\n");
      printf("*** to point to the correct copy of gts-config, and remove the file config.cache\n");
      printf("*** before re-running configure\n");
    } 
  else if ((gts_major_version != GTS_MAJOR_VERSION) ||
	   (gts_minor_version != GTS_MINOR_VERSION) ||
           (gts_micro_version != GTS_MICRO_VERSION))
    {
      printf("*** GTS header files (version %d.%d.%d) do not match\n",
	     GTS_MAJOR_VERSION, GTS_MINOR_VERSION, GTS_MICRO_VERSION);
      printf("*** library (version %d.%d.%d)\n",
	     gts_major_version, gts_minor_version, gts_micro_version);
    }
  else
    {
      if ((gts_major_version > major) ||
        ((gts_major_version == major) && (gts_minor_version > minor)) ||
        ((gts_major_version == major) && (gts_minor_version == minor) && (gts_micro_version >= micro)))
      {
        return 0;
       }
     else
      {
        printf("\n*** An old version of GTS (%d.%d.%d) was found.\n",
               gts_major_version, gts_minor_version, gts_micro_version);
        printf("*** You need a version of GTS newer than %d.%d.%d. The latest version of\n",
	       major, minor, micro);
        printf("*** GTS is always available from http://gts.sourceforge.net.\n");
        printf("***\n");
        printf("*** If you have already installed a sufficiently new version, this error\n");
        printf("*** probably means that the wrong copy of the gts-config shell script is\n");
        printf("*** being found. The easiest way to fix this is to remove the old version\n");
        printf("*** of GTS, but you can also set the GTS_CONFIG environment to point to the\n");
        printf("*** correct copy of gts-config. (In this case, you will have to\n");
        printf("*** modify your LD_LIBRARY_PATH enviroment variable, or edit /etc/ld.so.conf\n");
        printf("*** so that the correct libraries are found at run-time))\n");
      }
    }
  return 1;
}
],, no_gts=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CFLAGS="$ac_save_CFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_gts" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$GTS_CONFIG" = "no" ; then
       echo "*** The gts-config script installed by GTS could not be found"
       echo "*** If GTS was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the GTS_CONFIG environment variable to the"
       echo "*** full path to gts-config."
     else
       if test -f conf.gtstest ; then
        :
       else
          echo "*** Could not run GTS test program, checking why..."
          CFLAGS="$CFLAGS $GTS_CFLAGS"
          LIBS="$LIBS $GTS_LIBS"
          AC_TRY_LINK([
#include <gts.h>
#include <stdio.h>
],      [ return ((gts_major_version) || (gts_minor_version) || (gts_micro_version)); ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding GTS or finding the wrong"
          echo "*** version of GTS. If it is not finding GTS, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"
          echo "***"],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means GTS was incorrectly installed"
          echo "*** or that you have moved GTS since it was installed. In the latter case, you"
          echo "*** may want to edit the gts-config script: $GTS_CONFIG" ])
          CFLAGS="$ac_save_CFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
     GTS_CFLAGS=""
     GTS_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(GTS_CFLAGS)
  AC_SUBST(GTS_LIBS)
  rm -f conf.gtstest
])

dnl  DX_INSTALL_PATH
dnl  Tries to find the location where dx is installed if it
dnl  can not then it defaults to /usr/local/dx
dnl  --------------------------------------------------------
AC_DEFUN([DX_INSTALL_PATH],
[
AC_CACHE_CHECK([for dx install path], ac_cv_dx_install_path,
[
AC_MSG_RESULT(locating)
DX_DEFAULT_INST=/usr/local/dx
AC_CHECK_PROGS( DX, dx )

DX_PATH=""
if test -n "$DX" ; then
  AC_MSG_CHECKING([for path via "dx -whereami"])
  DX_PATH=`$DX -whereami | grep "installed in" | sed -e "s/installed in //" -e "s'^\(.*\)\/$'\1'"`
  if test -z "$DX_PATH" ; then
        AC_MSG_RESULT([warning: old version of dx script in path])
  elif test "x$ARCH" = "xintelnt" ; then
	DX_PATH=`cygpath -w -s "$DX_PATH"`
  fi
fi

if test -z "$DX_PATH" ; then
  AC_MSG_CHECKING([for /usr/local/bin/dx])
  if test -x "/usr/local/bin/dx" ; then
     DX_PATH=`/usr/local/bin/dx -whereami | grep "installed in" | sed -e "s/installed in //"`
  fi

  if test -z "$DX_PATH" ; then
        AC_MSG_WARN([Missing dx script--please install OpenDX first.])
  elif test "x$ARCH" = "xintelnt" ; then
	DX_PATH=`cygpath -w -s "$DX_PATH"`
  fi
fi
ac_cv_dx_install_path=$DX_PATH
])
DX_PATH=$ac_cv_dx_install_path
])
# Configure path for the GNU Scientific Library
# Christopher R. Gabriel <cgabriel@linux.it>, April 2000


AC_DEFUN([AM_PATH_GSL],
[
AC_ARG_WITH(gsl-prefix,[  --with-gsl-prefix=PFX   Prefix where GSL is installed (optional)],
            gsl_prefix="$withval", gsl_prefix="")
AC_ARG_WITH(gsl-exec-prefix,[  --with-gsl-exec-prefix=PFX Exec prefix where GSL is installed (optional)],
            gsl_exec_prefix="$withval", gsl_exec_prefix="")
AC_ARG_ENABLE(gsltest, [  --disable-gsltest       Do not try to compile and run a test GSL program],
		    , enable_gsltest=yes)

  if test "x${GSL_CONFIG+set}" != xset ; then
     if test "x$gsl_prefix" != x ; then
         GSL_CONFIG="$gsl_prefix/bin/gsl-config"
     fi
     if test "x$gsl_exec_prefix" != x ; then
        GSL_CONFIG="$gsl_exec_prefix/bin/gsl-config"
     fi
  fi

  AC_PATH_PROG(GSL_CONFIG, gsl-config, no)
  min_gsl_version=ifelse([$1], ,0.2.5,$1)
  AC_MSG_CHECKING(for GSL - version >= $min_gsl_version)
  no_gsl=""
  if test "$GSL_CONFIG" = "no" ; then
    no_gsl=yes
  else
    GSL_CFLAGS=`$GSL_CONFIG --cflags`
    GSL_LIBS=`$GSL_CONFIG --libs`

    gsl_major_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${gsl_major_version}" = "x" ; then
       gsl_major_version=0
    fi

    gsl_minor_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${gsl_minor_version}" = "x" ; then
       gsl_minor_version=0
    fi

    gsl_micro_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${gsl_micro_version}" = "x" ; then
       gsl_micro_version=0
    fi

    if test "x$enable_gsltest" = "xyes" ; then
      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
      CFLAGS="$CFLAGS $GSL_CFLAGS"
      LIBS="$LIBS $GSL_LIBS"

      rm -f conf.gsltest
      AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* my_strdup (const char *str);

char*
my_strdup (const char *str)
{
  char *new_str;
  
  if (str)
    {
      new_str = (char *)malloc ((strlen (str) + 1) * sizeof(char));
      strcpy (new_str, str);
    }
  else
    new_str = NULL;
  
  return new_str;
}

int main (void)
{
  int major = 0, minor = 0, micro = 0;
  int n;
  char *tmp_version;

  system ("touch conf.gsltest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = my_strdup("$min_gsl_version");

  n = sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) ;

  if (n != 2 && n != 3) {
     printf("%s, bad version string\n", "$min_gsl_version");
     exit(1);
   }

   if (($gsl_major_version > major) ||
      (($gsl_major_version == major) && ($gsl_minor_version > minor)) ||
      (($gsl_major_version == major) && ($gsl_minor_version == minor) && ($gsl_micro_version >= micro)))
    {
      exit(0);
    }
  else
    {
      printf("\n*** 'gsl-config --version' returned %d.%d.%d, but the minimum version\n", $gsl_major_version, $gsl_minor_version, $gsl_micro_version);
      printf("*** of GSL required is %d.%d.%d. If gsl-config is correct, then it is\n", major, minor, micro);
      printf("*** best to upgrade to the required version.\n");
      printf("*** If gsl-config was wrong, set the environment variable GSL_CONFIG\n");
      printf("*** to point to the correct copy of gsl-config, and remove the file\n");
      printf("*** config.cache before re-running configure\n");
      exit(1);
    }
}

],, no_gsl=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CFLAGS="$ac_save_CFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_gsl" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$GSL_CONFIG" = "no" ; then
       echo "*** The gsl-config script installed by GSL could not be found"
       echo "*** If GSL was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the GSL_CONFIG environment variable to the"
       echo "*** full path to gsl-config."
     else
       if test -f conf.gsltest ; then
        :
       else
          echo "*** Could not run GSL test program, checking why..."
          CFLAGS="$CFLAGS $GSL_CFLAGS"
          LIBS="$LIBS $GSL_LIBS"
          AC_TRY_LINK([
#include <stdio.h>
],      [ return 0; ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding GSL or finding the wrong"
          echo "*** version of GSL. If it is not finding GSL, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means GSL was incorrectly installed"
          echo "*** or that you have moved GSL since it was installed. In the latter case, you"
          echo "*** may want to edit the gsl-config script: $GSL_CONFIG" ])
          CFLAGS="$ac_save_CFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
#     GSL_CFLAGS=""
#     GSL_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(GSL_CFLAGS)
  AC_SUBST(GSL_LIBS)
  rm -f conf.gsltest
])


