AC_MSG_CHECKING([Checking for GMP])
SAVE_LIBS="$LIBS"
LIBS="$LIBS -lgmpxx -lgmp"
AC_TRY_LINK([#include <gmp.h>], 
        [__gmpz_init],
        [has_gmp=1],
        [has_gmp=0])
if test $has_gmp = 0; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
fi
LIBS="$SAVE_LIBS"