AC_MSG_CHECKING([Checking for PPL])
SAVE_LIBS="$LIBS"
LIBS="$LIBS -lppl"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([[#include <ppl.hh>]],[[Parma_Polyhedra_Library::version]])],
  [has_ppl=1],
  [has_ppl=0])
if test $has_ppl = 0; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
fi
LIBS="$SAVE_LIBS"