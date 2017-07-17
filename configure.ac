AC_INIT([dyngb], [2.0], [john.perry@usm.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
AC_PROG_CXX
# autoconf devs don't seem to have C++ configuration figured out...
# search_libs fails to include headers, while it does for C!?!?
# so pretty much everything from here on is hacked together
AC_LANG([C++])
#m4_include([ppl.m4])
#AM_PATH_PPL(1.1)
AC_CHECK_LIB([glpk], [glp_create_prob])
SAVE_LIBS="$LIBS"
LIBS="$LIBS -lppl"
m4_include([ac_checking_ppl.m4])
m4_include([ac_checking_gmp.m4])
m4_include([ax_require_defined.m4])
m4_include([ax_cxx_compile_stdcxx.m4])
AX_CXX_COMPILE_STDCXX(11)
# no way to check gmpxx yet...
#AC_SEARCH_LIBS([__gmpz_init], [gmp])
#AC_SEARCH_LIBS([version], [ppl])
#AC_CHECK_LIB([gmp])
#AC_CHECK_HEADER([gmpxx.h])
#AC_SEARCH_LIBS([__gmpz_init], [gmpxx])
AC_OUTPUT

