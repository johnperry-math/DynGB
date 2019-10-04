# DOXYGEN SUPPORT
include $(top_srcdir)/aminclude.am

# ensure the distribution of the doxygen configuration file
EXTRA_DIST = doxygen.cfg

#bin_PROGRAMS = test_cyclicn_mpi
bin_PROGRAMS = test_monomials test_hilbert_functions test_incremental_betti test_glpk_solver test_gmp test_cyclicn test_dynamic test_initial_analysis user_interface test_4by4 test_cab_es1 test_cab_es2 test_cab_es4 test_cab_es5 test_cab_es6 test_cab_es9 test_f4 test_f4_dynamic test_hash

DYNAMIC_LIB = algorithm_buchberger_dynamic.cpp dynamic_engine.cpp
F4_LIB = f4_dynamic.cpp f4_hash.cpp monomial_ideal.cpp betti.cpp
HILBERT_LIB = hilbert_functions.cpp dense_univariate_rational_poly.cpp dense_univariate_integer_poly.cpp
LP_SOLVERS_LIB = lp_solver.cpp skeleton.cpp glpk_solver.cpp ppl_solver.cpp
MONOMIAL_LIB = monomial.cpp indeterminate.cpp
ORDERING_LIB = monomial_ordering.cpp particular_orderings.cpp
POLY_LIB = fields.cpp polynomial_ring.cpp polynomial.cpp polynomial_array.cpp polynomial_double_buffered.cpp polynomial_geobucket.cpp polynomial_linked_list.cpp
REDUCTION_LIB = reduction_support.cpp
STATIC_LIB = algorithm_buchberger_basic.cpp critical_pair.cpp
STRAT_LIB = strategies.cpp sugar_strategy.cpp normal_strategy.cpp weighted_sugar_strategy.cpp

ALL_DYNAMIC_LIBS = $(DYNAMIC_LIB) $(F4_LIB) $(HILBERT_LIB) $(LP_SOLVERS_LIB) $(ORDERING_LIB) $(MONOMIAL_LIB) $(POLY_LIB) $(REDUCTION_LIB) $(STATIC_LIB) $(STRAT_LIB)

GMP_LINK = -lgmpxx -lgmp
DYNAMIC_LINK = $(GMP_LINK) -lppl -lglpk
PARALLEL_LINK = -lpthread -pthread

test_monomials_SOURCES = test_monomials.cpp polynomial_ring.cpp fields.cpp $(MONOMIAL_LIB) $(ORDERING_LIB)

test_hilbert_functions_SOURCES = test_hilbert_functions.cpp $(HILBERT_LIB) $(MONOMIAL_LIB) $(ORDERING_LIB) polynomial_linked_list.cpp polynomial.cpp fields.cpp monomial_ideal.cpp
test_hilbert_functions_LDFLAGS = $(GMP_LINK)

test_incremental_betti_SOURCES = test_incremental_betti.cpp $(MONOMIAL_LIB) $(ORDERING_LIB) betti.cpp

test_glpk_solver_SOURCES = test_glpk_solver.cpp $(LP_SOLVERS_LIB) $(MONOMIAL_LIB) $(ORDERING_LIB)
test_glpk_solver_LDFLAGS = $(DYNAMIC_LINK)

test_gmp_SOURCES = test_gmp.cpp
test_gmp_LDFLAGS = $(GMP_LINK)

test_cyclicn_SOURCES = test_cyclicn.cpp $(HILBERT_LIB) $(ORDERING_LIB) $(MONOMIAL_LIB) $(POLY_LIB) $(REDUCTION_LIB) $(STATIC_LIB) $(STRAT_LIB) algorithm_buchberger_explorer_serial.cpp f4_reduction.cpp f4_hash.cpp monomial_ideal.cpp
test_cyclicn_LDFLAGS = $(GMP_LINK) $(PARALLEL_LINK)

test_dynamic_SOURCES = test_dynamic.cpp $(DYNAMIC_LIB) $(HILBERT_LIB) $(LP_SOLVERS_LIB) $(ORDERING_LIB) $(MONOMIAL_LIB) $(POLY_LIB) $(REDUCTION_LIB) $(STATIC_LIB) $(STRAT_LIB) monomial_ideal.cpp betti.cpp
test_dynamic_LDFLAGS = $(DYNAMIC_LINK)

test_initial_analysis_SOURCES = test_initial_analysis.cpp $(ALL_DYNAMIC_LIBS)
test_initial_analysis_LDFLAGS = $(DYNAMIC_LINK)

user_interface_SOURCES = user_interface.cpp $(ALL_DYNAMIC_LIBS)
user_interface_LDFLAGS = $(DYNAMIC_LINK) $(PARALLEL_LINK)

test_4by4_SOURCES = test_4by4.cpp $(ALL_DYNAMIC_LIBS)
test_4by4_LDFLAGS = $(DYNAMIC_LINK)

test_cab_es1_SOURCES = test_cab_es1.cpp $(ALL_DYNAMIC_LIBS)
test_cab_es1_LDFLAGS = $(DYNAMIC_LINK)

test_cab_es2_SOURCES = test_cab_es2.cpp $(ALL_DYNAMIC_LIBS)
test_cab_es2_LDFLAGS = $(DYNAMIC_LINK)

test_cab_es4_SOURCES = test_cab_es4.cpp $(ALL_DYNAMIC_LIBS)
test_cab_es4_LDFLAGS = $(DYNAMIC_LINK)

test_cab_es5_SOURCES = test_cab_es5.cpp $(ALL_DYNAMIC_LIBS)
test_cab_es5_LDFLAGS = $(DYNAMIC_LINK)

test_cab_es6_SOURCES = test_cab_es6.cpp $(ALL_DYNAMIC_LIBS)
test_cab_es6_LDFLAGS = $(DYNAMIC_LINK)

test_cab_es9_SOURCES = test_cab_es9.cpp $(ALL_DYNAMIC_LIBS)
test_cab_es9_LDFLAGS = $(DYNAMIC_LINK)

test_hash_SOURCES = $(MONOMIAL_LIB) $(ORDERING_LIB) f4_hash.cpp test_hash.cpp

test_f4_SOURCES = test_f4.cpp $(HILBERT_LIB) $(ORDERING_LIB) $(MONOMIAL_LIB) $(POLY_LIB) $(REDUCTION_LIB) $(STATIC_LIB) $(STRAT_LIB) f4_reduction.cpp f4_hash.cpp monomial_ideal.cpp
test_f4_LDFLAGS = $(GMP_LINK) $(PARALLEL_LINK)

test_f4_dynamic_SOURCES = test_f4_dynamic.cpp  $(DYNAMIC_LIB) $(HILBERT_LIB) $(LP_SOLVERS_LIB) $(ORDERING_LIB) $(MONOMIAL_LIB) $(POLY_LIB) $(REDUCTION_LIB) $(STATIC_LIB) $(STRAT_LIB) $(F4_LIB)
test_f4_dynamic_LDFLAGS = $(DYNAMIC_LINK) $(GMP_LINK) $(PARALLEL_LINK)
