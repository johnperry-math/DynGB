# USAGE: make BUILD=[build | debug | profile]

#CXX = clang++
CXX = g++
#CXX = g++-mp-5
cxxflags.build = -Ofast -std=c++11 -I/usr/local/include -L/usr/local/lib
cxxflags.debug = -g -O0 -std=c++11 -I/usr/local/include -L/usr/local/lib
cxxflags.profile = -g -fno-inline -pg -O2 -std=c++11 -I/usr/local/include -L/usr/local/lib
CXXFLAGS = $(cxxflags.${BUILD})

MONOMIALS = monomial.o monomial_ordering.o particular_orderings.o monomial_ideal.o
POLYNOMIALS = indeterminate.o polynomial_ring.o polynomial.o
POLYNOMIAL_REPS = polynomial_array.o polynomial_double_buffered.o polynomial_geobucket.o polynomial_linked_list.o dense_univariate_integer_poly.o dense_univariate_rational_poly.o
STRATEGIES = strategies.o normal_strategy.o sugar_strategy.o wsugar_strategy.o
GBSUPPORT = critical_pair.o reduction_support.o $(STRATEGIES)
DYNAMIC = lp_solver.o skeleton.o glpk_solver.o dynamic_buchberger.o dynamic_engine.o ppl_solver.o
COMMALG = hilbert_functions.o betti.o
F4 = f4.o
NONDYNAMIC_OS = fields.o buchberger.o $(MONOMIALS) $(POLYNOMIALS) $(POLYNOMIAL_REPS) $(GBSUPPORT) $(COMMALG) $(F4)
ALL_OS = fields.o $(MONOMIALS) $(POLYNOMIALS) $(POLYNOMIAL_REPS) $(GBSUPPORT) $(COMMALG) $(DYNAMIC) $(F4)
EXPLORER = buchberger_explorer.o $(NONDYNAMIC_OS)

default: test_cyclicn test_hilbert test_incremental_betti test_dynamic test_monomials test_4by4 cab_es test_cyclicn_mpi test_glpk_solver test_initial_analysis user_interface test_f4

debug:
	make clean
	make BUILD=debug

build:
	make clean
	make BUILD=build

clean:
	rm *.o *.gch*
	rm -rf *.dSYM
	rm test_glpk_solver test_monomials test_hilbert test_incremental_betti
	rm test_cyclicn test_cyclicn_mpi test_4by4 test_dynamic test_cab_es? test_f4
	rm user_interface test_initial_analysis

cyclic4: $(ALL_OS)
	$(CXX) $(CXXFLAGS) -o test_cyclic4 $(NONDYNAMIC_OS) test_cyclic4.cpp

test_cyclicn: $(ALL_OS) buchberger.o dynamic_buchberger.o dynamic_engine.o skeleton.o buchberger_explorer_serial.o
	$(CXX) $(CXXFLAGS) -o test_cyclicn -lgmpxx -lgmp $(NONDYNAMIC_OS) buchberger_explorer_serial.o test_cyclicn.cpp

test_cyclicn_mpi: $(ALL_OS) buchberger_explorer.o
	#mpicxx $(CXXFLAGS) -o test_cyclicn_mpi $(EXPLORER) buchberger.o -lmpiP -lm -lbfd -liberty -lunwind test_cyclicn_mpi.cpp
	mpicxx $(CXXFLAGS) -o test_cyclicn_mpi -lgmp -lgmpxx $(EXPLORER) test_cyclicn_mpi.cpp

test_4by4: $(ALL_OS) buchberger.o dynamic_buchberger.o dynamic_engine.o skeleton.o test_4by4.cpp
	$(CXX) $(CXXFLAGS) -o test_4by4 -lglpk -lm -lgmpxx -lgmp -lppl $(ALL_OS) test_4by4.cpp buchberger.o

cab_es: $(ALL_OS) buchberger.o dynamic_buchberger.o dynamic_engine.o skeleton.o test_cab_es1 test_cab_es2 test_cab_es4 test_cab_es5 test_cab_es6 test_cab_es9

test_cab_es1: buchberger.o dynamic_buchberger.o dynamic_engine.o skeleton.o test_cab_es1.cpp
	$(CXX) $(CXXFLAGS) -o test_cab_es1 -lglpk -lm -lppl -lgmp -lgmpxx $(ALL_OS) test_cab_es1.cpp buchberger.o

test_cab_es2: buchberger.o dynamic_buchberger.o dynamic_engine.o skeleton.o test_cab_es2.cpp                                              
	$(CXX) $(CXXFLAGS) -o test_cab_es2 -lglpk -lm -lppl -lgmp -lgmpxx $(ALL_OS) test_cab_es2.cpp buchberger.o

test_cab_es4: buchberger.o dynamic_buchberger.o dynamic_engine.o skeleton.o test_cab_es4.cpp
	$(CXX) $(CXXFLAGS) -o test_cab_es4 -lglpk -lm -lppl -lgmp -lgmpxx $(ALL_OS) test_cab_es4.cpp buchberger.o

test_cab_es5: buchberger.o dynamic_buchberger.o dynamic_engine.o skeleton.o test_cab_es5.cpp
	$(CXX) $(CXXFLAGS) -o test_cab_es5 -lglpk -lm -lppl -lgmp -lgmpxx $(ALL_OS) test_cab_es5.cpp buchberger.o
                                     
test_cab_es6: buchberger.o dynamic_buchberger.o dynamic_engine.o skeleton.o test_cab_es6.cpp
	$(CXX) $(CXXFLAGS) -o test_cab_es6 -lglpk -lm -lppl -lgmp -lgmpxx $(ALL_OS) test_cab_es6.cpp buchberger.o

test_cab_es9: buchberger.o dynamic_buchberger.o dynamic_engine.o skeleton.o test_cab_es9.cpp
	$(CXX) $(CXXFLAGS) -o test_cab_es9 -lglpk -lm -lppl -lgmp -lgmpxx $(ALL_OS) buchberger.o test_cab_es9.cpp

buchberger.o: system_constants.hpp algorithm_buchberger_basic.hpp algorithm_buchberger_basic.cpp monomial.hpp strategies.hpp sugar_strategy.hpp weighted_sugar_strategy.hpp reduction_support.o
	$(CXX) $(CXXFLAGS) -c algorithm_buchberger_basic.hpp
	$(CXX) $(CXXFLAGS) -o buchberger.o -c algorithm_buchberger_basic.cpp

dynamic_engine.o: system_constants.hpp dynamic_engine.hpp dynamic_engine.cpp monomial.hpp monomial_ideal.hpp hilbert_functions.hpp betti.hpp polynomial.hpp lp_solver.o skeleton.o glpk_solver.o ppl_solver.o
	$(CXX) $(CXXFLAGS) -c dynamic_engine.hpp
	$(CXX) $(CXXFLAGS) -c dynamic_engine.cpp

dynamic_buchberger.o: system_constants.hpp algorithm_buchberger_dynamic.hpp algorithm_buchberger_dynamic.cpp monomial.o dynamic_engine.o skeleton.o hilbert_functions.o betti.o critical_pair.hpp sugar_strategy.hpp weighted_sugar_strategy.hpp reduction_support.o
	$(CXX) $(CXXFLAGS) -c algorithm_buchberger_dynamic.hpp
	$(CXX) $(CXXFLAGS) -o dynamic_buchberger.o -c algorithm_buchberger_dynamic.cpp

buchberger_explorer.o: system_constants.hpp algorithm_buchberger_basic.hpp algorithm_buchberger_basic.cpp monomial.hpp strategies.hpp algorithm_buchberger_explorer.hpp algorithm_buchberger_explorer.cpp hilbert_functions.hpp betti.o sugar_strategy.hpp weighted_sugar_strategy.hpp reduction_support.o
	$(CXX) $(CXXFLAGS) -c algorithm_buchberger_explorer.hpp
	mpicxx $(CXXFLAGS) -o buchberger_explorer.o -c algorithm_buchberger_explorer.cpp

buchberger_explorer_serial.o: system_constants.hpp algorithm_buchberger_basic.hpp algorithm_buchberger_basic.cpp monomial.hpp strategies.hpp algorithm_buchberger_explorer_serial.hpp algorithm_buchberger_explorer_serial.cpp hilbert_functions.hpp betti.o sugar_strategy.hpp weighted_sugar_strategy.hpp reduction_support.o
	$(CXX) $(CXXFLAGS) -c algorithm_buchberger_explorer_serial.hpp
	mpicxx $(CXXFLAGS) -o buchberger_explorer_serial.o -c algorithm_buchberger_explorer_serial.cpp

f4.o: f4_reduction.hpp f4_reduction.cpp system_constants.hpp fields.hpp monomial.hpp monomial_ordering.hpp polynomial.hpp critical_pair.hpp algorithm_buchberger_basic.hpp
	$(CXX) $(CXXFLAGS) -c f4_reduction.hpp
	$(CXX) $(CXXFLAGS) -o f4.o -c f4_reduction.cpp

fields.o: system_constants.hpp fields.hpp fields.cpp
	$(CXX) $(CXXFLAGS) -c fields.hpp
	$(CXX) $(CXXFLAGS) -c fields.cpp

indeterminate.o: polynomial_ring.hpp indeterminate.hpp indeterminate.cpp monomial.hpp
	$(CXX) $(CXXFLAGS) -c indeterminate.hpp
	$(CXX) $(CXXFLAGS) -c indeterminate.cpp

monomial_ordering.o: system_constants.hpp monomial_ordering.hpp monomial_ordering.cpp monomial.hpp
	$(CXX) $(CXXFLAGS) -c monomial_ordering.hpp
	$(CXX) $(CXXFLAGS) -c monomial_ordering.cpp

particular_orderings.o: system_constants.hpp monomial_ordering.o particular_orderings.hpp particular_orderings.cpp monomial.hpp goda.hpp indeterminate.o
	$(CXX) $(CXXFLAGS) -c goda.hpp
	$(CXX) $(CXXFLAGS) -c particular_orderings.hpp
	$(CXX) $(CXXFLAGS) -c particular_orderings.cpp

monomial.o: system_constants.hpp strategies.hpp monomial_ordering.o particular_orderings.o monomial.hpp monomial.cpp
	$(CXX) $(CXXFLAGS) -c monomial.hpp
	$(CXX) $(CXXFLAGS) -c monomial.cpp

polynomial_ring.o: system_constants.hpp polynomial_ring.hpp polynomial_ring.cpp indeterminate.o
	$(CXX) $(CXXFLAGS) -c polynomial_ring.hpp
	$(CXX) $(CXXFLAGS) -c polynomial_ring.cpp

polynomial.o: polynomial_ring.o system_constants.hpp polynomial.hpp polynomial.cpp monomial.hpp
	$(CXX) $(CXXFLAGS) -c polynomial.hpp
	$(CXX) $(CXXFLAGS) -c polynomial.cpp

polynomial_array.o: system_constants.hpp polynomial.o polynomial_array.hpp polynomial_array.cpp monomial.hpp
	$(CXX) $(CXXFLAGS) -c polynomial_array.hpp
	$(CXX) $(CXXFLAGS) -c polynomial_array.cpp

polynomial_double_buffered.o: system_constants.hpp polynomial.o polynomial_double_buffered.hpp polynomial_double_buffered.cpp monomial.hpp
	$(CXX) $(CXXFLAGS) -c polynomial_double_buffered.hpp
	$(CXX) $(CXXFLAGS) -c polynomial_double_buffered.cpp

polynomial_geobucket.o: system_constants.hpp polynomial.o polynomial_array.o polynomial_linked_list.o polynomial_geobucket.hpp polynomial_geobucket.cpp monomial.hpp
	$(CXX) $(CXXFLAGS) -c polynomial_geobucket.hpp
	$(CXX) $(CXXFLAGS) -c polynomial_geobucket.cpp

polynomial_linked_list.o: system_constants.hpp polynomial.o polynomial_linked_list.hpp polynomial_linked_list.cpp monomial.hpp
	$(CXX) $(CXXFLAGS) -c polynomial_linked_list.hpp
	$(CXX) $(CXXFLAGS) -c polynomial_linked_list.cpp

dense_univariate_integer_poly.o: system_constants.hpp dense_univariate_integer_poly.hpp dense_univariate_integer_poly.cpp
	$(CXX) $(CXXFLAGS) -c dense_univariate_integer_poly.hpp
	$(CXX) $(CXXFLAGS) -c dense_univariate_integer_poly.cpp

dense_univariate_rational_poly.o: system_constants.hpp dense_univariate_rational_poly.hpp dense_univariate_rational_poly.cpp
	$(CXX) $(CXXFLAGS) -c dense_univariate_rational_poly.hpp
	$(CXX) $(CXXFLAGS) -c dense_univariate_rational_poly.cpp

strategies.o: system_constants.hpp monomial.hpp strategies.hpp strategies.cpp
	$(CXX) $(CXXFLAGS) -c strategies.hpp
	$(CXX) $(CXXFLAGS) -c strategies.cpp

normal_strategy.o: system_constants.hpp strategies.o critical_pair.hpp critical_pair.cpp normal_strategy.hpp normal_strategy.cpp
	$(CXX) $(CXXFLAGS) -c normal_strategy.hpp
	$(CXX) $(CXXFLAGS) -c normal_strategy.cpp

lp_solver.o: system_constants.hpp monomial.hpp lp_solver.hpp lp_solver.cpp
	$(CXX) $(CXXFLAGS) -c lp_solver.hpp
	$(CXX) $(CXXFLAGS) -c lp_solver.cpp

skeleton.o: system_constants.hpp monomial.hpp monomial.cpp polynomial.hpp polynomial.cpp skeleton.hpp skeleton.cpp lp_solver.hpp
	$(CXX) $(CXXFLAGS) -c skeleton.hpp
	$(CXX) $(CXXFLAGS) -c skeleton.cpp

glpk_solver.o: system_constants.hpp monomial.hpp monomial.cpp polynomial.hpp polynomial.cpp lp_solver.hpp glpk_solver.hpp glpk_solver.cpp
	$(CXX) $(CXXFLAGS) -c glpk_solver.hpp
	$(CXX) $(CXXFLAGS) -c glpk_solver.cpp

ppl_solver.o: system_constants.hpp monomial.hpp monomial.cpp polynomial.hpp polynomial.cpp lp_solver.hpp ppl_solver.hpp ppl_solver.cpp
	$(CXX) $(CXXFLAGS) -c ppl_solver.hpp
	$(CXX) $(CXXFLAGS) -c ppl_solver.cpp

sugar_strategy.o: system_constants.hpp strategies.o critical_pair.hpp critical_pair.cpp sugar_strategy.hpp sugar_strategy.cpp
	$(CXX) $(CXXFLAGS) -c sugar_strategy.hpp
	$(CXX) $(CXXFLAGS) -c sugar_strategy.cpp

wsugar_strategy.o: system_constants.hpp strategies.o critical_pair.hpp critical_pair.cpp weighted_sugar_strategy.hpp weighted_sugar_strategy.cpp
	$(CXX) $(CXXFLAGS) -c weighted_sugar_strategy.hpp
	$(CXX) $(CXXFLAGS) -c -o wsugar_strategy.o weighted_sugar_strategy.cpp

critical_pair.o: system_constants.hpp strategies.o normal_strategy.o wsugar_strategy.o sugar_strategy.o polynomial.hpp polynomial_linked_list.hpp polynomial_array.hpp polynomial_geobucket.hpp polynomial_double_buffered.hpp polynomial.cpp polynomial_linked_list.cpp polynomial_array.cpp polynomial_double_buffered.cpp polynomial_geobucket.cpp
	$(CXX) $(CXXFLAGS) -c critical_pair.hpp
	$(CXX) $(CXXFLAGS) -c critical_pair.cpp

reduction_support.o: system_constants.hpp strategies.o polynomial.hpp polynomial.cpp reduction_support.hpp reduction_support.cpp
	$(CXX) $(CXXFLAGS) -c reduction_support.hpp
	$(CXX) $(CXXFLAGS) -c reduction_support.cpp

monomial_ideal.o: system_constants.hpp monomial.hpp monomial_ideal.cpp hilbert_functions.hpp betti.o polynomial.hpp skeleton.hpp
	$(CXX) $(CXXFLAGS) -c monomial_ideal.hpp
	$(CXX) $(CXXFLAGS) -c monomial_ideal.cpp

hilbert_functions.o: system_constants.hpp hilbert_functions.hpp hilbert_functions.cpp monomial.o fields.o polynomial_linked_list.o monomial_ideal.hpp dense_univariate_integer_poly.o dense_univariate_rational_poly.o
	$(CXX) $(CXXFLAGS) -c hilbert_functions.hpp
	$(CXX) $(CXXFLAGS) -c hilbert_functions.cpp

betti.o: system_constants.hpp betti.hpp betti.cpp monomial.o
	$(CXX) $(CXXFLAGS) -c betti.hpp
	$(CXX) $(CXXFLAGS) -c betti.cpp

test_hilbert: test_hilbert_functions.cpp hilbert_functions.o particular_orderings.o monomial.o monomial_ordering.o dense_univariate_integer_poly.o dense_univariate_rational_poly.o monomial_ideal.o
	$(CXX) $(CXXFLAGS) -o test_hilbert -lgmpxx -lgmp test_hilbert_functions.cpp $(NONDYNAMIC_OS)

test_incremental_betti: test_incremental_betti.cpp betti.o particular_orderings.o monomial.o monomial_ordering.o monomial_ideal.o hilbert_functions.o
	$(CXX) $(CXXFLAGS) -o test_incremental_betti -lgmp -lgmpxx test_incremental_betti.cpp $(NONDYNAMIC_OS)

test_dynamic: $(ALL_OS) monomial_ideal.o dynamic_engine.o dynamic_buchberger.o test_dynamic.cpp
	$(CXX) $(CXXFLAGS) -o test_dynamic -lglpk -lm -lgmp -lgmpxx -lppl test_dynamic.cpp  $(ALL_OS) buchberger.o

test_f4: $(ALL_OS) test_f4.cpp f4.o
	$(CXX) $(CXXFLAGS) -o test_f4 test_f4.cpp critical_pair.o fields.o strategies.o polynomial.o polynomial_array.o polynomial_linked_list.o particular_orderings.o monomial.o monomial_ordering.o polynomial_ring.o sugar_strategy.o buchberger.o wsugar_strategy.o reduction_support.o indeterminate.o polynomial_geobucket.o normal_strategy.o polynomial_double_buffered.o

test_monomials: monomial.o monomial_ordering.o particular_orderings.o test_monomials.cpp polynomial_ring.o
	$(CXX) $(CXXFLAGS) -o test_monomials test_monomials.cpp monomial.o monomial_ordering.o particular_orderings.o fields.o polynomial_ring.o indeterminate.o

test_glpk_solver: test_glpk_solver.cxx lp_solver.o glpk_solver.o monomial.o
	$(CXX) $(CXXFLAGS) -o test_glpk_solver -lglpk -lm -lgmp -lgmpxx -lppl lp_solver.o glpk_solver.o monomial.o particular_orderings.o monomial_ordering.o test_glpk_solver.cxx

test_initial_analysis: test_initial_analysis.cpp dynamic_buchberger.o monomial.o
	$(CXX) $(CXXFLAGS) -o test_initial_analysis -lglpk -lm -lgmp -lgmpxx -lppl buchberger.o lp_solver.o glpk_solver.o monomial.o particular_orderings.o monomial_ordering.o dynamic_buchberger.o hilbert_functions.o dynamic_engine.o polynomial_ring.o reduction_support.o ppl_solver.o betti.o fields.o strategies.o sugar_strategy.o monomial_ideal.o critical_pair.o polynomial_array.o indeterminate.o normal_strategy.o wsugar_strategy.o polynomial.o polynomial_geobucket.o dense_univariate_integer_poly.o dense_univariate_rational_poly.o polynomial_double_buffered.o polynomial_linked_list.o skeleton.o test_initial_analysis.cpp

user_interface: $(ALL_OS) user_interface.cpp
	$(CXX) $(CXXFLAGS) -o user_interface -g -std=c++11 -lppl -lgmp -lgmpxx -lglpk fields.o polynomial_ring.o monomial.o polynomial.o polynomial_array.o strategies.o indeterminate.o reduction_support.o lp_solver.o ppl_solver.o glpk_solver.o skeleton.o dynamic_engine.o hilbert_functions.o sugar_strategy.o monomial_ordering.o particular_orderings.o betti.o dense_univariate_rational_poly.o dense_univariate_integer_poly.o monomial_ideal.o wsugar_strategy.o critical_pair.o normal_strategy.o polynomial_geobucket.o polynomial_double_buffered.o polynomial_linked_list.o dynamic_buchberger.o user_interface.cpp buchberger.o