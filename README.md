# DynGB
Testbed for exploring dynamic algorithms to compute a Groebner basis

## Contents

   - [`Building`](#building)
   - [`Usage`](#usage)
   - [`Deciphering the output`](#deciphering-the-output)
   - [`Defining new input systems`](#defining-new-input-systems)
   - [`Some results`](#some-results)
   - [`Additional details`](#additional-details)

## Building

Depending on what I'm working on at the moment, some older stuff may not actually work. I do try to keep things in order.

### Prerequisites

You will need

   - a C++ compiler for C++14, with standard libraries
   - the [Gnu Multiple Precision Arithmetic library](https://gmplib.org)
   - the [Gnu Linear Programming Kit](https://www.gnu.org/software/glpk/)
   - the [Parma Polyhedra Library](https://www.bugseng.com/parma-polyhedra-library)

### Most systems

To build this, you will need the `autotools` libraries. I build this on several different Unix-based workstations using something along these lines.

    autoreconf
    automake --add-missing
    autoconf
    mkdir debug # unless it already exists
    cd debug
    ../configure
    make -j 6 some_executable

Details on `some_executable` are below, in the section [on usage](#usage).

This assumes that you have a `debug` folder; if not, you have to make it. The `autoreconf`, `automake`, and `autoconf` commands may not be necessary, depending on how agreeable your system finds the last one I pushed from. (I'm afraid I haven't quite figured out autotools.)

### macOS

If you are on macOS and installed PPL, GLPK, etc. via ports or brew or some other such system, you may need to specify the include and link directories:

    ../configure LDFLAGS='-L/opt/local/lib' CPPFLAGS='-I/opt/local/include'

Go ahead and specify compiler flags such as `-Ofast` etc. with `CPPFLAGS` and `CXXFLAGS`. For instance, I typically type out

    ../configure LDFLAGS='-L/opt/local/lib' CPPFLAGS='-I/opt/local/include -Ofast' CXXFLAGS='-I/opt/local/include -Ofast'

unless I want to debug, in which case I type out

    ../configure LDFLAGS='-L/opt/local/lib' CPPFLAGS='-I/opt/local/include -O0 -g' CXXFLAGS='-I/opt/local/include -O0 -g'

Now (re-)invoke `make`.

## Usage

### Executables that should build

For `some_executable` I generally make sure the following build:

  - [`test_cyclicn`](#test_cyclicn): an interface to a static Buchberger algorithm that runs on the Cyclic-n systems
  - [`test_dynamic`](#test_dynamic): an interface to a dynamic Buchberger algorithm that runs on the Cyclic-n systems
  - [`test_f4_dynamic`](#test_f4_dynamic): an interface to an F4 algorithm, either static or dynamic, that runs on the Cyclic-n systems
  - [`user_interface`](#user_interface): a program that runs any of the above

Details on using each follows.

### `test_cyclicn`

Options are:

  - `m=...` or `mod=...` or `modulus=...` where `...` is a prime integer; this is the modulus, or ring characteristic. Defaults to 43.
  - `n=...` or `num=...` or `numvars=...` where `...` is a positive integer; this is the number of variables.
  - `hom` or `homog` or `homogeneous` indicates that the system should be homogeneous. Leave it out, and it won't be.
  - `r=...` or `repr=...` or `representation=...` where `...` indicates what data structure an s-polynomial should have:
    - `1` indicates linked lists (very slow)
    - `2` indicates geobuckets (the default, slow)
    - `3` indicates [double-buffered polynomials](#double-buffered?) (pretty slow)
  - `ord=...` or `order=...` or `ordering=...` where `...` is
    - `generic` indicates a generic graded reverse lexicographic order (the default)
    - `grevlex` indicates the graded reverse lexicographic order
    - `lex` indicates the lexicographic ordering
    - `wgrevlex` indicates a graded lexicographic ordering with a weight vector; follow this by `n` positive integers
  - `strat=...` or `strategy=...` where `...` is
    - `normal` or `norm` indicates the normal strategy
    - `sugar` or `sug` indicates a sugar strategy (the default)
    - `wsugar` `wsug` indicates a weighted sugar strategy; follow this by `n` positive integers

The wisest options to use are the defaults. Other options might not work if I haven't tested them in a while, and the question(s) for which I implemented them are of the distant past.

#### Double-buffered?

The origin of all this code is in a class project from some years back. I looked into how "double-buffered arrays" might work as a representation for s-polynomials. The short answer is that they don't work much better than linked lists. I've left that in here, but it's probably not worth playing with.

#### Example

To compute a Groebner basis of the homogenized Cyclic-7 system over a finite field of characteristic 32003, use

    test_cyclicn m=32003 n=7 homog

### `test_dynamic`

Options are:

  - Accepts `test_cyclicn`'s options, except for the ordering (`ord`), which it shouldn't anyway.
  - `heur=...` or `heuristic=...`, where `...` chooses the leading term; this can be
    - `1` indicates the (standard) Hilbert heuristic, with ties broken by lex order on the monomials
    - `2` indicates the (standard) Hilbert heuristic, with ties broken by the monomial degree (the default)
    - `3` indicates the graded Hilbert heuristic, with ties broken by the monomial degree
    - `4` indicates the smallest number of critical pairs heuristic
    - `5` indicates the smallest number of critical pairs heuristic, but graded by the ordering
    - `6` indicates the Betti heuristic, with ties broken by the Hilbert function, and ties broken by the monomial degree
    - `7` indicates the "Big Betti" heuristic
    - `8` indicates the "Big Betti" heuristic, with ties broken by the Hilbert function, and ties broken by the monomial degree
    - `9` indicates the "smoothest degrees" heuristic
    - `10` indicates the "largest maximum component" heuristic
    - `11` indicates the heuristic of degree first, then the graded Hilbert function
    - `12` indicates the heuristic of degree first, then the (standard) Hilbert function
    - `13` indicates a random choice of monomial
  - `solver=...`, where `...` is
    - `skel` or `skeleton` indicates an implementation of Zolotykh's Skeleton algorithm (the default)
    - `glpk` or `GLPK` indicates the Gnu Linear Programming Kit
    - `ppl` or `PPL` indicates the Parma Polyhedra Library

The wisest option is *usually* the default, but PPL should work. Be advised that I haven't used GLPK in a while, and a quick test (April 2020) suggests that it's broken on some systems.

#### Example

To compute a Groebner basis of the homogenized Cyclic-7 system over a finite field of characteristic 32003, using the Parma Polyhedra Library to find border vectors, use

    test_dynamic m=32003 n=7 homog solver=ppl

### `test_f4_dynamic`

(Despite the name, this can also run a static F4 algorithm.)

Options are:

  - `m=...` or `mod=...` or `modulus=...` where `...` is a prime integer; this is the modulus, or the ring characteristic. Defaults to 43.
  - `n=...` or `num=...` or `numvars=...` where `...` is a positive integer; this is the number of variables.
  - `hom` or `homog` or `homogeneous` indicates that the systems should be homogenized. Default is `false`.
  - `static` directs the program to compute with the static F4 algorithm. Default is `false`.
  - `analyze ...` applies only in a dynamic context, where `...` is one of
    - `row` (the default) analyzes only one row of the F4 matrix before selecting a leading term
    - `matrix` analyzes the entire F4 matrix before selecting a leading term
  - `r=...` or `refinements=...` where `...` is a positive integer; this limits the number of refinements that the matrix analysis can perform before giving up on a matrix and completing it statically. The default is `0`, which means no limit.

#### Example

To compute a Groebner basis of the homogenized Cyclic-7 system over a finite field of characteristic 32003, use

    test_f4_dynamic m=32003 n=7 homog


### `user_interface`

This prompts for data interactively. I generally prepare the examples in a text file and redirect the input to the text file. ([Directions to do that are available below.](#defining-new-input-systems)) You can find a number of example systems in the directory `examples_for_user_interface`. The suffixes on each system (before `.txt`) indicate which algorithm they use:

  - no suffix uses the static Buchberger algorithm
  - `_f4` uses the static F4 algorithm
  - `_df4m` uses the dynamic F4 algorithm with whole-matrix analysis
  - `_df4r` uses the dynamic F4 algorithm with sequential-row analysis

#### Example

To compute a Groebner basis for the Katsura-10 system using whole-matrix analysis, use

    ./user_interface < ../example_systems_for_user_interface/katsura10_df4m.txt

## Defining new input systems

To specify your own input systems, go through [`user_interface`](#user_interface). It's possible to use the program interactively, but I find it easier to prepare a text file and pipe it to the program. If you define  `a_very_interesting_system.txt`, you can just run

    user_interface < .../a_very_interesting_system.txt

where `...` specifies the path to `a_very_interesting_system.txt`.

What you place in the file depends on whether you want a static, dynamic, Buchberger-style, or F4-style algorithm. Use the files in `examples_for_user_interface` as a guide.

   - All formats:
     1. For all of them, the first line indicates the reduction style: `buchberger` or `f4`.
     2. The second line indicates the field size. This must be a prime number; I have not implemented the rational field or finite fields of non-prime size.
     3. The third line indicates the number of indeterminates. This should be a positive number greater than 1. (1 might work, but why are you using Groebner bases?)
     4. The fourth line indicates whether you want to define the indeterminates' names. Specify `y` or `n`. 
     5. **If you specify yes**, list the indeterminates on the next line, separated by spaces.
     **Otherwise,** the program will assign indeterminates on its own (`x0`, `x1`, etc.) so you can move to the next step.
     7. The fifth or sixth line indicates the number of the ideal's generators. We'll call this number `m`; it must be a positive number greater than 1. (1 might work, but why are you using Groebner bases?)
     8. The next `m` lines are the ideal's generators. Specify each as a sum of terms, separating spaces (or asterisks) between factors, using carets to indicate exponents.
     **Tip:** The code is not made to accept decimals. If the generators contain decimals, multiply each by the power of 10 necessary to eliminate decimals. For example, `x + 1.1y` becomes `10x+11y`.
     9. The next line expects you to indicate whether you want a static algorithm or a dynamic algorithm. Indicate `s` or `d`.
   - The input for a static algorithm is now complete; the remaining options apply only to dynamic algorithms. You now have two routes, depending on whether your first line selected Buchberger- or F4-style.
     - For F4-style, indicate whether you want to use the DF4M variant (`matrix`) or the DF4R variant (`row`).
       - The input for DF4R is complete.
       - For DF4M, you may set a maximum number of refinements per any matrix. The default is 0, which allows any matrix to re-examine *after each row reduction* all the rows that have not yet been used in refinement. **This can take a very long time**, so setting this to a positive number `r` limits any matrix to `r` refinements, after which it completes the matrix reduction by static computation. The next matrix will renew the attempt to compute dynamically, so it still takes a pretty long time, just not quite as long in some systems. I don't found this option useful.
     - For Buchberger-style, the next lines indicate:
       1. The solver to use. Your choices are `skel` (an implementation of Zolotykh's Skeleton algorithm), `ppl` (Parma Polyhedra Library), or `glpk` (Gnu Linear Programming Kit). Recently I have relied on `ppl`, which computes exact cones somewhat more quickly than `skel`. *April 2020:* It seems that `glpk` does not work at the present time on some systems.
       1. The heuristic to use: (explained in the 2017 ISSAC paper)
          - `h` is Hilbert heuristic with standard grading, ties broken by degree
          - `gh` is Hilbert heuristic with grading weighted by the ordering, ties broken by degree
          - `b` is Betti heuristic, with ties broken by Hilbert, then by degree (I believe this is broken at the moment, or perhaps just slow)
          - `bb` is Big Betti heuristic, with ties broken by Hilbert, then by degree (I believe this is broken at the moment, or perhaps just slow)
          - `gb` is graded Betti heuristic, with ties broken by Hilbert, then by degree (I believe this is broken at the moment, or perhaps just slow)
          - `c` is minimal critical pairs heuristic
          - `d` is degree heuristic, with ties broken by Hilbert heuristic
          - `r` is the Evil Random heuristic, which chooses a compatible leading term at random and is a terrible choice except in systems with very small Groebner fans (where it can sometimes match or beat the others!)
       1. Whether to analyze the inputs first (`y` or `n`) to start with a better ordering than would an incremental computation that considers the polynomials by increasing degree. This might not work right now, but it did once.

## Deciphering the output

   - Upon launch:
     - The program will first print the system it's computing a Groebner basis for.
       **Tip:** Make sure it's the system you expect. I have sometimes discovered something doesn't work as expected because I didn't specify the input correctly, or there was a bug in the way `user_interface` read it.
   - When running:
     - All algorithms regularly report the number of critical pairs remaining, as well as the leading monomials of the polynomials generating each pair, and the degree of the current pair.
     - F4-style algorithms will report the matrix size and which row has most recently completed reduction.
     - Dynamic algorithms will report a change in ordering and which leading monomial is selected.
     - Dynamic algorithms also check to make sure that certain conditions are met.
       - The leading monomial selected is checked when the linear program is returned and solved. If this fails, an error message will appear, along the lines of `ERROR HERE` or `Did not order polynomial`. If you see this message, **this is a bug**; please report it, and include output if you can. This should not occur now, because I disabled an optimization in DF4M that causes occasional, random failure in some systems until I investigate and fix it.
       - After every new selection, the algorithm ensures that previous choices of leading terms remain constant. This will print some output along the lines of `is consistent?`, `fails again`, `recovered with...` **This is normal:** a desired leading term may or may not work on account of past choices. This sequence of questions and answers reflects the algorithm's investigation. (See the 2014 AAECC paper's discussion on adding constraints.)
     - When the system detects multiple cores, it will use them via the C++ `threads` library, and report the launching and joining of threads.
   - Upon termination:
     - The system will report basis size, leading terms of each polynomial, number of s-polynomials computed, number of monomials in basis (probably *not* distinct monomials), and internally-tracked execution time for several subprograms.
     - F4-style algorithms will report the number of row reductions performed, the time spent creating matrices, and the time spent row-reducing the matrices. They will also report the time spent computing critical pairs.
     - Dynamic algorithms will report the amount of time spent in dynamic overhead, which consists of analyzing terms and selecting a leading term. If additional divisibility criteria to detect useless pairs are enabled there will also be a non-zero report of those criteria. (As of mid-April 2020, they are disabled because they're slow.)

## Some results

The `results` directory contains the output of some larger and/or more interesting systems I have used. At the present time, interpret "more interesting" as, "non-deterministic results with DF4M." Some highlights: (all characteristic 32003 unless stated otherwise)

date | computation type | output size | time | mem | remarks
-|-|-|-|-|-|-
**Cyclic-9 homogeneous**|
14 April | DF4R | 1991 | 3975sec (1h 5m) | 2.6GB | 
14 April | DF4M | 2005 | 6050sec (1h 41m) | 4.0GB | 
14 April | Static | 5602 | 12757sec (3h 32m) | 3.6GB | 
**Cyclic-10 homogeneous** |
27 December | DF4R | 14200 | 12d 3h 58m | 166GB | characteristic 43, USM Magnolia cluster
27 December | Static | ? | ? | ? | characteristic 43, USM magnolia cluster, not reported in paper ,static still incomplete after 3+ months
**Eco-8** | | | | | example of nondeterminism in DF4M
8 April | DF4M | 33 (x4) | ~1sec | 93MB |
8 April | DF4M | 43 | ~1sec | 93MB |
8 April | DF4M | 38 | ~1sec | 93MB |
**Katsura-10** | | | | | example of nondeterminism in DF4M
9 April | DF4M | 423 | 1646sec (27m) | 345MB |
9 April | DF4M | 423 | 1399sec (23m) | 358MB |
9 April | DF4M | 423 | 419sec (7m) | 334MB |
**Kotsireas**  | | | | | example of nondeterminism in DF4M
8 April | DF4M | 17 | 3sec | 91MB |
8 April | DF4M | 17 | 3sec | 92MB |
8 April | DF4M | 18 | 2sec | 89MB |
8 April | DF4M | 17 | 4sec | 97MB |
8 April | DF4M | 17 | 3sec | 110MB |
8 April | DF4M | 17 | 4sec | 95MB |

## Additional details

The system includes documentation; I tried to be very generous with the comments, but didn't always succeed. You can extract it if you have installed [`doxygen`](http://www.doxygen.org/). Generate the documentation by running the script `make_documentation`.

In theory, `automake` should create the documentation for you, but for some reason the auto tools script I downloaded doesn't seem to work, which probably means I didn't read the fine print carefully enough. Nevertheless, installing doxygen and running `make_documentation` work for me.
