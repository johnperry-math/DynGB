/**
  \mainpage DynGB: Dynamic Gr&ouml;bner basis project files

  \author John Perry, john.perry@usm.edu
  \date 2014-present
  \copyright <a href="../../COPYING.txt">GNU Public License</a>; see \ref Copyright below

  \tableofcontents

  \section Overview Overview

  This project is designed to be a testbed/reference implementation
  for dynamic Gr&ouml;bner basis computation,
  using the algorithms described in \cite CaboaraDynAlg and \cite CaboaraPerry,
  along with some newer ideas.

  No claim to high efficiency or exemplary programming is implied.
  I wrote this to be relatively usable (compared to an earlier Sage implementation)
  and <i>easy to modify</i>, especially as regards modularity, polymorphism,
  and getting detailed data.

  \section Install Installation and dependencies

  For @em this program, a simple `make BUILD=build` should do.
  However, there are several prerequisite programs you need to install first:
    - a C++ compiler that understands C++11; 
    - <a href="https://gmplib.org/" target="_blank">GMP</a>,
      the Gnu Multi-Precision library \cite gmp;
    - <a href="https://www.gnu.org/software/glpk/" target="_blank">GLPK</a>,
      the Gnu Linear Programming Kit \cite glpk;
    - <a href="http://bugseng.com/products/ppl" target="_blank">PPL</a>,
      the Parma Polyhedra Library \cite BagnaraHZ08SCP.

  On a Linux system you can install these very easily (in Fedora I used `Apper`).
  All three should build without difficulty on a Macintosh.
  I have not tried to build on Windows, but I don&rsquo;t use hardware magic,
  so it @em ought to build.

  \section Usage Usage

  There are two ways to use the system.
    -# Write a program that builds a polynomial system and accesses the library
       directly. Most of the <a href="examples.html">examples</a> provided
       illustrate that approach. This is annoying and a rather difficult,
       though not as difficult as it was before the Indeterminate class.
    -# Use the `user_interface` file. A number of systems are defined in the
       directory `examples_for_user_interface`. The format is defined in the
       documentation to user_interface(). This is still annoying,
       but not quite so difficult.

  \section Status Current status

  As of January 2017:
  \li The code works consistently on many different examples.
    However, it is slow: I am not trying to reach
    <span style="font-variant:small-caps;">Singular</span>-level optimization
    (not at the current time, anyway). Typically, this code is one to two orders
    of magnitude slower than <span style="font-variant:small-caps;">Singular</span>.
  \li Nevertheless, it outperforms <span style="font-variant:small-caps;">Singular</span>
    on at least one system: Caboara&rsquo;s Example 2.
  \li Unless I&rsquo;m doing something very stupid, the weighted sugar strategy
    is an unmitigated disaster and should be avoided.
  \li The code is slow, though in the last week of July the implementation of an
    \f$O(1)\f$ memory manager cut the time required for the dynamic
    implementation by nearly 2/3. A very simple optimization of assigning an
    object&rsquo;s array to a local variable before entering a loop cut the time
    required for both dynamic and static by roughly 40%.

  @note Exponent-packing is currently turned off,
    since it doesn&rsquo;t seem to help enough to make it worth the trouble.

  @warning
    I implemented the exponent-packing in rather shoddy fashion:
    the first 8 exponents are cast to \c uint8_t , then shifted appropriately.
    Comparisons for \em other are made explicitly. Trouble will arise when one of
    the first 8 exponents exceeds \f$2^8-1\f$, though practically that
    hasn&rsquo;t been a problem until now. This may be easy to fix: if the
    exponent comparison passes, test all the variables explicitly, not just those
    after the first 8. But there are issues with arithmetic operations that
    could arise, as well; the multiplication and division operators are
    implemented semi-intelligently; that is, intelligently under the assumption
    that the exponents are valid. So problems could arise there even in the case
    where we fix the equality comparison.

  \section Todo To-do list
  In no particular order, aside from the indicated priority.
  See the <a href="todo.html">to-do page</a> for a full list
  (I may have missed some things).

  \subsection hipri Higher priority
  \todo These are the highest priority items:
  things I think are needed before I&rsquo;d call it &ldquo;ready.&rdquo;
  - Organize files into directories.
  - General improvements to efficiency based on profiling. (ongoing)
  - Implement simplex solver as oracle for DDM, compare with DDM
    (idea due to D. Lichtblau).
  - Optimize length() in Polynomial_Linked_List.
  - Add Fukuda and Prodon&rsquo;s cdd as an LP_Solver. \cite Fukuda_DoubleDescriptionRevisited
  - Bring polynomial iterators in line with C++ convention.
  - Implement other C++11 modernizations.
  - Generalize/improve the memory manager.
  - <span style="text-decoration:line-through;">Add PPL as an LP_Solver. \cite BagnaraHZ08SCP</span>
  - <span style="text-decoration:line-through;">Implement Caboara&rsquo;s examples.</span>
  - <span style="text-decoration:line-through;">Implement graded Hilbert numerators.</span>
  - <span style="text-decoration:line-through;">Implement or link to a simplex solver, compare with DDM.</span>
  - <span style="text-decoration:line-through;">Determine what's wrong with the \f$4\times4\f$ system. (Turns out nothing was wrong: the system is simply not amenable to polyhedra.)</span>
  - <span style="text-decoration:line-through;">Implement a global analysis at the beginning of the algorithm.</span>
  - <span style="text-decoration:line-through;">Implement Hilbert polynomials using multiple-precision arithmetic. (Denominators get too large for `long long`!!!)</span>

  \subsection mdpri Medium priority

  \todo These items would be nice, but aren&rsquo;t a big deal for me at present.
  - Improve rings and fields:
    - Create a general \c Ring class.
      - Build polynomial rings off rings, not off fields.
        This could be difficult, since we typically
        want polynomials to have invertible coefficients.
        It doesn&rsquo;t seem strictly necessary, though:
        S-polynomials and top-reductions, for instance, can be computed
        by multiplying by the leading coefficient of the other polynomial,
        rather than by dividing by one&rsquo;s own coefficient.
      - Implement Dense_Univariate_Integer_Polynomial as a proper \c Polynomial representation.
      - Create a general \c Euclidean_Ring class. Add to it the
        divide_by_common_term() function.
    - Create a general \c Field class.
      - Subclass Prime_Field to be \c Field.
      - Implement Dense_Univariate_Rational_Polynomial as a proper \c Polynomial representation.
      - Implement \f$\mathbb Q\f$ as a field.
  - Reimplement Double_Buffered_Polynomial so its arrays contain
    pointers to Monomial, rather than an expanded Monomial.
    See if that changes things.
  - Re-examine what&rsquo;s going on with masks, since the plus to efficiency
    doesn&rsquo;t seem worth the effort.
  - Implement marked polynomials with a dynamic algorithm that works
    practically in the grevlex order, with the marked term being the true leading
    monomial. This may be very inefficient to reduce.
  - Implement a \c Dictionary_Linked_Polynomial class, where any term points
    to one unique instance of a monomial, rather than having many copies of
    monomials in different polynomials. Upside is that equality test during
    canonicalization is instantaneous (compare pointers).
    Downsides may include finding/sorting the monomials, indirection.
  - Detach monomial ordering from monomials,
    since caching ordering data doesn&rsquo;t seem to help much?
  - <span style="text-decoration:line-through;">Implement a \c Polynomial_Builder
    class to help build polynomials more easily
    by reading from an input file. That way we don&rsquo;t have to write a fresh
    control program for each example system.</span> (see user_interface())
  - <span style="text-decoration:line-through;">Implement an Indeterminate class
    and a \c Polynomial_Builder class to help build polynomials more easily.</span>

  \subsection lopri Lower priority

  \todo I&rsquo;m not sure these are worth doing.
  - Most skeleton code seems to have little overhead, so most
    &ldquo;improvements&rdquo; related to that falls here:
    - Implement DDM with the Fukuda-Prodon criterion, compare to Zolotykh&rsquo;s.
    - Implement Roune&rsquo;s algorithms for Hilbert functions.
    - Compare each potential PP with all other potential PP&rsquo;s,
      reducing the number of false positives. [This does not seem to be necessary
      at the moment, as the overhead is quite small, but it is still a thought.]
    - Add a hash mechanism to the \c constraint class to help avoid redundnacy.
  - Create a Matrix_Ordering_Data class as a subset of Monomial_Order_Data.
  - Add an \c insert() function to Monomial_Node to insert another
    Polynomial_Linked_List, subsequently to be destroyed.
  - Think about computing all inverses of a small prime field immediately at startup.
  - Test matrix orderings more thoroughly.

  \section Apologia Apologia pro labora sua

  This probably has bugs I haven&rsquo;t yet worked out,
  but I have done a lot of bug-fixing, including the use of `valgrind`
  to identify and repair memory leaks.

  I originally implemented this in Sage.
  That was purely a proof-of-concept product;
  it was very slow, and I discovered later that it had a number of bugs.

  \li It is not easy to use the guts of
  <span style="font-variant:small-caps;">Singular</span> from Sage.
  In particular, the geobuckets.
  But even if that were possible&hellip;
  \li In all the computer algebra systems I&rsquo;ve looked at,
  a monomial ordering is part of the ring structure.
  At least in <span style="font-variant:small-caps;">Singular</span>,
  a &ldquo;wgrevlex&rdquo; ordering received a different structure than a
  &ldquo;grevlex&rdquo; ordering, in particular a disadvantageous structure.
  So my preliminary implementation in Singular worked, but tended to be a lot
  slower than <tt>std()</tt> <i>even though it did less work.</i>
  In addition, the implementation crashed often, for reasons I
  wasn&rsquo;t able to sort out, even with help from the developers.

  That is when I decided to develop this code.
  As it turns out, that was a good thing, because the original Sage version
  had a number of bugs that I discovered only while developing later versions.

  Although I wrote it from scratch, without a doubt it reflects what I saw
  in CoCoA and <span style="font-variant:small-caps;">Singular</span>.
  No great claim to originality or even usability is implied.
  The intent of this software is not to compete with those,
  but to provide a more robust launchpad to implement the algorithm there
  than I had before.

  \section Copyright Copyright details

  This file is part of DynGB.

  DynGB is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  Foobar is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with DynGB. If not, see <http://www.gnu.org/licenses/>.
*/

/**
  \example test_cyclic4.cpp
  This illustrates how to compute a Gr&ouml;bner basis of the Cyclic-4 system
  \f[
    x_1 + x_2 + x_3 + x_4,\\ x_1 x_2 + x_2 x_3 + x_3 x_4 + x_4 x_1,\\
    x_1 x_2 x_3 + x_2 x_3 x_4 + x_3 x_4 x_1 + x_4 x_1 x_2,\\ x_1 x_2 x_3 x_4 - 1
  \f]
  using this package.
*/

/**
  \example test_cyclicn.cpp
  This illustrates how to compute a Gr&ouml;bner basis of the Cyclic-\f$n\f$
  system \f[
    x_1 + \cdots + x_n,\\
    x_1 x_2 + x_2 x_3 + \cdots + x_n x_1,\\
    x_1 x_2 x_3 + x_2 x_3 x_4 + \cdots + x_n x_1 x_2,\\
    \vdots\\
    x_1 \cdots x_{n-1} + x_2 \cdots x_n + \cdots + x_n x_1 \cdots x_{n-2},\\
    x_1 \cdots x_n - 1
  \f]
  using this package. This version uses the @b static Buchberger algorithm.
  See test_dynamic.cpp for the dynamic version.
*/

/**
  \example test_dynamic.cpp
  This illustrates how to compute a Gr&ouml;bner basis of the Cyclic-\f$n\f$
  system \f[
    x_1 + \cdots + x_n,\\
    x_1 x_2 + x_2 x_3 + \cdots + x_n x_1,\\
    x_1 x_2 x_3 + x_2 x_3 x_4 + \cdots + x_n x_1 x_2,\\
    \vdots\\
    x_1 \cdots x_{n-1} + x_2 \cdots x_n + \cdots + x_n x_1 \cdots x_{n-2},\\
    x_1 \cdots x_n - 1
  \f]
  using this package. This version uses the @b dynamic Buchberger algorithm.
  See test_cyclicn.cpp for the dynamic version.
*/

/**
  \example test_4by4.cpp
  This illustrates how to compute a Gr&ouml;bner basis of the \f$4\times4\f$
  system described \cite YanGeobuckets. The \f$4\times4\f$ system consists
  of the entries of the equation
  \f[AB-BA,\f]
  where each entry in \f$A\f$ and \f$B\f$ is a unique variable.
  For instance, the first entry in \f$AB-BA\f$ is
  \f[-x_{12}x_{19} + x_1x_{20} + x_2x_{24} + x_{28}x_3 - x_{17}x_4 - x_{18}x_8.\f]
  This system is notable because our code outperforms
  <span style="font-variant:small-caps;">Singular</span> on this system,
  and because the Hilbert heuristic will not work here despite the fact that
  coefficients were @c uint64_t.
*/

/**
  \example test_cab_es1.cpp
  This illustrates how to compute a Gr&ouml;bner basis of the first example
  in @cite CaboaraDynAlg, \f[
    t^4 z b  + x^3 y a ,\\
    t x^8 y z  + 32002 a b^4 c d e ,\\
    x y^2 z^2 d  + z c^2 e^2 ,\\
    t x^2 y^3 z^4  + a b^2 c^3 e^2
  \f]
*/

/**
  \example test_cab_es2.cpp
  This illustrates how to compute a Gr&ouml;bner basis of the second example
  in @cite CaboaraDynAlg, \f[
    32002 y^82 a  + x^32 z^32 ,\\
    x^45  + 32002 y^13 z^21 b ,\\
    32002 y^33 z^12  + x^41 c ,\\
    32002 y^33 z^12 d  + x^22 ,\\
    x^5 y^17 z^22 e  + 32002 ,\\
    t x y z  + 32002
  \f] This system is notable because it was the first time our code outperformed
  <span style="font-variant:small-caps;">Singular</span> on a polynomial system.
*/

/**
  \example test_cab_es4.cpp
  This illustrates how to compute a Gr&ouml;bner basis of the second example
  in @cite CaboaraDynAlg, \f[
    31838 B  + 45 P  + 35 S  + 31967 ,\\
    35 P  + 31976 S  + 25 T  + 40 Z ,\\
    31838 B^2  + 25 P S  + 31985 T  + 15 W  + 30 Z ,\\
    15 P T  + 20 S Z  + 31994 W ,\\
    31992 B^3  + P W  + 2 T Z ,\\
    3 B^2  + 31992 B S  + 99 W
  \f]
*/

/**
  \example test_cab_es5.cpp
  This illustrates how to compute a Gr&ouml;bner basis of the second example
  in @cite CaboaraDynAlg, \f[
    b x  + 32002 a y ,
    d x  + 32002 c y  + 32002 d  + y ,\\
    a^2  + b^2  + 32002 r^2 ,\\
    c^2  + d^2  + 32002 s^2  + 32001 c  + 1 ,\\
    a^2  + b^2  + 32001 a c  + c^2  + 32001 b d  + d^2  + 32002 t^2
  \f]
*/

/**
  \example test_cab_es6.cpp
  This illustrates how to compute a Gr&ouml;bner basis of the second example
  in @cite CaboaraDynAlg, \f[
    u  + v  + y  + 32002 ,\\
    t  + 2 u  + z  + -3 ,\\
    t  + 2 v  + y  + 32002 ,\\
    32002 t  + 32002 u  + 32002 v  + x  + 32002 y  + 32002 z ,\\
    t u x^2  + 26650 y z^3 ,\\
    28222 t y  + v z
  \f]
*/

/**
  \example test_cab_es9.cpp
  This illustrates how to compute a Gr&ouml;bner basis of the second example
  in @cite CaboaraDynAlg, \f[
    32002 x^2 y z^4  + t ,\\
    32002 x^5 y^7  + u z^2 ,\\
    v x^3 z  + 32002 y^2 ,\\
    w z^5  + 32002 x y^3
  \f]
*/
