# DynGB
Testbed for exploring dynamic algorithms to compute a Groebner basis

To build this, you will need the autotools libraries.
I build this on several different workstations using something along these lines:

    autoreconf
    automake --add-missing
    autoconf
    cd debug
    ../configure
    make -j 6 some_executable

This assumes that you have a `debug` folder; if not, you have to make it.
The `autoreconf`, `automake`, and `autoconf` commands may not be necessary,
depending on how agreeable your system finds the last one I pushed from.
(I'm afraid I haven't quite figured out autotools.)

If you are on macOS (as I sometimes am) and have installed PPL, GLPK, etc. via ports
or brew or some other such system, you will need to specify the include and link
directories using something along these lines:

   ../configure LDFLAGS='-L/opt/local/lib' CPPFLAGS='-I/opt/local/include'

If you want compiler flags such as `-Ofast` etc., go ahead and specify them with
`CPPFLAGS` while you're at it.

Depending on what I'm working on at the moment, some older stuff may not actually work.
I do try to keep things in order.

For more details, you really do have to make and read the documentation.
For that, you will need to install the doxygen system.
You can generate the documentation by running the script `make_documentation`.
In theory, automake should create this for you, but for some reason
the directions I followed don't seem to work, which probably means
I didn't follow them correctly.
Nevertheless, installing doxygen and running the script work for me.
