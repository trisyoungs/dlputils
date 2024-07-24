DISCLAIMER
----------

These programs are provided in good faith, but are not guaranteed to work. 

Some are written for specific (molecular dynamics) systems and may not be applicable
to others.

If these are used in a research application, it is the responsibility of the user
to ensure there are no bugs and to validate that the codes provide correct, publishable
answers.

COMPILATION
-----------

dlputils has both GNU Make and CMake build files available.

** Configure/Make

Configure the build by first running:

./configure

Assuming no errors occur, all the serial utilities can be made by simply running:

./make

For specific MPI-dependent programs, these must be compiled individually. For instance:

./make sq

** CMake

An out-of-tree build is recommended. From a build directory of your choice, run:

cmake /path/to/dlputils-rNNN

..where NNN is the revision number. On Windows you will need to specify the generator
to use - for the MinGW toolchain, use:

cmake -G "MinGW Makefiles" /path/to/dlputils-rNNN

USAGE
-----

Brief descriptions of individual program usage can be found by running the code in question
with no arguments. See http://trisyoungs.github.io/dlputils for detailed instructions.
