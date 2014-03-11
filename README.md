# gubbins


gubbins is a free open-source quantum-mechanics simulator based on the
solution of the Schroedinger-Poisson equations.

At present, one-dimensional quantum heterostructures can be simulated
in a self-consistent as well as a non-selfconsistent manner. Future
support for two- and three-dimensional heterostructure is planned.

### What is needed to compile gubbins?

The libraries 

   * deal.II ( http://dealii.org/ )
   * PETSc ( http://www.mcs.anl.gov/petsc/ ), and 
   * SLEPc ( http://www.grycap.upv.es/slepc/ )

are required. That is, you need to successfully have installed deal.II
with PETSc and SLEPc on your machine before gubbins will configure.

### How do I compile gubbins?

The safe way to compile gubbins is "out of source". 

   * Make a build directory parallel to the gubbins root directory 
    
        $ mkdir gubbins-build;

   *  run cmake

        $ cmake ../gubbins;

   * finally build the library with

        $ make

      The shared library is built in gubbins/lib.

### How do I use gubbins?

To test gubbins, or to see what he can do, see the examples
directory. Compiling an example is the same as compiling the library
(except that we do not care about out of source builds).

Here's one way to compile step-0:

    $ cd examples/step-0
    $ cmake .
    $ make

With this the example "step-0" is built and linked against gubbins.

