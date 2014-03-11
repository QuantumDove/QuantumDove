# QuantumDove


QuantumDove is a free open-source quantum-mechanics simulator based on
the solution of the Schroedinger-Poisson equations.

At present, one-dimensional quantum heterostructures can be simulated
in a self-consistent as well as a non-selfconsistent manner. Future
support for two- and three-dimensional heterostructure is planned.

### What is needed to compile QuantumDove?

The libraries 

   * deal.II ( http://dealii.org/ )
   * PETSc ( http://www.mcs.anl.gov/petsc/ ), and 
   * SLEPc ( http://www.grycap.upv.es/slepc/ )

are required. That is, you need to successfully have installed deal.II
with PETSc and SLEPc as dependencies on your machine before
QuantumDove will configure.

### How do I compile QuantumDove?

The safe way to compile QuantumDove is "out of source". 

   * Make a build directory within the QuantumDove root directory 
    
        $ mkdir build;

   *  run cmake

        $ cmake ../;

   * finally build the library with

        $ make

The shared library is built in QuantumDove/lib as libqdove.so.a.b.c.

### How do I use QuantumDove?

To test QuantumDove, or to see what he can do, see the examples
directory. Compiling an example is the same as compiling the library
(except that we do not care about out of source builds).

Here's one way to compile step-0:

    $ cd examples/step-0
    $ cmake .
    $ make

With this the example "step-0" is built and linked against
QuantumDove.

