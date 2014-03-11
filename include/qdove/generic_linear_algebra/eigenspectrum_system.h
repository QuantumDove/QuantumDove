/* 
   Copyright (C) 2013 the gubbins authors
   
   Permission is hereby granted, free of charge, to any person
   obtaining a copy of this software and associated documentation
   files (the "Software"), to deal in the Software without
   restriction, including without limitation the rights to use, copy,
   modify, merge, publish, distribute, sublicense, and/or sell copies
   of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:
   
   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.
   
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/

#ifndef __gubbins_eigenspectrum_system_h
#define __gubbins_eigenspectrum_system_h

/* Types of eigenspectrum system */
#include <qdove/generic_linear_algebra/system_type.h>

/* Objects needed for a sparse-matrix eigenspectrum system */
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

namespace gubbins
{
  /**
   * A class denoting the objects that make up a generalised
   * eigenspectrum system \[(A-B\lambda)x=0\]
   */
  class EigenspectrumSystem
  {
  public:

    /**
       Constructor 
    */
    EigenspectrumSystem ();

    /**
       Destructor 
    */
    ~EigenspectrumSystem () {};
    
    /**
       Objects that make up a generalised eigenspectrum problem.
    */
    struct System
    {
      System() : n_eigenpairs(1) {}

      /**
	 System matrix  
       */
      dealii::PETScWrappers::SparseMatrix        A;

      /**
	 Mass matrix (also called the overlap matrix in some of the
	 literature).
      */
      dealii::PETScWrappers::SparseMatrix        B;

      /**
	 Eigenfunctions 
      */
      std::vector<dealii::PETScWrappers::Vector> x;

      /**
	 Eigenvalues 
      */
      std::vector<double>                        lambda;

      /**
	 Number of eigenpairs - deafults to 1 
      */
      unsigned int                               n_eigenpairs;
    };

    /**
	Set the system type 
    */
    void set_system_type (gubbins::SystemType &eigenspectrum_system_type);

    /**
       Reinitialise matrices, vectors, and values of the eigenspectrum
       system
    */
    void reinit ();
    
  private:

    /**
       Internal reference to a generalised eigenspectrum system
    */
    System system;
    
    /**
       Internal reference to the system type 
    */
    SystemType system_type;

  };
  
}

#endif // __gubbins_eigenspectrum_system_h

