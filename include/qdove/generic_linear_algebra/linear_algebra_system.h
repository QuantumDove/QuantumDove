/* 
   Copyright (C) 2013 the QuantumDove authors
   
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

#ifndef __gubbins_linear_algebra_system_h
#define __gubbins_linear_algebra_system_h

/* Types of linear algebra system */
#include <qdove/generic_linear_algebra/system_type.h>

/* Types of test space */
#include <qdove/base/test_space.h>

/* Objects needed for a sparse-matrix eigenspectrum system */
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

namespace gubbins
{
  /**
   * A class denoting the objects that make up a 
   * linear algebra system \[Ax=b\]
   */
  class LinearAlgebraSystem
  {
  public:

    /**
       Constructor 
    */
    LinearAlgebraSystem ();

    /**
       Destructor 
    */
    ~LinearAlgebraSystem () {};
    
    /**
       Objects that make up a generalised eigenspectrum problem. 
    */
    struct System
    {
      /**
	 System matrix  
      */
      dealii::PETScWrappers::SparseMatrix A;

      /**
	 System (right-hand-side) vector 
      */
      dealii::PETScWrappers::Vector       b;

      /**
	 Solution vector 
      */
      dealii::PETScWrappers::Vector       x;
    };

    /**
       Set the system type 
    */
    void set_system_type (gubbins::SystemType &linear_algebra_system_type);

    /**
       Reinitialise matrices and vectors of the linear algebra
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

#endif // __gubbins_linear_algebra_system_h

