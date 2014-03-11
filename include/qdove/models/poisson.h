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

#ifndef __gubbins_poisson_h
#define __gubbins_poisson_h

#include <qdove/base/trial_space.h>
#include <qdove/base/test_space.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

namespace gubbins
{
  
  namespace Poisson
  {

    /**
       \brief An implementation of Poisson's problem. 
       
       @author Toby D. Young 2013.
    */
    template <int dim>
      class Problem
      {
      public:
	
	/**
	   Constructor 
	*/
	Problem ();
	
	/**
	   Constructor. Generate the problem using these trial and
	   test spaces.
	*/
	Problem (gubbins::TrialSpace<dim> &trial,
		 gubbins::TestSpace<dim>  &test);
	
	/**
	   Destructor 
	*/
	~Problem ();
	
	/**
	   Reinitialise matrices and vectors. 
	*/
	void reinit ();
	
	/**
	   Assemble matrices and vectors. 
	*/
	void assemble (const dealii::PETScWrappers::Vector &rhs_function);
	
	/**
	   Solve the system. 
	*/
	unsigned int solve ();
	
	/**
	   Get the solution vector. 
	*/
	void get_solution_vector (dealii::PETScWrappers::Vector &vector);

      private:
	
	/**
	   Pointer to trial space.
	*/
	gubbins::TrialSpace<dim> *trial_space;

	/**
	   Pointer to test space. 
	*/
	gubbins::TestSpace<dim>  *test_space;
	
	/**
	   Flag indicating if the problem has been initialised. 
	*/
	bool init;
	
	/**
	   System matrix to the linear algebra equation set.
	*/
	dealii::PETScWrappers::SparseMatrix system_matrix;

	/**
	   System vector (or rhs vector) to the linear algebra equation set.
	*/
	dealii::PETScWrappers::Vector       system_vector;

	/**
	   Solution vector to the linear algebra equation set.
	*/
	dealii::PETScWrappers::Vector       solution_vector;
	
	/**
	   A matrix defining the row/column positions of constraints.
	*/
	dealii::ConstraintMatrix            constraints;

      };
    
  } // namespace Poisson
  
} // namespace gubbins

#endif // poisson_h
