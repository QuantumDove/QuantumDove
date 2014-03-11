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

#ifndef __gubbins_schroedinger_h
#define __gubbins_schroedinger_h

#include <gubbins/base/trial_space.h>
#include <gubbins/base/test_space.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/slepc_solver.h>

namespace gubbins
{

  namespace Schroedinger
  {

    /**
       \brief An implementation of Schroedinger's problem.
       
       @author Toby D. Young 2013.
    */
    template <int dim>
    class Problem
    {
    public:

      /**
          Constructor
      */
      Problem (unsigned int eigenpairs = 1);

      /**
	 Constructor. Generate the problem using these trial and
	 test spaces.
      */
      Problem (gubbins::TrialSpace<dim> &trial,
               gubbins::TestSpace<dim>  &test,
               unsigned int eigenpairs = 1);

      /**
         Destructor
      */
      ~Problem ();

      /**
	 Reinitialise matrices and vectors.
      */
      void reinit  ();

      /**
         Assemble the Schroedinger problem.
      */
      void assemble (const dealii::PETScWrappers::Vector &ke_function,
                     const dealii::PETScWrappers::Vector &pe_function);

      /**
          Solve the system.
      */
      unsigned int solve ();

      /**
         Get the solution eigenpairs.
      */
      void get_solution_eigenpairs (std::vector<double>                        &values,
                                    std::vector<dealii::PETScWrappers::Vector> &vectors);

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
	 System matrix (or hamiltonian matrix) to the eigenspectrum
	 problem.
      */
      dealii::PETScWrappers::SparseMatrix        system_matrix;

      /**
	 Overlap matrix (or mass matrix) to the eigenspectrum problem.
      */
      dealii::PETScWrappers::SparseMatrix        overlap_matrix;

      /**
	 Solution eigenvalues to the eigenspectrum problem.
      */
      std::vector<double>                        solution_values;

      /**
	 Solution eigenvectors to the eigenspectrum problem.
      */
      std::vector<dealii::PETScWrappers::Vector> solution_vectors;

      /**
	 A matrix defining the row/column positions of constraints.
      */
      dealii::ConstraintMatrix                   constraints;

      /**
         number of eigenpairs to solve for
      */
      const unsigned int n_eigenpairs;

      /**
         Flag indicating if the problem has been initialised.
      */
      bool init;
    };

  } // namespace Schroedinger

} // namespace gubbins

#endif // __gubbins_schroedinger_h
