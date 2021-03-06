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

#ifndef __qdove_schroedinger_h
#define __qdove_schroedinger_h

#include <qdove/base/trial_space.h>
#include <qdove/base/test_space.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/slepc_solver.h>

namespace qdove
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
      Problem (qdove::TrialSpace<dim> &trial,
               qdove::TestSpace<dim>  &test,
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
      qdove::TrialSpace<dim> *trial_space;

      /**
         Pointer to test space.
      */
      qdove::TestSpace<dim>  *test_space;

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

} // namespace qdove

#endif // __qdove_schroedinger_h
