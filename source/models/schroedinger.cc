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

#include <qdove/models/schroedinger.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace gubbins
{

  namespace Schroedinger
  {

    template <int dim>
    Problem<dim>::Problem (unsigned int eigenpairs)
      :
      n_eigenpairs (eigenpairs),
      init (false)
    {}

    template <int dim>
    Problem<dim>::Problem (gubbins::TrialSpace<dim> &trial,
                           gubbins::TestSpace<dim>  &test,
                           unsigned int eigenpairs)
      :
      trial_space (&trial),
      test_space (&test),
      n_eigenpairs (eigenpairs),
      init (false)
    {}

    template <int dim>
    Problem<dim>::~Problem ()
    {}

    // Setup the matrices and vectors
    template <int dim>
    void
    Problem<dim>::reinit ()
    {
      // distribute degrees of freedom on the finite element space
      test_space->dofs ().distribute_dofs (test_space->fe ());

      // Initialise system matrices and vectors.
      system_matrix.reinit (test_space->n_dofs (),
			    test_space->n_dofs (),
			    test_space->max_couplings_between_dofs ());
      
      overlap_matrix.reinit (test_space->n_dofs (),
			     test_space->n_dofs (),
			     test_space->max_couplings_between_dofs ());
      
      solution_vectors.resize (n_eigenpairs);
      for (unsigned int i=0; i<n_eigenpairs; ++i)
        solution_vectors[i].reinit (test_space->n_dofs ());
      
      solution_values.resize (n_eigenpairs);
      for (unsigned int i=0; i<n_eigenpairs; ++i)
        solution_values[i] = 0.;
      
      // Initialise boundary constraints
      constraints.clear ();
      dealii::DoFTools::make_zero_boundary_constraints (test_space->dofs (), constraints);
      constraints.close ();
      
      init = true;
    }

    // Assembly
    template <int dim>
    void
    Problem<dim>::assemble (const dealii::PETScWrappers::Vector &ke_function,
			    const dealii::PETScWrappers::Vector &pe_function)
    {
      // assert (false && "Pure virtual function called...");
      assert (init==true && "Problem has not been initialised");
      assert ((ke_function.size ()==test_space->n_dofs ()) && "Incompatible vector sizes.");
      assert ((pe_function.size ()==test_space->n_dofs ()) && "Incompatible vector sizes.");

      dealii::QGauss<dim> quadrature_formula (2);
      dealii::FEValues<dim> fe_values (test_space->fe (), quadrature_formula,
               dealii::update_values            |
               dealii::update_gradients         |
               dealii::update_quadrature_points |
               dealii::update_JxW_values);

      const unsigned int dofs_per_cell = test_space->n_dofs_per_cell ();
      const unsigned int n_q_points    = quadrature_formula.size ();

      std::vector<unsigned int> local_dof_indices (dofs_per_cell);

      // cell-wise representation of eigenspectrum problem
      dealii::FullMatrix<double> cell_system (dofs_per_cell, dofs_per_cell);
      dealii::FullMatrix<double> cell_overlap (dofs_per_cell, dofs_per_cell);

      // cell-wise representation of the function
      std::vector<double> cell_ke_function (n_q_points);
      std::vector<double> cell_pe_function (n_q_points);

      typename dealii::DoFHandler<dim>::active_cell_iterator
        cell = test_space->dofs ().begin_active (),
        endc = test_space->dofs ().end();

      // assemble matrices cell-wise
      for (; cell!=endc; ++cell)
      {
        cell_system  = 0;
        cell_overlap = 0;
        fe_values.reinit (cell);

        // get the representation of the function on this cell
        fe_values.get_function_values (ke_function, cell_ke_function);
        fe_values.get_function_values (pe_function, cell_pe_function);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
	      {
		
		// assemble local matrix
		cell_system(i,j)
		  +=
		  cell_ke_function[q_point]        *
		  fe_values.shape_grad (i,q_point) *
		  fe_values.shape_grad (j,q_point) *
		  fe_values.JxW(q_point)
		  +
		  cell_pe_function[q_point]         *
		  fe_values.shape_value (i,q_point) *
		  fe_values.shape_value (j,q_point) *
		  fe_values.JxW (q_point)
		  ;
		
		cell_overlap(i,j)
		  +=
		  fe_values.shape_value (i,q_point) *
		  fe_values.shape_value (j,q_point) *
		  fe_values.JxW (q_point);
	      }
	
        // Apply constraints and distribute local objects to global
        // objects.
        cell->get_dof_indices (local_dof_indices);
	
        constraints.
          distribute_local_to_global (cell_system,
				      local_dof_indices,
				      system_matrix);
	
        constraints.
          distribute_local_to_global (cell_overlap,
				      local_dof_indices,
				      overlap_matrix);
	
      } // cell
      
      system_matrix.compress (dealii::VectorOperation::add);
      overlap_matrix.compress (dealii::VectorOperation::add);
    }

    // Simple solver. \todo Figure out why scaling the matrices did
    // not help the iterative solver to converge. (Increasing the
    // convergence limit of the solver to 1e-24 did help however...).
    template <int dim>
    unsigned int
    Problem<dim>::solve ()
    {
      // const double factor = overlap_matrix (0,0);
      // assert ((factor!=0) && "A highly improbable internal error has occured in the schroedinger solver routine.");
      // system_matrix  /= factor;
      // overlap_matrix /= factor;

      dealii::SolverControl solver_control (n_eigenpairs*system_matrix.m (), 1e-24);
      dealii::SLEPcWrappers::SolverLAPACK lapack (solver_control);

      lapack.set_which_eigenpairs (EPS_SMALLEST_REAL);
      lapack.solve (system_matrix, overlap_matrix, solution_values, solution_vectors, n_eigenpairs);

      for (unsigned int i=0; i<n_eigenpairs; ++i)
      {
        constraints.distribute (solution_vectors[i]);
        const double overlap_norm_square = overlap_matrix.matrix_norm_square (solution_vectors[i]);
        solution_vectors[i] /= sqrt (overlap_norm_square);
      }

      return solver_control.last_step();
    }

    // Return the eigenpairs
    template <int dim>
    void
    Problem<dim>::get_solution_eigenpairs (std::vector<double>                        &values,
					   std::vector<dealii::PETScWrappers::Vector> &vectors)
    {
      vectors.resize (n_eigenpairs);
      for (unsigned int i=0; i<n_eigenpairs; ++i)
	{
	  vectors[i].reinit (test_space->n_dofs ());
	  vectors[i] = solution_vectors[i];
	}
      
      values.resize (n_eigenpairs);
      for (unsigned int i=0; i<n_eigenpairs; ++i)
	{
	  values[i] = solution_values[i];
	}
    }
    
    
  } // namespace Schroedinger

} // namespace gubbins

#include "schroedinger.inst"
