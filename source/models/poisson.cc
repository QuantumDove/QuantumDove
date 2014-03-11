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

#include <qdove/models/poisson.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace gubbins
{

  namespace Poisson
  {
    
    template <int dim>
    Problem<dim>::Problem ()
      :
      init (false)
    {}
    
    template <int dim>
    Problem<dim>::Problem (gubbins::TrialSpace<dim> &trial,
			   gubbins::TestSpace<dim>  &test)
      :
      trial_space (&trial),
      test_space (&test),
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
      
      system_vector.reinit (test_space->n_dofs ());
      
      solution_vector.reinit (test_space->n_dofs ());
      
      // Initialise boundary constraints
      constraints.clear ();
      dealii::DoFTools::make_zero_boundary_constraints (test_space->dofs (), constraints);
      constraints.close ();
      
      init = true;
    }
    
    // Simple CG solver
    template <int dim>
    void 
      Problem<dim>::assemble (const dealii::PETScWrappers::Vector &rhs_function)
    {
      // assert (false && "Pure virtual function called...   :-|");
      assert (init==true && "Problem has not been initialised");
      assert ((rhs_function.size ()==test_space->n_dofs ()) && "Incompatible vector sizes.");
      
      dealii::QGauss<dim> quadrature_formula (2);
      dealii::FEValues<dim> fe_values (test_space->fe (), quadrature_formula,
				       dealii::update_values            |
				       dealii::update_gradients         |
				       dealii::update_quadrature_points |
				       dealii::update_JxW_values);
      
      const unsigned int dofs_per_cell = test_space->n_dofs_per_cell ();
      const unsigned int n_q_points    = quadrature_formula.size ();

      std::vector<unsigned int> local_dof_indices (dofs_per_cell);

      // cell-wise representation of a system of equations
      dealii::FullMatrix<double> cell_system (dofs_per_cell, dofs_per_cell);
      dealii::Vector<double>     cell_rhs (dofs_per_cell);

      // cell-wise representation of the function
      std::vector<double> cell_rhs_function (n_q_points);

      typename dealii::DoFHandler<dim>::active_cell_iterator
       	cell = test_space->dofs ().begin_active (),
       	endc = test_space->dofs ().end();

      // assemble matrices cell-wise
      for (; cell!=endc; ++cell)
	{
	  cell_system = 0;
	  cell_rhs    = 0;
	  fe_values.reinit (cell);
	  
	  // get the representation of the function on this cell
	  fe_values.get_function_values (rhs_function, cell_rhs_function);

	  for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      {
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  {
		    
		    cell_system(i,j)
		      +=
		      fe_values.shape_grad (i,q_point) *
		      fe_values.shape_grad (j,q_point) *
		      fe_values.JxW(q_point);
		  } // i

		cell_rhs(j)
		  +=
		  cell_rhs_function[q_point]        *
		  fe_values.shape_value (j,q_point) *
		  fe_values.JxW (q_point);

	      } // j

	  // Apply constraints and distribute local objects to global
	  // objects.
	  cell->get_dof_indices (local_dof_indices);
	  
	  constraints.
	    distribute_local_to_global (cell_system,
					local_dof_indices,
					system_matrix);

	  constraints.
	    distribute_local_to_global (cell_rhs,
					local_dof_indices,
					system_vector);
	} // cell

      system_matrix.compress (dealii::VectorOperation::add);
      system_vector.compress (dealii::VectorOperation::add);
    }
    
    // Simple CG solver
    template <int dim>
    unsigned int 
      Problem<dim>::solve ()
    {
      dealii::SolverControl solver_control (solution_vector.size (),
					    1e-8*system_vector.l2_norm ());
      
      dealii::PETScWrappers::SolverCG cg (solver_control);
      dealii::PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);
      cg.solve (system_matrix, solution_vector, system_vector, preconditioner);
      
      return solver_control.last_step();
    }



    // Return the eigenpairs
    template <int dim>
    void 
    Problem<dim>::get_solution_vector (dealii::PETScWrappers::Vector &vector)
    {
       vector.reinit (test_space->n_dofs ());
       vector = solution_vector;
    }


    
  } // namespace Poisson

} // namespace gubbins

#include "poisson.inst"
