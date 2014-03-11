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

#ifndef __gubbins_test_space_h
#define __gubbins_test_space_h

#include <qdove/base/trial_space.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

namespace gubbins
{
  /**
     Implementation of a finite element test space.
  */
  template<int dim>
    class TestSpace
    {
    public:
      /**
	 Constructor 
      */
      TestSpace ();

      /**
	 Constructor 
      */
      TestSpace (gubbins::TrialSpace<dim> &trial_space);

      /**
	 Destructor 
      */
      ~TestSpace () 
	{
	  dof_handler.clear ();
	}

      /**
	 Return the number of degrees of freedom in this space 
      */
      unsigned int n_dofs ();

      /**
	 Return the number of degrees of freedom per cell 
      */
      unsigned int n_dofs_per_cell ();

      /**
	 Return the maximum number of couplings in this space 
      */
      unsigned int max_couplings_between_dofs ();

      /**
	 Return a reference to the dof handler used by this space 
      */
      dealii::DoFHandler<dim, dim> &dofs ();

      /**
	 Return a reference to the finite element used by this space
      */
      dealii::FESystem<dim, dim> &fe ();
      
    private:
      
      /**
	 The trial space associated with this test space.
      */
      gubbins::TrialSpace<dim> trial_space; 
      
      /**
	 Handler of degrees of freedom. 
      */
      dealii::DoFHandler<dim, dim> dof_handler;
      
      /**
	 Finite element. 
      */
      dealii::FESystem<dim, dim> finite_element;  
    };
}

#endif // __gubbins_test_space_h
