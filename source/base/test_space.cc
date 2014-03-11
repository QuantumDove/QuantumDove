/* 
   Copyright (C) 2013 the QDove authors
   
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

#include <qdove/base/test_space.h>

namespace gubbins
{
  template<int dim>
  TestSpace<dim>::TestSpace ()
    :
    finite_element (dealii::FE_Q<dim, dim> (1), 1)
  {}

  template<int dim>
  TestSpace<dim>::TestSpace (gubbins::TrialSpace<dim> &trial_space)
    :
    dof_handler (*(trial_space.triangulation ())),
    finite_element (dealii::FE_Q<dim, dim> (1), 1)
  {
    dof_handler.distribute_dofs (finite_element);
  }

  template<int dim>
  dealii::DoFHandler<dim, dim>&
  TestSpace<dim>::dofs () 
  {
    return this->dof_handler;
  }

  template<int dim>
  unsigned int
  TestSpace<dim>::max_couplings_between_dofs ()
  {
    return this->dof_handler.max_couplings_between_dofs ();
  }

  template<int dim>
  unsigned int 
  TestSpace<dim>::n_dofs ()
  {
    return this->dof_handler.n_dofs ();
  }

  template<int dim>
  unsigned int 
  TestSpace<dim>::n_dofs_per_cell ()
  {
    return this->finite_element.dofs_per_cell;
  }

  template<int dim>
  dealii::FESystem<dim, dim>&
  TestSpace<dim>::fe ()
  {
    return this->finite_element;
  }

} // namespace TestSpace

#include "test_space.inst"

