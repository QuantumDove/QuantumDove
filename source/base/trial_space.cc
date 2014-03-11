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

#include <qdove/base/trial_space.h>

namespace gubbins
{
  template<int dim>
  TrialSpace<dim>::TrialSpace ()
  {}

  template<int dim>
  TrialSpace<dim>::TrialSpace (dealii::Triangulation<dim> &triangulation)
    :
    triangulation_description (&triangulation),
    constraints_type (gubbins::Constraints::DirchletZero)
  {}

  template<int dim>
  dealii::Triangulation<dim>*
  TrialSpace<dim>::triangulation ()
  {
    return (this->triangulation_description);
  }

  template<int dim>
  void
  TrialSpace<dim>::make_unit_grid ()
  {}

}

#include "trial_space.inst"
