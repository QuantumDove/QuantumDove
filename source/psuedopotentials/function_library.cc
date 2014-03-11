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

#include <qdove/psuedopotentials/function_library.h>

namespace gubbins
{
  template<int dim>
  double
  IsotropicOscillator<dim>::value (const dealii::Point<dim> &p,
				   const unsigned int) const
  {
    return p.square ();
  }

  template<int dim>
  void
  IsotropicOscillator<dim>::value_list (const std::vector<dealii::Point<dim> > &points,
					std::vector<double>                    &values,
					const unsigned int) const
  {
    Assert (values.size () == points.size (),
            dealii::ExcDimensionMismatch(values.size (), points.size ()));
    
    for (unsigned int i=0; i<points.size(); ++i)
      {
        values[i] = points[i].square();
      }
  }
}

namespace gubbins
{
  template<int dim>
  double
  AnisotropicOscillator<dim>::value (const dealii::Point<dim> &p,
				     const unsigned int) const
  {
    // cdot this with spring constant first
    return p.square ();
  }

  template<int dim>
  void
  AnisotropicOscillator<dim>::value_list (const std::vector<dealii::Point<dim> > &points,
					  std::vector<double>                    &values,
					  const unsigned int                      component) const
  {
    Assert (values.size () == points.size (),
            dealii::ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0; i<points.size(); ++i)
      values[i] = gubbins::AnisotropicOscillator<dim>::value (points[i], component);
  }
}


#include "function_library.inst"
