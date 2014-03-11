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

#ifndef __qdove_function_library_h
#define __qdove_function_library_h

#include <deal.II/base/function.h>

namespace qdove
{
  /**
    Anisotropic harmonic oscillator.  \[f(x):=k^2x^2\]
   */
  template <int dim>
    class IsotropicOscillator
    :
    public dealii::Function<dim>
  {
  public:
    virtual double value (const dealii::Point<dim> &p,
                          const unsigned int        component = 0) const;

    virtual void value_list (const std::vector<dealii::Point<dim> > &points,
			     std::vector<double>                    &values,
			     const unsigned int                      component = 0) const;
  };

  /**
    Anisotropic harmonic oscillator.  \[f(x_i):=k_i^2x_i^2\]
   */
  template <int dim>
    class AnisotropicOscillator
    :
    public dealii::Function<dim>
  {
  public:
    virtual double value (const dealii::Point<dim> &p,
                          const unsigned int        component = 0) const;

    virtual void value_list (const std::vector<dealii::Point<dim> > &points,
                             std::vector<double>                    &values,
                             const unsigned int                      component = 0) const;
  };
}

#endif // __qdove_function_library_h
