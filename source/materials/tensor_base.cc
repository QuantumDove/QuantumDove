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

#include <cmath>
#include <cstdlib>
#include <iostream>

#include <qdove/materials/tensor_base.h>

namespace gubbins
{  
  
  template <int rank>
  TensorBase<rank>::~TensorBase ()
  {
    // wipe the underlying tensor
    (*this).clear ();
  }

  template <int rank>
  unsigned int
  TensorBase<rank>::n_components () const
  {
    // = dim^rank
    return static_cast<unsigned int> (std::pow (3, rank));
  }

  template <int rank>
  SymmetryFlag
  TensorBase<rank>::symmetry_flag () const
  {
    return this->symmetry;
  }
  
  template <int rank>
  void
  TensorBase<rank>::distribute () 
  {
    // Pure virtual function called.
    assert (false && "Virtual function called."); 
  }
  
} // namespace gubbins


#include "tensor_base.inst"
