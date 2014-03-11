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

#include <gubbins/materials/dielectric_field.h>

namespace gubbins
{  

  void DielectricTensor::distribute ()
  {
    // kill off old tensor representation.
    (*this).clear ();
    
    switch (this->symmetry)
      {
      case SymmetryFlag::trigonal:
	this->distribute_to_trigonal ();
	break;
	
      case SymmetryFlag::hexagonal:
	this->distribute_to_hexagonal ();
	break;
	
      case SymmetryFlag::cubic:
	this->distribute_to_cubic ();
	break;
	
      default:
	assert (false && "Not implemented yet.");
      }
  }
  
  void DielectricTensor::distribute_to_hexagonal ()
  {
    assert ((this->symmetry == SymmetryFlag::hexagonal) && "Internal error.");
    assert ((this->n_constants == 1) && "Internal error: The number of moduli neq 2");
    
    // Indices of moduli *must* come in the order: 11, 33. There is
    // no way to check this!
    
    // 11
    (*this)[0][0] = this->constants[0];
    
    // 22
    (*this)[1][1] = this->constants[0];
    
    // 33
    (*this)[2][2] = this->constants[1];
  }
  
  void DielectricTensor::distribute_to_trigonal ()
  {
    assert ((this->symmetry == SymmetryFlag::trigonal) && "Internal error.");
    assert ((this->n_constants == 2) && "Internal error: The number of moduli neq 2");
    
    // Indices of moduli *must* come in
    // the order: 11, 33. There is no way
    // to check this!
    
    // 11
    (*this)[0][0] = this->constants[0];
    
    // 22
    (*this)[1][1] = this->constants[0];
    
    // 33
    (*this)[2][2] = this->constants[1];
  }
    
  void DielectricTensor::distribute_to_cubic ()
  {
    assert ((this->symmetry == SymmetryFlag::cubic) && "Internal error.");
    assert ((this->n_constants == 1) && "Internal error: The number of moduli neq 1");
    
    // Indices of moduli *must* come in the order: 11. There is no
    // way to check this!
    
    // 11
    (*this)[0][0] = this->constants[0];
    
    // 22
    (*this)[1][1] = this->constants[0];
    
    // 33
    (*this)[2][2] = this->constants[0];
  }    
  
} // namespace gubbins



