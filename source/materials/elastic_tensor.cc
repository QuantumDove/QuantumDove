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

#include <qdove/materials/elastic_tensor.h>

namespace gubbins
{  

  void ElasticTensor::distribute ()
  {
    // kill off old tensor representation.
    (*this).clear ();
    
    switch (this->symmetry)
      {
      case SymmetryFlag::hexagonal:
	this->distribute_to_hexagonal ();
	break;
	
      case SymmetryFlag::trigonal:
	this->distribute_to_trigonal ();
	break;
	
      case SymmetryFlag::cubic:
	this->distribute_to_cubic ();
	break;
	
      default:
	assert (false && "Not implemented yet.");
      }
  }
  
  
  void ElasticTensor::distribute_to_hexagonal ()
  {
    assert ((this->symmetry == SymmetryFlag::hexagonal) && "Internal error.");
    assert ((this->n_constants == 5) && "Internal error: The number of moduli neq 5");
    
    // Indices of moduli *must* come in the order: 11, 12, 13, 33,
    // 44. There is no way to check this!
    
    // 11
    (*this)[0][0][0][0] = this->constants[0];
    
    // 22
    (*this)[1][1][1][1] = this->constants[0];
    
    // 12
    (*this)[0][0][1][1] = (*this)[1][1][0][0] = this->constants[1];
    
    // 13
    (*this)[0][0][2][2] = (*this)[2][2][0][0] = this->constants[2];
    
    // 23
    (*this)[1][1][2][2] = (*this)[2][2][1][1] = this->constants[2];
    
    // 33
    (*this)[2][2][2][2] = this->constants[3];
    
    // 44
    (*this)[1][2][1][2] = (*this)[1][2][2][1] = (*this)[2][1][1][2] = (*this)[2][1][2][1] = this->constants[4];
    
    // 55
    (*this)[0][2][0][2] = (*this)[0][2][2][0] = (*this)[2][0][0][2] = (*this)[2][0][2][0] = this->constants[4];
    
    // 66
    (*this)[0][1][0][1] = (*this)[0][1][1][0] = (*this)[1][0][0][1] = (*this)[1][0][1][0] = 0.5 * (this->constants[0]-this->constants[1]);
  }

  void ElasticTensor::distribute_to_trigonal ()
  {
    assert ((this->symmetry == SymmetryFlag::trigonal) && "Internal error.");
    assert ((this->n_constants == 6) && "Internal error: The number of moduli neq 6");
    
    // Indices of moduli *must* come in the order: 11, 12, 13, 14,
    // 33, 44. There is no way to check this!
    
    // trigonal has five subgroups
    // case group 32 = D_3
    
    // 11
    (*this)[0][0][0][0] = this->constants[0];
    
    // 12
    (*this)[0][0][1][1] = (*this)[1][1][0][0] = this->constants[1];
    
    // 13
    (*this)[0][0][2][2] = (*this)[2][2][0][0] = this->constants[2];
    
    // 14
    (*this)[0][0][1][2] = (*this)[0][0][2][1] = (*this)[1][2][0][0] = (*this)[2][1][0][0] = this->constants[3];
    
    // 22 = 11
    (*this)[1][1][1][1] = this->constants[0];
    
    // 23 = 13
    (*this)[1][1][2][2] = (*this)[2][2][1][1] = this->constants[2];
    
    // 24 = -14
    (*this)[1][1][1][2] = (*this)[1][1][2][1] = (*this)[1][2][1][1] = (*this)[2][1][1][1] = -1.0 * this->constants[3];
    
    // 33
    (*this)[2][2][2][2] = this->constants[4];
    
    // 44
    (*this)[1][2][1][2] = (*this)[1][2][2][1] = (*this)[2][1][1][2] = (*this)[2][1][2][1] = this->constants[5];
    
    // 55
    (*this)[0][2][0][2] = (*this)[0][2][2][0] = (*this)[2][0][0][2] = (*this)[2][0][2][0] = this->constants[5];
    
    // 56 = 14 = -24
    (*this)[0][2][0][1] = (*this)[0][2][1][0] = (*this)[2][0][0][1] = (*this)[2][0][1][0] = this->constants[3];
    (*this)[0][1][0][2] = (*this)[0][1][2][0] = (*this)[1][0][0][2] = (*this)[1][0][2][0] = this->constants[3];
    
    // 66
    (*this)[0][1][0][1] = (*this)[0][1][1][0] = (*this)[1][0][0][1] = (*this)[1][0][1][0] = 0.5 * (this->constants[0]-this->constants[1]);
  }   
  
  void ElasticTensor::distribute_to_cubic ()
  {
    assert ((this->symmetry == SymmetryFlag::cubic) && "Internal error.");
    assert ((this->n_constants == 3) && "Internal error: The number of moduli neq 3");
    
    // Indices of moduli *must* come in the order: 11][ 12][ 44. There
    // is no way to check this!
    
    // 11
    (*this)[0][0][0][0] = this->constants[0];
    
    // 22
    (*this)[1][1][1][1] = this->constants[0];
    
    // 33
    (*this)[2][2][2][2] = this->constants[0];
    
    // 12
    (*this)[0][0][1][1] = (*this)[1][1][0][0] = this->constants[0];
    
    // 13
    (*this)[0][0][2][2] = (*this)[2][2][0][0] = this->constants[1];
    
    // 23
    (*this)[1][1][2][2] = (*this)[2][2][1][1] = this->constants[1];
    
    // 44
    (*this)[1][2][1][2] = (*this)[1][2][2][1] = (*this)[2][1][1][2] = (*this)[2][1][2][1] = this->constants[2];
    
    // 55
    (*this)[0][2][0][2] = (*this)[0][2][2][0] = (*this)[2][0][0][2] = (*this)[2][0][2][0] = this->constants[2];
    
    // 66
    (*this)[0][1][0][1] = (*this)[0][1][1][0] = (*this)[1][0][0][1] = (*this)[1][0][1][0] = this->constants[2];
  }
  
} // namespace gubbins




