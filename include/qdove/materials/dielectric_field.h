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

#ifndef __qdove_dielectric_tensor_h
#define __qdove_dielectric_tensor_h

#include <qdove/materials/tensor_base.h>

namespace qdove
{
  
  /**
     The dielectric tensor of moduli.  
  */
  class DielectricTensor
    : 
    qdove::TensorBase<2>
  {
    
    DielectricTensor (const std::initializer_list<double> &moduli) 
      :  
    qdove::TensorBase<2> (moduli) 
      {}
    
    /**
       Virtual destructor.  
    */
    virtual ~DielectricTensor ();
    
    /**
       Distribute moduli onto this tensor field according to rules
       governing the crystal group symmetry.
    */
    void distribute (); 
    
    /**
       Access operator to the underlying tensor field base.
    */
    inline  
      qdove::TensorBase<2> operator* () const  
    {  
      return (*this);
    }
    
  private:

    /**
       Distribute dielectric moduli onto this tensor according to the
       rules governed by the hexagonal symmetry group. 
    */ 
    void distribute_to_hexagonal ();

    /**
       Distribute dielectric moduli onto this tensor according to the
       rules governed by the trigonal symmetry group.
    */
    void distribute_to_trigonal ();
    
    /**
       Distribute dielectric moduli onto this tensor according to the
       rules governed by the cubic symmetry group.
    */
    void distribute_to_cubic ();
    
  }; /* DielectricField */

} /* namespace qdove */

#endif /* __qdove_dielectric_tensor_h */


