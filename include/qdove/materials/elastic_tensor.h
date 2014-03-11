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

#include <qdove/materials/tensor_base.h>

#ifndef __gubbins_elastic_tensor_h
#define __gubbins_elastic_tensor_h

namespace gubbins
{  
  /**
     The elastic tensor of moduli.  
  */
  class ElasticTensor
    :  
    public TensorBase<4> 
  {
  public:
    
    /**
       Constructor.  
    */
    ElasticTensor (const std::initializer_list<double> &moduli) 
      :  
    gubbins::TensorBase<4> (moduli) 
      {}
    
    /**
       Virtual destructor.  
    */
    virtual ~ElasticTensor ();
    
    /**
       Distribute moduli onto this tensor field according to rules
       governing the crystal group symmetry.
    */
    void distribute (); 
    
    /**
       Access operator to the underlying tensor field base.
    */
    inline  
      gubbins::TensorBase<4> operator* () const  
    {  
      return (*this);
    }
    
  private:
    
    /**
       Distribute elastic moduli onto this tensor according to the
       rules governed by the hexagonal symmetry group.
    */ 
    void distribute_to_hexagonal ();
    
    /**
       Distribute elastic moduli onto this tensor according to the
       rules governed by the cubic symmetry group.
    */
    void distribute_to_cubic ();
    
    /**
       Distribute elastic moduli onto this tensor according to the
       rules governed by the trigonal symmetry group. 
    */
    void distribute_to_trigonal ();
    
  }; /* ElasticTensor */
  
} /* namespace gubbins */

#endif /* __gubbins_elastic_tensor_h */
