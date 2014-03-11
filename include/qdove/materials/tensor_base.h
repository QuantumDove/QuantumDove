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

#ifndef __gubbins_tensor_base_h
#define __gubbins_tensor_base_h

#include <qdove/materials/material_symmetry.h>

#include <deal.II/base/tensor.h>

#include <cassert>
#include <map>
#include <array>

namespace gubbins
{
  /**
     Base class for crystal tensors of empirical moduli. 
  */
  template <int rank>
    class TensorBase
    :
    public dealii::Tensor<rank, 3, double>
  {
  protected:
    
    /**
       Original symmetry flags handed to the constructor. 
    */
    SymmetryFlag symmetry;

  public:
      
    /**
       Constructor.  
    */
    TensorBase (const std::initializer_list<double> &moduli) 
      : 
    n_constants (moduli.size ()),
      constants (new double[n_constants]),
      tensor (dealii::Tensor<rank, 3, double> ())
	{
	  // This makes no sense
	  assert ((n_constants!=0) && "In valid state: A list of moduli can not be empty.");
	  std::copy (moduli.begin(), moduli.end(), this->constants);
	}
    
    /**
       Virtual destructor.  
    */
    virtual ~TensorBase ();
    
    /**
       Return the symmetry flag associated with this tensor field.
    */
    SymmetryFlag symmetry_flag () const; 
    
    /**
       Return the number of components in this tensor field.
    */
    unsigned int n_components () const; 
    
    /**
       Write gnuplot. 
    */
    void write_gnuplot (std::ostream &out = std::cout) const;
    
    /**
       Distribute moduli onto this tensor field according to rules
       governing the crystal group symmetry.
    */
    virtual 
      void distribute ();
    
    /**
       Assignment operator. Initialise this tensor with the T.
    */
    TensorBase<rank> &operator = (const TensorBase<rank> &T);
    
    /**
       Access operator to the underlying tensor associated with this
       tensor field.
    */
    inline
      dealii::Tensor<rank, 3, double> operator* () const
    {
      return (*this); 
    }
    
  
    /**
       The number of constants in this tensor. 
    */
    const unsigned int n_constants;
    
    /**
       An array of constants (handed to the constructor). 
    */
    double *constants; 

  private:
    
    /**
       Underlying tensor. 
    */
    dealii::Tensor<rank, 3, double> tensor;
    
  }; /* TensorBase */
  
} /* namespace gubbins */

#endif /* __gubbins_tensor_base_h */

