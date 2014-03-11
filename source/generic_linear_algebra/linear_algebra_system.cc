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

#include <cassert>
#include <gubbins/generic_linear_algebra/linear_algebra_system.h>

namespace gubbins
{
  /* Constructor */
  LinearAlgebraSystem::LinearAlgebraSystem () 
    :
    // These defaults are temporary
    system_type (StandardHermitian)
  {}

  void
  LinearAlgebraSystem::set_system_type (gubbins::SystemType &linear_algebra_system_type)
  {}

  void
  LinearAlgebraSystem::reinit ()
  {
    // conditional reflecting the symmetry of the system; the defaule is always Hermitian
    bool is_symmetric = true;

    // Cases where matrices are *not* symmetric
    if (system_type==StandardNonHermitian)
      is_symmetric = false;

    // Actually initialise the matrices 
    // (generalised_system.A).reinit (dof_handler.n_dofs (),
    // 			  dof_handler.n_dofs (),
    // 			  dof_handler.max_couplings_between_dofs (),
    //                    is_symmetric);
    
    // actually initialise vectors
    // (system.b).reinit (dof_handler.n_dofs ());    
    // (system.x).reinit (dof_handler.n_dofs ());    

    (void) is_symmetric;
  }
}
