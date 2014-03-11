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

#include <cassert>
#include <qdove/generic_linear_algebra/eigenspectrum_system.h>

namespace qdove
{
  /* Constructor */
  EigenspectrumSystem::EigenspectrumSystem () 
    :
    // These defaults are temporary
    system_type (GeneralisedHermitian)
  {}

  void
  EigenspectrumSystem::set_system_type (qdove::SystemType &eigenspectrum_system_type)
  {
    // Switch an enumerator to get problem type (generalised hermitian
    // or standard nonhermitian, or...)
    system_type = eigenspectrum_system_type; 
  }

  void
  EigenspectrumSystem::reinit ()
  {
    // conditional reflecting the symmetry of the system; the defaule is always Hermitian
    bool is_symmetric = true;

    // Cases where matrices are *not* symmetric
    if (system_type==GeneralisedNonHermitian ||
	system_type==StandardNonHermitian    )
      is_symmetric = false;

    // Actually initialise the matrices 
    if (is_symmetric)
      {
	// (generalised_system.A).reinit (dof_handler.n_dofs (),
	// 			  dof_handler.n_dofs (),
	// 			  dof_handler.max_couplings_between_dofs (),
	//                    true);
	
	// (generalised_system.B).reinit (dof_handler.n_dofs (),
	// 			  dof_handler.n_dofs (),
	// 			  dof_handler.max_couplings_between_dofs (),
	//                    true);
      }
    else
      {
	// (generalised_system.A).reinit (dof_handler.n_dofs (),
	// 			  dof_handler.n_dofs (),
	// 			  dof_handler.max_couplings_between_dofs (),
	//                    false);
	
	// (generalised_system.B).reinit (dof_handler.n_dofs (),
	// 			  dof_handler.n_dofs (),
	// 			  dof_handler.max_couplings_between_dofs (),
	//                    false);
      }

    // actually initialise eigen space (vectors and values).
    (system.x).resize (system.n_eigenpairs);
    // for (unsigned int i=0; i<system.n_eigenpairs; ++i)
    // (system.x)[i].reinit (dof_handler.n_dofs ());
    
    (system.lambda).resize (system.n_eigenpairs);
    for (unsigned int i=0; i<system.n_eigenpairs; ++i)
      (system.lambda)[i] = (0.);
  }
}
