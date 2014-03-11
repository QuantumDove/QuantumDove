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

#ifndef __gubbins_generic_eigenspectrum_solver_h
#define __gubbins_generic_eigenspectrum_solver_h

#include <qdove/generic_linear_algebra/eigenspectrum_system.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/slepc_spectral_transformation.h>

namespace gubbins
{
  class EigenspectrumSolver
  {
  public:

    /**
       Construcor 
    */
    EigenspectrumSolver () {};

    /**
       Destructor 
    */
    ~EigenspectrumSolver () {};

    /**
       Solve an eigenspectrum system 
    */
    void solve (EigenspectrumSystem::System &eigenspectrum_system);

  private:
  };
}

#endif // __gubbins_generic_eigenspectrum_solver_h
