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

#include <qdove/generic_linear_algebra/generic_eigenspectrum_solver.h>

namespace qdove
{
  void solve (EigenspectrumSystem::System &eigenspectrum_system)
  {
    // Sub-standard setup for solving linear algebra system
    dealii::SolverControl solver_control (1000, 1e-03);
    dealii::SLEPcWrappers::SolverKrylovSchur solver_ks (solver_control);

    // Instruction to solve
    solver_ks.solve (eigenspectrum_system.A,
		     eigenspectrum_system.B,
		     eigenspectrum_system.lambda,
		     eigenspectrum_system.x,
		     eigenspectrum_system.n_eigenpairs);
  }
}
