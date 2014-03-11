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

#ifndef __qdove_statistics_h
#define __qdove_statistics_h

#include <deal.II/lac/petsc_vector.h>

namespace qdove
{

  /**
     \brief Functions for working with Fermi-Dirac statistics.

     @author Toby D. Young and Karl Rupp 2013
  */
  namespace FermiDirac
  {
    
    /**
       Compute the density function for a given wavefunction such
       that: \f$\rho=\|\psi\|^2\f$.
    */
    void compute_density (const dealii::PETScWrappers::Vector &wavefunction,
			  dealii::PETScWrappers::Vector       &density);
    
    /**
       Compute the number density function from a given wavefunction
       and density of states.
       
       \todo Currently this only works at default temperature
       \f$300\,\f$K.
    */
    void compute_number_density (const dealii::PETScWrappers::Vector &wavefunction,
				 const double                        &energy_value,
				 const double                        &fermi_energy_value,
				 const dealii::PETScWrappers::Vector &density_of_states,
				 dealii::PETScWrappers::Vector       &density);
    
    /**
       Compute the number density function from a set of wavefunctions
       and densities of states.
       
       \todo Currently this only works at default temperature
       \f$300\,\f$K.
    */
    void compute_number_density (const std::vector<dealii::PETScWrappers::Vector> &wavefunction,
				 const std::vector<double>                        &energy_value,
				 const double                                     &fermi_energy_value,
				 const dealii::PETScWrappers::Vector              &density_of_states,
				 dealii::PETScWrappers::Vector                    &density);
    
    /**
       Compute the density of states from a given effective mass
       function.
    */
    void compute_density_of_states (const dealii::PETScWrappers::Vector &effective_mass_function,
				    dealii::PETScWrappers::Vector       &density_of_states);
    
  } // namespace FermiDirac
  
  
  /**
     \brief Functions for working with Bose-Einstein statistics. 

     @note This is an empty namespace.
  */
  namespace BoseEinstein
  {
    
  } // namespace BoseEinstein

} // namespace qdove

#endif // __qdove_statistics_h
