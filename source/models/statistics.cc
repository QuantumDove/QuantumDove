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

#include <gubbins/models/statistics.h>
#include <gubbins/materials/constants.h>


namespace gubbins
{

  namespace FermiDirac
  {

    void compute_density (const dealii::PETScWrappers::Vector &wavefunction,
			  dealii::PETScWrappers::Vector       &density)
    {
      // Assume that the input wavefunction is correct
      density.reinit (wavefunction.size ());

      // Do a straight point-wise multiplication
      for (unsigned int i=0; i<wavefunction.size (); ++i)
	density[i] = wavefunction[i] * wavefunction[i];
    }

    void compute_number_density (const dealii::PETScWrappers::Vector &wavefunction,
				 const double                        &energy_value,
				 const double                        &fermi_energy_value,
				 const dealii::PETScWrappers::Vector &density_of_states,
				 dealii::PETScWrappers::Vector       &density)
    {
      // Assume that the input wavefunction is correct
      assert (wavefunction.size ()==density_of_states.size () && "Incompatible vector sizes.");
      density.reinit (wavefunction.size ());

      // Compute the occupancy level from integrated Fermi-Dirac
      // statistics. Make use of fixed temperature 300K for now.
      const double kbt = gubbins::KB * 300;
      const double occupancy = kbt * std::log (1.+std::exp ( (fermi_energy_value-energy_value) / kbt) );

      // Do a straight point-wise multiplication
      for (unsigned int i=0; i<wavefunction.size (); ++i)
	density[i] = density_of_states[i] * wavefunction[i] * wavefunction[i] * occupancy;
    }

    void compute_number_density (const std::vector<dealii::PETScWrappers::Vector> &wavefunction,
				 const std::vector<double>                        &energy_values,
				 const double                                     &fermi_energy_value,
				 const dealii::PETScWrappers::Vector              &density_of_states,
				 dealii::PETScWrappers::Vector                    &density)
    {
      // Assume that the input wavefunction is correct
      assert (wavefunction.size ()==energy_values.size () && "Incompatible vector sizes.");
      assert (wavefunction[0].size ()==density_of_states.size () && "Incompatible vector sizes.");
      density.reinit (wavefunction[0].size ());

      // Make use of a temporary vector
      dealii::PETScWrappers::Vector tmp (density.size ());

      // Do a straight point-wise multiplication
      for (unsigned int i=0; i<wavefunction.size (); ++i)
         {
	   FermiDirac::compute_number_density (wavefunction[i], 
					       energy_values[i], 
					       fermi_energy_value, 
					       density_of_states, 
					       tmp);
	   density.add (tmp);
	 }
    }

    void compute_density_of_states (const dealii::PETScWrappers::Vector &effective_mass_function,
				    dealii::PETScWrappers::Vector       &density_of_states)
    {
      // Assume that the input wavefunction is correct
      density_of_states.reinit (effective_mass_function.size ());

      // Do a straight point-wise multiplication
      for (unsigned int i=0; i<density_of_states.size (); ++i)
	density_of_states[i] = (effective_mass_function[i]*gubbins::M0) / (gubbins::HBAR*gubbins::HBAR*gubbins::PI);
    }

    
  } // namespace FermiDirac

  namespace BoseEinstein
  {
    
  } // namespace FermiDirac

} // namespace gubbins

// #include "fermi_dirac.inst"
