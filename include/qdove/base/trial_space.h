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

#ifndef __gubbins_trial_space_h
#define __gubbins_trial_space_h

#include <deal.II/grid/tria.h>
#include <deal.II/lac/constraint_matrix.h>

namespace gubbins
{
  enum Constraints
  {
    DirchletZero
  }; // enum Constraints

  /**
     The trial space is defined by the grid and boundary constraints.
  */
  template<int dim>
    class TrialSpace
    {
    public:
      /**
	 Constructor 
      */
      TrialSpace ();

      /**
	 Constructor: Get a representation of the triangulation and boundary constraints 
      */
      TrialSpace (dealii::Triangulation<dim> &triangulation);

      /**
	 Destructor 
      */
      ~TrialSpace () {};

      /**
	 Return pointer to triangulation. 
      */
      dealii::Triangulation<dim>* triangulation ();

      /**
	 Make a unit grid. 
      */
      void make_unit_grid ();

    private:

      /**
	 Grid. 
      */
      dealii::Triangulation<dim> *triangulation_description;

      /**
	 Boundary constraints 
      */
      gubbins::Constraints constraints_type;

      /**
	 Matrix implementation of boundary constraints 
      */
      dealii::ConstraintMatrix constraint_matrix;
  };
}

#endif // __gubbins_trial_space_h
