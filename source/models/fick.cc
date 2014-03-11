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

#include <gubbins/models/fick.h>

namespace gubbins
{
  namespace Fick
  {

    // Solution

    template<int dim>
    Solution<dim>::Solution ()
      :
      dealii::Function<dim> (),
      profile_is_symmetric(true),
      init (false)
    {}

    template<int dim>
    Solution<dim>::Solution (gubbins::TrialSpace<dim> &trial,
           gubbins::TestSpace<dim>  &test)
      :
      dealii::Function<dim> (),
      trial_space (&trial),
      test_space (&test),
      profile_is_symmetric(true),
      init (false)
    {}

    template<int dim>
    Solution<dim>::~Solution ()
    {}

    template<int dim>
    void
    Solution<dim>::reinit ()
    {
      // set initialisation flag to true
      init = true;
    }

    template <int dim>
    void
    Solution<dim>::set_initial_length (const double initial_length)
    {
      assert ((initial_length!=0) &&
	      "Invalid number state: the initial length is equal to zero.");
      length = initial_length;
    }

    template <int dim>
    void
    Solution<dim>::set_initial_height (const double initial_height)
    {
      assert ((initial_height>0) && (initial_height<=1) &&
        "Invalid number state (out of bounds): The initial value must be in bounds 0<initial_value<1.");
      height = initial_height;
    }

    template <int dim>
    void
    Solution<dim>::set_initial_rate (const double initial_rate)
    {
      assert ((initial_rate>0) &&
        "Invalid number state: The initial rate must be greater than zero.");
      rate = initial_rate;
    }

    template <int dim>
    void
    Solution<dim>::set_symmetric_profile(const bool b) 
    { 
      profile_is_symmetric = b; 
    }

    template <int dim>
    inline
    double
    Solution<dim>::value (const dealii::Point<dim> &point,
                          const unsigned int        component) const
    {
      switch (dim)
      {
      case 1:
        // Simplified solution to Fick's laws 1d.
        if (profile_is_symmetric)
          return height * 0.5 * ( std::erf ((length-point[0]) / rate) + std::erf ((length+point[0]) / rate) );
        else
          return height * 0.5 * ( 1.0 + std::erf ((point[0]-length) / rate) );
        break;

      default:
        assert (false && "Not implemented yet");
      }
    }

    template <int dim>
    inline
    void
    Solution<dim>::value_list (const std::vector<dealii::Point<dim> > &points,
             std::vector<double>                    &value_list) const
    {
      assert ((points.size ()==value_list.size ()) && "Internal error.");

      // loop the above function over all points i
      for (unsigned int i=0; i<points.size (); ++i)
  value_list[i] = value (points[i]);
    }

    template<int dim>
    void
    Solution<dim>::interpolate_analytic_solution (dealii::PETScWrappers::Vector &interpolated_solution)
    {
      assert (init && "Problem is in a wrong state: init (must be called first)");

      // make sure the correct dof_handler is used
      interpolated_solution.reinit (test_space->n_dofs ());

      // and interpolate solution
      dealii::VectorTools::interpolate (test_space->dofs (),
                (*this),
                interpolated_solution);
    }

    // Problem

    template<int dim>
    Problem<dim>::Problem (gubbins::TrialSpace<dim> &trial,
         gubbins::TestSpace<dim>  &test)
      :
      trial_space (&trial),
      test_space (&test),
      init (false)
    {}

    template<int dim>
    Problem<dim>::~Problem ()
    {}

  } // namespace Fick

} // namespace gubbins

// Explicit instantiations
#include "fick.inst"
