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

#ifndef __qdove_fick_h
#define __qdove_fick_h

#include <qdove/base/trial_space.h>
#include <qdove/base/test_space.h>

#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/petsc_vector.h>

#include <fstream>
#include <iostream>
#include <sstream>

namespace qdove
{
  namespace Fick
  {

    /**
	\brief An implementation of Fick's problem. 

	@author Toby D. Young and Karl Rupp 2013
    */
    template <int dim>
    class Problem
    {
    public:

      /**
	 Constructor 
      */
      Problem ();

      /**
	 Constructor. Generate the problem using these trial and
	 test spaces.
      */
      Problem (qdove::TrialSpace<dim> &trial,
	       qdove::TestSpace<dim>  &test);
      
      /**
	 Destructor 
      */
      ~Problem ();

    private:
    
      /**
         Pointer to trial space.
      */  
      qdove::TrialSpace<dim> *trial_space;

      /**
         Pointer to test space.
      */
      qdove::TestSpace<dim>  *test_space;

      /**
	 Flag indicating if the problem has been initialised. 
      */
      bool init;

    }; // Problem


    /**
       \brief An implementation of Fick's problem. 

       This is an analytical solution to Fick's Laws, where the
       initial condition is either a rectangle function \f$\Omega \in\f$
       span\f$[-x:x]\f$, or a single transition at \f$x\f$.

       @author Toby D. Young and Karl Rupp 2013
    */
    template <int dim>
      class Solution
      :
      public dealii::Function<dim>
      {
      public:
	
	/**
	   Constructor.
	*/
	Solution ();
	
	/**
	   Constructor. Generate the problem using these trial and
	   test spaces.
	*/
	Solution (qdove::TrialSpace<dim> &trial,
		  qdove::TestSpace<dim>  &test);
	
	/**
	   Destructor 
	*/
	~Solution ();
	
	/**
	 Set initial lenght of diffusion. 
	*/
	void set_initial_length (const double initial_length);
	
	/**
	   Set initial value of diffusion. 
	*/
	void set_initial_height (const double initial_height);
	
	/**
	 Set "rate" of diffusion. 
	*/
	void set_initial_rate (const double initial_rate);
	
	/**
	   Set symmetry flag. If true, the returned solution is
	   symmetric around zero.  Otherwise, a transition from 0 to 1
	   at the initial length is returned.
	*/
	void set_symmetric_profile (const bool b);
	
	/**
	   Return the analytic solution to Fick's equation at this
	   point.
	*/
	virtual
	  double value (const dealii::Point<dim> &point,
			const unsigned int        component = 0) const;
	
	/**
	   Return a list of analytic solutions to Fick's equation at
	   these points. 
	*/
	virtual
	  void value_list (const std::vector<dealii::Point<dim> > &points,
			   std::vector<double>                    &value_list) const;
	
	/**
	   Interpolate the analytic solution to Fick's equation onto
	   this vector with this dof_handler.
	*/
	void interpolate_analytic_solution (dealii::PETScWrappers::Vector &interpolated_solution);
	
	/**
	   Reinitialise problem to default values.
	*/
	void reinit ();
	
      private:
	
	/**
	   Pointers to trial space
	*/
	qdove::TrialSpace<dim> *trial_space;
	
	/**
	   Pointers to test space
	*/
	qdove::TestSpace<dim>  *test_space;
	
	/**
	   Flag indicating whether a symmetric or non-symmetric
	   transition should be returned
	*/
	bool profile_is_symmetric;
	
	/**
	   Flag indicating if the problem has been initialised.
	*/
	bool init;

	/**
	   Half-length of the initial rectangle function (centredd on
	   zero).
	*/
	double length;
	
	/**
	   Initial value of the "concentration" (or height of the
	   rectangle function).
	*/
	double height;
	
	/**
	   Rate of "diffusion". 
	*/
	double rate;
	
      }; // Solution
    
  } // namespace Fick
  
} // namespace qdove

#endif // __qdove_fick_h
