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

#ifndef __gubbins_constants_h
#define __gubbins_constants_h


/** 
    A bunch of constants taken from Appendix A of Paul Harrisson's
    book, "Quantum Wells, Wires and Dots".
*/

namespace gubbins
{
  // Universal

  /**
     circumference over diameter 
  */
  const double PI      = 3.1415926535897932384;  

  /**
     Planck's constant over 2*pi 
  */
  const double HBAR    = 1.054588757e-34;        

  /**
     electronic charge           
  */
  const double E0      = 1.602189246e-19;        

  /**
     electron rest mass          
  */
  const double M0      = 9.109534e-31;           

  /**
     permittivity of free space  
  */
  const double EPSILON = 8.854187827e-12;        

  /**
     Boltzmann constant 
  */
  const double KB      = 1.3806488e-23;          

  // CdTe
  const double mstar_p_CdTe =  0.11;
  const double mstar_h_CdTe =  0.60;

  /**
     według nextnanoa' 
  */
  const double epsilon_CdTe = 10.6;              
  const double bandgap_CdTe =  1.606;

  // GaAs
  const double mstar_GaAs        =  0.067;
  const double permittivity_GaAs = 12.9 * EPSILON;

} // namespace gubbins

#endif // __gubbins_constants_h
