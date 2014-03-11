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

#ifndef __qdove_material_symmetry_h
#define __qdove_material_symmetry_h

namespace qdove
{
  /* TODO: Investigate how a union can constructed to access
     combinations of flags using group number. */
  
  /**
     Table of symmetry flags - quick reference. 
  */     
  enum SymmetryFlag
  {
    
    /**
       A null symmetry used do denote a unknown, unused, or no
       existing symmetry group.
    */
    null         = 0,     
    
    /**
       Actual symmetry groups 
    */
    triclinic    = 0x0001,
    monoclinic   = 0x0002,
    orthorhombic = 0x0003, 
    tetragonal   = 0x0004, 
    trigonal     = 0x0005,     
    hexagonal    = 0x0006,     
    cubic        = 0x0007
  }; 
  
  /**
     Output stream operator that outputs the name of the symmetry
     flag.
  */
  template <class STREAM>
    inline
    STREAM &operator << (STREAM &stream, const SymmetryFlag symmetry_flag)
    {
      stream << "   SymmetryFlag::Crystal ";
      
      if (symmetry_flag &SymmetryFlag::null)         stream << "null"; 
      if (symmetry_flag &SymmetryFlag::triclinic)    stream << "triclinic";
      if (symmetry_flag &SymmetryFlag::monoclinic)   stream << "monoclinic";
      if (symmetry_flag &SymmetryFlag::orthorhombic) stream << "orthorhombic";
      if (symmetry_flag &SymmetryFlag::tetragonal)   stream << "tetragonal";
      if (symmetry_flag &SymmetryFlag::trigonal)     stream << "trigonal";
      if (symmetry_flag &SymmetryFlag::hexagonal)    stream << "hexagonal";
      if (symmetry_flag &SymmetryFlag::cubic)        stream << "cubic";
      
      return symmetry_flag;
    }
} // qdove

#endif // __qdove_material_symmetry_h
