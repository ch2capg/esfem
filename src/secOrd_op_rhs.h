/*! \file secOrd_op_rhs.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 28.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef SECORD_OP_RHS_H
#define SECORD_OP_RHS_H 

#include <memory>
#include "esfem_fwd.h"

//! The evolving surface finite element method
namespace Esfem{
  //! Parabolic and elliptic second order operators
  namespace SecOrd_op{
    //! Scalar valued right-hand side of a PDE
    class Rhs{
    public:
      //! Get space-time and finite element space
      explicit Rhs(const Grid::Grid_and_time&);
      //! Needed for the pointer to data implementation technique
      ~Rhs();

      //! Assemble the load vector and add with some scale to the input.      
      void assemble_and_addScaled_to(Grid::Scal_FEfun&);
    private:
      struct Data;
      //! Data pointer
      std::unique_ptr<Data> d_ptr;
    };

    //! Vector valued right-hand side of a PDE
    class Vec_rhs{
    public:
      //! \copydoc Rhs::Rhs()
      explicit Vec_rhs(const Grid::Grid_and_time&);
      //! \copydoc Rhs:~Rhs()
      ~Vec_rhs();

      //! \copydoc Rhs::assemble_and_addScaled_to()
      void assemble_and_addScaled_to(Grid::Vec_FEfun&);
    private:
      struct Data;
      //! Data pointer
      std::unique_ptr<Data> d_ptr;
    };
  }
}

#endif // SECORD_OP_RHS_H
