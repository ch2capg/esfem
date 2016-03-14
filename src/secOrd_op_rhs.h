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

#include "esfem_fwd.h"

namespace Esfem{
  namespace SecOrd_op{
    class Rhs{
    public:
      explicit Rhs(const Grid::Grid_and_time&);
      ~Rhs();
      Rhs(const Rhs&) = delete;
      Rhs& operator=(const Rhs&) = delete;

      void assemble_and_addScaled_to(Grid::Scal_FEfun&);
    private:
      struct Data;
      Data* d_ptr {nullptr};
    };
  }
}

#endif // SECORD_OP_RHS_H

/*! Log:
 */
