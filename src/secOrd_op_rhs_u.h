/*! \file secOrd_op_rhs_u.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 04.02.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef SECORD_OP_RHS_U_H
#define SECORD_OP_RHS_U_H 

#include "esfem_fwd.h"

namespace Esfem{
  namespace SecOrd_op{
    class Rhs_u{
    public:
      explicit Rhs_u(const Io::Parameter&, const Grid::Grid_and_time&);
      ~Rhs_u();
      Rhs_u(const Rhs_u&) = delete;
      Rhs_u& operator=(const Rhs_u&) = delete;

      void assemble_and_addScaled_to(Grid::Scal_FEfun&) const;
      void assemble(Grid::Scal_FEfun&) const;
    private:
      struct Data;
      Data* d_ptr {nullptr};
    };
  }
}

#endif // SECORD_OP_RHS_U_H

/*! Log:
 */
