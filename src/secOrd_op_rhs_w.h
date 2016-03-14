/*! \file secOrd_op_rhs.h

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

#ifndef SECORD_OP_RHS_W_H
#define SECORD_OP_RHS_W_H 

#include "esfem_fwd.h"

namespace Esfem{
  namespace SecOrd_op{
    class Rhs_w{
    public:
      explicit Rhs_w(const Io::Parameter&,
		     const Grid::Grid_and_time&);
      ~Rhs_w();
      Rhs_w(const Rhs_w&) = delete;
      Rhs_w& operator=(const Rhs_w&) = delete;

      void assemble_and_addScaled_to(Grid::Scal_FEfun&) const;
    private:
      struct Data;
      Data* d_ptr {nullptr};
    };
  }
}

#endif // SECORD_OP_RHS_W_H

/*! Log:
 */
