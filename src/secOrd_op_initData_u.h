/*! \file secOrd_op_initData_u.h

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

#ifndef SECORD_OP_INITDATA_U_H
#define SECORD_OP_INITDATA_U_H 

#include "esfem_fwd.h"

namespace Esfem{
  namespace SecOrd_op{
    class Init_data_u{
    public:
      explicit Init_data_u(const Grid::Grid_and_time&);
      ~Init_data_u();
      Init_data_u(const Init_data_u&) = delete;
      Init_data_u& operator=(const Init_data_u&) = delete;
      
      void interpolate(Grid::Scal_FEfun&) const;
    private:
      struct Data;
      Data* d_ptr {nullptr};
    };
  }
}

#endif // SECORD_OP_INITDATA_U_H

/*! Log:
 */
