/*! \file secOrd_op_initData.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 30.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef SECORD_OP_INITDATA_H
#define SECORD_OP_INITDATA_H 

#include "esfem_fwd.h"

namespace Esfem{
  namespace SecOrd_op{
    class Init_data{
    public:
      explicit Init_data(const Grid::Grid_and_time&);
      ~Init_data();
      Init_data(const Init_data&) = delete;
      Init_data& operator=(const Init_data&) = delete;
      
      void interpolate(Grid::Scal_FEfun&);
    private:
      struct Data;
      Data* d_ptr {nullptr};
    };
  }
}

#endif // SECORD_OP_INITDATA_H

/*! Log:
 */
