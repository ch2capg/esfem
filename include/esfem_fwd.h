/*! \file esfem_fwd.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 26. Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 26.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef ESFEM_FWD_H
#define ESFEM_FWD_H 

namespace Esfem{
  namespace Io{
    class Parameter;
    class L2H1_calculator;
    class Error_stream;
    class Paraview;
    namespace Dgf{
      class Handler;
    }
  }
  namespace Grid{
    class Deformation;
    class Grid_and_time;
    class Scal_FEfun;
    class Vec_FEfun;
  }
  namespace SecOrd_op{
    class Init_data;
    class Init_data_u;
    class Init_data_w;
    class Rhs;
    class Rhs_u;
    class Rhs_w;
    class Linear_heat;
    class Brusselator;
  }
  namespace Impl{
    class Evolving_grid;    
  }
  enum class Growth {promoting, inhibiting};
}

#endif // ESFEM_FWD_H

/*! Log:
 */
