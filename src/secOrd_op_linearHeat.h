/*! \file secOrd_op_linearHeat.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 28. Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 28.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef SECORD_OP_LINEARHEAT_H
#define SECORD_OP_LINEARHEAT_H 

#include <memory>
#include "esfem_fwd.h"

namespace Esfem{
  namespace SecOrd_op{
    class Linear_heat{
    public:
      explicit Linear_heat(const Io::Parameter&, const Grid::Grid_and_time&);
      ~Linear_heat();
      Linear_heat(const Linear_heat&) = delete;
      Linear_heat& operator=(const Linear_heat&) = delete;
      
      void solve(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const;
      void apply_massMatrix_to(Grid::Scal_FEfun&) const;
      void mass_matrix(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const;
    private:
      struct Data;
      std::unique_ptr<Data> d_ptr;
    };    
  }
}

#endif // SECORD_OP_LINEARHEAT_H

/*! Log:
 */
