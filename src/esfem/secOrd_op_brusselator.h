/*! \file secOrd_op_brusselator.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Februar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 01.02.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef SECORD_OP_BRUSSELATOR_H
#define SECORD_OP_BRUSSELATOR_H 

#include "esfem_fwd.h"

namespace Esfem{
  enum class Growth {promoting, inhibiting};
  
  namespace SecOrd_op{
    class Brusselator{
    public:
      explicit Brusselator(const Io::Parameter&, const Grid::Grid_and_time&,
			   const Growth);
      explicit Brusselator(const Io::Parameter&, const Grid::Grid_and_time&,
			   const Growth, const Grid::Scal_FEfun&,
			   const Grid::Scal_FEfun&);

      ~Brusselator();
      Brusselator(const Brusselator&) = delete;
      Brusselator& operator=(const Brusselator&) = delete;

      void solve(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const;
      void mass_matrix(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const;
      void massMatrix_constOne(Grid::Scal_FEfun&) const;
      void assign_firstArg_quadMassMatrix(const Grid::Scal_FEfun&);
      void assign_secondArg_quadMassMatrix(const Grid::Scal_FEfun&);
    private:
      struct Data;
      Data* d_ptr {nullptr};
    };
  }
}


#endif // SECORD_OP_BRUSSELATOR_H

/*! Log:
 */
