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

#include <memory>
#include "esfem_fwd.h"

namespace Esfem{
  namespace SecOrd_op{
    class Brusselator{
    public:
      // This constructor creats two internal FE functions.
      explicit Brusselator(const Io::Parameter&, const Grid::Grid_and_time&,
			   const Growth);
      // This constructor keeps references of two external FE functions.
      explicit Brusselator(const Io::Parameter&, const Grid::Grid_and_time&,
			   const Growth, const Grid::Scal_FEfun&,
			   const Grid::Scal_FEfun&);

      ~Brusselator();
      Brusselator(const Brusselator&) = delete;
      Brusselator& operator=(const Brusselator&) = delete;

      void solve(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const;
      void mass_matrix(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const;
      void massMatrix_constOne(Grid::Scal_FEfun&) const;
      void add_massMatrixConstOne_to(Grid::Scal_FEfun&) const;
      
      // The following tow methods throw exceptions if you have constructed
      // this object with the second constructor.  
      void assign_firstArg_quadMassMatrix(const Grid::Scal_FEfun&);
      void assign_secondArg_quadMassMatrix(const Grid::Scal_FEfun&);

      // testing
      void operator()(const Grid::Scal_FEfun&, Grid::Scal_FEfun&) const;
    private:
      struct Data;
      std::unique_ptr<Data> d_ptr;
    };
  }
}


#endif // SECORD_OP_BRUSSELATOR_H

/*! Log:
 */
