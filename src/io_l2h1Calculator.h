/*! \file io_l2h1Calculator.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 27.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef IO_L2H1CALCULATOR_H
#define IO_L2H1CALCULATOR_H 

#include <memory>
#include "esfem_fwd.h"

namespace Esfem{
  namespace Io{
    class L2H1_calculator{
    public:
      explicit L2H1_calculator(const Grid::Grid_and_time&,
			       const Grid::Scal_FEfun& exact_solution,
			       const Grid::Scal_FEfun& numerical_solution);
      ~L2H1_calculator();
      L2H1_calculator(const L2H1_calculator&) = delete;
      L2H1_calculator& operator=(const L2H1_calculator&) = delete;
      
      double l2_err() const;
      double h1_err() const;
    private:
      struct Data;
      std::unique_ptr<Data> d_ptr;
    };
  }
}

#endif // IO_L2H1CALCULATOR_H

/*! Log:
 */
