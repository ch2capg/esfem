/*! \file esfem_operator.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 24. Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 24.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef ESFEM_OPERATOR_H
#define ESFEM_OPERATOR_H 

#include "esfem.h"

namespace Model{
  class Initial_data
    : public Dune::Fem::Function<FunctionSpace, Impl>;
  class PDErhs
    : public Dune::Fem::Function<FunctionSpace, Impl>;
  // InitialDataType initialData {timeProvider};
  // RHSFunctionType f {timeProvider};
  // class Model;	// FunctionSpace?
  // Dune::Fem::FunctionSpace<double, double, GridType::dimensionworld, 1>
  // Discrete_function::Function_space;
  // NonlinearModel model {timeProvider};
}

NonlinearOperator ellipticOp {xi, model, bdf_alpha_coeff.back()};


#endif // ESFEM_OPERATOR_H

/*! Log:
 */
