/*! \file operator.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 25. Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 25.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef OPERATOR_H
#define OPERATOR_H 

#include <config.h>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/operator.hh>
#include "discrete_function.h"
#include "parameter.h"

namespace Operator{
  // using Solver = Dune::Fem::CGInverseOperator<FE_function>;
  void interpolate_initData(const Dune::Fem::TimeProviderBase&,
			    Discrete_function::FE_function&);
  void assemble_RHS(const Dune::Fem::TimeProviderBase&,
		    Discrete_function::FE_function&);
  
  class Linear_heat :
    public Dune::Fem::Operator<Discrete_function::FE_function>
  {
  public:
    Linear_heat() = delete;
    explicit Linear_heat(const Parameter::PDE_data&,
			 const Dune::Fem::TimeProviderBase&);
    Linear_heat(const Linear_heat&) = delete;
    Linear_heat(Linear_heat&&) = delete;
    Linear_heat& operator=(const Linear_heat&) = delete;
    Linear_heat& operator=(Linear_heat&&) = delete;
    ~Linear_heat();

    void operator()(const Discrete_function::FE_function&,
		    Discrete_function::FE_function&) const override;
    void mass_matrix(const Discrete_function::FE_function&,
		     Discrete_function::FE_function&) const;
  private:
    struct Data;
    Data* d_ptr;
  };
}

#endif // OPERATOR_H

/*! Log:
 */
