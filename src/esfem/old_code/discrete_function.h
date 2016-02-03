/*! \file discrete_function.h

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

#ifndef DISCRETE_FUNCTION_H
#define DISCRETE_FUNCTION_H 

#include <config.h>
#include <string>
#include <dune/fem/space/lagrange.hh>	// discrete function space
#include <dune/fem/function/adaptivefunction.hh>	// discrete function
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/space/common/functionspace.hh>
#include "parameter.h"
#include "fe_grid.h"

namespace Discrete_function{
  using Function_space = Dune::Fem::
    FunctionSpace<double, double, FE_grid::Grid::dimensionworld, 1>;
  using FE_space = Dune::Fem::
    LagrangeDiscreteFunctionSpace<Function_space, FE_grid::Grid_part, POLORDER>;
  using FE_function = Dune::Fem::AdaptiveDiscreteFunction<FE_space>;
  using L2_norm = Dune::Fem::L2Norm<FE_grid::Grid_part>;
  using H1_norm = Dune::Fem::H1Norm<FE_grid::Grid_part>;

  void init_fef_and_norms(const Parameter::PDE_data&);
  
  // pseudo external variables
  L2_norm& l2norm(FE_grid::Grid_part* const = nullptr);
  H1_norm& h1norm(FE_grid::Grid_part* const = nullptr);
  FE_space& fe_space(FE_grid::Grid_part* const = nullptr);
  FE_function& numerical_solution(const std::string* const = nullptr,
				  FE_space* const = nullptr);
  FE_function& exact_solution(const std::string* const = nullptr,
				  FE_space* const = nullptr);
  FE_function& rhs(const std::string* const = nullptr,
		   FE_space* const = nullptr);
  // Right hand side for the LES
  FE_function& fef_feOperator(const std::string* const = nullptr,
			      FE_space* const = nullptr);
  // For the lineary implicit BDF method
  FE_function& load_vector(const std::string* const = nullptr,
			   FE_space* const = nullptr);
}

#endif // DISCRETE_FUNCTION_H

/*! Log:
 */
