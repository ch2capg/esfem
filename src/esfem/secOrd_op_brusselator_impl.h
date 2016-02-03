/*! \file secOrd_op_brusselator_impl.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Februar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 02.02.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef SECORD_OP_BRUSSELATOR_IMPL_H
#define SECORD_OP_BRUSSELATOR_IMPL_H 

#include <config.h>
#include <dune/fem/operator/common/operator.hh>
#include "esfem_fwd.h"
#include "secOrd_op_brusselator.h"
#include "grid.h"

using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using Entity = FE_function::DiscreteFunctionSpaceType::IteratorType::Entity;
using Local_function = FE_function::LocalFunctionType;

class Brusselator_op : public Dune::Fem::Operator<FE_function>{
public:
  explicit Brusselator_op(const Esfem::Io::Parameter&,
			  const Esfem::Grid::Grid_and_time&,
			  const Esfem::Growth,
			  const FE_function& triMassMatrix_firstArg,
			  const FE_function& triMassMatrix_secondArg);
  Brusselator_op(const Brusselator_op&) = delete;
  Brusselator_op& operator=(const Brusselator_op&) = delete;

  void operator()(const FE_function& rhs, FE_function& lhs) const override;
  void mass_matrix(const FE_function& rhs, FE_function& lhs) const;
  void massMatrix_constOne(FE_function&) const;
private:
  void heat_model(const Entity&, const Local_function& rhs_loc,
		  Local_function& lhs_loc) const;
  void quad_massMatrix_model(const Entity&, const Local_function& rhs_loc,
			     Local_function& lhs_loc) const;
  
  const FE_function& first_arg;
  const FE_function& second_arg;
  double massMatrix_lhs {0.};
  double stiffnessMatrix_lhs {0.};
  double quadMassMatrix_lhs {0.};
  double massMatrix_rhs {0.};
};

#endif // SECORD_OP_BRUSSELATOR_IMPL_H

/*! Log:
 */
