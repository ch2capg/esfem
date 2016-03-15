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
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlsolver.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>
#include "esfem_fwd.h"
#include "grid.h"

using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using Entity = FE_function::DiscreteFunctionSpaceType::IteratorType::Entity;
using Local_function = FE_function::LocalFunctionType;

using Linear_operator = Dune::Fem::ISTLLinearOperator<FE_function, FE_function>;
using Local_matrix = Linear_operator::LocalMatrixType;

class Brusselator_op
  : // public Dune::Fem::Operator<FE_function>,
    public Dune::Fem::DifferentiableOperator<Linear_operator>
{
public:
  explicit Brusselator_op(const Esfem::Io::Parameter&,
			  const Esfem::Grid::Grid_and_time&,
			  const Esfem::Growth,
			  const FE_function& triMassMatrix_firstArg,
			  const FE_function& triMassMatrix_secondArg);
  ~Brusselator_op();
  Brusselator_op(const Brusselator_op&) = delete;
  Brusselator_op& operator=(const Brusselator_op&) = delete;

  void operator()(const FE_function& rhs, FE_function& lhs) const override;
  void mass_matrix(const FE_function& rhs, FE_function& lhs) const;
  void massMatrix_constOne(FE_function&) const;

  void jacobian(const FE_function&, Linear_operator&) const;
private:
  struct Data;
  std::unique_ptr<Data> d_ptr;
  
  void heat_model(const Entity&, const Local_function& rhs_loc,
		  Local_function& lhs_loc) const;
  void quad_massMatrix_model(const Entity&, const Local_function& rhs_loc,
			     Local_function& lhs_loc) const;
  void jacobian_matrix_heat(const Entity&,
			    const Local_function&, Local_matrix&) const;
  void jacobian_matrix_quadMass(const Entity&,
				const Local_function&, Local_matrix&) const;
};

#endif // SECORD_OP_BRUSSELATOR_IMPL_H

/*! Log:
 */
