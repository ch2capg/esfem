/*! \file secOrd_op_brusselator_impl.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Februar 2016

     Implementation details for secOrd_op_brusselator_impl.h
     Created by Christian Power on 02.02.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <config.h>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include "secOrd_op_brusselator_impl.h"
#include "io_parameter.h"
#include "grid.h"

using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using Entity = FE_function::DiscreteFunctionSpaceType::IteratorType::Entity;
using Geometry = Entity::Geometry;
using Grid_part = FE_function::GridPartType;
using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;
using Local_function = FE_function::LocalFunctionType;
using Domain = FE_function::LocalFunctionType::DomainType;
using Range = FE_function::LocalFunctionType::RangeType;
using Jacobian_range = FE_function::LocalFunctionType::JacobianRangeType;

inline Range massMatrix_loc(const std::size_t pt, const Quadrature& q,
			    const Local_function& cf){
  Range u;
  cf.evaluate(q[pt], u);
  return u;  
}
inline Jacobian_range stiffnessMatrix_loc(const std::size_t pt, const Quadrature& q,
					  const Local_function& cf){
  Jacobian_range nabla_u;
  cf.jacobian(q[pt], nabla_u);
  return nabla_u;
}
inline Range quad_massMatrix_loc(const std::size_t pt, const Quadrature& q,
				 const Local_function& rhs_loc,
				 const Local_function& arg1_loc,
				 const Local_function& arg2_loc){
  Range u;
  Range tmp;
  rhs_loc.evaluate(q[pt], u);
  arg1_loc.evaluate(q[pt], tmp);
  u *= tmp;
  arg2_loc.evaluate(q[pt], tmp);
  u *= tmp;
  return u;
}


// ----------------------------------------------------------------------
// Implementation of secOrd_op_brusselator_impl.h

Brusselator_op::Brusselator_op(const Esfem::Io::Parameter& p,
			       const Esfem::Grid::Grid_and_time& gt,
			       const Esfem::Growth type,
			       const FE_function& quadMassMatrix_firstArg,
			       const FE_function& quadMassMatrix_secondArg)
  : first_arg {quadMassMatrix_firstArg}, second_arg {quadMassMatrix_secondArg}
{
  const auto dT = gt.time_provider().deltaT();
  switch(type){
  case Esfem::Growth::promoting:
    massMatrix_lhs = 1 + dT * p.tg_gamma();
    stiffnessMatrix_lhs = dT;
    quadMassMatrix_lhs = (-1) * dT * p.tg_gamma();
    massMatrix_rhs = dT * p.tg_gamma() * p.tg_a(); 
    break;
  case Esfem::Growth::inhibiting:
    massMatrix_lhs = 1;
    stiffnessMatrix_lhs = dT * p.tg_Dc();
    quadMassMatrix_lhs = dT * p.tg_gamma();
    massMatrix_rhs = dT * p.tg_gamma() * p.tg_b(); 
    break;
  default:
    throw std::logic_error {"Error in constructor of Data.  "
	"Unkown growth type."};
    break;
  };
}
void Brusselator_op::heat_model(const Entity& e, const Local_function& rhs_loc,
				Local_function& lhs_loc) const{
  const auto& g = e.geometry();
  Quadrature q {e, lhs_loc.order() + rhs_loc.order()};

  for(std::size_t pt = 0; pt < q.nop(); ++pt){
    // Lu = (M + tau A)u
    const auto& x = q.point(pt);
    const auto weight = q.weight(pt) * g.integrationElement(x);

    auto Mu = massMatrix_loc(pt, q, rhs_loc);
    auto Au = stiffnessMatrix_loc(pt, q, rhs_loc);

    Mu *= massMatrix_lhs * weight;
    Au *= stiffnessMatrix_lhs * weight;
    
    lhs_loc.axpy(q[pt], Mu, Au);
  }
}
void Brusselator_op::
quad_massMatrix_model(const Entity& e, const Local_function& rhs_loc,
		      Local_function& lhs_loc) const{
  const auto& g = e.geometry();
  const auto& arg1_loc = first_arg.localFunction(e);
  const auto& arg2_loc = second_arg.localFunction(e);
  Quadrature q
  {e, lhs_loc.order() + arg1_loc.order() + arg2_loc.order() + rhs_loc.order()};
  
  for(std::size_t pt = 0; pt < q.nop(); ++pt){
    const auto& x = q.point(pt);
    const auto weight = q.weight(pt) * g.integrationElement(x);

    auto Mu = quad_massMatrix_loc(pt, q, rhs_loc, arg1_loc, arg2_loc);

    Mu *= quadMassMatrix_lhs * weight;
    
    lhs_loc.axpy(q[pt], Mu);
  }
}
void Brusselator_op::operator()(const FE_function& rhs, FE_function& lhs) const{
  lhs.clear();
  const auto& df_space = lhs.space();
  for(const auto& entity : df_space){
    const auto& rhs_loc = rhs.localFunction(entity);
    auto lhs_loc = lhs.localFunction(entity);
    heat_model(entity, rhs_loc, lhs_loc);
    quad_massMatrix_model(entity, rhs_loc, lhs_loc);
  }  
}
void Brusselator_op::mass_matrix(const FE_function& rhs, FE_function& lhs) const{
  lhs.clear();
  const auto& df_space = lhs.space();
  for(const auto& entity : df_space){
    const auto& rhs_loc = rhs.localFunction(entity);
    auto lhs_loc = lhs.localFunction(entity);
    const auto& g = entity.geometry();
    Quadrature q {entity, lhs_loc.order() + rhs_loc.order()};
    for(std::size_t pt = 0; pt < q.nop(); ++pt){
      const auto& x = q.point(pt);
      const auto weight = q.weight(pt) * g.integrationElement(x);
      auto Mu = massMatrix_loc(pt, q, lhs_loc);

      Mu *= weight;
      lhs_loc.axpy(q[pt], Mu);
    }
  } 
}
void Brusselator_op::massMatrix_constOne(FE_function& lhs) const{
  lhs.clear();
  const auto& df_space = lhs.space();
  for(const auto& entity : df_space){
    auto lhs_loc = lhs.localFunction(entity);
    Quadrature q {entity, lhs_loc.order()};
    const auto& g = entity.geometry();
    for(std::size_t pt = 0; pt < q.nop(); ++pt){
      const auto& x = q.point(pt);
      const auto weight = q.weight(pt) * g.integrationElement(x);

      Range Mu {massMatrix_rhs * weight};
      lhs_loc.axpy(q[pt], Mu);
    } 
  }
}
// ----------------------------------------------------------------------
// Internal Implementation


/*! Log:
 */
