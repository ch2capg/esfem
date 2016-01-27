/*! \file operator_linearHeat.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 25. Januar 2016

     Implementation details for operator_linearHeat.h
     Created by Christian Power on 25.01.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */
#include "operator.h"
#include <dune/fem/quadrature/cachingquadrature.hh>

#ifdef DEBUG
#include <iostream>
#endif

using Discrete_function::FE_function;
using Geometry
= FE_function::DiscreteFunctionSpaceType::IteratorType::Entity::Geometry;
using Grid_part = FE_function::GridPartType;
using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;
using Domain = FE_function::LocalFunctionType::DomainType;
using Range = FE_function::LocalFunctionType::RangeType;
using Jacobian_range = FE_function::LocalFunctionType::JacobianRangeType;

void matrixFree_assembly(const double dT, const Geometry&, const Quadrature&,
			 const FE_function::LocalFunctionType&,
			 FE_function::LocalFunctionType&);
void massMatrixFree_assembly(const Geometry&, const Quadrature&,
			     const FE_function::LocalFunctionType&,
			     FE_function::LocalFunctionType&);
inline
Range mass_matrix(const std::size_t pt, const Quadrature& q,
		  const FE_function::LocalFunctionType& cf){
  Range u;
  cf.evaluate(q[pt], u);
  return u;
}
inline
Jacobian_range stiffness_matrix(const std::size_t pt, const Quadrature& q,
				const FE_function::LocalFunctionType& cf){
  Jacobian_range nabla_u;
  cf.jacobian(q[pt], nabla_u);
  return nabla_u;
}

// ----------------------------------------------------------------------
// Implementaion operator.h

struct Operator::Linear_heat::Data{
  double bdf_alpha_lead_coeff {1.};
  const Dune::Fem::TimeProviderBase& time_provider;
};

Operator::Linear_heat::Linear_heat(const Parameter::PDE_data& d,
				   const Dune::Fem::TimeProviderBase& tp)
  : d_ptr {new Data {d.bdf_alphas().back(), tp}}
    {}
Operator::Linear_heat::~Linear_heat(){
  delete d_ptr;
#ifdef DEBUG 
  std::cerr << "~Linear_heat(): deleted d_ptr\n";
#endif
}
void Operator::Linear_heat::operator()(const Discrete_function::FE_function& cfef,
				       Discrete_function::FE_function& fef) const{
  fef.clear();
  const auto& df_space = fef.space();
  for(const auto& entity : df_space){
    const auto& geometry = entity.geometry();
    const auto cfef_loc = cfef.localFunction(entity);
    auto fef_loc = fef.localFunction(entity);
    Quadrature quad {entity, cfef_loc.order() + fef_loc.order()};
    matrixFree_assembly(d_ptr -> time_provider.deltaT(), geometry, quad,
			cfef_loc, fef_loc);
  }
}
void Operator::Linear_heat::mass_matrix(const Discrete_function::FE_function& cfef,
					Discrete_function::FE_function& fef) const{
  fef.clear();
  const auto& df_space = fef.space();
  for(const auto& entity : df_space){
    const auto& geometry = entity.geometry();
    const auto cfef_loc = cfef.localFunction(entity);
    auto fef_loc = fef.localFunction(entity);
    Quadrature quad {entity, cfef_loc.order() + fef_loc.order()};
    massMatrixFree_assembly(geometry, quad, cfef_loc, fef_loc);
  }
}

// ----------------------------------------------------------------------
// Internal implementation

void matrixFree_assembly(const double dT, const Geometry& g,
			 const Quadrature& q,
			 const FE_function::LocalFunctionType& cf,
			 FE_function::LocalFunctionType& f){
  static_assert(Domain::dimension == 3, "Bad domain dimension.");
  static_assert(Range::dimension == 1, "Bad range dimension.");
  
  for(std::size_t pt = 0; pt < q.nop(); ++pt){
    // Lu = (M + tau A)u
    auto Mu = mass_matrix(pt, q, cf);
    auto Au = stiffness_matrix(pt, q, cf);

    const auto& x = q.point(pt);
    const auto weight = q.weight(pt) * g.integrationElement(x);

    Mu *= weight;
    Au *= dT * weight;
    f.axpy(q[pt], Mu, Au);
  }
}
void massMatrixFree_assembly(const Geometry& g,
			     const Quadrature& q,
			     const FE_function::LocalFunctionType& cf,
			 FE_function::LocalFunctionType& f){
  static_assert(Domain::dimension == 3, "Bad domain dimension.");
  static_assert(Range::dimension == 1, "Bad range dimension.");
  
  for(std::size_t pt = 0; pt < q.nop(); ++pt){
    // Lu = Mu
    auto Mu = mass_matrix(pt, q, cf);

    const auto& x = q.point(pt);
    const auto weight = q.weight(pt) * g.integrationElement(x);

    Mu *= weight;
    f.axpy(q[pt], Mu, 0);
  }
}
/*! Log:
 */
