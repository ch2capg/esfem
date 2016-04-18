/*! \file secOrd_op_rhs_impl.cpp
    \brief Implementation for `secOrd_op_rhs_impl.h`

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) April 2016

     Idea
     --------------------------------------------------

     Implementation of class `MyClass`

     \author Christian Power 
     \date 16. April 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#include <cmath>
#include "secOrd_op_rhs_impl.h"

//! Scalar valued dune finite element function
using Scal_fef = Esfem::Grid::Scal_FEfun::Dune_FEfun;
//! Vector valued dune finite element function
using Vec_fef = Esfem::Grid::Vec_FEfun::Dune_FEfun;
//! Geometry for an element
using Geometry
= Scal_fef::DiscreteFunctionSpaceType::IteratorType::Entity::Geometry;
//! Technicality for quadrature 
using Grid_part = Scal_fef::GridPartType;
//! Quadrature on a straight element
using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;
//! Implementing this
using Esfem::Impl::Rhs_fun;
//! Implementing this
using Esfem::Impl::Vec_rhs_fun;

// ----------------------------------------------------------------------
// Implementation of Rhs_fun

Rhs_fun::Rhs_fun(const Dune::Fem::TimeProviderBase& tpb)
  : tp {tpb}
{}
void Rhs_fun::evaluate(const Domain& d, Range& r) const{
  const double x = d[0];
  const double y = d[1];
  // const double z = d[2];
  const double t = tp.time();
  r = std::exp(-6.*t)*x*y;
}
Rhs_fun::Range Rhs_fun::operator()(const Domain& d) const{
  Range r {0};
  evaluate(d,r);
  return r;
}

// ----------------------------------------------------------------------
// Implementation of Vec_rhs_fun

Vec_rhs_fun::Vec_rhs_fun(const Dune::Fem::TimeProviderBase& tpb)
  : tp {tpb}
{}

void Vec_rhs_fun::evaluate(const Domain& d, Range& r) const{
  const double x = d[0];
  const double y = d[1];
  // const double z = d[2];
  const double t = tp.time();
  r[0] = std::exp(-6.*t)*x*y;
  r[1] = 0;
  r[2] = 0;
}
Vec_rhs_fun::Range Vec_rhs_fun::operator()(const Domain& d) const{
  Range r {0};
  evaluate(d,r);
  return r;
}

// ----------------------------------------------------------------------
// Helper functions

// void Esfem::Impl::assemble_RHS(const Rhs_fun& rhs_fun, Scal_fef& fef){
//   fef.clear();
//   const auto& df_space = fef.space();
//   for(const auto& entity : df_space){
//     const auto& geometry = entity.geometry();
//     const Quadrature quad {entity, 2 * df_space.order() + 1};
//     auto fef_local = fef.localFunction(entity);
//     for(std::size_t pt = 0; pt < quad.nop(); ++pt){
//       const auto& x = quad.point(pt);
//       Rhs_fun::Range fx {rhs_fun(geometry.global(x))};
//       // rhs_fun.evaluate(geometry.global(x), fx);
//       fx *= quad.weight(pt) * geometry.integrationElement(x);
//       fef_local.axpy(quad[pt], fx);
//     }  
//   }
//   fef.communicate();    
// }
// void Esfem::Impl::assemble_RHS(const Vec_rhs_fun& rhs_fun, Vec_fef& vfef){
//   vfef.clear();
//   const auto& df_space = vfef.space();
//   for(const auto& entity: df_space){
//     const auto& geometry = entity.geometry();
//     const Quadrature quad {entity, 2 * df_space.order() + 1};
//     auto vfef_local = vfef.localFunction(entity);
//     for(std::size_t pt = 0; pt < quad.nop(); ++pt){
//       const auto& x = quad.point(pt);
//       Vec_rhs_fun::Range fx {rhs_fun(geometry.global(x))};
//       fx *= quad.weight(pt) * geometry.integrationElement(x);
//       vfef_local.axpy(quad[pt], fx);
//     }
//   }
//   vfef.communicate();
// }
