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

#include "secOrd_op_rhs_impl.h"

using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using Geometry
= FE_function::DiscreteFunctionSpaceType::IteratorType::Entity::Geometry;
using Grid_part = FE_function::GridPartType;
using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;

// ----------------------------------------------------------------------
// Implementation of RHS_data

RHS_data::RHS_data(const Dune::Fem::TimeProviderBase& tpb)
  : tp {tpb}
{}
void RHS_data::evaluate(const Domain& d, Range& r) const{
  static_assert(Domain::dimension == 3, "Bad domain dimension.");
  static_assert(Range::dimension == 1, "Bad range dimension.");
  const double x = d[0];
  const double y = d[1];
  const double z = d[2];
  const double t = tp.time();
  r = std::exp(-6.*t)*x*y;
  // r =
    // #include "u_rhs.txt"
    // #include "/Users/christianpower/cpp/syntax/data/harmonicSphere_cpp.txt"
    ;
}
RHS_data::Range RHS_data::operator()(const Domain& d) const{
  Range r {0};
  evaluate(d,r);
  return r;
}

// ----------------------------------------------------------------------
// Helper functions

void Esfem::Impl::assemble_RHS(const RHS_data& rhs_fun, FE_function& fef){
  fef.clear();
  const auto& df_space = fef.space();
  for(const auto& entity : df_space){
    const auto& geometry = entity.geometry();
    const Quadrature quad {entity, 2 * df_space.order() + 1};
    auto fef_local = fef.localFunction(entity);
    massMatrix_for_entity(geometry, quad, rhs_fun, fef_local);
  }
  fef.communicate();    
}

void Esfem::Impl::massMatrix_for_entity(const Geometry& g, const Quadrature& q,
					const RHS_data& rhs_fun,
					FE_function::LocalFunctionType& f_loc){
  for(std::size_t pt = 0; pt < q.nop(); ++pt){
    const auto& x = q.point(pt);
    RHS_data::Range fx {rhs_fun(g.global(x))};
    // rhs_fun.evaluate(g.global(x), fx);
    fx *= q.weight(pt) * g.integrationElement(x);
    f_loc.axpy(q[pt], fx);
  }  
}
