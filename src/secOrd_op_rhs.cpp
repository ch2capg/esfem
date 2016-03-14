/*! \file secOrd_op_rhs.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 30. Januar 2016

     Implementation details for secOrd_op_rhs.h
     Created by Christian Power on 30.01.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <cmath>
#include <config.h>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include "secOrd_op_rhs.h"
#include "grid.h"

#ifdef DEBUG
#include <iostream>
#endif

using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using Geometry
= FE_function::DiscreteFunctionSpaceType::IteratorType::Entity::Geometry;
using Grid_part = FE_function::GridPartType;
using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;

class RHS_data
  : public Dune::Fem::Function<Esfem::Grid::Grid_and_time::Function_space, RHS_data>
{
public:
  using Base = Esfem::Grid::Grid_and_time::Function_space;
  using Domain = Base::DomainType;
  using Range = Base::RangeType;
  
  RHS_data() = delete;
  explicit RHS_data(const Dune::Fem::TimeProviderBase&);
  RHS_data(const RHS_data&) = delete;
  RHS_data& operator=(const RHS_data&) = delete;

  void evaluate(const Domain&, Range&) const;
  Range operator()(const Domain&) const;
private:
  const Dune::Fem::TimeProviderBase& tp;
};

void assemble_RHS(const RHS_data&, FE_function&);
void matrixFree_assembly(const Dune::Fem::TimeProviderBase&, const Geometry&,
			 const Quadrature&, FE_function::LocalFunctionType&);
void massMatrix_for_entity(const Geometry&, const Quadrature&, const RHS_data&,
			   FE_function::LocalFunctionType&);

// ----------------------------------------------------------------------
// Implementation esfem.h

struct Esfem::SecOrd_op::Rhs::Data{
  RHS_data rhs;
  const Dune::Fem::TimeProviderBase& tp;
  const Grid::Grid_and_time::FE_space& fe_space;
  Data(const Grid::Grid_and_time& gt)
    : rhs {gt.time_provider()}, tp {gt.time_provider()},
      fe_space {gt.fe_space()}
  {}
};

Esfem::SecOrd_op::Rhs::Rhs(const Grid::Grid_and_time& gt){
  try{
    d_ptr = new Data {gt};
  }
  catch(std::exception&){
    std::throw_with_nested(std::logic_error
			   {"Error in constructor of Grid_and_time."});
  }
  catch(...){
    throw std::logic_error {"Unknown error in constructor of Grid_and_time."};
  }
}
Esfem::SecOrd_op::Rhs::~Rhs(){
  delete d_ptr;
  d_ptr = nullptr;
#ifdef DEBUG
  std::cerr << "~Rhs(): delete d_ptr.\n";
#endif
}
void Esfem::SecOrd_op::Rhs::assemble_and_addScaled_to(Grid::Scal_FEfun& fef){
  static FE_function tmp {"local_tmp", d_ptr -> fe_space};
  assemble_RHS(d_ptr -> rhs, tmp);
  FE_function& dune_fef = fef;
  dune_fef.axpy(d_ptr -> tp.deltaT(), tmp); 
}

// ----------------------------------------------------------------------
// Internal implementation

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
  // r = std::exp(-6.*t)*x*y;
  r =
    // #include "u_rhs.txt"
#include "/Users/christianpower/cpp/syntax/data/harmonicSphere_cpp.txt"
    ;
}
RHS_data::Range RHS_data::operator()(const Domain& d) const{
  Range r {0};
  evaluate(d,r);
  return r;
}

void assemble_RHS(const RHS_data& rhs_fun, FE_function& fef){
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

void massMatrix_for_entity(const Geometry& g, const Quadrature& q,
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


/*! Log:
 */
