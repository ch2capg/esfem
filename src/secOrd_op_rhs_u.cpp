/*! \file secOrd_op_rhs_u.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Implementation details for secOrd_op_rhs_u.h
     Created by Christian Power on 04.02.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <cmath>
#include <config.h>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include "secOrd_op_rhs_u.h"
#include "grid.h"
#include "io_parameter.h"

#ifdef DEBUG
#include <iostream>
#endif

using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using Geometry
= FE_function::DiscreteFunctionSpaceType::IteratorType::Entity::Geometry;
using Grid_part = FE_function::GridPartType;
using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;

class RHS_data_u
  : public Dune::Fem::Function<Esfem::Grid::Grid_and_time::Function_space, RHS_data_u>
{
public:
  using Base = Esfem::Grid::Grid_and_time::Function_space;
  using Domain = Base::DomainType;
  using Range = Base::RangeType;
  
  RHS_data_u() = delete;
  explicit RHS_data_u(const Dune::Fem::TimeProviderBase&,
		    const Esfem::Io::Parameter&);
  RHS_data_u(const RHS_data_u&) = delete;
  RHS_data_u& operator=(const RHS_data_u&) = delete;

  void evaluate(const Domain&, Range&) const;
  Range operator()(const Domain&) const;
private:
  const Dune::Fem::TimeProviderBase& tp;
  const double a;
  const double gamma;
};

void assemble_RHS(const RHS_data_u&, FE_function&);
void matrixFree_assembly(const Dune::Fem::TimeProviderBase&, const Geometry&,
			 const Quadrature&, FE_function::LocalFunctionType&);
void massMatrix_for_entity(const Geometry&, const Quadrature&, const RHS_data_u&,
			   FE_function::LocalFunctionType&);

// ----------------------------------------------------------------------
// Implementation esfem.h

struct Esfem::SecOrd_op::Rhs_u::Data{
  RHS_data_u rhs;
  const Dune::Fem::TimeProviderBase& tp;
  const Grid::Grid_and_time::FE_space& fe_space;
  Data(const Io::Parameter& p, const Grid::Grid_and_time& gt)
    : rhs {gt.time_provider(), p}, tp {gt.time_provider()},
      fe_space {gt.fe_space()}
  {}
};

Esfem::SecOrd_op::Rhs_u::Rhs_u(const Io::Parameter& p, const Grid::Grid_and_time& gt){
  try{
    d_ptr = new Data {p, gt};
  }
  catch(std::exception&){
    std::throw_with_nested(std::logic_error
			   {"Error in constructor of Grid_and_time."});
  }
  catch(...){
    throw std::logic_error {"Unknown error in constructor of Grid_and_time."};
  }
}
Esfem::SecOrd_op::Rhs_u::~Rhs_u(){
  delete d_ptr;
  d_ptr = nullptr;
#ifdef DEBUG
  std::cerr << "~Rhs_u(): delete d_ptr.\n";
#endif
}
void Esfem::SecOrd_op::Rhs_u::assemble_and_addScaled_to(Grid::Scal_FEfun& fef) const{
  static FE_function tmp {"local_tmp", d_ptr -> fe_space};
  assemble_RHS(d_ptr -> rhs, tmp);
  FE_function& dune_fef = fef;
  dune_fef.axpy(d_ptr -> tp.deltaT(), tmp); 
}
void Esfem::SecOrd_op::Rhs_u::assemble(Grid::Scal_FEfun& fef) const{
  FE_function& dune_fef = fef;
  assemble_RHS(d_ptr -> rhs, dune_fef);
}

// ----------------------------------------------------------------------
// Internal implementation

RHS_data_u::RHS_data_u(const Dune::Fem::TimeProviderBase& tpb,
		   const Esfem::Io::Parameter& p)
  : tp {tpb}, a {p.tg_a()}, gamma {p.tg_gamma()}
{}
void RHS_data_u::evaluate(const Domain& d, Range& r) const{
  static_assert(Domain::dimension == 3, "Bad domain dimension.");
  static_assert(Range::dimension == 1, "Bad range dimension.");
  const double x = d[0];
  const double y = d[1];
  const double z = d[2];
  const double t = tp.time();
  const double g = gamma;
  const double u = std::exp(-6. * t) * x * y;
  // const double w = std::exp(-6. * t) * y * z;
  r = u
  //r =
    // - g * a + g * u	// rhs for w == 0.
    // tp.deltaT() * g * a 
    // u
    // tp.deltaT() * g * u * u * w
    // #include "/Users/christianpower/cpp/syntax/data/brusselator_u_cpp.txt"
    // #include "/Users/christianpower/cpp/syntax/data/brusselator_u_v2_cpp.txt"
    ;
}
RHS_data_u::Range RHS_data_u::operator()(const Domain& d) const{
  Range r {0};
  evaluate(d,r);
  return r;
}

void assemble_RHS(const RHS_data_u& rhs_fun, FE_function& fef){
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
			   const RHS_data_u& rhs_fun,
			   FE_function::LocalFunctionType& f_loc){
  for(std::size_t pt = 0; pt < q.nop(); ++pt){
    const auto& x = q.point(pt);
    RHS_data_u::Range fx {rhs_fun(g.global(x))};
    fx *= q.weight(pt) * g.integrationElement(x);
    f_loc.axpy(q[pt], fx);
  }  
}


/*! Log:
 */
