/*! \file secOrd_op_linearHeat.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 01. Februar 2016

     Implementation details for secOrd_op_linearHeat.h
     Created by Christian Power on 01.02.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <stdexcept>
#include <config.h>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/oemsolver.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include "secOrd_op_linearHeat.h"
#include "io_parameter.h"
#include "grid.h"

using namespace std;
using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using Solver = Dune::Fem::CGInverseOperator<FE_function>;
using Geometry
= FE_function::DiscreteFunctionSpaceType::IteratorType::Entity::Geometry;
using Grid_part = FE_function::GridPartType;
using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;
using Domain = FE_function::LocalFunctionType::DomainType;
using Range = FE_function::LocalFunctionType::RangeType;
using Jacobian_range = FE_function::LocalFunctionType::JacobianRangeType;

class Linear_heat_op : public Dune::Fem::Operator<FE_function>{
public:
  explicit Linear_heat_op(const Esfem::Io::Parameter&,
			  const Esfem::Grid::Grid_and_time&);
  Linear_heat_op(const Linear_heat_op&) = delete;
  Linear_heat_op& operator=(const Linear_heat_op&) = delete;

  void operator()(const FE_function& rhs, FE_function& lhs) const override;
  void mass_matrix(FE_function&);
  void mass_matrix(const FE_function& rhs, FE_function& lhs) const;
private:
  double bdf_alpha_lead_coeff {1.};
  const Dune::Fem::TimeProviderBase& tp;
  FE_function tmp_fef;
};

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
// Implementation secOrd_op_linearHeat.h

struct Esfem::SecOrd_op::Linear_heat::Data{  
  Linear_heat_op heat_op;
  Solver heat_solver;
  Data(const Io::Parameter& p, const Grid::Grid_and_time& gt)
    : heat_op {p, gt}, heat_solver {heat_op, p.eps(), p.eps()}
  {}
};
Esfem::SecOrd_op::Linear_heat::Linear_heat(const Io::Parameter& p,
					   const Grid::Grid_and_time& gt)
try : d_ptr {make_unique<Data>(p, gt)}
{}
catch(const std::exception&){
  std::throw_with_nested(std::logic_error{"Error in constructor of Linear_heat."});
 }
 catch(...){
   throw std::logic_error{"Unkown error in constructor of Linear_heat."};
 }

Esfem::SecOrd_op::Linear_heat::~Linear_heat() = default;
// {
//   delete d_ptr;
//   d_ptr = nullptr;
// #ifdef DEBUG
//   std::cerr << "~Linear_heat(): delete d_ptr\n";
// #endif
// }
void Esfem::SecOrd_op::Linear_heat::
solve(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const{
  const FE_function& fef1 = rhs;
  FE_function& fef2 = lhs;
  d_ptr -> heat_solver(fef1, fef2);
}
void Esfem::SecOrd_op::Linear_heat::
apply_massMatrix_to(Grid::Scal_FEfun& sfef) const{
  FE_function& fef = sfef;
  d_ptr -> heat_op.mass_matrix(fef);
}
void Esfem::SecOrd_op::Linear_heat::
mass_matrix(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const{
  const FE_function& rhs_ref = rhs;
  FE_function& lhs_ref = lhs;
  d_ptr -> heat_op.mass_matrix(rhs_ref, lhs_ref);
}
// ----------------------------------------------------------------------
// Internal Implementation

Linear_heat_op::Linear_heat_op(const Esfem::Io::Parameter& p,
			       const Esfem::Grid::Grid_and_time& gt)
  : bdf_alpha_lead_coeff {p.bdf_alphas().back()},
    tp {gt.time_provider()}, tmp_fef {"tmp_fef", gt.fe_space()}
{}
void Linear_heat_op::operator()(const FE_function& cfef, FE_function& fef) const{
  fef.clear();
  const auto& df_space = fef.space();
  for(const auto& entity : df_space){
    const auto& geometry = entity.geometry();
    const auto cfef_loc = cfef.localFunction(entity);
    auto fef_loc = fef.localFunction(entity);
    Quadrature quad {entity, cfef_loc.order() + fef_loc.order()};
    matrixFree_assembly(tp.deltaT(), geometry, quad,
			cfef_loc, fef_loc);
  }
}
void Linear_heat_op::mass_matrix(FE_function& fef){
  tmp_fef.assign(fef);
  try{
    fef.clear();
    const auto& df_space = fef.space();
    for(const auto& entity : df_space){
      const auto& geometry = entity.geometry();
      const auto& cfef_loc = tmp_fef.localFunction(entity);
      auto fef_loc = fef.localFunction(entity);
      Quadrature quad {entity, fef_loc.order()};
      massMatrixFree_assembly(geometry, quad, cfef_loc, fef_loc);
    }
  }
  catch(...){
    fef.assign(tmp_fef);
    throw std::runtime_error{"Error in Linear_heat::mass_matrix()"};
  }
}
void Linear_heat_op::mass_matrix(const FE_function& rhs, FE_function& lhs) const{
  lhs.clear();
  const auto& df_space = lhs.space();
  for(const auto& entity : df_space){
    const auto& geometry = entity.geometry();
    const auto& cfef_loc = rhs.localFunction(entity);
    auto fef_loc = lhs.localFunction(entity);
    Quadrature quad {entity, fef_loc.order() + cfef_loc.order()};
    massMatrixFree_assembly(geometry, quad, cfef_loc, fef_loc);
  }
}

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
    f.axpy(q[pt], Mu);
  }
}

/*! Log:
 */
