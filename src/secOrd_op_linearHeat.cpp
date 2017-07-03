/*! \file secOrd_op_linearHeat.cpp
    \brief Some linear elliptic operators 

     Revision history
     --------------------------------------------------

          Revised by Christian Power July 2017
          Originally written by Christian Power
               (power22c@gmail.com) 01. Februar 2016

    \author Christian Power 
    \date 03. Juli 2017
    \copyright Copyright (c) 2017 Christian Power.  All rights reserved.
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
#include "alePaper_impl.h"

using namespace std;
using namespace Esfem;
namespace eo = Esfem::elliptic_operators;
namespace ale = alePaper;
using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using Solver = Dune::Fem::CGInverseOperator<FE_function>;
using Entity = FE_function::DiscreteFunctionSpaceType::IteratorType::Entity;
using Geometry = Entity::Geometry;
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
  fef.communicate();
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

//! Implementing ale_movement()
namespace ale_movement_impl{
  //! Operator for the dune solver
  class Linear_ale_op : public Dune::Fem::Operator<FE_function>{
  public:
    explicit Linear_ale_op(const Esfem::Io::Parameter& p,
			   const Esfem::Grid::Grid_and_time& gt)
      :bdf_alpha_lead_coeff {p.bdf_alphas().back()},
       tp {gt.time_provider()}, tmp_fef {"tmp_fef", gt.fe_space()}
    {}
    Linear_ale_op(const Linear_ale_op&) = delete;
    Linear_ale_op& operator=(const Linear_ale_op&) = delete;

    void operator()(const FE_function& rhs, FE_function& lhs) const override;
    void mass_matrix(const FE_function& rhs, FE_function& lhs) const;
  private:
    double bdf_alpha_lead_coeff {1.};
    const Dune::Fem::TimeProviderBase& tp;
    FE_function tmp_fef;
  };

  //! Cf. Brusselator_scheme::ale_aleMovement()
  class  ale_operator : public eo::elliptic_op{
    Linear_ale_op heat_op;
    Solver heat_solver;
  public:
    ale_operator(const Io::Parameter& p,
		 const Grid::Grid_and_time& gt)
      :heat_op {p, gt}, heat_solver {heat_op, p.eps(), p.eps()}
    {}
    void solve(const Grid::Scal_FEfun& rhs, 
	       Grid::Scal_FEfun& lhs) const override{
      const FE_function& fef1 = rhs;
      FE_function& fef2 = lhs;
      heat_solver(fef1, fef2);
    }
    void mass_matrix(const Grid::Scal_FEfun& rhs, 
		     Grid::Scal_FEfun& lhs) const override{
      const FE_function& rhs_ref = rhs;
      FE_function& lhs_ref = lhs;
      heat_op.mass_matrix(rhs_ref, lhs_ref);
    }
  };  
  
  // Implementing Linear_ale_op

  //! Assemble the ALE matrix
  /*! Formula should be: ∫ Uʰ (Wʰ - Vʰ) · ∇ϕʰ  ∀ϕʰ FE-functions
   */
  inline Jacobian_range ale_matrix(const Domain& x, const double t){
    constexpr int dim = 3;
    auto L = [](const double t){ return 1. + .2 * sin(4.*M_PI*t); };
    auto K = [](const double t){ return .1 + .05 * sin(2.*M_PI*t); };
    auto diff = [dT = 1e-6](auto fun, const double t){ 
      return (fun(t + dT) - fun(t))/dT;
    };
    
    Jacobian_range wmv; // w - v = v_\ale - v_\usual
    wmv[0][0] = x[0] * diff(K, t)/ K(t);
    wmv[0][1] = x[1] * diff(K, t)/ K(t);
    wmv[0][2] = x[2] * diff(L, t)/ L(t);

    ale::levelset_grad lsg;
    ale::vector_type v(3);
    for(int i{}; i < dim; ++i) v(i) = x[i];
    ale::vector_type res(3);
    lsg(v, res, t);

    for(int i{}; i < dim; ++i) wmv[0][i] -= res(i);
    return wmv;
  }
  
  void ale_assembly(const Entity& e, const double t,
		    const double dT, const Geometry& g,
		    const Quadrature& q,
		    const FE_function::LocalFunctionType& cf,
		    FE_function::LocalFunctionType& f){
    static_assert(Domain::dimension == 3, "Bad domain dimension.");
    static_assert(Range::dimension == 1, "Bad range dimension.");
  
    for(std::size_t pt = 0; pt < q.nop(); ++pt){
      // Lu = (M + tau A)u
      auto Mu = mass_matrix(pt, q, cf);
      auto Au = stiffness_matrix(pt, q, cf);
      Domain x_global = g.global(Dune::coordinate(q[pt]));
      auto Bu = ale_matrix(x_global, t);
      Au += Bu; // dT comes later
      
      const auto& x = q.point(pt);
      const auto weight = q.weight(pt) * g.integrationElement(x);

      Mu *= weight;
      Au *= dT * weight;
      f.axpy(q[pt], Mu, Au);
    }
  }
  
  void Linear_ale_op::operator()(const FE_function& cfef, FE_function& fef) const{
    fef.clear();
    const auto& df_space = fef.space();
    for(const auto& entity : df_space){
      const auto& geometry = entity.geometry();
      const auto cfef_loc = cfef.localFunction(entity);
      auto fef_loc = fef.localFunction(entity);
      Quadrature quad {entity, cfef_loc.order() + fef_loc.order()};
      ale_assembly(entity, tp.time(), tp.deltaT(), geometry, quad,
		   cfef_loc, fef_loc);
    }
    fef.communicate();
  }
  void Linear_ale_op::mass_matrix(const FE_function& rhs, FE_function& lhs) const{
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
}

std::unique_ptr<eo::elliptic_op> 
eo::ale_movement(const Io::Parameter& p,
		 const Grid::Grid_and_time& gt){
  return make_unique<ale_movement_impl::ale_operator>(p, gt);
}
