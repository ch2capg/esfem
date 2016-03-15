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

#ifdef DEBUG
#include <iostream>
#endif

using namespace std;
using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using FE_space = FE_function::DiscreteFunctionSpaceType;
using Entity = FE_space::IteratorType::Entity;
using Geometry = Entity::Geometry;
using Grid_part = FE_function::GridPartType;
using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;
using Local_function = FE_function::LocalFunctionType;
using Domain = Local_function::DomainType;
using Range = Local_function::RangeType;
using Jacobian_range = Local_function::JacobianRangeType;
using Range_field = Local_function::RangeFieldType;

static constexpr int dim_domain
= Esfem::Grid::Grid_and_time::Function_space::dimDomain;
static constexpr int dim_range
= Esfem::Grid::Grid_and_time::Function_space::dimRange;

using Linear_operator = Dune::Fem::ISTLLinearOperator<FE_function, FE_function>;
using Local_matrix = Linear_operator::LocalMatrixType;
using Diffusion_tensor = Dune::FieldMatrix<Range_field, dim_domain, dim_domain>;

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
inline Range jacobian_quadMass_loc(const std::size_t pt, const Quadrature& q,
				   const Local_function& arg1_loc,
				   const Local_function& arg2_loc){
  Range u;
  Range tmp;
  arg1_loc.evaluate(q[pt], u);
  arg2_loc.evaluate(q[pt], tmp);
  u *= tmp;
  return u;
}

void prepare_linearOperator_matrix(const FE_space&, Linear_operator&);
std::size_t calculate_matrix_row_size(const FE_function&);

// ----------------------------------------------------------------------
// Implementation of secOrd_op_brusselator_impl.h

struct Brusselator_op::Data{
  const FE_function& first_arg;
  const FE_function& second_arg;
  std::vector<Range> vec_range {};
  std::vector<Jacobian_range> vec_jacRange {};
  double massMatrix_lhs {0.};
  double stiffnessMatrix_lhs {0.};
  double quadMassMatrix_lhs {0.};
  double massMatrix_rhs {0.};
  Data(const Esfem::Io::Parameter&, const Esfem::Grid::Grid_and_time&,
       const Esfem::Growth, const FE_function&, const FE_function&);
};

Brusselator_op::Data::Data(const Esfem::Io::Parameter& p,
			   const Esfem::Grid::Grid_and_time& gt,
			   const Esfem::Growth type,
			   const FE_function& quadMassMatrix_firstArg,
			   const FE_function& quadMassMatrix_secondArg)
  : first_arg {quadMassMatrix_firstArg},
    second_arg {quadMassMatrix_secondArg},
    vec_range(calculate_matrix_row_size(first_arg)),
    vec_jacRange(calculate_matrix_row_size(first_arg))
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

Brusselator_op::Brusselator_op(const Esfem::Io::Parameter& p,
			       const Esfem::Grid::Grid_and_time& gt,
			       const Esfem::Growth type,
			       const FE_function& quadMassMatrix_firstArg,
			       const FE_function& quadMassMatrix_secondArg)
try : d_ptr {make_unique<Data>
	  (p, gt, type, 
	   quadMassMatrix_firstArg, quadMassMatrix_secondArg)}
{}
catch(const std::exception&){
   std::throw_with_nested(std::logic_error {"Error in constructor of "
  "Brusselator_op."});
}
catch(...){
  throw std::logic_error {"Unkown error in constructor of Brusselator_op."};
}
Brusselator_op::~Brusselator_op() = default;
// {
//   delete d_ptr;
//   d_ptr = nullptr;
// #ifdef DEBGUG
//   std::cerr << "~Brusselator_op(): delete d_ptr.\n";
// #endif
// }
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
      auto Mu = massMatrix_loc(pt, q, rhs_loc);

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

      Range Mu {d_ptr -> massMatrix_rhs * weight};
      lhs_loc.axpy(q[pt], Mu);
    } 
  }
}
void Brusselator_op::jacobian(const FE_function& fef, Linear_operator& matrix) const{
  const auto& df_space = fef.space();

  prepare_linearOperator_matrix(df_space, matrix);
  
  for(const auto& entity : df_space){
    const auto fef_loc = fef.localFunction(entity);
    auto matrix_loc = matrix.localMatrix(entity, entity);

    jacobian_matrix_heat(entity, fef_loc, matrix_loc);
    jacobian_matrix_quadMass(entity, fef_loc, matrix_loc);
  }
  matrix.communicate();
}

// ------------------------------------------------------------
// private members

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

    Mu *= d_ptr -> massMatrix_lhs * weight;
    Au *= d_ptr -> stiffnessMatrix_lhs * weight;
    
    lhs_loc.axpy(q[pt], Mu, Au);
  }
}
void Brusselator_op::
quad_massMatrix_model(const Entity& e, const Local_function& rhs_loc,
		      Local_function& lhs_loc) const{
  const auto& g = e.geometry();
  const auto& arg1_loc = d_ptr -> first_arg.localFunction(e);
  const auto& arg2_loc = d_ptr -> second_arg.localFunction(e);
  Quadrature q
  {e, lhs_loc.order() + arg1_loc.order() + arg2_loc.order() + rhs_loc.order()};
  
  for(std::size_t pt = 0; pt < q.nop(); ++pt){
    const auto& x = q.point(pt);    
    const auto integral_factor
      = q.weight(pt) * g.integrationElement(x)
      * d_ptr -> quadMassMatrix_lhs;
    
    auto Mu = quad_massMatrix_loc(pt, q, rhs_loc, arg1_loc, arg2_loc);
    Mu *= integral_factor;
    
    lhs_loc.axpy(q[pt], Mu);
  }
}
void Brusselator_op::
jacobian_matrix_heat(const Entity& entity,
		     const Local_function& fef_loc, Local_matrix& matrix_loc) const{
  auto& vec_range = d_ptr -> vec_range;
  auto& vec_jacRange = d_ptr -> vec_jacRange;
  
  const auto& geometry = entity.geometry();
  const auto& base_set = matrix_loc.domainBasisFunctionSet();
  
  Quadrature q {entity, 2 * fef_loc.order()};
  for(std::size_t pt = 0; pt < q.nop(); ++pt){
    const auto& x = q.point(pt);
    const auto integral_factor = q.weight(pt) * geometry.integrationElement(x);

    base_set.evaluateAll(q[pt], vec_range);
    base_set.jacobianAll(q[pt], vec_jacRange);

    // auto u0 = massMatrix_loc(pt, q, fef_loc);
    // auto jac_u0 = stiffnessMatrix_loc(pt, q, fef_loc);

    for(std::size_t localCol = 0; localCol < base_set.size(); ++localCol){
      Range basisFunction = vec_range[localCol];
      basisFunction *= d_ptr -> massMatrix_lhs;
      Jacobian_range grad_basisFunction = vec_jacRange[localCol];
      grad_basisFunction *= d_ptr -> stiffnessMatrix_lhs;
      matrix_loc
	.column(localCol)
	.axpy(vec_range, vec_jacRange,
	      basisFunction, grad_basisFunction,
	      integral_factor);
    }
  }
}
void Brusselator_op::
jacobian_matrix_quadMass(const Entity& entity,
			 const Local_function& fef_loc, Local_matrix& matrix_loc) const{
  auto& vec_range = d_ptr -> vec_range;
  const auto& arg1_loc = d_ptr -> first_arg.localFunction(entity);
  const auto& arg2_loc = d_ptr -> second_arg.localFunction(entity);

  const auto& geometry = entity.geometry();
  const auto& base_set = matrix_loc.domainBasisFunctionSet();
  
  Quadrature q {entity, 2 * fef_loc.order() + arg1_loc.order() + arg2_loc.order()};
  for(std::size_t pt = 0; pt < q.nop(); ++pt){
    const auto& x = q.point(pt);
    const auto integral_factor
      = q.weight(pt) * geometry.integrationElement(x)
      * d_ptr -> quadMassMatrix_lhs;

    base_set.evaluateAll(q[pt], vec_range);
    const auto u0 = jacobian_quadMass_loc(pt, q, arg1_loc, arg2_loc);

    for(std::size_t localCol = 0; localCol < base_set.size(); ++localCol){
      Range basisFunction = vec_range[localCol] * u0;
      matrix_loc
	.column(localCol)
	.axpy(vec_range, basisFunction, integral_factor);
    }
  }
}

// ----------------------------------------------------------------------
// Internal Implementation

std::size_t calculate_matrix_row_size(const FE_function& fef){
  const auto& df_space = fef.space();
  return df_space.blockMapper().maxNumDofs() * df_space.localBlockSize;
  // localBlockSize is equal to 1 for scalar functions
}
void prepare_linearOperator_matrix(const FE_space& df_space,
				   Linear_operator& matrix){
  Dune::Fem::DiagonalStencil<FE_space, FE_space> stencil {df_space, df_space};
  matrix.reserve(stencil);
  matrix.clear();
}
void jacobian_matrix_heat(const Local_function&, Local_matrix&, const Entity&,
			  std::vector<Range>&, std::vector<Jacobian_range>&){
  
}

/*! Log:
 */
