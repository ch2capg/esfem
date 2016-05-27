/*! \file secOrd_op_solutionDriven_impl.cpp
    \brief Implementation of `secOrd_op_solutionDriven_impl.h`

    Revision history
    --------------------------------------------------

         Revised by Christian Power April 2016
         Revised by Christian Power March 2016
         Originally written by Christian Power
              (power22c@gmail.com) March 2016

    Idea
    --------------------------------------------------

    Implementation of helper classes of `secOrd_op_solutionDriven.h`.


    \author Christian Power
    \date 26. April 2016
    \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
*/

#include <dassert.h>
#include "esfem_error.h"
#include "secOrd_op_solutionDriven_impl.h"

using std::size_t;
using Esfem::Impl::MCF_op; 
using Vector_fef = MCF_op::Vector_fef;
using Scalar_fef = MCF_op::Scalar_fef;
//! Vector valued finite element space
using Vec_fe_space = Vector_fef::DiscreteFunctionSpaceType;
template<typename T>
using Local_function = typename T::LocalFunctionType;
template<typename T>
using Domain = typename Local_function<T>::DomainType;
template<typename T>
using Range = typename Local_function<T>::RangeType;
using Geometry = MCF_op::Geometry;
using Grid_part = MCF_op::Grid_part;
using Quadrature = MCF_op::Quadrature;
template<typename T>
using Jacobian_range = typename Local_function<T>::JacobianRangeType;
// Forward typedef
// using Linear_operator = MCF_op::Linear_operator;



// static void prepare_linearOperator_matrix(const Vec_fe_space&, Linear_operator&);
// static std::size_t calculate_matrix_row_size(const Vector_fef&);

// ----------------------------------------------------------------------
// Implementation of Esfem::Impl::MCF_op

// ------------------------------------------------------------
// Public interface

MCF_op::MCF_op(const Io::Parameter& p,
	       const Grid::Grid_and_time& g,
	       const Scalar_fef& u_input)
  :alpha {p.velocity_regularization()},
   delta {p.surface_growthFactor()},
   epsilon {p.mcf_regularization()},
   eps {p.eps()},
   tp {g.time_provider()},
   u {u_input}
{}

void MCF_op::operator()(const Vector_fef& rhs, Vector_fef& lhs) const{
  // (M + (alpha + epsilon * tau) A) X
  lhs.clear();
  const auto& df_space = lhs.space();
  for(const auto& entity : df_space){
    const auto& geometry = entity.geometry();
    const auto rhs_loc = rhs.localFunction(entity);
    auto lhs_loc = lhs.localFunction(entity);
    Quadrature quad {entity, rhs_loc.order() + lhs_loc.order()};
    mcf_lhs_matrixFree_assembly(geometry, quad, rhs_loc, lhs_loc);
  }
  lhs.communicate();
}

void MCF_op::brusselator_rhs(const Vector_fef& rhs, Vector_fef& lhs) const{
  // (M + alpha * A) X + tau * delta * M(u^n, surfaceNormal)
  lhs.clear();
  const auto& df_space = lhs.space();
  for(const auto& entity : df_space){
    const auto& geometry = entity.geometry();
    const auto rhs_loc = rhs.localFunction(entity);
    const auto u_loc = u.localFunction(entity);	// Scalar valued
    auto lhs_loc = lhs.localFunction(entity);
    Quadrature quad {entity, rhs_loc.order() + lhs_loc.order()};
    mcf_rhs_matrixFree_assembly(geometry, quad, rhs_loc, u_loc, lhs_loc);
  }
  lhs.communicate();
}

// void MCF_op::jacobian_matrix(const Vector_fef& fef, Linear_operator& matrix) const{
//   const auto& df_space = fef.space();
//   prepare_linearOperator_matrix(df_space, matrix);
//   for(const auto& entity : df_space){
//     const auto fef_loc = fef.localFunction(entity);
//     auto matrix_loc = matrix.localMatrix(entity, entity);
//     // jacobian_matrix_heat(entity, fef_loc, matrix_loc);
//     // jacobian_matrix_quadMass(entity, fef_loc, matrix_loc);
//   }
//   matrix.communicate();
// }

// ------------------------------------------------------------
// Private interface

void MCF_op::mcf_lhs_matrixFree_assembly(const Geometry& g,
					 const Quadrature& q,
					 const Local_function<Vector_fef>& cf,
					 Local_function<Vector_fef>& f) const{
  for(size_t pt = 0; pt < q.nop(); ++pt){
    // (M + (alpha + epsilon * tau) A) X
    const auto& x = q.point(pt);
    const auto integral_factor = q.weight(pt) * g.integrationElement(x);

    auto X_p = evaluate(pt, q, cf);
    X_p *= integral_factor;
    
    auto dX_p = jacobian(pt, q, cf);
    dX_p *= integral_factor * (alpha + epsilon * tp.deltaT());

    f.axpy(q[pt], X_p, dX_p);
  }
}
void MCF_op::mcf_rhs_matrixFree_assembly(const Geometry& g,
					 const Quadrature& q,
					 const Local_function<Vector_fef>& cf,
					 const Local_function<Scalar_fef>& u_loc,
					 Local_function<Vector_fef>& f) const{
  for(size_t pt = 0; pt < q.nop(); ++pt){
    // (M + alpha * A) X + tau * delta * M(u^n, surfaceNormal)
    const auto& x = q.point(pt);
    const auto integral_factor = q.weight(pt) * g.integrationElement(x);

    auto n_p = surface_normal(g);
      
    const auto u_p = evaluate(pt, q, u_loc);
    n_p *= u_p * tp.deltaT() * delta * integral_factor;

    auto X_p = evaluate(pt, q, cf);
    X_p *= integral_factor;
    X_p += n_p;
    
    auto dX_p = jacobian(pt, q, cf);
    dX_p *= alpha * integral_factor;

    f.axpy(q[pt], X_p, dX_p);	
  }
}

MCF_op::Range<Vector_fef>
MCF_op::surface_normal(const Geometry& g) const{
  static_assert(dim_vec_domain == 3, "Bad dimension");
  const auto basis = oriented_basis(g);
  auto normal = nonUnit_normal(basis);
  const auto norm = euclidean_norm(normal);
  Assert::dynamic<Assert::level(7), Esfem::SolutionDriven_error>
    (norm > eps,
     Assert::compose(__FILE__, __LINE__, "Norm of normal vector almost vanishes."));
  normal /= norm;
  return normal;
}

// ----------------------------------------------------------------------
// Implementation helper function

std::vector<MCF_op::Domain<Vector_fef> >
Esfem::Impl::oriented_basis(const MCF_op::Geometry& g){
  const Domain<Vector_fef> p0 = g.corner(0);
  const Domain<Vector_fef> p1 = g.corner(1);
  const Domain<Vector_fef> p2 = g.corner(2);

  // Oriented basis of tangent space
  // Draw a picture to understand this.
  const Domain<Vector_fef> v = p2 - p0;
  const Domain<Vector_fef> w = p1 - p0;
  return {v, w};
}

MCF_op::Range<Vector_fef> Esfem::Impl::nonUnit_normal
(const std::vector<MCF_op::Domain<Vector_fef> >& basis){
  Assert::dynamic<Assert::level(7), Esfem::SolutionDriven_error>
    (basis.size() == 2, "Dimension of basis bad.");
  const auto& v = basis[0];
  const auto& w = basis[1];

  // Cross product formula
  Range<Vector_fef> normal;
  normal[0] = v[1] * w[2] - v[2] * w[1];
  normal[1] = - v[0] * w[2] + v[2] * w[0];
  normal[2] = v[0] * w[1] - v[1] * w[0];

  return normal;
}

// ----------------------------------------------------------------------
// Internal Implementation

// std::size_t calculate_matrix_row_size(const Vector_fef& fef){
//   const auto& df_space = fef.space();
//   return df_space.blockMapper().maxNumDofs() * df_space.localBlockSize;
//   // localBlockSize is equal to 1 for scalar functions
// }

// void prepare_linearOperator_matrix(const Vec_fe_space& df_space,
// 				   Linear_operator& matrix){
//   Dune::Fem::DiagonalStencil<Vec_fe_space, Vec_fe_space>
//     stencil {df_space, df_space};
//   matrix.reserve(stencil);
//   matrix.clear();
// }
