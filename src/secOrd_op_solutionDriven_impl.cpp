/*! \file secOrd_op_solutionDriven_impl.cpp
    \author Christian Power
    \date 18. March 2016

    \brief Implementation of `secOrd_op_solutionDriven_impl.h`

     Revision history
     --------------------------------------------------

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Idea
     --------------------------------------------------

     Implementation of helper classes of `secOrd_op_solutionDriven.h`.

         Created by Christian Power on 18.03.2016
         Copyright (c) 2016 Christian Power.  All rights reserved.

*/

using std::size_t;
using Esfem::Impl::MCF_op; 
using Vector_fef = MCF_op::Vector_fef;
using Scalar_fef = MCF_op::Scalar_fef;
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

// ----------------------------------------------------------------------
// Implementation of Esfem::Impl::MCF_op

// ------------------------------------------------------------
// Public interface

MCF_op::MCF_op(const Io::Parameter& p,
	       const Grid::Grid_and_time& g,
	       const Scalar_fef& u_input)
  : alpha {p.velocity_regularization()},
  delta {p.surface_growthFactor()},
  epsilon {p.mcf_regularization()},
  eps {p.eps()},
  tp {g.time_provider()},
  u {u_input}
{}

void MCF_op::operator()(const Vector_fef& rhs, Vector_fef& lhs) const{
  // (M + (alpha + epsilon * tau) A) X
  lhs.clear();
  const auto& df_space = fef.space();
  for(const auto& entity : df_space){
    const auto& geometry = entity.geometry();
    const auto rhs_loc = rhs.localFunction(entity);
    auto lhs_loc = lhs.localFunction(entity);
    Quadrature quad {entity, rhs_loc.order() + lhs_loc.order()};
    mcf_lhs_matrixFree_assembly(geometry, quad, rhs_loc, lhs_loc);
  }
  lhs.communicate();
}

void MCF_op::rhs(const Vector_fef& rhs, Vector_fef& lhs){
  // (M + alpha * A) X + tau * delta * M(u^n, surfaceNormal)
  lhs.clear();
  const auto& df_space = fef.space();
  for(const auto& entity : df_space){
    const auto& geometry = entity.geometry();
    const auto& rhs_loc = rhs.localFunction(entity);
    const auto& u_loc = u.localFunction(entity);	// Scalar valued 
    auto lhs_loc = lhs.localFunction(entity);
    Quadrature quad {entity, fef_loc.order()};
    mcf_rhs_matrixFree_assembly(geometry, quad, rhs_loc, u_loc, lhs_loc);
  }
  lhs.communicate();
}

// ------------------------------------------------------------
// Private interface

void MCF_op::mcf_lhs_matrixFree_assembly(const Geometry& g,
					 const Quadrature& q,
					 const Local_function<Vector_fef>& cf,
					 Local_function<Vector_fef>& f){
  for(size_t pt = 0; pt < q.nop(); ++pt){
    // (M + (alpha + epsilon * tau) A) X
    const auto& x = q.point(pt);
    const auto integral_factor = q.weight(pt) * g.integrationElement(x);

    auto MX = integral_factor * mass_matrix(pt, q, cf);
    auto AX = integral_factor * (alpha + epsilon * tp.deltaT()) 
      * stiffness_matrix(pt, q, cf);

    f.axpy(q[pt], MX, AX);
  }
}
void mcf_rhs_matrixFree_assembly(const Geometry& g,
				 const Quadrature& q,
				 const Local_function<Vector_fef>& cf,
				 const Local_function<Scalar_fef>& u_loc,
				 Local_function<Vector_fef>& f){
  for(size_t pt = 0; pt < q.nop(); ++pt){
    // (M + alpha * A) X + tau * delta * M(u^n, surfaceNormal)
    const auto& x = q.point(pt);
    const auto integral_factor = q.weight(pt) * g.integrationElement(x);

    auto MX = integral_factor * mass_matrix(pt, q, cf);
    MX += integral_factor * tp.deltaT() * delta 
      * massMatrix_surfaceNormal(pt, q, u_loc);
    auto AX = integral_factor * alpha * stiffness_matrix(pt, q, cf) ;

    f.axpy(q[pt], MX, AX);	
  }
}

void Esfem::Impl::
massMatrixFree_assembly(const Geometry& g,
			const Quadrature& q,
			const Local_function<Vector_fef>& cf,
			Local_function<Vector_fef>& f){
  for(size_t pt = 0; pt < q.nop(); ++pt){
    // MX
    auto Mu = mass_matrix(pt, q, cf);

    const auto& x = q.point(pt);
    const auto integral_factor = q.weight(pt) * g.integrationElement(x);

    Mu *= integral_factor;
    f.axpy(q[pt], Mu);	// maybe error
  }
}

// ----------------------------------------------------------------------
// Implementation of helper functions


