/*! \file secOrd_op_solutionDriven_impl.h
    \author Christian Power
    \date 18. March 2016

    \brief Helper classes for `Esfem::SecOrd_op::Solution_driven`

     Revision history
     --------------------------------------------------

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Idea
     --------------------------------------------------

     Implementation of `secOrd_op_solutionDriven.h`

         Created by Christian Power on 18.03.2016
         Copyright (c) 2016 Christian Power.  All rights reserved.

*/

#ifndef SECORD_OP_SOLUTIONDRIVEN_IMPL_H
#define SECORD_OP_SOLUTIONDRIVEN_IMPL_H

#include <cmath>
#include <vector>
#include <config.h>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include "io_parameter.h"
#include "grid.h"

namespace Esfem{
  namespace Impl{
    class MCF_op 
      : public Dune::Fem::Operator<Esfem::Grid::Vec_FEfun::Dune_FEfun>
    {
    public:
      using Vector_fef = Esfem::Grid::Vec_FEfun::Dune_FEfun; 
      using Scalar_fef = Esfem::Grid::Scal_FEfun::Dune_FEfun;
      template<typename T>
      using Local_function = typename T::LocalFunctionType;
      template<typename T>
      using Domain = typename Local_function<T>::DomainType;
      template<typename T>
      using Range = typename Local_function<T>::RangeType;
      using RangeField = typename Local_function<Vector_fef>::RangeFieldType;

      MCF_op(const Io::Parameter&,
	     const Grid::Grid_and_time&,
	     const Scalar_fef& u_input);
      /*!< \brief C.f. constructor of `Esfem::SecOrd_op::Solution_driven`
       */
      MCF_op(const MCF_op&) = delete;
      MCF_op& operator=(const MCF_op&) = delete;
      
      void operator()(const Vector_fef& rhs, Vector_fef& lhs) const override;
      /*!< \brief Applying the mean curvature operator.  
	          The cg solver uses this.

	The precise formulation reads as
        \f{equation*}{
	  \parentheses[\big]{M_3^n + (\alpha + \varepsilon\tau) A_3^n} 
	  \nodalValue{X}^{n+1} 
	\f}
       */
      void rhs(const Vector_fef& rhs, Vector_fef& lhs);
      /*!< \brief Generates rhs for the linear system.

	The new value of the finite element function will be
        \f{equation*}{
	  (M_3^n + \alpha A_3^n) \nodalValue{X}^n 
	  + \tau \delta M_3^n(\nodalValue{u}^n, 
	  \nodalValue{\surfaceNormal}).
	\f}
       */

      /*! \name Variable containing dimensions
       */
      //@{
      static constexpr int dim_vec_domain = Domain<Vector_fef>::dimension;
      static constexpr int dim_vec_range = Range<Vector_fef>::dimension;
      static constexpr int dim_scalar_domain = Domain<Scalar_fef>::dimension;
      static constexpr int dim_scalar_range = Range<Scalar_fef>::dimension;
      //@}

      // ------------------------------------------------------------
      // Typedef for helper functions
      using Geometry 
      = Vector_fef::DiscreteFunctionSpaceType::IteratorType::Entity::Geometry;
      using Grid_part = Vector_fef::GridPartType;
      using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;
      template<typename T>
      using Jacobian_range = typename Local_function<T>::JacobianRangeType;

      static_assert(dim_vec_domain == dim_vec_range, "Bad dimension");
      static_assert(dim_scalar_range == 1, "Bad dimension");
    private:
      void mcf_lhs_matrixFree_assembly(const Geometry&,
				       const Quadrature&,
				       const Local_function<Vector_fef>&,
				       Local_function<Vector_fef>&) const;
      /*!< \brief Used by `MCF_op::operator()`
	
	###Pseudo code

	`(M + (alpha + epsilon * tau) A) X`.
       */
      void mcf_rhs_matrixFree_assembly(const Geometry&,
				       const Quadrature&,
				       const Local_function<Vector_fef>&,
				       const Local_function<Scalar_fef>& u_loc,
				       Local_function<Vector_fef>&) const;
      /*!< \brief Used by `MCF_op::rhs`

	###Pseudo code

	`(M + alpha * A) X + tau * delta * M(u^n, surfaceNormal)`
       */      
      Range<Vector_fef> surface_normal(const Geometry&) const;
      /*!< \brief Calculates the outwards pointing normal vector.
      */
      
      const double alpha; 
      /*!< Velocity Laplace regularization parameter by Lubich */
      const double delta; /*!< Tumor growth parameter */
      const double epsilon; 
      /*!< Mean curvature regularization parameter by Elliott */
      const double eps; /*!< Generic tolerance */
      const Dune::Fem::TimeProviderBase& tp; /*!< Time step provider */
      const Scalar_fef& u; 
      /*!< Concentration of the growth promoting chemical substance */
    };
    /*!< \brief C.f. `Esfem::SecOrd_op::Solution_driven`
     */

    /*! \name Helper functions for `MCF_op`
     */
    //@{
    std::vector<MCF_op::Domain<MCF_op::Vector_fef> >
    oriented_basis(const MCF_op::Geometry&);
    /*!< \brief Calculates an oriented basis for the tangent space.
                Assumes `grid dimension == 2` and
                `world dimension == 3`.
    */
    MCF_op::Range<MCF_op::Vector_fef> nonUnit_normal
    (const std::vector<MCF_op::Domain<MCF_op::Vector_fef> >& basis);
    /* \brief Calculates an non-normalized normal vector via
              the cross product formula.  Assumes `basis.size() == 2`
	      and that the basis is correctly oriented.
     */
    inline MCF_op::RangeField
    euclidean_norm(const MCF_op::Range<MCF_op::Vector_fef>& v){
      static_assert(MCF_op::dim_vec_range == 3, "Bad dimension");
      return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    /*!< \brief Euclidean norm for a 3-dimensional vector. */
    inline MCF_op::Range<MCF_op::Vector_fef>
    evaluate(const std::size_t pt,
	     const MCF_op::Quadrature& q,
	     const MCF_op::Local_function<MCF_op::Vector_fef>& cf){
      MCF_op::Range<MCF_op::Vector_fef> X;
      cf.evaluate(q[pt], X);
      return X;
    }
    inline MCF_op::Range<MCF_op::Scalar_fef>
    evaluate(const std::size_t pt,
	     const MCF_op::Quadrature& q,
	     const MCF_op::Local_function<MCF_op::Scalar_fef>& cf){
      MCF_op::Range<MCF_op::Scalar_fef> u;
      cf.evaluate(q[pt], u);
      return u;
    }
    /*!< \brief Returns `X(p)` */
    inline MCF_op::Jacobian_range<MCF_op::Vector_fef>
    jacobian(const std::size_t pt,
	     const MCF_op::Quadrature& q,
	     const MCF_op::Local_function<MCF_op::Vector_fef>& cf){
      MCF_op::Jacobian_range<MCF_op::Vector_fef> nabla_X;
      cf.jacobian(q[pt], nabla_X);
      return nabla_X;
    }
    /*!< \brief Returns `dX(p)`*/
    //@}

  }	// namespace Impl
}	// namespace Esfem

#endif // SECORD_OP_SOLUTIONDRIVEN_IMPL_H
