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
#include <config.h>
#include <dune/fem/operator/common/operator.hh>
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

      /*! \name Functions returning dimensions
       */
      //@{
      constexpr int dim_vec_domain() { return Domain<Vector_fef>::dimension; }
      constexpr int dim_vec_range() { return Range<Vector_fef>::dimension; }
      constexpr int dim_scalar_domain() { return Domain<Scalar_fef>::dimension; }
      constexpr int dim_scalar_range() { return Range<Scalar_fef>::dimension; }
      //@}

      // ------------------------------------------------------------
      // Typedef for helper functions
      using Geometry 
      = Vector_fef::DiscreteFunctionSpaceType::IteratorType::Entity::Geometry;
      using Grid_part = Vector_fef::GridPartType;
      using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;
      template<typename T>
      using Jacobian_range = typename Local_function<T>::JacobianRangeType;

      static_assert(dim_vec_domain() == dim_vec_range(), "Bad dimension");
      static_assert(dim_scalar_range() == 1, "Bad dimension");
    private:
      void mcf_lhs_matrixFree_assembly(const Geometry&,
				       const Quadrature&,
				       const Local_function<Vector_fef>&,
				       Local_function<Vector_fef>&);
      /*!< \brief Used by `MCF_op::operator()`
	
	###Pseudo code

	`(M + (alpha + epsilon * tau) A) X`.
       */
      void mcf_rhs_matrixFree_assembly(const Geometry&,
				       const Quadrature&,
				       const Local_function<Vector_fef>&,
				       const Local_function<Scalar_fef>& u_loc,
				       Local_function<Vector_fef>&);
      /*!< \brief Used by `MCF_op::rhs`

	###Pseudo code

	`(M + alpha * A) X + tau * delta * M(u^n, surfaceNormal)`
       */
      RangeField euclidean_norm(const Range<Vector_fef>&);
      /*!< \brief Euclidean norm for dimension 3 vectors. */
      
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
    MCF_op::Range<Vector_fef>
    massMatrix_surfaceNormal(const std::size_t pt,
			     const MCF_op::Quadrature&,
			     const MCF_op::Local_function<Scalar_fef>& u_loc);
    inline MCF_op::Range<Vector_fef>
    mass_matrix(const std::size_t pt,
		const MCF_op::Quadrature& q,
		const MCF_op::Local_function<Vector_fef>& cf){
      MCF_op::Range<Vector_fef> X;
      cf.evaluate(q[pt], X);
      return X;
    }
    /*!< \brief Calculates `MX` */
    inline MCF_op::Jacobian_range<Vector_fef>
    stiffness_matrix(const std::size_t pt,
		     const MCF_op::Quadrature& q,
		     const MCF_op::Local_function<Vector_fef>& cf){
      MCF_op::Jacobian_range<Vector_fef> nabla_X;
      cf.jacobian(q[pt], nabla_X);
      return nabla_X;
    }
    /*!< \brief Calculates `AX`*/
    //@}

    // ----------------------------------------------------------------------
    // Implementation class inline function

    inline RangeField MCF_op::euclidian_norm(const Range<Vector_fef>& v){
      static_assert( dim_vec_range() == 3, "Bad dimension");
      const auto rv = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
      Assert::dynamic<Assert::level(7), Esfem::SolutionDriven_error>
	(rv > eps, __FILE__, __LINE__, "Norm is too small.");
      return rv;
    }
}

#endif // SECORD_OP_SOLUTIONDRIVEN_IMPL_H
