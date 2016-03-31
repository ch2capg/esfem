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
      using Scalar_fef = Esfem::Grid::Scal_FEfun::Dune_FEfun;
      /*!< \brief \f$ f\colon \R^3 \to \R \f$ */
      using Vector_fef = Esfem::Grid::Vec_FEfun::Dune_FEfun;
      /*!< \brief \f$ f\colon \R^3 \to \R^3 \f$ */
      template<typename T>
      using Local_function = typename T::LocalFunctionType;
      /*!< \brief Finite element function restricted to a triangle. */
      template<typename T>
      using Domain = typename Local_function<T>::DomainType;
      template<typename T>
      using Range = typename Local_function<T>::RangeType;
      using RangeField = typename Local_function<Vector_fef>::RangeFieldType;
      /*!< \brief Simply \f$ \R \f$ */

      MCF_op(const Io::Parameter&,
	     const Grid::Grid_and_time&,
	     const Scalar_fef& u_input);
      /*!< \sa Constructor of `Esfem::SecOrd_op::Solution_driven`
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
      void rhs(const Vector_fef& rhs, Vector_fef& lhs) const;
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
      /*!< \brief Geometry means access to the finite element. */
      using Grid_part = Vector_fef::GridPartType;
      /*!< \brief Auxiliary template argument */
      using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;
      template<typename T>
      using Jacobian_range = typename Local_function<T>::JacobianRangeType;
      /*!< \brief Jacobian of the local finite element function
	\warning `Jacobian_range` is used like this: `jac[0][1]`.
      */
      
      static_assert(dim_vec_domain == dim_vec_range, "Bad dimension");
      static_assert(dim_scalar_range == 1, "Bad dimension");
    private:
      void mcf_lhs_matrixFree_assembly(const Geometry&,
				       const Quadrature&,
				       const Local_function<Vector_fef>&,
				       Local_function<Vector_fef>&) const;
      /*!< \brief Used by operator()
	
	The matrix vector formulation reads as
	\f{equation*}{
	  \parentheses[\big]{M + (\alpha + \varepsilon \tau) A}
	  \nodalValue{X}
	\f}
       */
      void mcf_rhs_matrixFree_assembly(const Geometry&,
				       const Quadrature&,
				       const Local_function<Vector_fef>&,
				       const Local_function<Scalar_fef>& u_loc,
				       Local_function<Vector_fef>&) const;
      /*!< \brief Used by rhs()

	The matrix vector formulation reads as
	\f{equation*}{
	  (M + \alpha A) \nodalValue{X} + \tau \delta 
	  M(\nodalValue{u}, \surfaceNormal)
	\f}
       */      
      Range<Vector_fef> surface_normal(const Geometry&) const;
      /*!< \brief Calculates the outwards pointing normal vector. */
      
      const double alpha; 
      /*!< \brief Velocity Laplace regularization parameter by Lubich */
      const double delta; /*!< \brief Tumor growth parameter */
      const double epsilon; 
      /*!< \brief Mean curvature regularization parameter by Elliott */
      const double eps; /*!< \brief Generic tolerance */
      const Dune::Fem::TimeProviderBase& tp; /*!< \brief Time step provider */
      const Scalar_fef& u; 
      /*!< \brief Concentration of the growth promoting chemical substance */
    };
    /*!< \sa `Esfem::SecOrd_op::Solution_driven` */

    std::vector<MCF_op::Domain<MCF_op::Vector_fef> >
    oriented_basis(const MCF_op::Geometry&);
    /*!< \brief Calculates an oriented basis for the tangent space.
                Assumes `grid dimension == 2` and
                `world dimension == 3`.
	 \return Two (numerical) vectors which represent an
	         oriented basis of the tangent space.
    */
    MCF_op::Range<MCF_op::Vector_fef> nonUnit_normal
    (const std::vector<MCF_op::Domain<MCF_op::Vector_fef> >& basis);
    /*!< \brief Calculates an non-normalized normal vector via
                the cross product formula.
         \warning Assumes `basis.size() == 2`
	           and that the basis is correctly oriented.
	 \return Outward pointing non unit normal vector.
     */
    inline MCF_op::RangeField
    euclidean_norm(const MCF_op::Range<MCF_op::Vector_fef>& v){
      static_assert(MCF_op::dim_vec_range == 3, "Bad dimension");
      return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    /*!< \brief Euclidean norm for a 3-dimensional vector. */

    /*! \name Evaluation helper functions
        \param pt Quadrature point obtained by `MCF_op::Quadrature`.
        \param q The quadrature object itself
        \param cf The local finit element function

	Evaluates the local finite element function
	on the quadrature points.  I have to pass the point and
	the quadrature
	object itself since I do not know what its precise type is.
	I do not consider using template a good practise for this
	problem.

	\todo Figure out the type of `q[pt]`.
     */
    //@{
    inline MCF_op::Range<MCF_op::Vector_fef>
    evaluate(const std::size_t pt,
	     const MCF_op::Quadrature& q,
	     const MCF_op::Local_function<MCF_op::Vector_fef>& cf){
      MCF_op::Range<MCF_op::Vector_fef> X;
      cf.evaluate(q[pt], X);
      return X;
    }
    /*! \brief Returns `X(p)` */
    inline MCF_op::Range<MCF_op::Scalar_fef>
    evaluate(const std::size_t pt,
	     const MCF_op::Quadrature& q,
	     const MCF_op::Local_function<MCF_op::Scalar_fef>& cf){
      MCF_op::Range<MCF_op::Scalar_fef> u;
      cf.evaluate(q[pt], u);
      return u;
    }
    /*!< \brief Returns `u(p)` */
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
