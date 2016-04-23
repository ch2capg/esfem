/*! \file secOrd_op_solutionDriven.h
    \brief Surface operator for the solution driven paper

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Idea
     --------------------------------------------------

     Provides a class that perfoms the Dziuk mean curvature flow
     scheme with an additional right-hand side.


    \author Christian Power
    \date 23. April 2016
    \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
*/

#ifndef SECORD_OP_SOLUTIONDRIVEN_H
#define SECORD_OP_SOLUTIONDRIVEN_H

#include <memory>
#include "esfem_fwd.h"

namespace Esfem{
  namespace SecOrd_op{
    class Solution_driven{
    public:
      //! Get PDE parameter, time and growth promoting function \f$u\f$.
      /*! Extracts from `Io::Parameter` the parameter \f$ \alpha \f$ 
	and \f$ \delta \f$.
	`Grid::Grid_and_time` provides dynamical time steps via `time_provider`.
	More precisely we need the method `velocity_regularization`,
	`surface_growthFactor` and `mcf_regularization` from `Io::Parameter`.
	\post Grid_and_time and `u_wrapper` outlive this object. 
       */
      Solution_driven(const Io::Parameter&, const Grid::Grid_and_time&,
		      const Grid::Scal_FEfun& u_wrapper);
      //! Required for the pointer to implementation technique.
      ~Solution_driven();

      //! Solve regularized mean curvature equation
      /*! \param rhs \f$\nodalValue{Y}\f$
	\retval lhs Solution of
	\f$
	\parentheses[\big]{M_3^n + (\alpha + \varepsilon\tau) A_3^n}
	\nodalValue{X}^{n+1} = \nodalValue{Y}
	\f$
	\remark The solver is the conjugated gradient method.
	\pre Parameter `rhs` must be assembled.
	\sa brusselator_rhs(), Vec_rhs
       */      
      void solve(const Grid::Vec_FEfun& rhs, Grid::Vec_FEfun& lhs) const;
      //! Generates right-hand side without load vector for the linear system.
      /*! \param rhs \f$\nodalValue{X}^n\f$
	\retval lhs 
	\f$ 
	(M_3^n + \alpha A_3^n) \nodalValue{X}^n
	+ \tau \delta M_3^n(\nodalValue{u}^n, \nodalValue{\surfaceNormal})
	\f$
	\pre \f$u\f$, which is given in the constructor, is still valid.  
	The nodal values of `rhs` should be the nodes itself.
	\sa Identity
      */
      void brusselator_rhs(const Grid::Vec_FEfun& rhs, Grid::Vec_FEfun& lhs) const;
    private:
      struct Data;
      //! Pointer to data members
      std::unique_ptr<Data> d_ptr;
    };
    /*!< \brief Solve the surface partial differential equation for
                solution driven paper.

     Partial differential equation
     ==================================================

     Parameter
     --------------------------------------------------

     \f$\varepsilon, \delta, a \in \R\f$ for the surface equation, where 
     \f$\varepsilon\f$ and \f$\alpha\f$ are small regularization parameter.

     Smooth problem
     --------------------------------------------------

     For given \f$u\colon \surface \to \R\f$ search for
     \f$X\colon \surface_{0} \times [0,T] \to \R^{m+1}\f$ such that
     \f{align*}{
       v - \alpha \laplaceBeltrami v = {} &  
       \parentheses[\big]{\varepsilon (-\meanCurvature) + \delta u} \surfaceNormal
       = \varepsilon \laplaceBeltrami X + \delta u \surfaceNormal, \\
       \dell_{t} X = {} & v(X).
     \f}

     Finite element discretization
     --------------------------------------------------

     For given \f$\nodalValue{u}\colon I \to \R^{N}\f$ search for
     \f$\nodalValue{X}\colon I \to \R^{3N}\f$ (surface nodal values) such that
     \f{equation*}{
       \parentheses[\big]{M(X) + \alpha A(X)} \dell_t X = 
       \varepsilon A(X)X + \delta M(u,\nodalValue{\surfaceNormal})
     \f}

     Full discretization
     --------------------------------------------------

     For given \f$\nodalValue{X}^n\f$ and \f$\nodalValue{u}^n\f$ solve for 
     \f$\nodalValue{X}^{n+1}\f$
        \f{equation*}{
	  \parentheses[\big]{M_3^n + (\alpha + \varepsilon\tau) A_3^n} 
	  \nodalValue{X}^{n+1} 
	  =  (M_3^n + \alpha A_3^n) \nodalValue{X}^n 
	  + \tau \delta M_3^n(\nodalValue{u}^n, 
	  \nodalValue{\surfaceNormal}),
	\f}
     where \f$\nodalValue{\surfaceNormal}^n\f$ is elementwise normal.  
     */
  } 
} // Esfem

#endif // SECORD_OP_SOLUTIONDRIVEN_H
