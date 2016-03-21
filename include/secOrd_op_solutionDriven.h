/*! \file secOrd_op_solutionDriven.h
    \author Christian Power
    \date 17. March 2016

    \brief Surface operator for the solution driven paper

     Revision history
     --------------------------------------------------

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Idea
     --------------------------------------------------

     Provides a class that perfoms the Dziuk mean curvature flow
     scheme with an additional right-hand side.


         Created by Christian Power on 17.03.2016
         Copyright (c) 2016 Christian Power.  All rights reserved.

*/

#ifndef SECORD_OP_SOLUTIONDRIVEN_H
#define SECORD_OP_SOLUTIONDRIVEN_H

#include <memory>
#include "esfem_fwd.h"

namespace Esfem{
  namespace SecOrd_op{
    class Solution_driven{
    public:
      Solution_driven(const Io::Parameter&, const Grid::Grid_and_time&,
		      const Grid::Scal_FEfun& u_wrapper);
      /*!< Extracts from `Io::Parameter` the parameter \f$ \alpha \f$ 
	and \f$ \delta \f$.
	`Grid::Grid_and_time` provides dynamical time steps via `time_provider`.
	`u` is consistent with the \f$ u \f$ in the description above.
	More precisely we need the method `velocity_regularization`,
	`surface_growthFactor` and `mcf_regularization` from `Io::Parameter`.
       */
      ~Solution_driven();

      void solve(const Grid::Vec_FEfun& rhs, Grid::Vec_FEfun& lhs) const;
      /*!< \brief Solve `lsh` in the lineare system `A * lsh = rhs`.
	
	We solve precicely the equation
	    (M_3^n + (\alpha + \varepsilon\tau) A_3^n)*lhs = rhs,
	via conjugated gradient method.
       */
      /* 
	\f{equation*}{
	  \parentheses[\big]{M_3^n + (\alpha + \varepsilon\tau) A_3^n} lhs = rhs,
	\f}
       */
      void rhs(const Grid::Vec_FEfun& rhs, Grid::Vec_FEfun& lhs) const;
      /*!< \brief Generates rhs for the linear system.

	The new value of the finite element function will be
       */
      /*
	\f{equation*}{
	  (M_3^n + \alpha A_3^n) \nodalValue{X}^n
	  + \tau \delta M_3^n(\nodalValue{u}^n, \nodalValue{\surfaceNormal}).
	\f}
       */
    private:
      struct Data;
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
       \dell_t \parentheses[\big]{M(X) w} + D_c A(X) w
       = \gamma \parentheses[\big]{b M(X) \nodalValue{1} 
       - M(X; u,u) w}.
     \f}

     Full discretization
     --------------------------------------------------

     For given \f$\nodalValue{X}^n\f$ and \f$\nodalValue{u}^n\f$ solve for 
     \f$\nodalValue{X}^{n+1}\f$
        \f{equation*}{
	  \parentheses[\big]{M_3^n + (\alpha + \varepsilon\tau) A_3^n} \nodalValue{X}^{n+1} 
	  =  (M_3^n + \alpha A_3^n) \nodalValue{X}^n + \tau \delta M_3^n(\nodalValue{u}^n, 
	  \nodalValue{\surfaceNormal}),
	\f}
     where \f$\nodalValue{\surfaceNormal}^n\f$ is elementwise normal.  
     */
  } 
} // Esfem

#endif // SECORD_OP_SOLUTIONDRIVEN_H
