/*! \file secOrd_op_linearHeat.h
    \author Christian Power
    \date 3. July 2017

    \brief Standard heat resp. diffusion equation for an 
           evolving or stationary surface problem

     Revision history
     --------------------------------------------------

          Revised by Christian Power July 2017
          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) 28. Januar 2016

     Idea
     --------------------------------------------------

     Provides a class that performs the standard Dziuk Elliott evolving surface
     finite element discretization with implicit euler time discretization.

     Partial differential equation
     ==================================================

     Parameter
     --------------------------------------------------

     - Right-hand side function f 
     - Function u_0 for start time t_0 
     - Surface velocity v
     
     Smooth problem
     --------------------------------------------------
     
     Search for \f$u\colon \surface \to \R\f$ for
     \f{equation*}{
       \matd u + u \diver(v) - \laplaceBeltrami u = f
     \f}
     
     Finite element discretization
     --------------------------------------------------

     Search for \f$\nodalValue{u}\colon I \to \R^N \f$ for
     \f{equation*}{
       \dell_t \parentheses[\big]{M(t) \nodalValue{u} } + A(t) \nodalValue{u}
       = M(t)\nodalValue{Pf},
     \f}
     where \f$ \nodalValue{Pf} \f$ are the nodal values of the 
     \f$ L^2 \f$-projection of \f$ f \f$.

     Full discretization
     --------------------------------------------------
     
     Given \f$ \nodalValue{u}^n \f$ solve for \f$ \nodalValue{u}^{n+1} \f$
     \f{equation*}{
       (M\nodalValue{u})^{n+1} + \tau (A \nodalValue{u})^{n+1}
       = (M \nodalValue{w})^n + \tau (M \nodalValue{Pf})^{n+1} 
     \f}

     Created by Christian Power on 28.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef SECORD_OP_LINEARHEAT_H
#define SECORD_OP_LINEARHEAT_H 

#include <memory>
#include "esfem_fwd.h"

namespace Esfem{
  //! Some elliptic operators
  /*! They may be interpreted as pseudo elliptic operators, since they
      contain the time step as a parameter */
  namespace elliptic_operators{
    //! Generic elliptic operator
    struct elliptic_op{
      //! Must be done, since it is suppressed
      elliptic_op() = default;
      //! Abstract base class
      elliptic_op(const elliptic_op&) = delete;      
      //! Abstract base class
      elliptic_op& operator=(const elliptic_op&) = delete;
      //! Solve a linear elliptic system
      virtual void solve(const Grid::Scal_FEfun& rhs, 
			 Grid::Scal_FEfun& lhs) const = 0;
      //! Apply mass matrix to lhs
      /*! \post Values of lhs are overwritten. */
      virtual void mass_matrix(const Grid::Scal_FEfun& rhs, 
			       Grid::Scal_FEfun& lhs) const = 0;
      //! Abstract base class
      virtual ~elliptic_op(){}
    };
    std::unique_ptr<elliptic_op> ale_movement(const Io::Parameter&,
					      const Grid::Grid_and_time&);
  }
  namespace SecOrd_op{
    using namespace elliptic_operators;
    class Linear_heat{
    public:
      explicit Linear_heat(const Io::Parameter&, const Grid::Grid_and_time&);
      ~Linear_heat();
      
      void solve(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const;
      void apply_massMatrix_to(Grid::Scal_FEfun&) const;
      void mass_matrix(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const;
    private:
      struct Data;
      std::unique_ptr<Data> d_ptr;
    };
  }
}

#endif // SECORD_OP_LINEARHEAT_H
