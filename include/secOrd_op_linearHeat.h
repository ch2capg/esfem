/*! \file secOrd_op_linearHeat.h
    \author Christian Power
    \date 17. March 2016

    \brief Standard heat resp. diffusion equation for an 
           evolving or stationary surface problem

     Revision history
     --------------------------------------------------

          Revised by Christian Power dd.mm.yyyy
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

         This programm implements a basic expression calculator.
         Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 28.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef SECORD_OP_LINEARHEAT_H
#define SECORD_OP_LINEARHEAT_H 

#include <memory>
#include "esfem_fwd.h"

namespace Esfem{
  namespace SecOrd_op{
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
