/*! \file grid_ode.h
    \brief ODE integrators for Grid_and_time 

     Revision history
     --------------------------------------------------

          Revised by Christian Power Juli 2017
          Originally written by Christian Power
               (power22c@gmail.com) Juli 2017

    \author Christian Power 
    \date 02. Juli 2017
    \copyright Copyright (c) 2017 Christian Power.  All rights reserved.
 */

#ifndef GRID_ODE_GUARD
#define GRID_ODE_GUARD 1

#include "grid.h"
#include <memory>

namespace Esfem{
  namespace singleStep_integrator{
    //! Single step integrator
    class ss_int{
    public:
      //! \f$\R^3\f$
      using Domain = Esfem::Grid::Deformation::Domain;
      //! \f$\R^3\f$
      using Range = Esfem::Grid::Deformation::Range;
      ss_int() = default;      
      //! Do one integration step
      virtual void integrate(const Domain& x, Range& y, 
			     const double t, const double dT) const = 0;
      //! Abstract base class
      virtual ~ss_int(){}
    };    
    //! Factory function
    /*! My home made implicit Euler method */
    std::unique_ptr<ss_int> capg_ie();
  }
}

#endif // GRID_ODE_GUARD
