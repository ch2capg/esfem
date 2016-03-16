/*! \file secOrd_op_identity_impl.h
    \author Christian Power
    \date 9. March 2016

    \brief The actuall dune class and function that is behind secOrd_op_identity.h

     Revision history:
     --------------------------------------------------

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Idea
     --------------------------------------------------

     Here comes the dune implementation.

         Created by Christian Power on 16.03.2016
         Copyright (c) 2016 Christian Power.  All rights reserved.
     
*/

#ifndef SECORD_OP_IDENTITY_IMPL_H
#define SECORD_OP_IDENTITY_IMPL_H

#include <config.h>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include "grid.h"

namespace Esfem{
  namespace Impl{
    struct Identity_impl
      : Dune::Fem::Function<Esfem::Grid::Grid_and_time::Function_space,      
			    Identity_impl>
    {
      using FE_space = Esfem::Grid::Grid_and_time::Function_space;
      using Domain = FE_space::DomainType;
      using Range = FE_space::RangeType;
      static_assert(Domain::dimension == Range::dimension, "Bad dimension");

      void evaluate(const Domain& d, Range& r) const { r = d; }
    };
  }
}

#endif // SECORD_OP_IDENTITY_IMPL_H
