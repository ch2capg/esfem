/*! \file secOrd_op_rhs_impl.cpp
    \brief Implementation for `secOrd_op_rhs_impl.h`

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) April 2016

     Idea
     --------------------------------------------------

     Implementation of class `MyClass`

     \author Christian Power 
     \date 16. April 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#include <cmath>
#include "secOrd_op_rhs_impl.h"
#include "esfem_error.h"

// // Scalar valued dune finite element function
// using Scal_fef = Esfem::Grid::Scal_FEfun::Dune_FEfun;
// // Vector valued dune finite element function
// using Vec_fef = Esfem::Grid::Vec_FEfun::Dune_FEfun;
// // Geometry for an element
// using Geometry
// = Scal_fef::DiscreteFunctionSpaceType::IteratorType::Entity::Geometry;
// // Technicality for quadrature 
// using Grid_part = Scal_fef::GridPartType;
// // Quadrature on a straight element
// using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;
//! Implementing this
using Esfem::Impl::Rhs_fun;
//! Implementing this
using Esfem::Impl::Vec_rhs_fun;

// ----------------------------------------------------------------------
// Implementation of Rhs_fun

Rhs_fun::Rhs_fun(const Dune::Fem::TimeProviderBase& tpb, const Growth type)
  : tp {tpb}
{
  switch(type){
  case Growth::promoting:
    fun_impl = [tp_ptr = &tp](const Domain& d, Range& r){
      const double x = d[0];
      const double y = d[1];
      // const double z = d[2];
      const double t = tp_ptr -> time();
      r = std::exp(-6.*t) * x * y;
    };
    break;
  case Growth::inhibiting:
    fun_impl = [tp_ptr = &tp](const Domain& d, Range& r){
      // const double x = d[0];
      const double y = d[1];
      const double z = d[2];
      const double t = tp_ptr -> time();
      r = std::exp(-6.*t) * y * z;
    };
    break;
  default:
    throw Rhs_error {"Impossible error in Rhs_fun()"};
    break;
  };
}

// ----------------------------------------------------------------------
// Implementation of Vec_rhs_fun

Vec_rhs_fun::Vec_rhs_fun(const Dune::Fem::TimeProviderBase& tpb)
  : tp {tpb}
{}

void Vec_rhs_fun::evaluate(const Domain& d, Range& r) const{
  const double x = d[0];
  const double y = d[1];
  // const double z = d[2];
  const double t = tp.time();
  r[0] = std::exp(-6.*t)*x*y;
  r[1] = 0;
  r[2] = 0;
}
Vec_rhs_fun::Range Vec_rhs_fun::operator()(const Domain& d) const{
  Range r {0};
  evaluate(d,r);
  return r;
}
