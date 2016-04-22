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
     \date 22. April 2016
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
      const double z = d[2];
      const double t = tp_ptr -> time();
      r = -((pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(-0.500000000000000/(exp(-0.500000000000000*t) + 1) + 0.500000000000000) + (pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(-0.500000000000000/(exp(-0.500000000000000*t) + 1) + 0.500000000000000) + (pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(-0.500000000000000/(exp(-0.500000000000000*t) + 1) + 0.500000000000000))*x*y*exp(-6*t) + 4*pow(x,3)*y*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 4*x*pow(y,3)*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 2*x*y*(-0.500000000000000/(exp(-0.500000000000000*t) + 1) + 0.500000000000000)*exp(-6*t) - 6*x*y*exp(-6*t) + 2*(pow(x,3)*y*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + (pow(x,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x/(pow(x,2) + pow(y,2) + pow(z,2)))*y*exp(-6*t) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + (4*pow(x,2)*pow(y,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(x,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*x*y/(pow(x,2) + pow(y,2) + pow(z,2)) + (4*pow(x,2)*pow(y,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(y,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*x*y/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(x*pow(y,3)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + x*(pow(y,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - y/(pow(x,2) + pow(y,2) + pow(z,2)))*exp(-6*t) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + 2*(2*pow(x,2)*y*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - y*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*x*z/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(2*x*pow(y,2)*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*y*z/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(2*x*y*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1);
    };
    break;
  case Growth::inhibiting:
    fun_impl = [tp_ptr = &tp](const Domain& d, Range& r){
      const double x = d[0];
      const double y = d[1];
      const double z = d[2];
      const double t = tp_ptr -> time();
      r = -((pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(-0.500000000000000/(exp(-0.500000000000000*t) + 1) + 0.500000000000000) + (pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(-0.500000000000000/(exp(-0.500000000000000*t) + 1) + 0.500000000000000) + (pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(-0.500000000000000/(exp(-0.500000000000000*t) + 1) + 0.500000000000000))*y*z*exp(-6*t) + 4*pow(x,2)*pow(y,3)*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 4*pow(x,3)*y*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 4*pow(x,2)*y*pow(z,3)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 2*y*z*(-0.500000000000000/(exp(-0.500000000000000*t) + 1) + 0.500000000000000)*exp(-6*t) - 6*y*z*exp(-6*t) + 2*(pow(x,3)*y*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + (pow(x,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x/(pow(x,2) + pow(y,2) + pow(z,2)))*y*exp(-6*t) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + (4*pow(x,2)*pow(y,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(x,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*x*y/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(pow(y,3)*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + (pow(y,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - y/(pow(x,2) + pow(y,2) + pow(z,2)))*z*exp(-6*t) - y*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + (4*pow(y,2)*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(y,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*y*z/(pow(x,2) + pow(y,2) + pow(z,2)) + (4*pow(y,2)*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(z,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*y*z/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(y*pow(z,3)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + y*(pow(z,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - z/(pow(x,2) + pow(y,2) + pow(z,2)))*exp(-6*t) - y*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1);
    };
    break;
  default:
    throw Rhs_error {Assert::compose(__FILE__, __LINE__, "Rhs_fun()")};
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
  const double z = d[2];
  const double t = tp.time();
  r[0] = 2*(-0.400000000000000*x*y*exp(-6*t) + (2.00000000000000/(exp(-0.500000000000000*t) + 1) + 0.0100000000000000*exp(-0.500000000000000*t) + 0.0100000000000000)*(-0.500000000000000/(exp(-0.500000000000000*t) + 1) + 0.500000000000000) + 0.0100000000000000*exp(-0.500000000000000*t) + 0.0100000000000000)*x/sqrt(4*pow(x,2) + 4*pow(y,2) + 4*pow(z,2));
  r[1] = 2*(-0.400000000000000*x*y*exp(-6*t) + (2.00000000000000/(exp(-0.500000000000000*t) + 1) + 0.0100000000000000*exp(-0.500000000000000*t) + 0.0100000000000000)*(-0.500000000000000/(exp(-0.500000000000000*t) + 1) + 0.500000000000000) + 0.0100000000000000*exp(-0.500000000000000*t) + 0.0100000000000000)*y/sqrt(4*pow(x,2) + 4*pow(y,2) + 4*pow(z,2));
  r[2] = 2*(-0.400000000000000*x*y*exp(-6*t) + (2.00000000000000/(exp(-0.500000000000000*t) + 1) + 0.0100000000000000*exp(-0.500000000000000*t) + 0.0100000000000000)*(-0.500000000000000/(exp(-0.500000000000000*t) + 1) + 0.500000000000000) + 0.0100000000000000*exp(-0.500000000000000*t) + 0.0100000000000000)*z/sqrt(4*pow(x,2) + 4*pow(y,2) + 4*pow(z,2));
}
Vec_rhs_fun::Range Vec_rhs_fun::operator()(const Domain& d) const{
  Range r {0};
  evaluate(d,r);
  return r;
}
