/*! \file secOrd_op_rhs_impl.cpp
    \brief Implementation for `secOrd_op_rhs_impl.h`

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) April 2016

     Idea
     --------------------------------------------------

     Insert the sage-syntax solution of the right-hand side here.

     \author Christian Power 
     \date 26. April 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#include <cmath>
#include "config.h"
#include <dune/fem/io/parameter.hh>
#include <dassert.h>
#include "esfem_error.h"
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
using Esfem::Impl::Rhs_fun;
using Esfem::Impl::Vec_rhs_fun;
using Dune::Fem::Parameter;
// ----------------------------------------------------------------------
// Implementation of Rhs_fun

Rhs_fun::Rhs_fun(const Dune::Fem::TimeProviderBase& tpb, const Growth type)
  :tp {tpb},
   rE {Parameter::getValue<double>("logistic_growth.r_end", 2.)},
   r0 {Parameter::getValue<double>("logistic_growth.r_start", 1.)},
   k {Parameter::getValue<double>("logistic_growth.steepness", .5)},
   Dc {Parameter::getValue<double>("tumor_growth.heat.Dc", 10.)},
   ep {Parameter::getValue<double>("tumor_growth.surface.epsilon", .01)},
   alpha {Parameter::getValue<double>("tumor_growth.surface.alpha", .01)},
   delta {Parameter::getValue<double>("tumor_growth.surface.delta", .4)}    
{
  dassert(rE > r0, Assert::compose(__FILE__, __LINE__, "r_end <= r_start"));
  dassert(r0 > 0, Assert::compose(__FILE__, __LINE__, "r_start < 0"));
  dassert(k > 0, Assert::compose(__FILE__, __LINE__, "Steepness non-positive"));
  // Other parameter are tested in Io::Parameter
  
  switch(type){
  case Growth::promoting:
    fun_impl = [tp_ptr = &tp, rA = r0, k = k, rE = rE]
      (const Domain& d, Range& r){
      const double x = d[0];
      const double y = d[1];
      const double z = d[2];
      const double t = tp_ptr -> time();
      r = 2*k*x*y*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1)*exp(-6*t) - (k*(pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1) + k*(pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1) + k*(pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1))*x*y*exp(-6*t) + 4*pow(x,3)*y*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 4*x*pow(y,3)*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) - 6*x*y*exp(-6*t) + 2*(pow(x,3)*y*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + (pow(x,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x/(pow(x,2) + pow(y,2) + pow(z,2)))*y*exp(-6*t) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + (4*pow(x,2)*pow(y,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(x,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*x*y/(pow(x,2) + pow(y,2) + pow(z,2)) + (4*pow(x,2)*pow(y,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(y,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*x*y/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(x*pow(y,3)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + x*(pow(y,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - y/(pow(x,2) + pow(y,2) + pow(z,2)))*exp(-6*t) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + 2*(2*pow(x,2)*y*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - y*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*x*z/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(2*x*pow(y,2)*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*y*z/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(2*x*y*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1);
    };
    break;
  case Growth::inhibiting:
    fun_impl = [tp_ptr = &tp, rA = r0, k = k, rE = rE, al = alpha, Dc = Dc]
      (const Domain& d, Range& r){
      const double x = d[0];
      const double y = d[1];
      const double z = d[2];
      const double t = tp_ptr -> time();
      r = 2*k*x*y*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1)*exp(-6*t) - (k*(pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1) + k*(pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1) + k*(pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1))*x*y*exp(-6*t) - 6*x*y*exp(-6*t) + (4*pow(x,3)*y*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 4*x*pow(y,3)*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 2*(pow(x,3)*y*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + (pow(x,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x/(pow(x,2) + pow(y,2) + pow(z,2)))*y*exp(-6*t) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + (4*pow(x,2)*pow(y,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(x,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*x*y/(pow(x,2) + pow(y,2) + pow(z,2)) + (4*pow(x,2)*pow(y,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(y,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*x*y/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(x*pow(y,3)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + x*(pow(y,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - y/(pow(x,2) + pow(y,2) + pow(z,2)))*exp(-6*t) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + 2*(2*pow(x,2)*y*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - y*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*x*z/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(2*x*pow(y,2)*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*y*z/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(2*x*y*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1))*Dc;
      // r = 2*k*y*z*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1)*exp(-6*t) - (k*(pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1) + k*(pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1) + k*(pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1))*y*z*exp(-6*t) - 6*y*z*exp(-6*t) + (4*pow(x,2)*pow(y,3)*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 4*pow(x,3)*y*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 4*pow(x,2)*y*pow(z,3)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 2*(pow(x,3)*y*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + (pow(x,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x/(pow(x,2) + pow(y,2) + pow(z,2)))*y*exp(-6*t) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + (4*pow(x,2)*pow(y,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(x,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*x*y/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(pow(y,3)*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + (pow(y,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - y/(pow(x,2) + pow(y,2) + pow(z,2)))*z*exp(-6*t) - y*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + (4*pow(y,2)*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(y,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*y*z/(pow(x,2) + pow(y,2) + pow(z,2)) + (4*pow(y,2)*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(z,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*y*z/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(y*pow(z,3)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + y*(pow(z,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - z/(pow(x,2) + pow(y,2) + pow(z,2)))*exp(-6*t) - y*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1))*Dc;
    };
    break;
  default:
    throw Rhs_error {Assert::compose(__FILE__, __LINE__, "Rhs_fun()")};
    break;
  };
}

void Rhs_fun::dassert(const bool assertion, const std::string& msg){
  Assert::dynamic<Assert::level(Assert::default_level), Esfem::Parameter_error>
    (assertion, msg);
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
