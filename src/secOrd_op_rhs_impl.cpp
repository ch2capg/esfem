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
#include <numeric>
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
using Esfem::Impl::sls_rhs;
using Esfem::Impl::sd_rhs;
using Dune::Fem::Parameter;
using namespace std;

// ----------------------------------------------------------------------
// Implementation of Rhs_fun

Rhs_fun::Rhs_fun(const Dune::Fem::TimeProviderBase& tpb, const Growth type)
  :tp {tpb}
{
  const double rE {Parameter::getValue<double>("logistic_growth.r_end", 2.)};
  const double rA {Parameter::getValue<double>("logistic_growth.r_start", 1.)};
  const double k {Parameter::getValue<double>("logistic_growth.steepness", .5)};
  const double Dc {Parameter::getValue<double>("tumor_growth.heat.Dc", 10.)};
  const double ep {Parameter::getValue<double>("tumor_growth.surface.epsilon", .01)};
  const double al {Parameter::getValue<double>("tumor_growth.surface.alpha", 1e-3)};
  const double delta {Parameter::getValue<double>("tumor_growth.surface.delta", .4)};

  dassert(rE > rA, Assert::compose(__FILE__, __LINE__, "r_end <= r_start"));
  dassert(rA > 0, Assert::compose(__FILE__, __LINE__, "r_start < 0"));
  dassert(k > 0, Assert::compose(__FILE__, __LINE__, "Steepness non-positive"));
  // Other parameter are tested in Io::Parameter
  
  switch(type){
  case Growth::promoting:
    fun_impl = [&tp = tp, rA, k, rE]
      (const Domain& d, Range& r){
      const double x = d[0];
      const double y = d[1];
      const double z = d[2];
      const double t = tp.time();
      r = 2*k*x*y*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1)*exp(-6*t) - (k*(pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1) + k*(pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1) + k*(pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1))*x*y*exp(-6*t) + 4*pow(x,3)*y*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 4*x*pow(y,3)*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) - 6*x*y*exp(-6*t) + 2*(pow(x,3)*y*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + (pow(x,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x/(pow(x,2) + pow(y,2) + pow(z,2)))*y*exp(-6*t) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + (4*pow(x,2)*pow(y,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(x,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*x*y/(pow(x,2) + pow(y,2) + pow(z,2)) + (4*pow(x,2)*pow(y,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(y,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*x*y/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(x*pow(y,3)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + x*(pow(y,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - y/(pow(x,2) + pow(y,2) + pow(z,2)))*exp(-6*t) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + 2*(2*pow(x,2)*y*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - y*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*x*z/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(2*x*pow(y,2)*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*y*z/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(2*x*y*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1);
    };
    break;
  case Growth::inhibiting:
    fun_impl = [&tp = tp, rA, k, rE, al, Dc]
      (const Domain& d, Range& r){
      const double x = d[0];
      const double y = d[1];
      const double z = d[2];
      const double t = tp.time();
      r = 2*k*y*z*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1)*exp(-6*t) - (k*(pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1) + k*(pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1) + k*(pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*(rA/(rA*(exp(-k*t) - 1) - rE*exp(-k*t)) + 1))*y*z*exp(-6*t) - 6*y*z*exp(-6*t) + (4*pow(x,2)*pow(y,3)*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 4*pow(x,2)*y*pow(z,3)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),3) + 2*(2*pow(x,2)*y*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - y*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(x,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + 2*(2*x*pow(y,2)*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*x*y/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(pow(y,3)*z*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + (pow(y,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - y/(pow(x,2) + pow(y,2) + pow(z,2)))*z*exp(-6*t) - y*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1) + 2*(2*x*y*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - x*y*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*x*z/(pow(x,2) + pow(y,2) + pow(z,2)) + (4*pow(y,2)*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(y,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(y,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*y*z/(pow(x,2) + pow(y,2) + pow(z,2)) + (4*pow(y,2)*pow(z,2)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - pow(z,2)*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)) - (pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1)*exp(-6*t))*y*z/(pow(x,2) + pow(y,2) + pow(z,2)) + 2*(y*pow(z,3)*exp(-6*t)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) + y*(pow(z,3)/pow(pow(x,2) + pow(y,2) + pow(z,2),2) - z/(pow(x,2) + pow(y,2) + pow(z,2)))*exp(-6*t) - y*z*exp(-6*t)/(pow(x,2) + pow(y,2) + pow(z,2)))*(pow(z,2)/(pow(x,2) + pow(y,2) + pow(z,2)) - 1))*Dc;
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
  : tp {tpb},
  alpha {Parameter::getValue<double>("tumor_growth.surface.alpha", 1e-3)},
  epsilon {Parameter::getValue<double>("tumor_growth.surface.epsilon", .01)},
  r_start {Parameter::getValue<double>("logistic_growth.r_start", 1.)},
  r_end {Parameter::getValue<double>("logistic_growth.r_end", 2.)},
  k {Parameter::getValue<double>("logistic_growth.steepness", .5)},
  delta {Parameter::getValue<double>("tumor_growth.surface.delta", .4)}
{ 
  // no checking has to be done, since other classes already do this
  // cache[0] = tp.time() - 1; // cache[0] != tp.time()
  // update_cache();
}

void Vec_rhs_fun::update_cache() const{
  if(cache[0] == tp.time()) return;
  const double t = tp.time();
  const double e_kt = exp(-k * t);
  const double r_t = r_end * r_start / (r_end * e_kt + r_start * (1 - e_kt) );
  cache[0] = t;
  cache[1] = k * ( 1 - r_t / r_end);
  cache[2] = 2 * ( alpha * cache[1] + epsilon);
  cache[3] = delta * exp(-6 * t);
}

void Vec_rhs_fun::evaluate(const Domain& d, Range& r) const{
  r = 0;
  /*
  const double x = d[0];
  const double y = d[1];
  const double z = d[2];
  const double norm_square = x*x + y*y + z*z;
  const double factor = -(1. + (alpha + epsilon) * 2. / norm_square);
  r[0] = factor * d[0];
  r[1] = factor * d[1];
  r[2] = factor * d[2];
  */
  /*
  update_cache();
  const double factor = cache[2] * abs_d + cache[3] / abs_d - cache[4] * x * y;
  r[0] = factor * x / abs_d;
  r[1] = factor * y / abs_d;
  r[2] = factor * z / abs_d;
  */
}
Vec_rhs_fun::Range Vec_rhs_fun::operator()(const Domain& d) const{
  Range r {0};
  evaluate(d,r);
  return r;
}

sls_rhs::sls_rhs(const Grid::Grid_and_time& gt) 
  :tp {gt.time_provider()},
   lvec {"lvec", gt.vec_fe_space()},
   rA {Parameter::getValue<double>("logistic_growth.r_start", 1.)},
   rE {Parameter::getValue<double>("logistic_growth.r_end", 2.)},
   a {Parameter::getValue<double>("tumor_growth.surface.alpha", 1e-3)},
   e {Parameter::getValue<double>("tumor_growth.surface.epsilon", .01)},
   k {Parameter::getValue<double>("logistic_growth.steepness", .5)}
{}
auto sls_rhs::operator()(const dom& d) const -> ran{
  const auto 
    norm = sqrt(inner_product(&d[0], &d[0]+ dom::dimension, &d[0], 0.)),
    a_til = k * (1 - norm/rE),
    mc = dim / norm,
    factor = a_til + (a * a_til + e ) * mc / norm;
  ran r = d;
  r *= factor;
  return r;
}
void sls_rhs::addScaled_to(Grid::Vec_FEfun& rhs){
  assemble_RHS(*this, lvec);
  Grid::Vec_FEfun::Dune_FEfun& fef = rhs;
  fef.axpy(tp.deltaT(), lvec);
}

sd_rhs::sd_rhs(const Grid::Grid_and_time& gt)
  :tp {gt.time_provider()},
   lvec {"lvec", gt.vec_fe_space()},
   a {Parameter::getValue<double>("tumor_growth.surface.alpha", 1e-3)},
   e {Parameter::getValue<double>("tumor_growth.surface.epsilon", .01)}
{}
auto sd_rhs::operator()(const dom& d) const -> ran{
  const auto fac = (1. + (a+e)*dim*exp(-2. * tp.time()));
  ran r {d};
  r *= fac;
  return r;
}
void sd_rhs::addScaled_to(Grid::Vec_FEfun& rhs){
  assemble_RHS(*this, lvec);
  Grid::Vec_FEfun::Dune_FEfun& fef = rhs;
  fef.axpy(tp.deltaT(), lvec);
}
