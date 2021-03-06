/*! \file secOrd_op_initData_impl.cpp
    \brief Implementation of secOrd_op_initData_impl.h

     Revision history
     --------------------------------------------------

          Revised by Christian Power June 2016
          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) Februar 2016

     Idea
     --------------------------------------------------

     Implementing `Explicit_initial_data` and `Random_initial_data`.

     \author Christian Power 
     \date 15. June 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <config.h>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include "secOrd_op_initData_impl.h"
#include "io_parameter.h"
#include "esfem_error.h"


using Esfem::Impl::Explicit_initial_data;
using Esfem::Impl::Random_initial_data;
using Esfem::Impl::Analytic_velocity;
using Esfem::Impl::sphere_1EF;
using Esfem::Impl::sphere_2EF;
using Esfem::Impl::sphere_3EF;
using Esfem::Impl::const_fct;
using Esfem::Impl::sphere_eigenFun;
using Esfem::Impl::sphere_mcf_sol;
using Esfem::Impl::sls_iData;
using Esfem::Impl::sls_v_iData;
using Esfem::Impl::sd_iData;
using Dune::Fem::Parameter;
using namespace std;
//! \f$ \R^3 \f$
using Vec_domain = Analytic_velocity::Domain;
//! \f$ \R^3 \f$
using Vec_range = Analytic_velocity::Range;

// ----------------------------------------------------------------------
// Some static inline analytic expression

//! Helper for eoc_velocity()
/*! \param t Current time \f$t\f$
  \returns \f$r(t)/r_{end} = \frac{r_0}{r_{end} e^{-kt} + r_0 (1-e^{-kt})}\f$
*/
static inline double r_div_rEnd(const double t){
  const double r0 = 1., r_end = 2., kt = .5 * t;
  return r0 / (r_end * exp(-kt) + r0 * (1 - exp(-kt)) );
}
//! Velocity for the solution driven 2016 paper
/*! Eoc means that this velocity was used to generate the eoc tables.
  \param t Current time \f$t\f$
  \param d Position  \f$x\f$ at time \f$t\f$
  \retval r \f$v(x,t) = k \Bigl(1 - \frac{r(t)}{r_{end}}\Bigr) x\f$
  \sa r_div_rEnd
 */
static inline void eoc_velocity(const double t, const Vec_domain& d, Vec_range& r){
  const double k = .5, // should be read from a parameter file
    r_d_rEnd = r_div_rEnd(t),
    factor =   k * (1 - r_d_rEnd);
  r[0] = d[0] * factor;
  r[1] = d[1] * factor;
  r[2] = d[2] * factor;
}

// ----------------------------------------------------------------------
// Implementation Explicit_initial_data

Explicit_initial_data::
Explicit_initial_data(const Esfem::Grid::Grid_and_time& gt,
		      const Esfem::Growth type)
  :tp {gt.time_provider()}
{
  switch(type){
  case Growth::promoting:
    fun_impl = // [tp_ptr = &tp]
      [&tp = tp](const Domain& d, Range& r){
      const double x = d[0];
      const double y = d[1];
      const double t = tp.time();
      r = x * y * std::exp(-6. * t);
    };
    break;
  case Growth::inhibiting:
    fun_impl = [&tp = tp](const Domain& d, Range& r){
      const double y = d[1];
      const double z = d[2];
      const double t = tp.time();
      r = y * z * std::exp(-6. * t);
    };
    break;
  default:
    throw InitData_error {Assert::compose(__FILE__, __LINE__, "Bad Growth type")};
    break;
  };
}

// ----------------------------------------------------------------------
// Implementation Random_initial_data

Random_initial_data::
Random_initial_data(const Esfem::Io::Parameter& p,
		    const Esfem::Growth type)
  :Random_initial_data {hom_value(p, type), pertubation(p, type)}
{
  std::cout << print_configuration(p, type) << std::endl;
}
Random_initial_data::
Random_initial_data(const double hom_value,
		    const double pertubation)
  :random_fun {std::bind(Random_dist {hom_value, hom_value + pertubation},
			 Random_engine {})}
{}

// ----------------------------------------------------------------------
// sphere_1EF

sphere_1EF::sphere_1EF(const Grid::Grid_and_time& gt)
  :tp {gt.time_provider()} {}
sphere_1EF* sphere_1EF::clone(){
  return new sphere_1EF {*this};
}
void sphere_1EF::interpolate(Grid::Scal_FEfun& fef) const{
  using Fef = Grid::Scal_FEfun::Dune_FEfun;
  Dune::LagrangeInterpolation<Fef>::interpolateFunction(*this, fef); 
}
void sphere_1EF::evaluate(const domain& x, range& y)const{
  y = x[0]*x[1]*exp(-6*tp.time());
}

// ----------------------------------------------------------------------
// const_fct

const_fct::const_fct(const Grid::Grid_and_time& gt)
  :tp {gt.time_provider()} {}
const_fct* const_fct::clone(){
  return new const_fct {*this};
}
void const_fct::interpolate(Grid::Scal_FEfun& fef) const{
  using Fef = Grid::Scal_FEfun::Dune_FEfun;
  Dune::LagrangeInterpolation<Fef>::interpolateFunction(*this, fef); 
}
void const_fct::evaluate(const domain&, range& y)const{
  y = 1.;
}

// ----------------------------------------------------------------------
// sphere_eigenFun

sphere_eigenFun::sphere_eigenFun(const Grid::Grid_and_time& gt)
  :tp {gt.time_provider()} {}
sphere_eigenFun* sphere_eigenFun::clone(){
  return new sphere_eigenFun {*this};
}
void sphere_eigenFun::interpolate(Grid::Vec_FEfun& vfef) const{
  using Vfef = Esfem::Grid::Vec_FEfun::Dune_FEfun;
  Dune::LagrangeInterpolation<Vfef>::interpolateFunction(*this, vfef); 
}
void sphere_eigenFun::evaluate(const Domain& x, Range& y) const{
  y[0] = x[0]*x[1]; // xy
  y[1] = x[1]*x[2]; // yz
  y[2] = x[0]*x[2]; // xz
  y *= exp(-6.*tp.time());
}

// ----------------------------------------------------------------------
// sphere_mcf_sol

sphere_mcf_sol::sphere_mcf_sol(const Grid::Grid_and_time& gt) 
  :tp {gt.time_provider()} {}
sphere_mcf_sol* sphere_mcf_sol::clone(){
  return new sphere_mcf_sol {*this};
}
void sphere_mcf_sol::interpolate(Grid::Vec_FEfun& rhs) const{
  using vfef = Esfem::Grid::Vec_FEfun::Dune_FEfun;
  Dune::LagrangeInterpolation<vfef>::interpolateFunction(*this, rhs);
}
void sphere_mcf_sol::evaluate(const Domain& x, Range& y) const{
  auto norm = 0.;
  for(int i = 0; i < Domain::dimension; ++i) norm += x[i]*x[i];
  norm = sqrt(norm);
  const auto rt = sqrt(1. - 2 * 2 * tp.time());
  y = x;
  y *= rt / norm;
}

// ----------------------------------------------------------------------
// sls_iData

sls_iData::sls_iData(const Grid::Grid_and_time& gt)
  :tp {gt.time_provider()},
   rA {Parameter::getValue<double>("logistic_growth.r_start", 1.)},
   rE {Parameter::getValue<double>("logistic_growth.r_end", 2.)},
   k {Parameter::getValue<double>("logistic_growth.steepness", .5)}
{}  
void sls_iData::interpolate(Grid::Vec_FEfun& rhs) const{
  using vfef = Esfem::Grid::Vec_FEfun::Dune_FEfun;
  Dune::LagrangeInterpolation<vfef>::interpolateFunction(*this, rhs);
}
void sls_iData::evaluate(const Domain& x, Range& y) const{
  const auto 
    norm = sqrt(inner_product(&x[0], &x[0]+Domain::dimension, &x[0], 0.)),
    ekt = exp(-k*tp.time()),
    lgf = rE*rA/(rE * ekt + rA * (1 - ekt));
  y = x;
  y *= lgf / norm;
}

sls_v_iData::sls_v_iData(const Grid::Grid_and_time& gt)
  :tp {gt.time_provider()},
   rA {Parameter::getValue<double>("logistic_growth.r_start", 1.)},
   rE {Parameter::getValue<double>("logistic_growth.r_end", 2.)},
   k {Parameter::getValue<double>("logistic_growth.steepness", .5)}
{}  
void sls_v_iData::interpolate(Grid::Vec_FEfun& rhs) const{
  using vfef = Esfem::Grid::Vec_FEfun::Dune_FEfun;
  Dune::LagrangeInterpolation<vfef>::interpolateFunction(*this, rhs);
}
void sls_v_iData::evaluate(const Domain& x, Range& y) const{
  const auto 
    ekt = exp(-k*tp.time()),
    rt = rE*rA/(rE * ekt + rA * (1. - ekt));
  y = x;
  y *= k * (1. - rt/ rE);
}

sd_iData::sd_iData(const Grid::Grid_and_time& gt) :tp {gt.time_provider()} {}
void sd_iData::interpolate(Grid::Vec_FEfun& rhs) const{
  using vfef = Esfem::Grid::Vec_FEfun::Dune_FEfun;
  Dune::LagrangeInterpolation<vfef>::interpolateFunction(*this, rhs);
}
void sd_iData::evaluate(const Domain& x, Range& y) const{
  const auto
    norm = sqrt(inner_product(&x[0], &x[0]+Domain::dimension, &x[0], 0.)),
    fac = exp(tp.time()) / norm;
  y = x;
  y *= fac;
}
// ----------------------------------------------------------------------
// Analytic_velocity

Analytic_velocity::Analytic_velocity(const Esfem::Grid::Grid_and_time& gt)
  :tp {gt.time_provider()}
{}

void Analytic_velocity::evaluate(const Domain& d, Range& r) const{
  const double t = tp.time();
  eoc_velocity(t, d, r);
}

// ----------------------------------------------------------------------
// helper functions

double Esfem::Impl::
hom_value(const Esfem::Io::Parameter& p, const Esfem::Growth type){
  double rv = 0.;
  switch(type){
  case Esfem::Growth::promoting:
    rv = p.u_hom_value();
    break;
  case Esfem::Growth::inhibiting:
    rv = p.w_hom_value();
    break;
  default:
    throw InitData_error{Assert::compose(__FILE__, __LINE__, "Bad Growth type")};
    break;
  };
  return rv;
}

double Esfem::Impl::
pertubation(const Esfem::Io::Parameter& p, const Esfem::Growth type){
  double rv = 0.;
  switch(type){
  case Esfem::Growth::promoting:
    rv = p.u_pertubation();
    break;
  case Esfem::Growth::inhibiting:
    rv = p.w_pertubation();
    break;
  default:
    throw InitData_error{Assert::compose(__FILE__, __LINE__, "Bad Growth type")};
    break;    
  };
  return rv;
}

std::string Esfem::Impl::print_configuration(const Esfem::Io::Parameter& p,
					     const Esfem::Growth type){
  std::ostringstream oss;
  oss << "Using Random_initial_data():\n"
    "Random distribution: uniform_real_distribution<> with\n"
      << "hom_value: " << hom_value(p, type) << '\n'
      << "pertubation: " << pertubation(p, type);
  return oss.str();
}

std::string Esfem::Impl::dof_filename(const Io::Parameter& p, const Growth type){
  std::string rv {};
  switch(type){
  case Growth::promoting:
    rv = p.u_init_dof();
    break;
  case Growth::inhibiting:
    rv = p.w_init_dof();
    break;
  default:
    throw InitData_error {Assert::compose(__FILE__, __LINE__, "Bad Growth type")};
    break;
  };
  return rv;
}

// ----------------------------------------------------------------------
// Implementation Init_data::Data

Esfem::SecOrd_op::Init_data::Data::Data(const Grid::Grid_and_time& gt,
					const Growth type)
  :eid_ptr {std::make_unique<Explicit_initial_data>(gt,type)}
{}

Esfem::SecOrd_op::Init_data::Data::Data(const Io::Parameter& p, const Growth type)
  :dof_io_filename {Impl::dof_filename(p, type)},
   rid_ptr {std::make_unique<Random_initial_data>(p, type)}
{}

// ----------------------------------------------------------------------
// Exact_velocity::Data

Esfem::SecOrd_op::Exact_velocity::Data::Data(const Grid::Grid_and_time& gt)
  :v_fun {gt}
{}
