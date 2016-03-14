/*! \file secOrd_op_initData_impl.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 21. Februar 2016

     Implementation details for secOrd_op_initData_impl.h
     Created by Christian Power on 21.02.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include "secOrd_op_initData_impl.h"
#include "io_parameter.h"

// ----------------------------------------------------------------------
// Implementation Explicit_initial_data

Esfem::Explicit_initial_data::
Explicit_initial_data(const Esfem::Grid::Grid_and_time& gt)
  : tp {gt.time_provider()}
{
  std::cout << "Using Explicit_initial_data():\n"
    "std::exp(-6.*t)*x*y" << std::endl;
}
void Esfem::Explicit_initial_data::
evaluate(const Domain& d, Range& r) const{
  static_assert(Domain::dimension == 3, "Bad domain dimension.");
  static_assert(Range::dimension == 1, "Bad range dimension.");
  const double x = d[0];
  const double y = d[1];
  const double t = tp.time();
  r = std::exp(-6.*t)*x*y;
}

// ----------------------------------------------------------------------
// Implementation Random_initial_data

Esfem::Random_initial_data::
Random_initial_data(const Esfem::Io::Parameter& p,
		    const Esfem::Growth type)
  : Random_initial_data {hom_value(p, type), pertubation(p, type)}
{
  std::cout << print_configuration(p, type) << std::endl;
}
Esfem::Random_initial_data::
Random_initial_data(const double hom_value,
		    const double pertubation)
  : random_fun (std::bind(Random_dist {hom_value, hom_value + pertubation},
			  Random_engine {}))
{}
void Esfem::Random_initial_data::
evaluate(const Domain&, Range& q) const{
  q = random_fun(); 
}  

// ----------------------------------------------------------------------
// helper functions

double Esfem::
hom_value(const Esfem::Io::Parameter& p, const Esfem::Growth type) try{
  double rv = 0.;
  switch(type){
  case Esfem::Growth::promoting:
    rv = p.u_hom_value();
    break;
  case Esfem::Growth::inhibiting:
    rv = p.w_hom_value();
    break;
  default:
    throw std::logic_error {"Invalid Growth type for hom_value()."};
    break;
  };
  return rv;
 }
 catch(const std::exception&){
   std::throw_with_nested(std::logic_error
			  {"Error in constructor of "
			      "Random_initial_data in hom_value()."});
 }
 catch(...){
   throw std::logic_error {"Unkown error in constructor of "
       "Random_initial_data in hom_value()"};
 }

double Esfem::
pertubation(const Esfem::Io::Parameter& p, const Esfem::Growth type) try{
  double rv = 0.;
  switch(type){
  case Esfem::Growth::promoting:
    rv = p.u_pertubation();
    break;
  case Esfem::Growth::inhibiting:
    rv = p.w_pertubation();
    break;
  default:
    throw std::logic_error {"Invalid Growth type for pertubation()."};
    break;    
  };
  return rv;
 }
 catch(const std::exception&){
   std::throw_with_nested(std::logic_error
			  {"Error in constructor of "
			      "Random_initial_data in pertubation()."});
 }
 catch(...){
   throw std::logic_error {"Unkown error in constructor of "
       "Random_initial_data in pertubation()"};
 }
std::string Esfem::print_configuration(const Esfem::Io::Parameter& p,
				       const Esfem::Growth type) try{
  std::ostringstream oss;
  oss << "Using Random_initial_data():\n"
    "Random distribution: uniform_real_distribution<> with\n"
      << "hom_value: " << hom_value(p, type) << '\n'
      << "pertubation: " << pertubation(p, type);
  return oss.str();
 }
 catch(const std::exception&){
   std::throw_with_nested(std::logic_error {"Error in print_configuration()."});
 }
std::string Esfem::dof_filename(const Io::Parameter& p, const Growth type) try{
  std::string rv {};
  switch(type){
  case Growth::promoting:
    rv = p.u_init_dof();
    break;
  case Growth::inhibiting:
    rv = p.w_init_dof();
    break;
  default:
    throw std::logic_error {"Unkown Growth type"};
    break;
  };
  return rv;
 }
 catch(const std::exception&){
   std::throw_with_nested(std::logic_error {"Error in dof_filename()."});
 }
/*! Log:
 */
