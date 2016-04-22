/*! \file secOrd_op_initData_impl.cpp
    \brief Implementation of secOrd_op_initData_impl.h

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) Februar 2016

     Idea
     --------------------------------------------------

     Implementing `Explicit_initial_data` and `Random_initial_data`.

     \author Christian Power 
     \date 22. April 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include "secOrd_op_initData_impl.h"
#include "io_parameter.h"
#include "esfem_error.h"

//! Implementing this
using Esfem::Impl::Explicit_initial_data;
//! Implementing this
using Esfem::Impl::Random_initial_data;


// ----------------------------------------------------------------------
// Implementation Explicit_initial_data

Explicit_initial_data::
Explicit_initial_data(const Esfem::Grid::Grid_and_time& gt,
		      const Esfem::Growth type)
  :tp {gt.time_provider()}
{
  switch(type){
  case Growth::promoting:
    fun_impl = [tp_ptr = &tp](const Domain& d, Range& r){
      const double x = d[0];
      const double y = d[1];
      const double t = tp_ptr -> time();
      r = x * y * std::exp(-6. * t);
    };
    break;
  case Growth::inhibiting:
    fun_impl = [tp_ptr = &tp](const Domain& d, Range& r){
      const double y = d[1];
      const double z = d[2];
      const double t = tp_ptr -> time();
      r = y * z * std::exp(-6. * t);
    };
    break;
  default:
    throw InitData_error {Assert::compose(__FILE__, __LINE__, "Bad Growth type")};
    break;
  };
  std::cout << "Using Explicit_initial_data():\n"
	    << "u = std::exp(-6.*t)*x*y\n"
	    << "w = std::exp(-6.*t)*y*z" << std::endl;
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
