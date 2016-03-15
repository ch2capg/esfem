/*! \file secOrd_op_initData.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     Implementation details for secOrd_op_initData.h
     Created by Christian Power on 30.01.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <config.h>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include "secOrd_op_initData.h"
#include "secOrd_op_initData_impl.h"
#include "io_dof.h"

using namespace std;
using Esfem::Explicit_initial_data;
using Esfem::Random_initial_data;

struct Esfem::SecOrd_op::Init_data::Data{  
  const std::string dof_io_filename {};
  unique_ptr<Explicit_initial_data> eid_ptr;
  unique_ptr<Random_initial_data> rid_ptr;
  Data(const Grid::Grid_and_time&);
  // eid_ptr constructor
  Data(const Io::Parameter&, const Growth);
  // rid_ptr constructor
  //  ~Data();
};

Esfem::SecOrd_op::Init_data::Init_data(const Grid::Grid_and_time& gt) 
try : d_ptr = make_unique<Data>(gt);
{}
catch(const std::exception&){
  std::throw_with_nested(logic_error {"Error in constructor of Init_data."});

catch(...){
  throw logic_error {"Unknown error in constructor of Init_data."};
}

Esfem::SecOrd_op::Init_data::Init_data(const Io::Parameter& p,
				       const Growth type) try{
  d_ptr = new Data {p, type};
}
catch(const std::exception&){
  std::throw_with_nested(std::logic_error{"Error in constructor of Init_data."});
}
catch(...){
  throw logic_error{"Unknown error in constructor of Init_data."};
}

 Esfem::SecOrd_op::Init_data::~Init_data() = default;
// {
//   delete d_ptr;
//   d_ptr = nullptr;
// #ifdef DEBUG
//   std::cerr << "~Init_data(): delete d_ptr.\n";
// #endif
// }
void Esfem::SecOrd_op::Init_data::interpolate(Grid::Scal_FEfun& fef) const try{
  using std::begin;
  using std::end;
  using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
  
  const auto eid_ptr = d_ptr -> eid_ptr;
  const auto rid_ptr = d_ptr -> rid_ptr;
  const auto& ofname = d_ptr -> dof_io_filename;
  
  if(eid_ptr)
    Dune::LagrangeInterpolation<FE_function>::
      interpolateFunction(*eid_ptr, fef);
  else if(rid_ptr){
    Dune::LagrangeInterpolation<FE_function>::
      interpolateFunction(*rid_ptr, fef);
    Io::dof_to_file(begin(fef), end(fef), ofname);
  }
  else
    throw logic_error {"No valid pointer"};
 }
 catch(const std::exception&){
   std::throw_with_nested(logic_error {"Error in Init_data::interpolate()."});
 }
// ----------------------------------------------------------------------
// Implementation Init_data::Data

Esfem::SecOrd_op::Init_data::Data::Data(const Grid::Grid_and_time& gt)
  try : eid_ptr {make_unique<Explicit_initial_data>(gt)}
 {}
 catch(const std::exception&){
   std::throw_with_nested(logic_error{"Error in constructor of "
	 "Init_data::Data."});
 }
 catch(...){
   throw logic_error {"Unkown error in constructor of Init_data::Data."};
 }

Esfem::SecOrd_op::Init_data::Data::Data(const Io::Parameter& p, const Growth type)
try : dof_io_filename {dof_filename(p, type)},
  rid_ptr {make_unique<Random_initial_data>(p, type)}
{}
 catch(const std::exception&){
   std::throw_with_nested(logic_error
			  {"Error in constructor of Init_data::Data."});
 }
 catch(...){
   throw logic_error {"Unkown error in constructor of Init_data::Data."};
 }

// Esfem::SecOrd_op::Init_data::Data::~Data(){
//   delete eid_ptr;
//   eid_ptr = nullptr;
//   delete rid_ptr;
//   rid_ptr = nullptr;
// #ifdef DEBUG
//   std::cerr << "~Init_data::Data(): delete eid_ptr, delete rid_ptr.\n";
// #endif 
// }

/*! Log:
 */
