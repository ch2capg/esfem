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
#include "esfem_error.h"

//! For convenience 
using namespace std;
//! Implementing this
using Esfem::SecOrd_op::Init_data;

Esfem::SecOrd_op::Init_data::
Init_data(const Grid::Grid_and_time& gt, const Growth type) 
  :d_ptr {make_unique<Data>(gt, type)}
{}

Esfem::SecOrd_op::Init_data::Init_data(const Io::Parameter& p,
				       const Growth type) 
  :d_ptr {make_unique<Data>(p, type)}
{}

 Esfem::SecOrd_op::Init_data::~Init_data() = default;

void Esfem::SecOrd_op::Init_data::interpolate(Grid::Scal_FEfun& fef) const{
  using std::begin;
  using std::end;
  using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
  
  const auto& eid_ptr = d_ptr -> eid_ptr;
  const auto& rid_ptr = d_ptr -> rid_ptr;
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
    throw InitData_error {Assert::compose(__FILE__, __LINE__, "Null pointer")};
}
