/*! \file secOrd_op_initData_u.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Implementation details for secOrd_op_initData_u.h
     Created by Christian Power on 04.02.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <stdexcept>
#include <cmath>
#include <config.h>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include "secOrd_op_initData_u.h"
#include "grid.h"

#ifdef DEBUG
#include <iostream>
#endif

class Initial_data_u
  : public Dune::Fem::Function<Esfem::Grid::Grid_and_time::Function_space,
			       Initial_data_u>
{
public:
  using Base = Esfem::Grid::Grid_and_time::Function_space;
  using Domain = Base::DomainType;
  using Range = Base::RangeType;
  
  explicit Initial_data_u(const Dune::Fem::TimeProviderBase&);
  Initial_data_u(const Initial_data_u&) = delete;
  Initial_data_u& operator=(const Initial_data_u&) = delete;

  void evaluate(const Domain&, Range&) const;
private:
  const Dune::Fem::TimeProviderBase& tp;
};

// ----------------------------------------------------------------------
// Implementation esfem.h

struct Esfem::SecOrd_op::Init_data_u::Data{
  Initial_data_u u0;
  Data(const Grid::Grid_and_time& gt)
    : u0 {gt.time_provider()}
  {}
};

Esfem::SecOrd_op::Init_data_u::Init_data_u(const Grid::Grid_and_time& gt){
  try{
    d_ptr = new Data {gt};
  }
  catch(const std::exception&){
    std::throw_with_nested(std::logic_error{"Error in constructor of Init_data_u."});
  }
  catch(...){
    throw std::logic_error{"Unknown error in constructor of Init_data_u."};
  }
}
Esfem::SecOrd_op::Init_data_u::~Init_data_u(){
  delete d_ptr;
  d_ptr = nullptr;
#ifdef DEBUG
  std::cerr << "~Init_data_u(): delete d_ptr.\n";
#endif
}
void Esfem::SecOrd_op::Init_data_u::interpolate(Grid::Scal_FEfun& fef) const{
  using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
  Dune::LagrangeInterpolation<FE_function>::interpolateFunction(d_ptr -> u0, fef);  
}

// ----------------------------------------------------------------------
// Implementation Initial_data_u

Initial_data_u::Initial_data_u(const Dune::Fem::TimeProviderBase& tpb)
  : tp {tpb}
{}
void Initial_data_u::evaluate(const Domain& d, Range& r) const{
  static_assert(Domain::dimension == 3, "Bad domain dimension.");
  static_assert(Range::dimension == 1, "Bad range dimension.");
  const double x = d[0];
  const double y = d[1];
  const double t = tp.time();
  r = std::exp(-6.*t)*x*y;
}

/*! Log:
 */
