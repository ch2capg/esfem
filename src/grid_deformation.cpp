/*! \file grid_deformation.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     Implementation details for grid_deformation.h
     Created by Christian Power on 27.01.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include "grid.h"
#include "grid_GridAndTime_impl.h"

#ifdef DEBUG
#include <iostream>
#endif

using Domain = Esfem::Grid::Deformation::Domain;
using Range = Esfem::Grid::Deformation::Range;
static_assert(Esfem::Grid::Deformation::Domain::dimension == 3,
	      "Bad domain dimension.");
static_assert(Esfem::Grid::Deformation::Range::dimension == 3,
	      "Bad range dimension.");

inline void identity(const Domain& x, Range& y) noexcept{  
  y[0] = x[0]; 
  y[1] = x[1]; 
  y[2] = x[2]; 
}

// ----------------------------------------------------------------------
// Implementaion of Deformation

struct Esfem::Grid::Deformation::Data{
  const Impl::Evolving_grid* eg_ptr {nullptr};
  const Dune::Fem::TimeProviderBase* tp_ptr {nullptr};
  Data() = default;
  Data(const Impl::Evolving_grid& eg) : eg_ptr {&eg} {}
};

Esfem::Grid::Deformation::Deformation()
  :  d_ptr {std::make_unique<Data>()}
{}

Esfem::Grid::Deformation::Deformation(const Impl::Evolving_grid& eg)
  : d_ptr {std::make_unique<Data>(eg)}
{}

Esfem::Grid::Deformation::~Deformation() = default;

void Esfem::Grid::Deformation::
set_timeProvider(const Dune::Fem::TimeProviderBase& tp){
  d_ptr -> tp_ptr = &tp;
}
void Esfem::Grid::Deformation::evaluate(const Domain& x, Range& y) const{
  // double t = d_ptr -> tp_ptr->time();
  identity(x,y);
  // const auto eg = *(d_ptr -> eg_ptr);
  // y = eg[x];
}

/*! Log:
 */
