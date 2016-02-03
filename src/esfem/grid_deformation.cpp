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

#ifdef DEBUG
#include <iostream>
#endif

using Domain = Esfem::Grid::Deformation::Domain;
using Range = Esfem::Grid::Deformation::Range;

inline void identity(const Domain& x, Range& y){
  static_assert(Domain::dimension == 3, "Bad domain dimension.");
  static_assert(Range::dimension == 3, "Bad range dimension.");
  
  y[0] = x[0];
  y[1] = x[1];
  y[2] = x[2];		  
}

// ----------------------------------------------------------------------
// Implementaion of Deformation

struct Esfem::Grid::Deformation::Data{
  const Dune::Fem::TimeProviderBase* tp_ptr {nullptr};
};

Esfem::Grid::Deformation::Deformation()
  : d_ptr {new Data{} }
    {};
// Esfem::Grid::Deformation::Deformation(const Io::Parameter& p)
// {}
Esfem::Grid::Deformation::~Deformation(){
  delete d_ptr;
  d_ptr = nullptr;
#ifdef DEBUG
  std::cerr << "~Deformation(): delete d_ptr.\n";
#endif
}
void Esfem::Grid::Deformation::
set_timeProvider(const Dune::Fem::TimeProviderBase& tp){
  d_ptr -> tp_ptr = &tp;
}
void Esfem::Grid::Deformation::evaluate(const Domain& x, Range& y) const{
  // double t = d_ptr -> tp_ptr->time();
  identity(x,y);
}

/*! Log:
 */
