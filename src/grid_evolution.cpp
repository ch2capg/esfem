/*! \file Grid_evolution.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 19. Januar 2016

     Implementation details for evolution.h
     Created by Christian Power on 19.01.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include "grid_evolution.h"
#ifdef DEBUG
#include <iostream>
#endif

Grid_evolution::DeformationCoordFunction::
DeformationCoordFunction(const Parameter::PDE_data& d)
{}
void Grid_evolution::DeformationCoordFunction::
set_timeProvider(Dune::Fem::TimeProviderBase& tp){
  tp_ptr = &tp;
}
void Grid_evolution::DeformationCoordFunction::
evaluate(const Domain& x, Range& y) const{  
  // double t = tp_ptr->time();
  y[0] = x[0];
  y[1] = x[1];
  y[2] = x[2];		
}
Grid_evolution::DeformationCoordFunction&
Grid_evolution::deformation(const Parameter::PDE_data& d){
  static Grid_evolution::DeformationCoordFunction deformation {};
  return deformation;
}
/*! Log:
 */
