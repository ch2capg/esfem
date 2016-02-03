/*! \file grid_GridAndTime.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     Implementation details for grid_GridAndTime.h
     Created by Christian Power on 27.01.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <config.h>
#include <stdexcept>
#include <dune/fem/solver/timeprovider.hh>
#include "grid.h"
#include "io_parameter.h"

#ifdef DEBUG
#include <iostream>
#endif

using Deformation = Esfem::Grid::Deformation;
using Host_grid = Esfem::Grid::Grid_and_time::Host_grid;
using HGrid_ptr = Dune::GridPtr<Host_grid>;
using Grid = Esfem::Grid::Grid_and_time::Grid;
using Grid_part = Esfem::Grid::Grid_and_time::Grid_part;
using Time_provider = Dune::Fem::GridTimeProvider<Grid>;
using FE_space = Esfem::Grid::Grid_and_time::FE_space;

Host_grid& loadBalance_and_deref(HGrid_ptr&);

// ----------------------------------------------------------------------
// Implemenation Esfem::Grid::Grid_and_time

struct Esfem::Grid::Grid_and_time::Data{
  Deformation d;
  HGrid_ptr hg_ptr;
  Grid g;
  Grid_part gp;
  Time_provider tp;
  FE_space fes;
  Data(const Io::Parameter& p);
};

Esfem::Grid::Grid_and_time::Grid_and_time(const Io::Parameter& p){
  try{
    d_ptr = new Data {p};
  }
  catch(const std::exception&){
    std::throw_with_nested(std::logic_error
			   {"Error in constructor of Grid_and_time."});
  }
  catch(...){
    throw std::logic_error{"Unknown error in constructor of Grid_and_time."};
  }
}
Esfem::Grid::Grid_and_time::~Grid_and_time(){
  delete d_ptr;
  d_ptr = nullptr;
#ifdef DEBUG
  std::cerr << "~Grid_and_time(): delete d_ptr.\n";
#endif
}
void Esfem::Grid::Grid_and_time::next_timeStep(const double dT){
  d_ptr -> tp.next(dT);
}
Dune::Fem::TimeProviderBase& Esfem::Grid::Grid_and_time::time_provider(){
  return d_ptr -> tp;
}
const Dune::Fem::TimeProviderBase& Esfem::Grid::Grid_and_time::time_provider() const{
  return d_ptr -> tp;
}
Grid& Esfem::Grid::Grid_and_time::grid() const{
  return d_ptr -> g;
}
Grid_part& Esfem::Grid::Grid_and_time::grid_part() const{
  return d_ptr -> gp;
}
FE_space& Esfem::Grid::Grid_and_time::fe_space() const{
  return d_ptr -> fes;
}

// ----------------------------------------------------------------------
// Internal implementation

Host_grid& loadBalance_and_deref(HGrid_ptr& hg_ptr){
  hg_ptr -> loadBalance();
  return *hg_ptr;
}

Esfem::Grid::Grid_and_time::Data::Data(const Io::Parameter& p)
  : d {},
    hg_ptr {p.grid()},
    g {loadBalance_and_deref(hg_ptr), d},
    gp {g},
    tp {p.start_time(), g},
    fes {gp}
{
  d.set_timeProvider(tp);
  tp.init(p.global_timeStep());
}

/*! Log:
 */
