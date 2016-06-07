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
#include <dune/common/exceptions.hh>
#include "grid.h"
#include "grid_GridAndTime_impl.h"
#include "io_parameter.h"
#include "esfem_error.h"

using namespace std;

using Esfem::Grid_error;
using Deformation = Esfem::Grid::Deformation;
using Host_grid = Esfem::Grid::Grid_and_time::Host_grid;
using HGrid_ptr = Dune::GridPtr<Host_grid>;
using Grid = Esfem::Grid::Grid_and_time::Grid;
using Grid_part = Esfem::Grid::Grid_and_time::Grid_part;
using Time_provider = Dune::Fem::GridTimeProvider<Grid>;
using FE_space = Esfem::Grid::Grid_and_time::FE_space;
using Vec_FE_space = Esfem::Grid::Grid_and_time::Vec_FE_space;

Host_grid& loadBalance_and_deref(HGrid_ptr&);

// ----------------------------------------------------------------------
// Implemenation Esfem::Grid::Grid_and_time::Data

struct Esfem::Grid::Grid_and_time::Data{
  // Impl::Evolving_grid eg;
  Deformation d;
  HGrid_ptr hg_ptr;
  Grid g;
  Grid_part gp;
  Time_provider tp;
  FE_space fes;
  Vec_FE_space vfes;
  Data(const Io::Parameter&);
  Data(const Io::Parameter&, const std::string& dgf_file,
       const double t0);
};

Esfem::Grid::Grid_and_time::Data::Data(const Io::Parameter& p)
  :// d {new Impl::Evolving_grid {eg}},
   d {p.grid()}, hg_ptr {p.grid()},
   g {loadBalance_and_deref(hg_ptr), d},
   gp {g}, tp {p.start_time(), g},
   fes {gp}, vfes {gp} {
  d.set_timeProvider(tp);
  tp.init(p.global_timeStep());
}

Esfem::Grid::Grid_and_time::Data::
Data(const Io::Parameter& p, const std::string& dgf_file,
     const double t0)
  :d {}, hg_ptr {dgf_file}, g {loadBalance_and_deref(hg_ptr), d},
   gp {g}, tp {t0, g}, fes {gp}, vfes {gp} {
  d.set_timeProvider(tp);
  tp.init(p.global_timeStep());
}

// ----------------------------------------------------------------------
// Implemenation Esfem::Grid::Grid_and_time


Esfem::Grid::Grid_and_time::Grid_and_time(const Io::Parameter& p)
try :d_ptr {std::make_unique<Data>(p)}
{}
catch(const Dune::Exception& e){
  ostringstream oss;
  oss << "Constructor\n"
      << "Dune error: " << e << std::endl;
  throw Grid_error {oss.str()};
 }
 catch(...){
   throw_with_nested(Grid_error {"Constructor"});
 }

Esfem::Grid::Grid_and_time::
Grid_and_time(const Io::Parameter& p, const std::string& dgf_file, const double t0)
try : d_ptr {std::make_unique<Data>(p, dgf_file, t0)}
{}
catch(const std::exception&){
  throw_with_nested(Grid_error {"Constructor"});
 }
 catch(const Dune::Exception& e){
   ostringstream oss;
   oss << "Constructor\n"
       << "Dune error: " << e << std::endl;
   throw Grid_error {oss.str()};
 }
 catch(...){
   throw Grid_error {"Constructor, unknown error"};
 }

Esfem::Grid::Grid_and_time::~Grid_and_time() = default;

void Esfem::Grid::Grid_and_time::next_timeStep(const double dT){
  d_ptr -> tp.next(dT);
}
// void Esfem::Grid::Grid_and_time::new_nodes(const Vec_FEfun& vfef){
//   auto eg = d_ptr -> eg;
//   eg = vfef;
// }
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
Vec_FE_space& Esfem::Grid::Grid_and_time::vec_fe_space() const{
  return d_ptr -> vfes;
}

// ----------------------------------------------------------------------
// Internal implementation

Host_grid& loadBalance_and_deref(HGrid_ptr& hg_ptr){
  hg_ptr -> loadBalance();
  return *hg_ptr;
}

/*! Log:
 */
