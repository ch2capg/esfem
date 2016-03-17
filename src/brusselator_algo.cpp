/*! \file brusselator_algo.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Implementation details for brusselator_algo.h

     Created by Christian Power on 16.03.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
*/

#include <config.h>
#include "brusselator_algo.h"
#include "brusselator_algo_impl.h"

using namespace Esfem;
using Esfem::Brusselator_scheme;
using Esfem::SecOrd_op::Identity;
using Scal_FEfun_set = Esfem::FEfun_set<Esfem::Grid::Scal_FEfun>;
using Vec_FEfun_set = Esfem::FEfun_set<Esfem::Grid::Vec_FEfun>;

void Esfem::brusselator_algo(int argc, char** argv){
  Dune::Fem::MPIManager::initialize(argc, argv);

  constexpr auto parameter_file =
    "/home/power/cpp/DISS_surfaces/data/tumor_parameter.txt";  

  Brusselator_scheme fem {argc, argv, parameter_file};

  for(long it = 0; it < fem.prePattern_timeSteps(); ++it){
    fem.pre_pattern_action();
    fem.next_timeStep();
  }
  // intermediate_action(fem); // Do we need this?
  // fem.next_timeStep(); // Do we need this?
  for(long it = 0; it < fem.pattern_timeSteps(); ++it){
    fem.pattern_action();
    fem.next_timeStep();
  }
  fem.final_action();
}

// ----------------------------------------------------------------------
// Implementation of Brusselator_scheme::Data
 
struct Brusselator_scheme::Data{
  Identity identity {};
  Esfem::Io::Parameter data;
  const Esfem::Io::Dgf::Handler dgf_handler;
  Err_stream estream;
  Esfem::Grid::Grid_and_time fix_grid;
  Scal_FEfun_set u;
  Scal_FEfun_set w;
  Vec_FEfun_set surface;
  Data(int argc, char** argv, const std::string& parameter_fname);
};

Brusselator_scheme::Data::Data(int argc, char** argv, const std::string& parameter_fname)
  : 
  data {argc, argv, parameter_fname},
  dgf_handler {data.grid()},
  estream {data},
  fix_grid {data},
  u {"u", fix_grid},
  w {"w", fix_grid},
  surface {"surface", fix_grid}
{}

// ----------------------------------------------------------------------
// Brusselator_scheme implemenation

Brusselator_scheme::~Brusselator_scheme() = default;

Brusselator_scheme::Brusselator_scheme(int argc, char** argv, const std::string& parameter_fname)
try : d_ptr {std::make_unique<Data>(argc, argv, parameter_fname)}
{
  pre_loop_action();
  next_timeStep();
}
catch(const std::exception&){
  throw_with_nested(std::runtime_error {"Error in constructor of Brusselator_scheme."});
 }
 catch(...){
   throw std::runtime_error {"Unknown error in constructor of Brusselator_scheme."};
 }

void Brusselator_scheme::next_timeStep(){
  d_ptr -> fix_grid.next_timeStep(d_ptr -> data.global_timeStep());
}
long Brusselator_scheme::prePattern_timeSteps() const{
  return d_ptr -> data.prePattern_timeSteps();
}
long Brusselator_scheme::pattern_timeSteps() const{
  return d_ptr -> data.pattern_timeSteps();
}

// ------------------------------------------------------------
// Brusselator_scheme loop action

void Brusselator_scheme::pre_loop_action(){
  const Init_data initData_loc {d_ptr -> data};
  const Err_cal errCal_loc {d_ptr -> fix_grid, d_ptr -> u, d_ptr -> w};
  Io::Paraview paraview_loc {d_ptr -> data, d_ptr -> fix_grid, 
      d_ptr -> u.fun, d_ptr -> w.fun};
  Solver solver {d_ptr -> data, d_ptr -> fix_grid, d_ptr -> u, d_ptr -> w};
  
  first_interpolate(d_ptr -> identity, initData_loc,
		    d_ptr -> u, d_ptr -> w, d_ptr -> surface);

  generate_header_line(d_ptr -> estream);
  write_error_line(d_ptr -> estream, d_ptr -> fix_grid, errCal_loc);
  
  paraview_loc.write();

  massMatrix_rhsLes(solver, d_ptr -> u, d_ptr -> w);

  d_ptr -> surface.write(d_ptr -> dgf_handler, "./");
}
void Brusselator_scheme::pre_pattern_action(){
  const Grid::Grid_and_time grid_loc 
  {d_ptr -> data, 
      compose_dgfName(d_ptr -> surface.fun.name()), 
      d_ptr -> fix_grid.time_provider().time()};
  Scal_FEfun_set u_loc {d_ptr -> u, grid_loc};
  Scal_FEfun_set w_loc {d_ptr -> w, grid_loc};
  const Err_cal errCal_loc {grid_loc, u_loc, w_loc};
  Solver solver_loc {d_ptr -> data, grid_loc, u_loc, w_loc};
  Io::Paraview paraview_loc 
  {d_ptr -> data, grid_loc, u_loc.fun, w_loc.fun};

  massMatrixConstOne_rhsLes(solver_loc, u_loc, w_loc);
  // assemble_and_addScaled_rhsLes(rhs, u, w);	// testing rhs
  solve_pde(solver_loc, u_loc, w_loc);
  massMatrix_rhsLes(solver_loc, u_loc, w_loc);
  // update_exactSolution(init_data, u, w);	// testing rhs

  write_error_line(d_ptr -> estream, d_ptr -> fix_grid, errCal_loc);
  paraview_loc.write();  

  d_ptr -> u = u_loc;
  d_ptr -> w = w_loc;
  // Perhaps do something with surface?
}
void Brusselator_scheme::pattern_action(){
}
void Brusselator_scheme::final_action(){
  d_ptr -> u >> d_ptr -> dgf_handler;
}
