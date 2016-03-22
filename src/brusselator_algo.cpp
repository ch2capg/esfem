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
#include "secOrd_op_solutionDriven.h"

using namespace Esfem;
using Esfem::Brusselator_scheme;
using Scal_FEfun_set = Brusselator_scheme::Scal_FEfun_set;
using Vec_FEfun_set = Brusselator_scheme::Vec_FEfun_set;

void Esfem::brusselator_algo(int argc, char** argv){
  Dune::Fem::MPIManager::initialize(argc, argv);

  constexpr auto parameter_file =
    "/home/power/cpp/DISS_surfaces/data/tumor_parameter.txt";  

  Brusselator_scheme fem {argc, argv, parameter_file};

  fem.prePattern_loop();
  fem.intermediate_action(); 
  fem.pattern_loop();
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
  Solver solver_prePattern;
  Io::Paraview paraview_prePattern;
  Err_cal errCal_prePattern;
  Data(int argc, char** argv, const std::string& parameter_fname);
};

Brusselator_scheme::Data::Data(int argc, char** argv, const std::string& parameter_fname)
  : data {argc, argv, parameter_fname},
  dgf_handler {data.grid()},
  estream {data},
  fix_grid {data},
  u {"u", fix_grid},
  w {"w", fix_grid},
  surface {"surface", fix_grid},
  solver_prePattern {data, fix_grid, u, w},
  paraview_prePattern {data, fix_grid, u.fun, w.fun},
  errCal_prePattern {fix_grid, u, w}
{}

// ----------------------------------------------------------------------
// Brusselator_scheme implemenation

Brusselator_scheme::Brusselator_scheme(int argc, char** argv, const std::string& parameter_fname)
try : data {argc, argv, parameter_fname},
  io {data},
  fix_grid {data},
  fef {fix_grid}
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
  fix_grid.next_timeStep(data.global_timeStep());
}
long Brusselator_scheme::prePattern_timeSteps() const{
  return data.prePattern_timeSteps();
}
long Brusselator_scheme::pattern_timeSteps() const{
  return data.pattern_timeSteps();
}

// ------------------------------------------------------------
// Brusselator_scheme loop action

void Brusselator_scheme::pre_loop_action(){
  PreLoop_helper helper {*this};
  helper.first_interpolate();
  helper.headLine_in_errFile();
  helper.plot_errors_in_errFile();
  helper.plot_paraview();
  helper.prepare_rhs();
  
  // const Init_data initData_loc {d_ptr -> data};
  // const Err_cal errCal_loc {d_ptr -> fix_grid, d_ptr -> u, d_ptr -> w};
  // Io::Paraview paraview_loc {d_ptr -> data, d_ptr -> fix_grid, 
  //     d_ptr -> u.fun, d_ptr -> w.fun};
  // Solver solver {d_ptr -> data, d_ptr -> fix_grid, d_ptr -> u, d_ptr -> w};
  // first_interpolate(d_ptr -> identity, initData_loc,
  // 		    d_ptr -> u, d_ptr -> w, d_ptr -> surface);
  // generate_header_line(d_ptr -> estream);
  // write_error_line(d_ptr -> estream, d_ptr -> fix_grid, errCal_loc);
  // paraview_loc.write();
  // massMatrix_rhsLes(solver, d_ptr -> u, d_ptr -> w);
}
void Brusselator_scheme::prePattern_loop(){
  PrePattern_helper helper {*this};
  for(long it = 0; it < prePattern_timeSteps(); ++it){
    helper.finalize_rhs();
    helper.solve_pde();
    helper.prepare_rhs();
    helper.plot_errors_in_errFile();
    helper.plot_paraview();
    next_timeStep();
  }
  // massMatrixConstOne_rhsLes(d_ptr -> solver_prePattern,
  // 			    d_ptr -> u, d_ptr -> w);
  // solve_pde(d_ptr -> solver_prePattern,
  // 	    d_ptr -> u, d_ptr -> w);
  // massMatrix_rhsLes(d_ptr -> solver_prePattern,
  // 		    d_ptr -> u, d_ptr -> w);
  // write_error_line(d_ptr -> estream, d_ptr -> fix_grid, errCal_loc);
  // paraview_loc.write();  
}
void Brusselator_scheme::intermediate_action(){
  fef.u.write(io.dgf_handler, "./intermediate_");
  fef.w.write(io.dgf_handler, "./intermediate_");
  fef.surface.write(io.dgf_handler, "./");
  intermediate_surface_rhs();
}
void Brusselator_scheme::pattern_loop(){
  for(long it = 0; it < pattern_timeSteps(); ++it){
    solve_surfacePDE();
    Pattern_helper helper {*this};
    helper.finalize_scalarPDE_rhs();
    helper.solve_scalarPDE();
    helper.prepare_all_rhs();
    helper.plot_errors_in_errFile();
    helper.plot_paraview();
    next_timeStep();
  }
  Helper_surface hs {d_ptr};
  hs.solve_for_surface();
  Helper_uw huw {d_ptr};
  huw.solve_pde();
  huw.write_error_line();
  huw.paraview_plot();
}
void Brusselator_scheme::final_action(){
  fef.u.write(io.dgf_handler, "./final_");
  fef.w.write(io.dgf_handler, "./final_");
  fef.surface.write(io.dgf_handler, "./final_"); 
}


void Brusselator_scheme::pre_pattern_action(){  
  // const Grid::Grid_and_time grid_loc 
  // {d_ptr -> data, 
  //     compose_dgfName(d_ptr -> surface.fun.name()), 
  //     d_ptr -> fix_grid.time_provider().time()};
  // Scal_FEfun_set u_loc {d_ptr -> u, grid_loc};
  // Scal_FEfun_set w_loc {d_ptr -> w, grid_loc};
  // const Err_cal errCal_loc {grid_loc, u_loc, w_loc};
  // Solver solver_loc {d_ptr -> data, grid_loc, u_loc, w_loc};
  // Io::Paraview paraview_loc 
  // {d_ptr -> data, grid_loc, u_loc.fun, w_loc.fun};

  // massMatrixConstOne_rhsLes(solver_loc, u_loc, w_loc);
  // assemble_and_addScaled_rhsLes(rhs, u, w);	// testing rhs
  // solve_pde(solver_loc, u_loc, w_loc);
  // massMatrix_rhsLes(solver_loc, u_loc, w_loc);
  // update_exactSolution(init_data, u, w);	// testing rhs

  massMatrixConstOne_rhsLes(d_ptr -> solver_prePattern,
			    d_ptr -> u, d_ptr -> w);
  solve_pde(d_ptr -> solver_prePattern,
	    d_ptr -> u, d_ptr -> w);
  massMatrix_rhsLes(d_ptr -> solver_prePattern,
		    d_ptr -> u, d_ptr -> w);
  
  write_error_line(d_ptr -> estream, d_ptr -> fix_grid, errCal_loc);
  paraview_loc.write();  

  d_ptr -> u = u_loc;
  d_ptr -> w = w_loc;
  // Perhaps do something with surface?
}
void Brusselator_scheme::pattern_action(){
  Helper_surface hs {d_ptr};
  hs.solve_for_surface();
  Helper_uw huw {d_ptr};
  huw.solve_pde();
  huw.write_error_line();
  huw.paraview_plot();
  
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
}

// ----------------------------------------------------------------------
// Implementation of Brusselator_scheme's helper struct's constructors

struct Brusselator_scheme::Fef::Fef(Grid::Grid_and_time& gt)
  : u {"u", gt},
  w {"w", gt},
  surface {"surface", gt}
{}

struct Brusselator_scheme::Io::Io(const Io::Parameter& p)
  : dgf_handler {p.grid()},
  estream {p}
{}
