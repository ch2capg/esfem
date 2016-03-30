/*! \file brusselator_algo.cpp
    \author Christian Power
    \date 30. March 2016

    \brief Implementation for `brusselator_algo.h`

     Revision history
     --------------------------------------------------

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Idea
     --------------------------------------------------
     
     Implementation of Esfem::brusselator_algo() and
     the class `Esfem::Brusselator_scheme`

         Created by Christian Power on 30.03.2016
         Copyright (c) 2016 Christian Power. All rights reserved.
*/

#include <config.h>
#include "brusselator_algo.h"
#include "brusselator_algo_impl.h"
#include "secOrd_op_solutionDriven.h"
#include "esfem_error.h"

#ifndef PFILE
#error Give a complete path for parameter file in macro variable PFILE.
#endif

using Bruss_error = Esfem::BrusselatorScheme_error;
/*!< \brief Shorter name for convenience */
using Esfem::Brusselator_scheme;
/*!< \brief We essentially implement this class in this source file. */
using Scal_FEfun_set = Esfem::Grid::Scal_FEfun_set;
/*!< \brief Four functions of type \f$ f\colon \R^3 \to \R \f$ */
using Vec_FEfun_set = Esfem::Grid::Vec_FEfun_set;
/*!< \brief Four functions of type \f$ f\colon \R^3 \to \R^3 \f$ */
using Vector_solver = Esfem::SecOrd_op::Solution_driven;
/*!< \brief Solver for the `X` respectively `surface` variable */

void Esfem::brusselator_algo(int argc, char** argv){
  Dune::Fem::MPIManager::initialize(argc, argv);
  constexpr auto parameter_file = PFILE;
  Brusselator_scheme fem {argc, argv, parameter_file};
  fem.prePattern_loop();
  fem.intermediate_action(); 
  fem.pattern_loop();
  fem.final_action();
}

// ----------------------------------------------------------------------
// Brusselator_scheme implementation

// ------------------------------------------------------------
// Brusselator_scheme public 

Brusselator_scheme::
Brusselator_scheme(int argc, char** argv,
		   const std::string& parameter_fname)
try : data {argc, argv, parameter_fname},
  io {data},
  fix_grid {data},
  fef {fix_grid}
{
  // TODO: test if folder for tmp file exists
  pre_loop_action();
}
catch(const std::exception&){
  throw_with_nested(Bruss_error {"Constructor."}); 
 }
 catch(...){
   throw Bruss_error {"Unknown error in constructor."};
 }

// --------------------------------------------------
// Brusselator_scheme loop action

void Brusselator_scheme::prePattern_loop(){
  std::cerr << "prePattern_timeSteps(): " << prePattern_timeSteps() << std::endl;
    
  PrePattern_helper helper {*this};
  for(long it = 0; it < prePattern_timeSteps(); ++it, next_timeStep()){
    // helper.finalize_rhs();
    helper.rhs();
    helper.solve_pde();
    // helper.prepare_rhs();
    // helper.plot_errors_in_errFile();
    helper.plot_paraview();
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
  switch(prePattern_timeSteps()){
  case 0: // heat.starttime == heat.pattern.endtime
    std::cerr << "heat.starttime == heat.pattern.endtime" << std::endl;
    fef.u.read(io.dgf_handler, "./intermediate_");
    fef.w.read(io.dgf_handler, "./intermediate_");
    break;
  default:
    std::cerr << "heat.starttime != heat.pattern.endtime" << std::endl;
    fef.u.write(io.dgf_handler, "./intermediate_");
    fef.w.write(io.dgf_handler, "./intermediate_");
    fef.surface.write(io.dgf_handler, "./");
    break;
  };
}
void Brusselator_scheme::pattern_loop(){
  std::cerr << "pattern_timeSteps(): " << pattern_timeSteps() << std::endl;
  
  for(long it = 0; it < pattern_timeSteps(); ++it, next_timeStep()){
    // solve_surfacePDE();
    rhs_and_solve_SPDE();
    Pattern_helper helper {*this};
    helper.finalize_scalarPDE_rhs();
    helper.solve_scalarPDE();
    // helper.prepare_all_rhs();
    // helper.plot_errors_in_errFile();
    helper.plot_paraview();
  }
  // Helper_surface hs {d_ptr};
  // hs.solve_for_surface();
  // Helper_uw huw {d_ptr};
  // huw.solve_pde();
  // huw.write_error_line();
  // huw.paraview_plot();
}
void Brusselator_scheme::final_action(){
  fef.u.write(io.dgf_handler, "./final_");
  fef.w.write(io.dgf_handler, "./final_");
  fef.surface.write(io.dgf_handler, "./final_"); 
}

// ------------------------------------------------------------
// Brusselator_scheme private

void Brusselator_scheme::pre_loop_action(){
  PreLoop_helper helper {*this};
  helper.first_interpolate();
  helper.headLine_in_errFile();
  // helper.plot_errors_in_errFile();
  helper.plot_paraview();
  // helper.prepare_rhs();
  next_timeStep();
  // ----------------------------------------------------------------------
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
void Brusselator_scheme::rhs_and_solve_SPDE(){
  std::cerr << "rhs_and_solve_SPDE()." << std::endl;
  RhsAndSolve_helper helper {*this};
  std::cerr << "RhsAndSolve_helper()." << std::endl;
  helper.scalar_massMatrix();
  std::cerr << "helper.scalar_massMatrix()." << std::endl;
  helper.solve_surface_and_save();
  std::cerr << "helper.solve_surface_and_save()." << std::endl;
  // // constructor 
  // const Grid::Grid_and_time grid
  // {data, compose_dgfName(fef.surface.fun.name()), 
  //    fix_grid.time_provider().time()};
  // Scal_FEfun_set u {fef.u, grid};
  // Scal_FEfun_set w {fef.w, grid};
  // Vec_FEfun_set X {fef.surface, grid};
  // Scalar_solver ss  {data, grid, u, w};
  // Vector_solver vs {data, grid, u.fun};
  // // scalar part
  // ss.u.mass_matrix(u.rhs_les);
  // ss.w.mass_matrix(w.rhs_les);
  // fef.u = u;
  // fef.w = w;
  // // vector part
  // vs.rhs(X.fun, X.rhs_les);
  // vs.solve(X.rhs_les, X.fun);
  // fef.surface.fun = X.fun;
  // fef.surface.write(io.dgf_handler, "./");   
}
// // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// // delete this
// void Brusselator_scheme::intermediate_surface_rhs(){
//   // const Grid::Grid_and_time grid_loc 
//   //{data, compose_dgfName(fef.surface.fun.name()), 
//   //    fix_grid.time_provider().time()};
//   // Scal_FEfun_set u_loc {fef.u, grid_loc};
//   // Vec_FEfun_set surface_loc {fef.surface, grid_loc};
//   const auto& u = fef.u;
//   const auto& w = fef.w;
//   Scalar_solver solver {data, fix_grid, u, w};
//   solver.u.mass_matrix(u.rhs_les);
//   solver.w.mass_matrix(w.rhs_les);
// }
// //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// // delete this
// void Brusselator_scheme::solve_surfacePDE(){
//   const Grid::Grid_and_time grid
//   {data, compose_dgfName(fef.surface.fun.name()), 
//      fix_grid.time_provider().time()};
//   const Scal_FEfun_set u {fef.u, grid};
//   Vec_FEfun_set X {fef.surface, grid};
//   SecOrd_op::Solution_driven solver {data, grid, u};
//   solver.rhs(X.fun, X.rhs);
//   solver.solve(X.rhs, X.fun);
//   fef.surface.fun = X.fun;
//   fef.surface.write(io.dgf_handler, "./"); 
// }
// void Brusselator_scheme::pre_pattern_action(){  
//   // const Grid::Grid_and_time grid_loc 
//   // {d_ptr -> data, 
//   //     compose_dgfName(d_ptr -> surface.fun.name()), 
//   //     d_ptr -> fix_grid.time_provider().time()};
//   // Scal_FEfun_set u_loc {d_ptr -> u, grid_loc};
//   // Scal_FEfun_set w_loc {d_ptr -> w, grid_loc};
//   // const Err_cal errCal_loc {grid_loc, u_loc, w_loc};
//   // Solver solver_loc {d_ptr -> data, grid_loc, u_loc, w_loc};
//   // Io::Paraview paraview_loc 
//   // {d_ptr -> data, grid_loc, u_loc.fun, w_loc.fun};

//   // massMatrixConstOne_rhsLes(solver_loc, u_loc, w_loc);
//   // assemble_and_addScaled_rhsLes(rhs, u, w);	// testing rhs
//   // solve_pde(solver_loc, u_loc, w_loc);
//   // massMatrix_rhsLes(solver_loc, u_loc, w_loc);
//   // update_exactSolution(init_data, u, w);	// testing rhs

//   massMatrixConstOne_rhsLes(d_ptr -> solver_prePattern,
// 			    d_ptr -> u, d_ptr -> w);
//   solve_pde(d_ptr -> solver_prePattern,
// 	    d_ptr -> u, d_ptr -> w);
//   massMatrix_rhsLes(d_ptr -> solver_prePattern,
// 		    d_ptr -> u, d_ptr -> w);

//   write_error_line(d_ptr -> estream, d_ptr -> fix_grid, errCal_loc);
//   paraview_loc.write();  

//   d_ptr -> u = u_loc;
//   d_ptr -> w = w_loc;
//   // Perhaps do something with surface?
// }

// !!!!!!!!!!!!!!!!!!!!
// TODO: delete this
// void Brusselator_scheme::pattern_action(){
//   Helper_surface hs {d_ptr};
//   hs.solve_for_surface();
//   Helper_uw huw {d_ptr};
//   huw.solve_pde();
//   huw.write_error_line();
//   huw.paraview_plot();
//   
//   const Grid::Grid_and_time grid_loc 
//   {d_ptr -> data, 
//       compose_dgfName(d_ptr -> surface.fun.name()), 
//       d_ptr -> fix_grid.time_provider().time()};
//   Scal_FEfun_set u_loc {d_ptr -> u, grid_loc};
//   Scal_FEfun_set w_loc {d_ptr -> w, grid_loc};
//   const Err_cal errCal_loc {grid_loc, u_loc, w_loc};
//   Solver solver_loc {d_ptr -> data, grid_loc, u_loc, w_loc};
//   Io::Paraview paraview_loc 
//   {d_ptr -> data, grid_loc, u_loc.fun, w_loc.fun};
// }

// --------------------------------------------------
// Flow control

long Brusselator_scheme::prePattern_timeSteps() const{
  return data.prePattern_timeSteps();
}
long Brusselator_scheme::pattern_timeSteps() const{
  return data.pattern_timeSteps();
}

// ----------------------------------------------------------------------
// Implementation of structs
// Brusselator_scheme::Fef and Brusselator_scheme::Io

Brusselator_scheme::Fef::Fef(const Esfem::Grid::Grid_and_time& gt)
  : u {"u", gt},
  w {"w", gt},
  surface {"surface", gt}
{}

Brusselator_scheme::Io::Io(const Esfem::Io::Parameter& p)
  : dgf_handler {p.grid()},
  u {p},
  w {p}
{}
