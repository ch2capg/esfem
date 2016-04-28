/*! \file brusselator_algo.cpp
    \brief Implementation for `brusselator_algo.h`

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Idea
     --------------------------------------------------
     
     Implementation of Esfem::brusselator_algo() and
     the class `Esfem::Brusselator_scheme`

     \author Christian Power 
     \date 28. April 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
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
  // fem.prePattern_loop();
  // fem.intermediate_action(); 
  fem.pattern_loop();
  fem.final_action();
}

// ----------------------------------------------------------------------
// Brusselator_scheme implementation

// ------------------------------------------------------------
// Brusselator_scheme public 

/*! \todo Test if folder for tmp file exists. */
Brusselator_scheme::
Brusselator_scheme(int argc, char** argv,
		   const std::string& parameter_fname)
try :data {argc, argv, parameter_fname},
  io {data},
  fix_grid {data},
  norm {fix_grid},  
  fef {fix_grid}, 
  // exact {fix_grid, data} // constructor for random initial data
  exact {fix_grid} // constructor for analytic initial data
{
  pre_loop_action(); // initialize member fef
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
  PrePattern_helper helper {*this};
  for(long it = 0; it < prePattern_timeSteps(); ++it, next_timeStep()){
    helper.rhs();
    helper.solve_pde();
    helper.plot_paraview();
  }
}
void Brusselator_scheme::intermediate_action(){
  const auto uw_path = fef.tmpFile_path + "intermediate_";
  switch(prePattern_timeSteps()){
  case 0: // heat.starttime == heat.pattern.endtime
    fef.u.read(io.dgf_handler, uw_path);
    fef.w.read(io.dgf_handler, uw_path);
    break;
  default:
    fef.u.write(io.dgf_handler, uw_path);
    fef.w.write(io.dgf_handler, uw_path);
    // fef.surface.write(io.dgf_handler, fef.tmpFile_path);
    break;
  };
}
void Brusselator_scheme::pattern_loop(){
  for(long it = 0; it < pattern_timeSteps(); ++it){
    update_surface();
    update_scalar_solution();
    rhs_and_solve_SPDE();
    update_velocity();    
    error_on_intSurface(); // Error on surface(t_n)
    next_timeStep();
    Pattern_helper helper {*this};
    helper.finalize_scalarPDE_rhs();
    helper.solve_scalarPDE();
    // helper.errors_on_numSurface();
    // helper.plot_paraview();
  }
}
void Brusselator_scheme::final_action(){
  fef.u.write(io.dgf_handler, "./final_");
  fef.w.write(io.dgf_handler, "./final_");
  fef.surface.write(io.dgf_handler, "./final_"); 
}

// ------------------------------------------------------------
// Brusselator_scheme private

void Brusselator_scheme::update_surface(){
  io.identity.interpolate(fef.surface.exact);

  // Save old surface for the velocity
  fef.surface.app = fef.surface.fun; 
}
void Brusselator_scheme::update_velocity(){
  exact.v.interpolate(fef.velocity.exact);

  // update fef.velocity.fun
  const auto dT = fix_grid.time_provider().deltaT();
  auto xNew_ptr = fef.surface.fun.cbegin();
  auto xOld_ptr = fef.surface.app.cbegin();
  for(auto v_ptr = fef.velocity.fun.begin();
      v_ptr != fef.velocity.fun.end();
      ++v_ptr, ++xNew_ptr, ++xOld_ptr)
    *v_ptr = (*xNew_ptr - *xOld_ptr)/ dT;
}
void Brusselator_scheme::update_scalar_solution(){
  exact.u.interpolate(fef.u.exact);
  exact.w.interpolate(fef.w.exact);
}

void Brusselator_scheme::error_on_intSurface(){
  const auto dT = fix_grid.time_provider().deltaT();
  // scalar_error
  io.u << dT << '\t'
       << norm.l2_err(fef.u.exact, fef.u.fun) << '\t'
       << norm.h1_err(fef.u.exact, fef.u.fun) << std::endl;
  io.w << dT << '\t'
       << norm.l2_err(fef.w.exact, fef.w.fun) << '\t'
       << norm.h1_err(fef.w.exact, fef.w.fun) << std::endl;
  // surface_error
  io.surface << dT << '\t'
	     << norm.l2_err(fef.surface.exact, fef.surface.app)
	     << '\t'
	     << norm.h1_err(fef.surface.exact, fef.surface.app)
	     << std::endl;
  // velocity error
  io.velocity << dT << '\t'
	      << norm.l2_err(fef.velocity.fun, fef.velocity.exact)
	      << '\t'
	      << norm.h1_err(fef.velocity.fun, fef.velocity.exact)
	      << std::endl;
}

void Brusselator_scheme::pre_loop_action(){
  PreLoop_helper helper {*this};
  // helper.random_initialValues();
  helper.analytic_initialValues();
  helper.headLine_in_errFile();
  helper.save_surface();
  // helper.plot_errors_in_errFile();
  // helper.plot_paraview();
  // next_timeStep();
}
void Brusselator_scheme::rhs_and_solve_SPDE(){
  RhsAndSolve_helper helper {*this};
  helper.scalar_massMatrix();
  helper.brusselator_rhs();
  helper.addScaled_surfaceLoadVector();
  helper.solve_surface_and_save();
}
// ----------------------------------------------------------------------
// Implementation of structs
// Brusselator_scheme::Fef and Brusselator_scheme::Io

Brusselator_scheme::Fef::Fef(const Esfem::Grid::Grid_and_time& gt)
  :u {"u", gt}, w {"w", gt},
   surface {"surface", gt}, velocity {"velocity", gt}
{}

Brusselator_scheme::Io::Io(const Esfem::Io::Parameter& p)
  :dgf_handler {p.grid()}, u {"_u", p}, w {"_w", p},
   surface {"_X", p}, velocity {"_v",p}
{}

Brusselator_scheme::Init_data::Init_data(const Esfem::Grid::Grid_and_time& gt)
  : u {gt, Growth::promoting}, w {gt, Growth::inhibiting}, v {gt}
{}
Brusselator_scheme::Init_data::Init_data(const Esfem::Grid::Grid_and_time& gt, const Esfem::Io::Parameter& p)
  :u {p, Growth::promoting}, w {p, Growth::inhibiting}, v {gt}
{}
