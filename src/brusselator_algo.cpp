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
     \date 23. April 2016
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
  fef {fix_grid},
  //  exact {data} // constructor for random initial data
  exact {fix_grid}
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
    fef.surface.write(io.dgf_handler, fef.tmpFile_path);
    break;
  };
}
void Brusselator_scheme::pattern_loop(){
  for(long it = 0; it < pattern_timeSteps(); ++it, next_timeStep()){
    rhs_and_solve_SPDE();
    Pattern_helper helper {*this};
    helper.finalize_scalarPDE_rhs();
    helper.solve_scalarPDE();
    helper.update_exactSolutions();
    helper.plot_errors_in_errFile();
    helper.plot_paraview();
  }
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
  // helper.random_initialValues();
  helper.analytic_initialValues();
  helper.headLine_in_errFile();
  // helper.plot_errors_in_errFile();
  helper.plot_paraview();
  // helper.prepare_rhs();
  next_timeStep();
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
  :u {"u", gt}, w {"w", gt}, surface {"surface", gt}
{}

Brusselator_scheme::Io::Io(const Esfem::Io::Parameter& p)
  :dgf_handler {p.grid()}, u {"_u", p}, w {"_w", p}
{}
Brusselator_scheme::Init_data::Init_data(const Esfem::Grid::Grid_and_time& gt)
  : u {gt, Growth::promoting}, w {gt, Growth::inhibiting}
{}
Brusselator_scheme::Init_data::Init_data(const Esfem::Io::Parameter& p)
  :u {p, Growth::promoting}, w {p, Growth::inhibiting} 
{}
