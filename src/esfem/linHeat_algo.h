/*! \file linHeat_algo.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Februar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 03.02.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef LINHEAT_ALGO_H
#define LINHEAT_ALGO_H 

#include <config.h>
// #include <dune/fem/solver/cginverseoperator.hh>
// #include <dune/fem/solver/oemsolver.hh>
// #include <dune/fem/quadrature/quadrature.hh>
// #include <dune/common/fmatrix.hh>
#include "esfem.h"

void linHeat_algo(int argc, char** argv){
  using namespace Esfem;
  
  Dune::Fem::MPIManager::initialize(argc, argv);

  const auto parameter_file =
    "/Users/christianpower/cpp/DISS_surfaces/data/tumor_parameter.txt";  
  Io::Parameter data {argc, argv, parameter_file};
#ifdef DEBUG
  std::clog << data << std::endl;
#endif

  Grid::Grid_and_time grid {data};
  Grid::Scal_FEfun exact_solution {"exact_solution", grid};
  Grid::Scal_FEfun numerical_solution {"numerical_solution", grid};
  Grid::Scal_FEfun tmp_fef {"tmp", grid};
  
  SecOrd_op::Init_data init_data {grid};
  SecOrd_op::Rhs rhs {grid};
  SecOrd_op::Linear_heat solver {data, grid};
  
  const Io::L2H1_calculator err_cal {grid, exact_solution, numerical_solution};
  Io::Error_stream err_log {data};
  err_log << std::scientific;
  Io::Paraview paraview_plot {data, grid, exact_solution, numerical_solution};
  
  init_data.interpolate(numerical_solution);
  init_data.interpolate(exact_solution);
  solver.mass_matrix(numerical_solution, tmp_fef);

  err_log << data.global_timeStep() << ' '
	  << err_cal.l2_err() << ' ' << err_cal.h1_err() << std::endl;
  // paraview_plot.write();
  
  grid.next_timeStep(data.global_timeStep());
  for(long it =0; it < data.max_timeSteps(); ++it){
    // Mu^n+1 + tau Au^n+1 = M^n + tau f^n+1
    
    rhs.assemble_and_addScaled_to(tmp_fef);
    solver.solve(tmp_fef, numerical_solution);
    solver.mass_matrix(numerical_solution, tmp_fef);

    init_data.interpolate(exact_solution);
    err_log << data.global_timeStep() << ' '
	    << err_cal.l2_err() << ' ' << err_cal.h1_err() << std::endl;
    // paraview_plot.write();
    
    grid.next_timeStep(data.global_timeStep());
  }
}

#endif // LINHEAT_ALGO_H

/*! Log:
 */
