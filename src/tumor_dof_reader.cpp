/*! \file tumor_dof_reader.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 18.12.2015

     Implementation details for tumor_dof_reader.h
     Created by Christian Power on 18.12.2015
     Copyright (c) 2015 Christian Power. All rights reserved.
 */

#include <config.h>
#include "tumor_dof_reader.h"

#include <dune/grid/geometrygrid.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/oemsolver.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/dataoutput.hh>

#include <algorithm>
#include <iostream>

#include "tg_dune_bdf.hpp"
#include "tg_deformation.hh"
#include "elliptic.hh"
#include "heat.hh"
#include "rhs.hh"
#include "dune_typedef_tg.hpp"
#include "/Users/christianpower/cpp/ODE_Solver/bdf.h"
#include "iodof.h"

struct DataOutputParameters:
  public Dune::Fem::
  LocalParameter<Dune::Fem::DataOutputParameters, DataOutputParameters> {
  DataOutputParameters(const std::string name, const int step)
    : name_(name), step_( step )
  {}
  DataOutputParameters(const DataOutputParameters& other)
    : step_( other.step_ )
  {}
  std::string prefix() const {
    std::stringstream s;
    s << name_ << step_ << "-";
    return s.str();
  }
private:
  std::string name_;
  int step_;
};

void tumorGrowth_read_VM_data(){
  
  double t_0 = 
    Dune::Fem::Parameter::getValue<double>("heat.starttime",0.0);
  double dT = 
    Dune::Fem::Parameter::getValue<double>("heat.timestep",0.1);
  double t_end = 
    Dune::Fem::Parameter::getValue<double>("heat.endtime",0.6);
  const int time_step_no_max = (t_end - t_0)/dT + .1;

  // prepare grid from DGF file
  const std::string gridkey =
    Dune::Fem::IOInterface::defaultGridKey( GridType::dimension );
  const std::string gridfile =
    Dune::Fem::Parameter::getValue< std::string >( gridkey );
  if( Dune::Fem::MPIManager::rank() == 0 )
    std::cout << "Loading macro grid: " << gridfile << std::endl;

  Dune::GridPtr<HostGridType> hostGrid {gridfile};
  hostGrid->loadBalance();

  // create grid
  BDF::EvoMapType evoMap {gridfile};
  evoMap.save_original_vertices();
  DeformationCoordFunction deformation {evoMap};
  GridType grid (*hostGrid, deformation);
  GridPartType gridPart (grid);
  DiscreteFunctionSpaceType dfSpace (gridPart);
  Vec_FE_Space r3dfSpace (gridPart);

  Dune::Fem::GridTimeProvider<GridType> timeProvider(t_0, grid);
  deformation.set_time_provider(timeProvider);

  // create vec_FE-functions
  Vec_FE_Fun next_surface {"next_surface", r3dfSpace};	// surface for next step
  Vec_FE_Fun surface_cg_rhs {"surface_cg_rhs", r3dfSpace};

  // create FE-functions
  DiscreteFunctionType u_n {"u_n", dfSpace};		// Numerical solution
  DiscreteFunctionType u_approx {"u_approx", dfSpace};	// Approximaion to u_n
  DiscreteFunctionType u_cg_rhs {"u_cg_rhs", dfSpace};
  DiscreteFunctionType w_n {"w_n", dfSpace};		// Numerical solution
  DiscreteFunctionType w_approx {"w_approx", dfSpace};	// Approximaion to u_n
  DiscreteFunctionType w_cg_rhs {"w_cg_rhs", dfSpace};
  DiscreteFunctionType tmp_fef {"tmp_fef", dfSpace};	// tmp finite element function

  using MyIOTuple = Dune::tuple<DiscreteFunctionType*, DiscreteFunctionType*>;
  using MyOutput = Dune::Fem::DataOutput<GridType, MyIOTuple>;
  MyIOTuple ioTuple(&u_n, &w_n);

  const std::string o_filename = Dune::Fem::Parameter::getValue<std::string>
			 ("tumor_growth.io.uw.output_name",
			  "output/tg_video-");
  MyOutput uw_dataOutput(grid, ioTuple, DataOutputParameters(o_filename, 0) );

  const std::string u_prefix = Dune::Fem::Parameter::
    getValue<std::string>("tumor_growth.io.u_prefix",
			  "/Users/christianpower/cpp/DISS_surfaces/"
			  "data/test_video/u_Sphere_4.txt");
  const std::string w_prefix = Dune::Fem::Parameter::
    getValue<std::string>("tumor_growth.io.w_prefix",
			  "/Users/christianpower/cpp/DISS_surfaces/"
			  "data/test_video/w_Sphere_4.txt");
  const std::string x_prefix = Dune::Fem::Parameter::
    getValue<std::string>("tumor_growth.io.X_prefix",
			  "/Users/christianpower/cpp/DISS_surfaces/"
			  "data/test_video/X_Sphere_4.txt");
  IO_dune_fem u_input(u_prefix, IO_direction::read);
  IO_dune_fem w_input(w_prefix, IO_direction::read);
  IO_dune_fem x_input(x_prefix, IO_direction::read);

  for(int time_step_no = 0; 
      time_step_no <= time_step_no_max; 
      timeProvider.next(dT), ++time_step_no){
    u_input(u_n);
    w_input(w_n);
    x_input(next_surface);
    evoMap.evolve(next_surface);
    
    uw_dataOutput.write(timeProvider);
  }
}

/*! Log:
 */
