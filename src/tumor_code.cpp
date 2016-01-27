/*! \file tumor_code.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 21.12.2015

     Implementation details for tumor_code.h
     Created by Christian Power on 21.12.2015
     Copyright (c) 2015 Christian Power. All rights reserved.
 */

#include "tumor_code.h"

#include <config.h>

#include <iostream>
#include <algorithm>

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

#include "tumor_bdf.h"
#include "tumor_deformation.h"
#include "tumor_growth.h"
#include "tumor_typedef.h"
#include "iodof.h"
// #include "/Users/christianpower/cpp/ODE_Solver/bdf.h"

class DataOutputParameters :
  public Dune::Fem::
  LocalParameter<Dune::Fem::DataOutputParameters, DataOutputParameters> {
public:
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

void tumor_growth_code(){
  // get time from parameter file
  const double t_0 = 
    Dune::Fem::Parameter::getValue<double>("heat.starttime",0.0);
  const double dT = 
    Dune::Fem::Parameter::getValue<double>("heat.timestep",0.1);
  const double t_end = 
    Dune::Fem::Parameter::getValue<double>("heat.endtime",0.6);
  const double t_sg =
    Dune::Fem::Parameter::getValue<double>("tumor_growth.surface.start_growth", 5.);
  const int time_step_no_max = (t_end - t_0)/dT + .1;
  const int pattern_time_step_no_max = (t_sg - t_0)/dT + .1;
  const int sg_time_step_no_max = (t_end - t_sg - t_0)/dT + .1;
  
  if((sg_time_step_no_max + pattern_time_step_no_max) != time_step_no_max)
    throw std::runtime_error {"Error in tumor_growth_code():\nSomething"
	" is wrong with the time step integers!"};
  
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

  R3Identity r3identity {timeProvider};
  
  // Surface operator
  const double solverEps 
    = Dune::Fem::Parameter::getValue<double>("heat.solvereps", 1e-8);
  
  // Tumor growth model & operator
  // Surface model
  const double alpha 
    = Dune::Fem::Parameter::getValue<double>("tumor_growth.surface.alpha", 1e-3);
  const double epsilon
    = Dune::Fem::Parameter::getValue<double>("tumor_growth.surface.epsilon", 1e-2);
  const double delta
    = Dune::Fem::Parameter::getValue<double>("tumor_growth.surface.delta", .4);

  const Tumor_surface_model surface_model {timeProvider, alpha, epsilon, delta};  
  const Tumor_surface_operator surface_operator {surface_model, u_approx};
  Vec_CG_Solver surface_solver {surface_operator, solverEps, solverEps};
  
  // Brusselator model
  const double a 
    = Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.a", .1); 
  const double b
    = Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.b", .9); 
  const double dc
    = Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.Dc", 10.); 
  const double gamma
    = Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.gamma", 100.);
  using Growth = Tumor_growth::Growth;

  const Tumor_Brusselator_model u_model {dT, gamma, a, Growth::promoting};
  const Tumor_Brusselator_model w_model {dT, gamma, b, Growth::inhibiting, dc};

  // Initial data
  const double hom_solution_u = a + b;
  if(hom_solution_u * hom_solution_u < solverEps)
    throw std::runtime_error {"Error in tumor_growth_code(): Bad \"a+b\"."};
  const double hom_solution_w = b / (hom_solution_u * hom_solution_u);
  const double pertu = Dune::Fem::Parameter::
    getValue<double>("tumor_growth.heat.init_data.pertubation", .01);

  const Initial_data_uniform initial_data_u {hom_solution_u, pertu};
  const Initial_data_uniform initial_data_w {hom_solution_w, pertu};

  // Brusselator operators
  const Tumor_Brusselator_operator u_operator {u_model, u_approx, w_approx};
  const Tumor_Brusselator_operator w_operator {w_model, u_approx, u_approx};
  LinearInverseOperatorType u_solver {u_operator, solverEps, solverEps};
  LinearInverseOperatorType w_solver {w_operator, solverEps, solverEps};

  // Norms and file output
  L2NormType l2norm(gridPart);
  H1NormType h1norm(gridPart);

  using TG_IOTuple = Dune::tuple<DiscreteFunctionType*, DiscreteFunctionType*>;
  TG_IOTuple uw_tuple(&u_n, &w_n);
  using TG_dataOutput = Dune::Fem::DataOutput<GridType, TG_IOTuple>;
  const std::string o_filename 
  {Dune::Fem::Parameter::getValue<std::string>("tumor_growth.io.uw.output_name",
					       "../output/ALE_LiteratureExample2-")};
  TG_dataOutput uw_dataOutput(grid, uw_tuple, DataOutputParameters{o_filename, 0});
  

  // starting calculation
  timeProvider.init(dT);     // Do your first action before you enter the for-loop.
  int time_step_no = 0;

  // evoMap.save_original_vertices(); already done
  InterpolationType::interpolateFunction(initial_data_u, u_n);
  InterpolationType::interpolateFunction(initial_data_w, w_n);

  u_approx.assign(u_n);
  w_approx.assign(w_n);
  // This space is useful for BDF codes
  uw_dataOutput.write(timeProvider);


  // pre-growth/ pattern stage

  std::cout << "Pattern stage" << std::endl;
  for(;
      time_step_no <= pattern_time_step_no_max;
      timeProvider.next(dT), ++time_step_no){

    u_operator.generate_rhs_old_surface(u_n, u_cg_rhs);
    u_operator.generate_rhs_new_surface(tmp_fef);
    u_cg_rhs.axpy(1., tmp_fef);
    u_solver(u_cg_rhs, u_n);
    u_approx.assign(u_n);

    w_operator.generate_rhs_old_surface(w_n, w_cg_rhs);
    w_operator.generate_rhs_new_surface(tmp_fef);	// this can be done faster
    w_cg_rhs.axpy(1., tmp_fef);
    w_solver(w_cg_rhs, w_n);
    w_approx.assign(w_n);
    
    uw_dataOutput.write(timeProvider);
  }

  std::cout << "Growth stage" << std::endl;
  for(; 
      time_step_no <= time_step_no_max; 
      timeProvider.next(dT), ++time_step_no){
  
    R3Interpolation::interpolateFunction(r3identity, next_surface); 
  
    surface_operator.generate_rhs(next_surface, surface_cg_rhs);
  
    u_operator.generate_rhs_old_surface(u_n, u_cg_rhs);
    w_operator.generate_rhs_old_surface(w_n, w_cg_rhs);
  
    surface_solver(surface_cg_rhs, next_surface);
  
    evoMap.evolve(next_surface);
  
    u_operator.generate_rhs_new_surface(tmp_fef);
    u_cg_rhs.axpy(1., tmp_fef);
    w_operator.generate_rhs_new_surface(tmp_fef);	// this can be done faster
    w_cg_rhs.axpy(1., tmp_fef);
  
    u_solver(u_cg_rhs, u_n);
    w_solver(w_cg_rhs, w_n);
  
    u_approx.assign(u_n);
    w_approx.assign(w_n);
    
    uw_dataOutput.write(timeProvider);
  }  
}

// oid tumor_growth_code_old(){ 
//  // get time from parameter file
//  double t_0 = 
//    Dune::Fem::Parameter::getValue<double>("heat.starttime",0.0);
//  double dT = 
//    Dune::Fem::Parameter::getValue<double>("heat.timestep",0.1);
//  double t_end = 
//    Dune::Fem::Parameter::getValue<double>("heat.endtime",0.6);
//  const int time_step_no_max = (t_end - t_0)/dT + .1;
// 
//  // prepare grid from DGF file
//  const std::string gridkey =
//    Dune::Fem::IOInterface::defaultGridKey( GridType::dimension );
//  const std::string gridfile =
//    Dune::Fem::Parameter::getValue< std::string >( gridkey );
//  if( Dune::Fem::MPIManager::rank() == 0 )
//    std::cout << "Loading macro grid: " << gridfile << std::endl;
// 
//  Dune::GridPtr<HostGridType> hostGrid {gridfile};
//  hostGrid->loadBalance();
// 
//  // create grid
//  BDF::EvoMapType evoMap {gridfile};
//  evoMap.save_original_vertices();
//  DeformationCoordFunction deformation {evoMap};
//  GridType grid (*hostGrid, deformation);
//  GridPartType gridPart (grid);
//  DiscreteFunctionSpaceType dfSpace (gridPart);
//  Vec_FE_Space r3dfSpace (gridPart);
// 
//  Dune::Fem::GridTimeProvider<GridType> timeProvider(t_0, grid);
//  deformation.set_time_provider(timeProvider);
// 
//  // create vec_FE-functions
//  Vec_FE_Fun next_surface {"next_surface", r3dfSpace};	// surface for next step
//  Vec_FE_Fun surface_cg_rhs {"surface_cg_rhs", r3dfSpace};
// 
//  // create FE-functions
//  DiscreteFunctionType u_n {"u_n", dfSpace};		// Numerical solution
//  DiscreteFunctionType u_approx {"u_approx", dfSpace};	// Approximaion to u_n
//  DiscreteFunctionType u_cg_rhs {"u_cg_rhs", dfSpace};
//  DiscreteFunctionType w_n {"w_n", dfSpace};		// Numerical solution
//  DiscreteFunctionType w_approx {"w_approx", dfSpace};	// Approximaion to u_n
//  DiscreteFunctionType w_cg_rhs {"w_cg_rhs", dfSpace};
//  DiscreteFunctionType tmp_fef {"tmp_fef", dfSpace};	// tmp finite element function
// 
//  // Abstract functions
//  const double mean_value_u
//    = Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.init"
// 					     "_data.mean_value.u", 1.);
//  const double variance_u
//    = Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.init"
// 					     "_data.variance.u", .5);
//  const double mean_value_w
//    = Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.init"
// 					     "_data.mean_value.w", .9);
//  const double variance_w
//    = Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.init"
// 					     "_data.variance.w", .5);
// 
//  const Initial_data_u initial_data_u {mean_value_u, variance_u};
//  const Initial_data_w initial_data_w {mean_value_w, variance_w};
//  R3Identity r3identity {timeProvider};
//  
//  // Tumor growth model & operator
//  // Surface model
//  const double alpha 
//    = Dune::Fem::Parameter::getValue<double>("tumor_growth.surface.alpha", 1e-3);
//  const double epsilon
//    = Dune::Fem::Parameter::getValue<double>("tumor_growth.surface.epsilon", 1e-2);
//  const double delta
//    = Dune::Fem::Parameter::getValue<double>("tumor_growth.surface.delta", .4);
// 
//  const Tumor_surface_model surface_model {timeProvider, alpha, epsilon, delta};
// 
//  // Surface operator
//  const double solverEps 
//    = Dune::Fem::Parameter::getValue<double>("heat.solvereps", 1e-8);
//  
//  const Tumor_surface_operator surface_operator {surface_model, u_approx};
//  Vec_CG_Solver surface_solver {surface_operator, solverEps, solverEps};
// 
//  // Brusselator model
//  const double a 
//    = Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.a", .1); 
//  const double b
//    = Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.b", .9); 
//  const double dc
//    = Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.Dc", 10.); 
//  const double gamma
//    = Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.gamma", 100.);
//  using Growth = Tumor_growth::Growth;
// 
//  const Tumor_Brusselator_model u_model {dT, gamma, a, Growth::promoting};
//  const Tumor_Brusselator_model w_model {dT, gamma, b, Growth::inhibiting, dc};
//  
//  // Brusselator operators
//  const Tumor_Brusselator_operator u_operator {u_model, u_approx, w_approx};
//  const Tumor_Brusselator_operator w_operator {w_model, u_n, u_n};
//  LinearInverseOperatorType u_solver {u_operator, solverEps, solverEps};
//  LinearInverseOperatorType w_solver {w_operator, solverEps, solverEps};
// 
//  // Norms and file output
//  L2NormType l2norm(gridPart);
//  H1NormType h1norm(gridPart);
//  // std::ofstream l2h1error_ofs {Dune::Fem::
//  //     Parameter::getValue<std::string>("fem.io.errorFile", "../output/l2h1error"),
//  //     std::ios_base::app};
// 
//  // Paraview output
//  IOTupleType ioTuple(&u_n, &w_n);
//  const int step = 0;	// is there any resonable choice, whitout adaptivity?
//  DataOutputType 
//    dataOutput(grid, ioTuple, 
// 		 DataOutputParameters(Dune::Fem::Parameter::getValue<std::string>
// 				      ("tumor_growth.io.u.output_name",
// 				       "../output/ALE_LiteratureExample-"),
// 				      step) );
//  // auto data_writer = [](DataOutputType& u_data, DataOutputType& w_data,
//  // 			const Dune::Fem::GridTimeProvider<GridType>& tp){
//  //   u_data.write(tp);
//  //   w_data.write(tp);
//  // };
//  // helper function for 'err_vec'
//  // auto calc_err_vec = [&u_np1, &exact_solution, &err_vec]{
//  //   err_vec.assign(u_np1);
//  //   err_vec.axpy((-1.), exact_solution);
//  // };
// 
//  const std::string u_prefix = Dune::Fem::Parameter::
//    getValue<std::string>("tumor_growth.io.u_prefix",
//  			  "/Users/christianpower/cpp/DISS_surfaces/"
//  			  "data/test_video/u_Sphere_4.txt");
//  const std::string w_prefix = Dune::Fem::Parameter::
//    getValue<std::string>("tumor_growth.io.w_prefix",
//  			  "/Users/christianpower/cpp/DISS_surfaces/"
//  			  "data/test_video/w_Sphere_4.txt");
//  const std::string x_prefix = Dune::Fem::Parameter::
//    getValue<std::string>("tumor_growth.io.X_prefix",
//  			  "/Users/christianpower/cpp/DISS_surfaces/"
//  			  "data/test_video/X_Sphere_4.txt");
//  IO_dune_fem u_input(u_prefix, IO_direction::write);
//  IO_dune_fem w_input(w_prefix, IO_direction::write);
//  IO_dune_fem x_input(x_prefix, IO_direction::write);
//  
//  // DiscreteFunctionType uN_read {"uN_read", dfSpace};		
//  // DiscreteFunctionType wN_read {"wN_read", dfSpace};		
//  // Vec_FE_Fun surface_read {"surface_read", r3dfSpace};
//  
//  // ----------------------------------------------------------------------
//  // starting calculation
//  timeProvider.init(dT);     // Do your first action before you enter the for-loop.
//  int time_step_no = 0;
// 
//  // evoMap.save_original_vertices(); already done
//  InterpolationType::interpolateFunction(initial_data_u, u_n);
//  InterpolationType::interpolateFunction(initial_data_w, w_n);
// 
//  u_approx.assign(u_n);
//  w_approx.assign(w_n);
// 
//  // This space is useful for BDF codes
// 
//  dataOutput.write(timeProvider);
//  // data_writer(u_dataOutput, w_dataOutput, timeProvider);
//  
//  for(; 
//      time_step_no <= time_step_no_max; 
//      timeProvider.next(dT), ++time_step_no){
//  
//    R3Interpolation::interpolateFunction(r3identity, next_surface); 
//  
//    surface_operator.generate_rhs(next_surface, surface_cg_rhs);
//  
//    u_operator.generate_rhs_old_surface(u_n, u_cg_rhs);
//    w_operator.generate_rhs_old_surface(w_n, w_cg_rhs);
// 
//    u_input(u_n);
//    w_input(w_n);
//    x_input(next_surface);
//    
//    surface_solver(surface_cg_rhs, next_surface);
//  
//    evoMap.evolve(next_surface);
//  
//    u_operator.generate_rhs_new_surface(tmp_fef);
//    u_cg_rhs.axpy(1., tmp_fef);
//    w_operator.generate_rhs_new_surface(tmp_fef);	// this can be done faster
//    w_cg_rhs.axpy(1., tmp_fef);
//  
//    u_solver(u_cg_rhs, u_n);
//    w_solver(w_cg_rhs, w_n);
//  
//    u_approx.assign(u_n);
//    w_approx.assign(w_n);
//    
//    dataOutput.write(timeProvider);
//    // data_writer(u_dataOutput, w_dataOutput, timeProvider);
//    
// / u_input(uN_read);
// / w_input(wN_read);
// / x_input(surface_read);
// / std::cout << std::scientific
// / 	      << "u_n - uN_read: " << l2norm.distance(u_n, uN_read) << '\n'
// /   	      << "w_n - wN_read: " << l2norm.distance(w_n, wN_read) << '\n'
// / 	      << "next_surface - surface_read: "
// / 	      << l2norm.distance(next_surface, surface_read) << '\n'
// / 	      << std::defaultfloat << std::endl;
//  }

/*! Log:
 */
