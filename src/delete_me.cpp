// begin 2014nonlinear experiment
void BDF::nonlinear_algorithm(){
  // get time from parameter file
  double t_0 = 
    Dune::Fem::Parameter::getValue<double>("heat.starttime",0.0);
  double dT = 
    Dune::Fem::Parameter::getValue<double>("heat.timestep",0.1);
  double	t_end = 
    Dune::Fem::Parameter::getValue<double>("heat.endtime",0.6);
  const int time_step_no_max = (t_end - t_0)/dT + .1;

  // setting up BDF coefficients
  const int bdf_no =
    Dune::Fem::Parameter::getValue<double>("heat.bdf", 1);
  const std::vector<double> bdf_alpha_coeff { NUMERIK::bdf_alpha_coeff(bdf_no) };
  const std::vector<double> bdf_gamma_coeff { NUMERIK::bdf_gamma_coeff(bdf_no) };
  // bdf_*_coeff.back() is the lead coefficient of the polynomial

  // debugging: print BDF coeffs
  // for(double d : bdf_gamma_coeff)
  // 	std::cout << d << ' ';
  // std::cout << std::endl;
  // for(double d : bdf_alpha_coeff)
  // 	std::cout << d << ' ';
  // std::cout << '\n' << bdf_alpha_coeff.back() << std::endl;

  // prepare grid from DGF file
  const std::string gridkey =
    Dune::Fem::IOInterface::defaultGridKey( GridType::dimension );
  std::cout << "gridkey: " << gridkey << std::endl;
  const std::string gridfile =
    Dune::Fem::Parameter::getValue< std::string >( gridkey );
  if( Dune::Fem::MPIManager::rank() == 0 )
    std::cout << "Loading macro grid: " << gridfile << std::endl;
  Dune::GridPtr< HostGridType > hostGrid {gridfile};
  hostGrid ->loadBalance();

  // create grid
  DeformationCoordFunction deformation {};
  GridType grid(*hostGrid, deformation);
  GridPartType gridPart(grid);
  DiscreteFunctionSpaceType dfSpace(gridPart);

  Dune::Fem::GridTimeProvider<GridType> timeProvider(t_0, grid);
  deformation.set_time_provider(timeProvider);

  // create FE-functions
  DiscreteFunctionType U_np1("U_np1", dfSpace);	// Numerical solution
  // In a math book U_np1 == U_n+1.  
  DiscreteFunctionType rhs("rhs", dfSpace);	// Right hand side for the LES
  DiscreteFunctionType load_vector("load_vector", dfSpace);
  DiscreteFunctionType xi("xi", dfSpace);		// For the lineary implicit BDF method
  DiscreteFunctionType exact_solution("exact_solution", dfSpace);	
  DiscreteFunctionType err_vec("dof_U_np1_minus_exact_solution", dfSpace);
  // Holds U_np1 - exact_solution
  std::vector<DiscreteFunctionType> prev_steps_U_nmk;	// All previous U_{n-k}
  for(int i = 0; i < bdf_no; ++i)
    prev_steps_U_nmk.push_back( 
			       DiscreteFunctionType {"U_nm" + std::to_string(bdf_no - (i+1)), dfSpace} );
  std::vector<DiscreteFunctionType> prev_steps_M_U_nmk;	// All previous M_U_{n-k}
  for(int i = 0; i < bdf_no; ++i)
    prev_steps_M_U_nmk.push_back(
				 DiscreteFunctionType {"M_U_nm" + std::to_string(bdf_no - (i+1)), dfSpace});

  // stupid check
  if(prev_steps_U_nmk.size() != prev_steps_M_U_nmk.size())
    throw std::runtime_error("ERROR in bdf_cycle().");

  // File output
  L2NormType l2norm(gridPart);
  H1NormType h1norm(gridPart);
  std::ofstream l2h1error_ofs
  {Dune::Fem::Parameter::getValue<std::string>("fem.io.errorFile", 
					       "output/l2h1error"),
      std::ios_base::app};
  // Paraview output
  IOTupleType ioTuple(&err_vec);
  const int step = 0;	// is there any resonable choice, whitout adaptivity?
  DataOutputType 
    dataOutput(grid, ioTuple, 
	       DataOutputParameters(Dune::Fem::Parameter::getValue<std::string>
				    ("fem.io.outputName",
				     "../output/ALE_LiteratureExample-"),
				    step) );
  // helper function for 'err_vec'
  auto calc_err_vec = [&]{
    err_vec.assign(U_np1);
    err_vec.axpy((-1.), exact_solution);
  };

  InitialDataType initialData {timeProvider};
  RHSFunctionType f {timeProvider};
  NonlinearModel model {timeProvider};
  NonlinearOperator ellipticOp {xi, model, bdf_alpha_coeff.back()};
  const double solverEps =
    Dune::Fem::Parameter::getValue<double>("heat.solvereps", 1e-8); 
  LinearInverseOperatorType solver(ellipticOp, solverEps, solverEps);	// CG-solver

  const std::string ofile_name = Dune::Fem::Parameter::
    getValue<std::string>("fem.io.errorFile", "output/l2h1error");
	
  IO_dune_fem hd_dof_com {ofile_name, IO_direction::write};
	
  auto write_error = [&](std::ostream& os)
  // helper function for 'l2h1error_ofs'
    {
      // os << std::defaultfloat << timeProvider.deltaT() << ' ' 
      // << std::scientific 
      // << l2norm.distance(exact_solution, U_np1) << ' '
      // << h1norm.distance(exact_solution, U_np1) << std::endl;
      hd_dof_com(U_np1);
    };

  // initial steps
  timeProvider.init(dT);     // Do your first action before you enter the for loop.
  int time_step_no = 0;

  // debugging
  // std::cout << "\n=========="  << std::endl;

  // initialize 'prev_steps_U_nmk' and 'prev_steps_M_U_mmk'
  for(;
      time_step_no < std::min(prev_steps_U_nmk.size(), size_t(time_step_no_max));
      ++time_step_no, timeProvider.next(dT)){
		
    // debugging
    // std::cout << std::defaultfloat << "time = " << timeProvider.time() 
    //  		  << std::scientific << std::endl;

    InterpolationType::interpolateFunction(initialData, exact_solution);
    prev_steps_U_nmk.at(time_step_no).assign(exact_solution);
    ellipticOp.mass_matrix(exact_solution, prev_steps_M_U_nmk.at(time_step_no));

    U_np1.assign(exact_solution);

    // output
    // calc_err_vec();
    // dataOutput.write(timeProvider);
    write_error(l2h1error_ofs);

    // debugging
    // std::cout << "exact_solution = " << *exact_solution.dbegin() << '\n'
    // 		  << "U_np1 = " << *(U_np1.dbegin()) << '\n'
    // 		  << "prev_steps_M_U_nmk.at(time_step_no) = " 
    // 		  << *prev_steps_M_U_nmk.at(time_step_no).dbegin()
    // 		  << std::endl;
  }

  // debugging
  // std::cout << "\n=========="  << std::endl;

  for(; 
      time_step_no <= time_step_no_max; 
      timeProvider.next(dT), ++time_step_no){

    // debugging
    // std::cout << std::defaultfloat << "time = " << timeProvider.time() 
    // 		  << std::scientific << std::endl;

    // 'rhs', i.e.
    // calculating: rhs = (-1) · (∑ᵢ₌₁ᵏ αᵢ₋₁ (Mu)ⁿ⁻ᵏ⁺ⁱ)
    //				xi = ∑ᵢ₌₁ᵏ γᵢ₋₁ uⁿ⁻ᵏ⁺ⁱ
    // NOTE for implicit euler: rhs = Mⁿ uⁿ and xi = uⁿ
    rhs.clear();	// set dof to zero
    // xi.clear();
    for(size_t i=0; i < prev_steps_U_nmk.size(); ++i){
      rhs.axpy(bdf_alpha_coeff.at(i), prev_steps_M_U_nmk.at(i));
      // xi.axpy(bdf_gamma_coeff.at(i), prev_steps_U_nmk.at(i));
      // std::cout << "rhs.axpy(bdf_alpha_coeff.at(i), "
      // 	"prev_steps_M_U_nmk.at(i)) = " << *rhs.dbegin() << std::endl;
    }
    // std::cout << "xi = " << *xi.dbegin() << std::endl;
    rhs *= (-1);
    // std::cout << "rhs *= (-1) : " << *rhs.dbegin() << std::endl;
    assembleRHS(f, load_vector); 
    // std::cout << "assembleRHS(f, load_vector) = " 
    // 		  << *load_vector.dbegin() << std::endl;
    rhs.axpy(timeProvider.deltaT(), load_vector);
    // std::cout << "rhs.axpy(timeProvider.deltaT(), load_vector) = " 
    // 		  << *rhs.dbegin() << std::endl;

    // 'solver'
    solver(rhs, U_np1);
    // std::cout << "solver(rhs, U_np1) = " << *U_np1.dbegin() << std::endl;

    // cycle bdf values
    for(size_t i=0; i < prev_steps_U_nmk.size() - 1; ++i){
      prev_steps_U_nmk.at(i).assign(prev_steps_U_nmk.at(i+1));
      prev_steps_M_U_nmk.at(i).assign(prev_steps_M_U_nmk.at(i+1));
    }
    prev_steps_U_nmk.back().assign(U_np1);
    ellipticOp.mass_matrix(U_np1, prev_steps_M_U_nmk.back());

    // output
    InterpolationType::interpolateFunction(initialData, exact_solution);		
    write_error(l2h1error_ofs);
    // calc_err_vec();
    // dataOutput.write(timeProvider);
    // if(time_step_no == time_step_no_max){
    // std::cout << std::defaultfloat
    // 		  << "Time = " << timeProvider.time() 
    // 		  << std::scientific << std::endl;
    // InterpolationType::interpolateFunction(initialData, exact_solution);		
    // write_error(l2h1error_ofs);
    // write_error(std::cout);
    // }

    // debugging
    // std::cout << "U_np1 = " << *U_np1.dbegin()
    // 		  << "\n----------" << std::endl;
  }
  l2h1error_ofs.close();
}
