/*********************************************************************
 *  dune_heat_algorithm.hpp                                          *
 *                                                                   *
 *  This header constains the algorithm and an auxiliar class which  *
 *  where original in the heat.cc file                               *
 *                                                                   *
 *  Revision history:                                                *
 *  none                                                             *
 *                                                                   *
 *                                                                   *
 *  Created by Christian Power on 19.06.14.                          *
 *  Copyright (c) 2014 Christian Power. All rights reserved.         *
 *                                                                   *
 *********************************************************************/

#ifndef DUNE_HEAT_ALGORITHM_HPP
#define DUNE_HEAT_ALGORITHM_HPP

// local includes
#include "dune_typedef_heat.hpp"
#include "/Users/christianpower/cpp/ODE_Solver/implicit_euler.h"
#include "/Users/christianpower/cpp/ODE_Solver/bdf.h"

// Remember: For every experiment modify
//	in 'heat.hh'
//	- RHSfunction
//	- InitialData
//	- HeatModel
// in 'elliptic.hh'
//	- EllitpicOperator
// in 'deformation.hh'
//	- DeformationCoordFunction

// DataOutputParameters
// --------------------

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
		// s << "ALE_LiteratureExample-" << step_ << "-";
		s << name_ << step_ << "-";
		return s.str();
	}
private:
	std::string name_;
	int step_;
};


void BDF::ie_heat_algorithm(){
	double t_0 = 
		Dune::Fem::Parameter::getValue<double>("heat.starttime",0.0);
	double dT = 
		Dune::Fem::Parameter::getValue<double>("heat.timestep",0.1);
	double	t_end = 
		Dune::Fem::Parameter::getValue<double>("heat.endtime",0.6);
	const int itno = (t_end - t_0)/dT + .1;
 
	// prepare grid from DGF file
	const std::string gridkey =
		Dune::Fem::IOInterface::defaultGridKey( GridType::dimension );
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

	L2NormType l2norm(gridPart);
	H1NormType h1norm(gridPart);
	std::ofstream l2h1error_ofs
	{Dune::Fem::Parameter::getValue<std::string>("fem.io.errorFile", 
												 "../output/l2h1error")};

	Dune::Fem::GridTimeProvider<GridType> timeProvider(t_0, grid);
	timeProvider.init(dT);     // Do your first action before you enter the for loop.

	deformation.set_time_provider(timeProvider);
	 // deformation is aware of timeProvider

	DiscreteFunctionType U_n("U_n",dfSpace), rhs("rhs", dfSpace), 
		load_vector("load_vector", dfSpace), exactSolution("exactSolution",dfSpace);

	InitialDataType initialData {timeProvider};
	InterpolationType::interpolateFunction( initialData, U_n );
	InterpolationType::interpolateFunction( initialData, exactSolution );
	// solutionContainer is used for higher order BDF methods,
	// U_n == uⁿ, rhs == tmp1, tmp == tmp2
	// (create discrete functions for intermediate functionals)

	// to get the degree of freedom vector
	// for(auto it = U_n.dbegin(); it != U_n.dend(); ++it)
	//	std::cout << *it << std::endl;
	
	RHSFunctionType f {timeProvider};
	// This expression is very long; it's literally the same f as in L(u) = f;
	// together with the function/method in rhs.hh, 

	IOTupleType ioTuple(&U_n);
	const int step = 0;	// is there any resonable choise, whitout adaptivity?
	DataOutputType 
		dataOutput(grid, ioTuple, 
				   DataOutputParameters(Dune::Fem::Parameter::getValue<std::string>
										("fem.io.outputName",
										 "../output/ALE_LiteratureExample-"),
										step) );
	// dataOutput.write(timeProvider);

	auto write_error = [&](std::ostream& os){
		os << std::defaultfloat << timeProvider.time() << ' ' 
		<< std::scientific 
		<< l2norm.distance(exactSolution, U_n) << ' '
		<< h1norm.distance(exactSolution, U_n) << std::endl;
	};

	// write_error(l2h1error_ofs);
	
	NonlinearModel model {timeProvider};
	NonlinearOperator ellipticOp {U_n, model};
	const double solverEps =
		Dune::Fem::Parameter::getValue<double>( "heat.solvereps", 1e-8 );
	LinearInverseOperatorType solver(ellipticOp, solverEps, solverEps);
	// Initializing CG-Solver
	// CAPG: it's very slow. Can't this be speed up somehow?

 	// dataOutput.write(timeProvider);

	ellipticOp.get_xi(U_n);
	ellipticOp.mass_matrix(U_n, rhs);
	// calculating: rhs = Mⁿ uⁿ  

	// debugging
	std::cout << std::scientific << "\n==========\n"
			  << "U_n.dbegin() = " << *(U_n.dbegin()) << std::endl;
	std::cout << "rhs.dbegin() = " << *(rhs.dbegin()) << std::endl;

	int iter = 0;
	for(timeProvider.next(dT); iter < itno; timeProvider.next(dT), ++iter)
		// put your loop-action here, but not the last action
	{
		std::cout << std::defaultfloat << "==========\n"
				  << "In for loop: iter = "<< iter << ", "
				  << "time = " << timeProvider.time() << std::endl;
		
		// f.setTime(timeProvider.time());	// obselete
		// chance the time for the (long) RHS function f to tⁿ⁺¹
		// ('f' from -∆u + \div(v) u + \matdot{u} = f)

		// initialData.setTime(timeProvider.time());	// obselete

		assembleRHS(f, load_vector);
		std::cout << "assembleRHS(f, load_vector) = " << *(load_vector.dbegin()) << std::endl;

		// assemly stiffness/load vector; the vector is called 'rhs'; 
		// in the sense of above it is fⁿ⁺¹

		rhs.axpy(timeProvider.deltaT(), load_vector);
		// rhs += Δt * load_vector
		// Just calculated: rhs == MⁿUⁿ + Δt fⁿ⁺¹
		
		std::cout << "rhs.axpy = " << *(rhs.dbegin()) << std::endl;

		//InterpolationType::interpolateFunction( initialData, U_n );
		solver(rhs, U_n);
		// Solve:  (Mⁿ⁺¹ + Δt A(Uⁿ)) uⁿ⁺¹ = rhs,
		// CAPG remark: this is very very slow.

		std::cout << "solver(rhs, U_n) = " << *(U_n.dbegin()) << std::endl;

		// dataOutput.write(timeProvider);

		InterpolationType::interpolateFunction(initialData, exactSolution);		
		// write_error(l2h1error_ofs);
		
		ellipticOp.mass_matrix(U_n, rhs);
		std::cout << "mass_matrix(U_n, rhs) = " << *(rhs.dbegin()) << std::endl;

		// calculating: rhs = Mⁿ uⁿ  
		ellipticOp.get_xi(U_n);

		// debugging
		std::cout << std::scientific
				  << "U_n.dbegin() = " << *(U_n.dbegin()) << std::endl;
		std::cout << "rhs.dbegin() = " << *(rhs.dbegin()) << std::endl;
		std::cout << "load_vector.dbegin() = " << *(load_vector.dbegin())
				  << std::endl;
	}
	l2h1error_ofs.close();
}

void BDF::bdf_algorithm(const int bdf_no){
	if( bdf_no < 1 || bdf_no > 6) 
		throw std::runtime_error("ERROR in BDF::bdf_algorithm():"
								 " Bad bdf_no = " + std::to_string(bdf_no) + "\n"
								 "Enter an integer between 1 and 6");

	double t_0 = 
		Dune::Fem::Parameter::getValue<double>("heat.starttime",0.0);
	double dT = 
		Dune::Fem::Parameter::getValue<double>("heat.timestep",0.1);
	double	t_end = 
		Dune::Fem::Parameter::getValue<double>("heat.endtime",0.6);
	const int time_step_no_max = (t_end - t_0)/dT + .1;
 
	// prepare grid from DGF file
	const std::string gridkey =
		Dune::Fem::IOInterface::defaultGridKey(GridType::dimension);
	const std::string gridfile =
		Dune::Fem::Parameter::getValue<std::string>(gridkey);
	if(Dune::Fem::MPIManager::rank() == 0)
		std::cout << "Loading macro grid: " << gridfile << std::endl;
	Dune::GridPtr<HostGridType> hostGrid {gridfile};
	hostGrid->loadBalance();

	// create grid
	DeformationCoordFunction deformation {};
	GridType grid(*hostGrid, deformation);
	GridPartType gridPart(grid);
	DiscreteFunctionSpaceType dfSpace(gridPart);

	L2NormType l2norm(gridPart);
	H1NormType h1norm(gridPart);
	std::ofstream l2h1error_ofs
	{Dune::Fem::Parameter::getValue<std::string>("fem.io.errorFile", 
												 "../output/l2h1error")};

	Dune::Fem::GridTimeProvider<GridType> timeProvider(t_0, grid);
	timeProvider.init(dT);     // Do your first action before you enter the for loop.

	deformation.set_time_provider(timeProvider);
	 // deformation is aware of timeProvider

	DiscreteFunctionType U_n("U_n",dfSpace), rhs("rhs", dfSpace), 
		load_vector("load_vector", dfSpace), exactSolution("exactSolution",dfSpace);
	DiscreteFunctionType ie_U_n("ie_U_n", dfSpace);
	DiscreteFunctionType ie_rhs("ie_rhs", dfSpace);

	std::vector<DiscreteFunctionType> prev_M_n_U_n;
	for(size_t i = 0; i < bdf_no; ++i)
		prev_M_n_U_n.push_back(
			DiscreteFunctionType {"U_n-" + std::to_string(i+1), dfSpace});
	
	
	IOTupleType ioTuple(&U_n);
	const int step = 0;	// is there any resonable choise, whitout adaptivity?
	DataOutputType 
		dataOutput(grid, ioTuple, 
				   DataOutputParameters(Dune::Fem::Parameter::getValue<std::string>
										("fem.io.outputName",
										 "../output/ALE_LiteratureExample-"),
										step) );
	// dataOutput.write(timeProvider);

	auto write_error = [&](std::ostream& os){
		os << std::defaultfloat << timeProvider.time() << ' ' 
		<< std::scientific 
		<< l2norm.distance(exactSolution, U_n) << ' '
		<< h1norm.distance(exactSolution, U_n) << std::endl;
	};

	InitialDataType initialData {timeProvider};
	RHSFunctionType f {timeProvider};
	NonlinearModel model {timeProvider};
	NonlinearOperator ellipticOp {U_n, model};
	const double solverEps =
		Dune::Fem::Parameter::getValue<double>( "heat.solvereps", 1e-8 );
	LinearInverseOperatorType solver(ellipticOp, solverEps, solverEps);

	// first time step
	int time_step_no = 0;
	// timeProvider.time() == t_0;

	InterpolationType::interpolateFunction(initialData, U_n);
	InterpolationType::interpolateFunction(initialData, exactSolution);
	InterpolationType::interpolateFunction(initialData, ie_U_n);

	// write_error(l2h1error_ofs);
 	// dataOutput.write(timeProvider);

	ellipticOp.get_xi(U_n);
	ellipticOp.mass_matrix(ie_U_n, ie_rhs);
	ellipticOp.mass_matrix(ie_U_n, prev_M_n_U_n.at(0) );
	// calculating: rhs = Mⁿ uⁿ  

	// debugging
	std::cout << std::scientific << "\n==========\n"
			  << "ie_U_n = " << *(ie_U_n.dbegin()) << '\t'
		      << "ie_rhs = " << *(ie_rhs.dbegin()) << '\n'
			  << "U_n = " << *(U_n.dbegin()) << '\t'
			  << "prev_M_n_U_n.at(0) = " << *(prev_M_n_U_n.at(0).dbegin()) << std::endl;


	for(timeProvider.next(dT); time_step_no < time_step_no_max; timeProvider.next(dT), ++time_step_no)
		// put your loop-action here, but not the last action
	{
		std::cout << std::defaultfloat << "==========\n"
				  << "In for loop: time_step_no = "<< time_step_no << ", "
				  << "time = " << timeProvider.time() << std::endl;
		
		// f.setTime(timeProvider.time());	// obselete
		// chance the time for the (long) RHS function f to tⁿ⁺¹
		// ('f' from -∆u + \div(v) u + \matdot{u} = f)

		// initialData.setTime(timeProvider.time());	// obselete

		assembleRHS(f, load_vector);
		std::cout << "assembleRHS(f, load_vector) = " << *(load_vector.dbegin()) << std::endl;

		// assemly stiffness/load vector; the vector is called 'rhs'; 
		// in the sense of above it is fⁿ⁺¹

		rhs.axpy(timeProvider.deltaT(), load_vector);
		// rhs += Δt * load_vector
		// Just calculated: rhs == MⁿUⁿ + Δt fⁿ⁺¹
		
		std::cout << "rhs.axpy = " << *(rhs.dbegin()) << std::endl;

		//InterpolationType::interpolateFunction( initialData, U_n );
		solver(rhs, U_n);
		// Solve:  (Mⁿ⁺¹ + Δt A(Uⁿ)) uⁿ⁺¹ = rhs,
		// CAPG remark: this is very very slow.

		std::cout << "solver(rhs, U_n) = " << *(U_n.dbegin()) << std::endl;

		// dataOutput.write(timeProvider);

		InterpolationType::interpolateFunction(initialData, exactSolution);		
		// write_error(l2h1error_ofs);
		
		ellipticOp.mass_matrix(U_n, rhs);
		std::cout << "mass_matrix(U_n, rhs) = " << *(rhs.dbegin()) << std::endl;

		// calculating: rhs = Mⁿ uⁿ  
		ellipticOp.get_xi(U_n);

		// debugging
		std::cout << std::scientific
				  << "U_n.dbegin() = " << *(U_n.dbegin()) << std::endl;
		std::cout << "rhs.dbegin() = " << *(rhs.dbegin()) << std::endl;
		std::cout << "load_vector.dbegin() = " << *(load_vector.dbegin())
				  << std::endl;
	}
	l2h1error_ofs.close();
}


// begin 2014nonlinear experiment
void BDF::nonlinear_algorithm(){
	double t_0 = 
		Dune::Fem::Parameter::getValue<double>("heat.starttime",0.0);
	double dT = 
		Dune::Fem::Parameter::getValue<double>("heat.timestep",0.1);
	double	t_end = 
		Dune::Fem::Parameter::getValue<double>("heat.endtime",0.6);
	const int itno = (t_end - t_0)/dT + .1;

	// setting up BDF coefficients
	const int bdf_no =
		Dune::Fem::Parameter::getValue<double>("heat.bdf", 1);
	const std::vector<double> bdf_alpha_coeff { NUMERIK::bdf_alpha_coeff(bdf_no) };
	const std::vector<double> bdf_gamma_coeff { NUMERIK::bdf_gamma_coeff(bdf_no) };
 
	// prepare grid from DGF file
	const std::string gridkey =
		Dune::Fem::IOInterface::defaultGridKey( GridType::dimension );
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

	L2NormType l2norm(gridPart);
	H1NormType h1norm(gridPart);
	std::ofstream l2h1error_ofs
	{Dune::Fem::Parameter::getValue<std::string>("fem.io.errorFile", 
												 "../output/l2h1error")};

	Dune::Fem::GridTimeProvider<GridType> timeProvider(t_0, grid);
	timeProvider.init(dT);     // Do your first action before you enter the for loop.

	deformation.set_time_provider(timeProvider);	
	 // deformation is aware of timeProvider

	// Setting up dof vectors
	DiscreteFunctionType U_n("U_n", dfSpace);	// Numerical solution
	DiscreteFunctionType rhs("rhs", dfSpace);	// Right hand side for the LES
	DiscreteFunctionType bdf_tmp("bdf_tmp", dfSpace);	
	 // Purpose of bdf_tmp:
	 // For the first k-time steps it is a tmp variable.
	 // Later it holds the sum from the previous BDF steps
	DiscreteFunctionType load_vector("load_vector", dfSpace);
	DiscreteFunctionType exact_solution("exact_solution", dfSpace);	
	DiscreteFunctionType err_vec("U_n_minus_exact_solution", dfSpace);
	 // Holds U_n - exact_solution

	auto calc_err_vec = [&]
	// helper function: caclulates err_vec - exact_solution
		{
			err_vec.assign(U_n);
			err_vec.axpy((-1.), exact_solution);
		};

	std::vector<DiscreteFunctionType> prev_steps_vec;
	for(int i = 0; i < bdf_no; ++i)
		prev_steps_vec.push_back( 
			DiscreteFunctionType {"U_n" + std::to_string(i+1), dfSpace} );

 	// debugging
	DiscreteFunctionType ie_U_n("ie_U_n", dfSpace);
	DiscreteFunctionType ie_rhs("ie_rhs", dfSpace);
	DiscreteFunctionType ie_load_vector("ie_load_vector", dfSpace);

	// Setting up PDE/ discrete operator
	InitialDataType initialData {timeProvider};
	InterpolationType::interpolateFunction( initialData, ie_U_n );
	
	RHSFunctionType f {timeProvider};
	NonlinearModel model {timeProvider};
	NonlinearOperator ellipticOp {ie_U_n, model};
	const double solverEps =
		Dune::Fem::Parameter::getValue<double>( "heat.solvereps", 1e-8 );
	LinearInverseOperatorType solver(ellipticOp, solverEps, solverEps);

	// Setting up paraview output
	IOTupleType ioTuple(&U_n);
	const int step = 0;	// is there any resonable choise, whitout adaptivity?
	DataOutputType 
		dataOutput(grid, ioTuple, 
				   DataOutputParameters(Dune::Fem::Parameter::getValue<std::string>
										("fem.io.outputName",
										 "../output/ALE_LiteratureExample-"),
										step) );
	// helper function
	auto write_error = [&](std::ostream& os){
		os << std::defaultfloat << timeProvider.time() << ' ' 
		<< std::scientific 
		<< l2norm.distance(exact_solution, U_n) << ' '
		<< h1norm.distance(exact_solution, U_n) << std::endl;
	};



	std::cout << '\n' << std::endl;
	// initialize prev_step_vec[i] dof
	int iter = (-1) * bdf_no;
	for(int i =0; 
		(i < prev_steps_vec.size()) && (iter < itno); 
		++i, ++iter, timeProvider.next(dT)){
		std::cout << std::defaultfloat << "==========\n"
				  << "In initial steps: iter = "<< iter << ", "
				  << "time = " << timeProvider.time() << std::endl;

		InterpolationType::interpolateFunction(initialData, bdf_tmp);
		ellipticOp.mass_matrix(bdf_tmp, prev_steps_vec.at(i));
		InterpolationType::interpolateFunction(initialData, exact_solution);
		U_n.assign(bdf_tmp);

		// debugging
		// if(i == 0){
		// 	InterpolationType::interpolateFunction(initialData, ie_U_n);
		// 	std::cout << "interpolate = " << *(ie_U_n.dbegin()) << std::endl;
		// 	ellipticOp.mass_matrix(ie_U_n, ie_rhs);
		// 	ellipticOp.get_xi(ie_U_n);
		// }
		// else{
		// 	assembleRHS(f, ie_load_vector);
		// 	ie_rhs.axpy(timeProvider.deltaT(), ie_load_vector);
		// 	solver(ie_rhs, ie_U_n);
		// 	ellipticOp.get_xi(ie_U_n);
		// }
	
		// calc_err_vec();
		// dataOutput.write(timeProvider);
		// write_error(l2h1error_ofs);

		std::cout << std::scientific
				  << "U_n.dbegin() = " << *(U_n.dbegin())  << '\t'
				  << "ie_U_n.dbegin() = " << *(ie_U_n.dbegin()) << '\n'
				  << "ie_rhs.dbegin() = " << *(ie_rhs.dbegin()) <<  std::endl;
	}

	ellipticOp.get_xi(U_n);
	ellipticOp.mass_matrix(U_n, rhs);

	// two helper functions

	auto bdf_rhs = [&] 
	// calculating: rhs = (-1) · (∑ᵢ₌₁ᵏ αᵢ₋₁ (Mu)ⁿ⁻ᵏ⁺ⁱ)
	//				bdf_tmp = ∑ᵢ₌₁ᵏ γᵢ₋₁ (Mu)ⁿ⁻ᵏ⁺ⁱ
	// NOTE for implicit euler: rhs = bdf_tmp = Mⁿ uⁿ
		{
			
			rhs.clear();	// set dof to zero
			bdf_tmp.clear();
			
			for(size_t i=0; i < prev_steps_vec.size(); ++i){
				rhs.axpy(bdf_alpha_coeff.at(i), prev_steps_vec.at(i));
				bdf_tmp.axpy(bdf_gamma_coeff.at(i), prev_steps_vec.at(i));
			}
			rhs *= (-1);
			ellipticOp.get_xi(bdf_tmp);
		};

	auto bdf_cycle = [&]
	// use this after calculating the numerical solution
		{
			// std::cout << "bdf_cycle at timeProvider.time() = "
			// << timeProvider.time() << std::endl;
			for(size_t i=0; i < prev_steps_vec.size() - 1; ++i)
				prev_steps_vec.at(i).assign(prev_steps_vec.at(i+1));
			ellipticOp.mass_matrix(U_n, prev_steps_vec.back());
		};
	

	for(; 
		iter < itno; 
		timeProvider.next(dT), ++iter)
		// put your loop-action here, but not the last action
	{
		std::cout << std::defaultfloat << "==========\n"
				  << "In for loop: iter = "<< iter << ", "
				  << "time = " << timeProvider.time() << std::endl;
		
		bdf_rhs();		
		assembleRHS(f, load_vector);
		// assemly stiffness/load vector; the vector is called 'load_vector'; 
		// in the sense of above it is fⁿ⁺¹

		rhs.axpy(timeProvider.deltaT(), load_vector);
		// rhs += Δt * load_vector
		// Just calculated: rhs == MⁿUⁿ + Δt fⁿ⁺¹

		rhs /= bdf_alpha_coeff.back();

		//InterpolationType::interpolateFunction( initialData, U_n );
		solver(rhs, U_n);
		// Solve:  (Mⁿ⁺¹ + Δt A(Uⁿ)) uⁿ⁺¹ = rhs,
		// CAPG remark: this is very very slow.

		bdf_cycle();
		InterpolationType::interpolateFunction(initialData, exact_solution);		
		calc_err_vec();

		// dataOutput.write(timeProvider);
		// if(timeProvider.time() > t_end - 3*dT)
			// write_error(l2h1error_ofs);	// write_error(std::cout);	

		assembleRHS(f, ie_load_vector);
		std::cout << "assembleRHS(f, ie_load_vector) = " << *(ie_load_vector.dbegin()) << std::endl;
		ie_rhs.axpy(timeProvider.deltaT(), ie_load_vector);
		std::cout << "ie_rhs.axpy = " << *(ie_rhs.dbegin()) << std::endl;
		solver(ie_rhs, ie_U_n);
		std::cout << "solver(ie_rhs, ie_U_n) = " << *(ie_U_n.dbegin()) << std::endl;
		ellipticOp.mass_matrix(ie_U_n, ie_rhs);
		std::cout << "mass_matrix(ie_U_n, ie_rhs) = " << *(ie_rhs.dbegin()) << std::endl;
		ellipticOp.get_xi(ie_U_n);

		// debugging
		std::cout << std::scientific
				  << "U_n.dbegin() = " << *(U_n.dbegin())  << "\t\t"
				  << "ie_U_n.dbegin() = " << *(ie_U_n.dbegin()) << std::endl;
		std::cout << "rhs.dbegin() = " << *(rhs.dbegin()) << "\t\t"
				  << "ie_rhs.dbegin() = " << *(ie_rhs.dbegin()) << std::endl;
		std::cout << "load_vector.dbegin() = " << *(load_vector.dbegin()) << '\t'
				  << "ie_load_vector.dbegin() = " << *(ie_load_vector.dbegin()) 
				  << std::endl;
		std::cout << "bdf_tmp.dbegin() = " << *(bdf_tmp.dbegin()) << std::endl;

		for(size_t i=0; i < prev_steps_vec.size(); ++i){
			std::cout << "prev_steps_vec.at(" << i << ").dbegin() = "
					  << *(prev_steps_vec.at(i).dbegin()) << std::endl;
		}

	}
	l2h1error_ofs.close();
}
// end 2014nonlinear experiment


// // begin paperALE experiment
// 
// // This has to be modified for the ODE solver
// struct EigenFullPivLuSolv{
// 	Eigen::Vector3d linear_solve(const Eigen::Matrix3d& A, 
// 								 const Eigen::Vector3d& b) const {
// 		return A.fullPivLu().solve(b);  
// 	}
// };
// struct EigenNorm{
// 	double norm(const Eigen::Vector3d& v) const {
// 		return v.norm();
// 	}
// };
// struct SurfaceVectorfield{
// 	void set_time(const double time) { t = time; }
// 	double get_time() { return t; }
// 	Eigen::Vector3d evaluate(const Eigen::Vector3d vec) const {
// 		double x = vec(0), y = vec(1), z = vec(2);
// 		return Eigen::Vector3d { 
// 			16.*M_PI*pow(x,3)*cos(2*M_PI*t)/( (pow(y,2) + pow(z,2) + 16.*pow(x,2)/pow(sin(2.*M_PI*t) + 4.,2) )*pow(sin(2.*M_PI*t) + 4.,3) ), 
// 				4.*M_PI*pow(x,2)*y*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16.*pow(x,2)/pow(sin(2.*M_PI*t) + 4.,2) )*pow(sin(2.*M_PI*t) + 4.,2)),
// 				4.*M_PI*pow(x,2)*z*cos(2.*M_PI*t)/((pow(y,2) + pow(z,2) + 16.*pow(x,2)/pow(sin(2*M_PI*t) + 4.,2) )*pow(sin(2.*M_PI*t) + 4,2) ) };
// 	}
// 	Eigen::Matrix3d jacobian(const Eigen::Vector3d vec) const {
// 		double x = vec(0), y = vec(1), z = vec(2);
// 		Eigen::Matrix3d jac_f_t_vec;
// 		jac_f_t_vec << 
// 			// first row
// 			48*M_PI*pow(x,2)*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,3)) - 512*M_PI*pow(x,4)*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,5)), 
// 			-32*M_PI*pow(x,3)*y*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,3)), 
// 			-32*M_PI*pow(x,3)*z*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,3)),
// 			// second row
// 			8*M_PI*x*y*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)) - 128*M_PI*pow(x,3)*y*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,4)), 
// 			-8*M_PI*pow(x,2)*pow(y,2)*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)) + 4*M_PI*pow(x,2)*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)), 
// 			-8*M_PI*pow(x,2)*y*z*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)), 
// 			// third row
// 			8*M_PI*x*z*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)) - 128*M_PI*pow(x,3)*z*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,4)), 
// 			-8*M_PI*pow(x,2)*y*z*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)), 
// 			-8*M_PI*pow(x,2)*pow(z,2)*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)) + 4*M_PI*pow(x,2)*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2));
// 		return jac_f_t_vec;
// 	}
// private:
// 	double t;
// };
// struct EigenIdMatrix3d{
// 	Eigen::Matrix3d identity_matrix() const { return Eigen::Matrix3d::Identity(); }
// };
// 
// 
// typedef NUMERIK::ImplicitEuler<EigenFullPivLuSolv, EigenNorm, SurfaceVectorfield, 
// 							   EigenIdMatrix3d, Eigen::Matrix3d, Eigen::Vector3d>						
// EvoMapIntegrator;
// 
// void BDF::algorithm()
// {
// 	/**********************************************************************
// 	 * Pseudocode:
// 	 * u⁰ = interpol(solution) >> plot
// 	 * container = M⁰u⁰ >> save
// 	 * for-loop n = 0 ↦ N with N·Δt + t_begin = t_end
// 	 *  t += ∆τ (set time)
// 	 *  tmp = container == Mⁿuⁿ
// 	 *  rhs = fⁿ⁺¹ == load Vector
// 	 *  tmp += Δτ · rhs (solution == Mⁿuⁿ + Δτ·rhs)
// 	 *  uⁿ⁺¹ == solution = (Mⁿ⁺¹ + Δτ Aⁿ⁺¹)⁻¹ tmp >> save
// 	 *  container = Mⁿ⁺¹ uⁿ⁺¹
// 	 * end-for-loop
// 	 *  
// 	 **********************************************************************/
//  
// 	// create grid from DGF file
// 	const std::string gridkey =
// 		Dune::Fem::IOInterface::defaultGridKey( GridType::dimension );
// 	const std::string gridfile =
// 		Dune::Fem::Parameter::getValue< std::string >( gridkey );
// 	if( Dune::Fem::MPIManager::rank() == 0 )
// 		std::cout << "Loading macro grid: " << gridfile << std::endl;
// 
// 	Dune::GridPtr< HostGridType > hostGrid {gridfile};
// 	hostGrid ->loadBalance();
// 
// 	BDF::EvoMapType evoMap {gridfile};
// 	DeformationCoordFunction deformation {evoMap};
// 	 // constructs an unsorted map with domain and range equal the vertexes from the
// 	 // initial grid
//   
// 	GridType grid(*hostGrid, deformation);
// 	GridPartType gridPart(grid);
// 	DiscreteFunctionSpaceType dfSpace(gridPart);
// 
// 	L2NormType l2norm(gridPart);
// 	H1NormType h1norm(gridPart);
// 	std::ofstream l2h1error_ofs
// 	{Dune::Fem::Parameter::getValue<std::string>("fem.io.errorFile", 
// 												 "../output/l2h1error")};
// 	//std::ios_base::app};
// 	//l2h1error_ofs << std::scientific;
// 	//l2h1error_ofs << "L2Error\tH1Error" << std::scientific << std::endl;
// 	//std::ostringstream l2h1error_helper_oss;
// 	//l2h1error_helper_oss << "time\ttimestep" << std::endl;
// 
// 	double t_0 = 
// 		Dune::Fem::Parameter::getValue<double>("heat.starttime",0.0);
// 	double dT = 
// 		Dune::Fem::Parameter::getValue<double>("heat.timestep",0.1);
// 	double	t_end = 
// 		Dune::Fem::Parameter::getValue<double>("heat.endtime",0.6);
// 	std::cout << (t_end - t_0)/dT + .1 << std::endl;
//     const int itno = (t_end - t_0)/dT + .1;
// 
// 	Dune::Fem::GridTimeProvider< GridType > timeProvider(t_0, grid);
// 	timeProvider.init(dT);     // Do your first action before you enter the for loop.
// 	// std::cout << timeProvider.time() << '\t' 
// 	// 		  << timeProvider.deltaT() << "\tstarting action."
// 	// 		  << std::endl;
// 	// l2h1error_helper_oss << timeProvider.time() << '\t' 
// 	//  					 << timeProvider.deltaT() << "\tstarting action."
// 	//  					 << std::endl;
// 
// 	DiscreteFunctionType solutionContainer( "U_nContainer", dfSpace ),
// 		U_n("U_n",dfSpace), rhs( "rhs", dfSpace ), tmp( "tmp", dfSpace ),
// 		exactSolution("exactSolution",dfSpace);
// 	InitialDataType initialData;
// 	InterpolationType::interpolateFunction( initialData, U_n );
// 	InterpolationType::interpolateFunction( initialData, exactSolution );
// 	// solutionContainer is used for higher order BDF methods,
// 	// U_n == uⁿ, rhs == tmp1, tmp == tmp2
// 	// (create discrete functions for intermediate functionals)
// 
// 	RHSFunctionType f;
// 	// This expression is very long; it's literally the same f as in L(u) = f;
// 	// together with the function/method in rhs.hh, 
// 
// 	IOTupleType ioTuple(&U_n);
// 	const int step = 0;	// is there any resonable choise, whitout adaptivity?
// 	DataOutputType 
// 		dataOutput(grid, ioTuple, 
// 				   DataOutputParameters(Dune::Fem::Parameter::getValue<std::string>
// 										("fem.io.outputName",
// 										 "../output/ALE_LiteratureExample-"),
// 										step) );
// 
// 	auto write_error = [&](std::ostream& os){
// 		os << std::defaultfloat << timeProvider.time() << ' '
// 		   << std::scientific 
// 		   << l2norm.distance(exactSolution, U_n) << ' '
// 		   << h1norm.distance(exactSolution, U_n) << std::endl;
// 	};
// 	write_error(l2h1error_ofs);
// 	
// 	EllipticOperatorType
// 		ellipticOp( HeatModelType {timeProvider} );
// 	
// 	const double solverEps =
// 		Dune::Fem::Parameter::getValue< double >( "heat.solvereps", 1e-8 );
// 	LinearInverseOperatorType solver( ellipticOp, solverEps, solverEps );
// 	// Initializing CG-Solver
// 	// CAPG: it's very slow. Can't this be speed up somehow?
// 
// 	dataOutput.write( timeProvider );
// 
// 	ellipticOp.mass_matrix(U_n, rhs);
// 	// calculating: rhs = Mⁿ uⁿ  (uⁿ == solutionContainer)
// 
// 	// std::vector<double> timecontainer, l2errorcontainer, h1errorcontainer;
// 	int iter = 0;
// 	for(timeProvider.next(dT); iter < itno; timeProvider.next(dT), ++iter)
// 		// put your loop-action here, but not the last action
// 	{
// 		// for debugging reasons
// 		// std::cout << timeProvider.time() << '\t' << timeProvider.deltaT() 
// 		// 		  << "\tfor loop action." << std::endl;
// 		// l2h1error_helper_oss << timeProvider.time() << '\t' << timeProvider.deltaT() 
// 		//  					 << "\tfor loop action." << std::endl;
// 
// 		f.setTime(timeProvider.time());
// 		// chance the time for the (long) RHS function f to tⁿ⁺¹
// 		// ('f' from -∆u + \div(v) u + \matdot{u} = f)
// 
// 		initialData.setTime(timeProvider.time());
// 
// 		std::array<double,2> time_array { 
// 			{timeProvider.time()- timeProvider.deltaT(), timeProvider.time()} 
// 		};
// 		EvoMapIntegrator impl_euler { NUMERIK::NewtonParameters<> {1e-10} };
// 		evoMap.evolve(time_array.begin(), time_array.end(), impl_euler);
// 
// 		assembleRHS( f, tmp );
// 		// assemly stiffness/load vector; the vector is called 'tmp'; 
// 		// in the sense of above it is fⁿ⁺¹
// 
// 		rhs.axpy( timeProvider.deltaT(), tmp );
// 		// rhs += Δt * tmp (i.e. rhs += \Delta t * tmp)
// 		// Just calculated: ( uⁿ + ∆t EllipticOperator(uⁿ)) + ∆t fⁿ
// 
// 		//InterpolationType::interpolateFunction( initialData, U_n );
// 		solver( rhs, U_n );
// 		// Solve:  (id + \Delta t EllipticOperator) u_{n+1} = rhs,
// 		// where rhs = ( u_n + \Delta t EllipticOperator(u_n)) + \Delta t f_n
// 		// CAPG remark: this is very very slow.
// 
// 		dataOutput.write( timeProvider );
// 
// 		InterpolationType::interpolateFunction( initialData, exactSolution );		
// 		write_error(l2h1error_ofs);
// 		
// 		ellipticOp.mass_matrix(U_n, rhs);
// 		// calculating: rhs = Mⁿ uⁿ  (uⁿ == solutionContainer)
// 
// 		// l2errorcontainer.push_back(l2norm.distance(exactSolution, U_n));
// 		// h1errorcontainer.push_back(h1norm.distance(exactSolution, U_n));
// 		// timecontainer.push_back(timeProvider.time());
// 	}
// 	// l2h1error_ofs << "#\n" << l2h1error_helper_oss.str() << '#' << std::endl;
// 
// 	// l2h1error_ofs << timeProvider.deltaT()   << ' '
// 	// 			  << timecontainer.back()    << ' '
// 	// 			  << l2errorcontainer.back() << ' '	
// 	// 			  << h1errorcontainer.back() << std::endl;
// 	l2h1error_ofs.close();
// }
// // end paperALE experiment

#endif // #ifndef DUNE_HEAT_ALGORITHM_HPP
