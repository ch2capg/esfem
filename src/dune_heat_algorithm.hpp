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

#include <algorithm>

// local includes
#include "dune_typedef_heat.hpp"
// #include "/Users/christianpower/cpp/ODE_Solver/implicit_euler.h"
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
	/**********************************************************************
	 * Pseudocode:
	 * u⁰ = interpol(solution) >> plot
	 * container = M⁰u⁰ >> save
	 * for-loop n = 0 ↦ N with N·Δt + t_begin = t_end
	 *  t += ∆τ (set time)
	 *  tmp = container == Mⁿuⁿ
	 *  rhs = fⁿ⁺¹ == load Vector
	 *  tmp += Δτ · rhs (solution == Mⁿuⁿ + Δτ·rhs)
	 *  uⁿ⁺¹ == solution = (Mⁿ⁺¹ + Δτ Aⁿ⁺¹)⁻¹ tmp >> save
	 *  container = Mⁿ⁺¹ uⁿ⁺¹
	 * end-for-loop
	 *  
	 **********************************************************************/

	// get time from parameter file
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
												 "../output/l2h1error"),
			std::ios_base::app};

	Dune::Fem::GridTimeProvider<GridType> timeProvider(t_0, grid);
	timeProvider.init(dT);     // Do your first action before you enter the for loop.
	deformation.set_time_provider(timeProvider);

	DiscreteFunctionType U_n {"U_n", dfSpace};
	DiscreteFunctionType rhs {"rhs", dfSpace}; 
	DiscreteFunctionType load_vector {"load_vector", dfSpace};
	DiscreteFunctionType exactSolution {"exactSolution", dfSpace};
	DiscreteFunctionType xi {"xi", dfSpace};

	InitialDataType initialData {timeProvider};

	// to get the degree of freedom vector
	// for(auto it = U_n.dbegin(); it != U_n.dend(); ++it)
	//	std::cout << *it << std::endl;
	
	RHSFunctionType f {timeProvider};
	// This expression is very long; it's literally the same f as in L(u) = f;
	// together with the function/method in rhs.hh, 

	IOTupleType ioTuple(&U_n);
	const int step = 0;	// is there any resonable choice, whitout adaptivity?
	DataOutputType 
		dataOutput(grid, ioTuple, 
				   DataOutputParameters(Dune::Fem::Parameter::getValue<std::string>
										("fem.io.outputName",
										 "../output/ALE_LiteratureExample-"),
										step) );
	// Usage: dataOutput.write(timeProvider);

	auto write_error = [&](std::ostream& os){
		os << std::defaultfloat << timeProvider.deltaT() << ' ' 
		<< std::scientific 
		<< l2norm.distance(exactSolution, U_n) << ' '
		<< h1norm.distance(exactSolution, U_n) << std::endl;
	};

	// Usage: write_error(l2h1error_ofs);
	
	NonlinearModel model {timeProvider};
	NonlinearOperator ellipticOp {xi, model};
	const double solverEps =
		Dune::Fem::Parameter::getValue<double>("heat.solvereps", 1e-8);
	LinearInverseOperatorType solver(ellipticOp, solverEps, solverEps);
	// Initializing CG-Solver
	// CAPG: it's very slow. Can't this be speed up somehow?

	InterpolationType::interpolateFunction(initialData, U_n);
	InterpolationType::interpolateFunction(initialData, exactSolution);

 	// dataOutput.write(timeProvider);

	// debugging
	// std::cout << "\n==========\n"
	// 		  << std::defaultfloat << "time = " << timeProvider.time() 
	// 		  << std::scientific << std::endl;
	// std::cout << "U_n = " << *(U_n.dbegin()) << std::endl;

	ellipticOp.mass_matrix(U_n, rhs);
	xi.assign(U_n);
	// calculating: rhs = Mⁿ uⁿ  
	// std::cout << "ellipticOp.mass_matrix(U_n, rhs) = " << *rhs.dbegin() << std::endl;


	int iter = 0;
	for(timeProvider.next(dT); iter < itno; timeProvider.next(dT), ++iter)
		// put your loop-action here, but not the last action
	{
		// std::cout << "----------\n"
		// 		  << std::defaultfloat<< "time = " << timeProvider.time()
		// 		  << std::scientific << std::endl;

		// f.setTime(timeProvider.time());	// obselete
		// chance the time for the (long) RHS function f to tⁿ⁺¹
		// ('f' from -∆u + \div(v) u + \matdot{u} = f)

		// initialData.setTime(timeProvider.time());	// obselete

		assembleRHS(f, load_vector);
		// std::cout << "assembleRHS(f, load_vector) = " << *load_vector.dbegin() 
		// 		  << std::endl;
		// assemly stiffness/load vector; the vector is called 'rhs'; 
		// in the sense of above it is fⁿ⁺¹

		rhs.axpy(timeProvider.deltaT(), load_vector);
		// std::cout << "rhs.axpy(timeProvider.deltaT(), load_vector) = "
		// 		  << *rhs.dbegin() << std::endl;
		// rhs += Δt * load_vector
		// Just calculated: rhs == MⁿUⁿ + Δt fⁿ⁺¹

		// InterpolationType::interpolateFunction(initialData, exactSolution);		
		// ellipticOp.get_xi(exactSolution);
		// std::cout << "ellipticOp.get_xi(exactSolution) = " << *exactSolution.dbegin()
		//			 << std::endl;

		//InterpolationType::interpolateFunction( initialData, U_n );
		solver(rhs, U_n);
		// std::cout << "solver(rhs, U_n) = " << *U_n.dbegin() << std::endl;
		// Solve:  (Mⁿ⁺¹ + Δt A(Uⁿ)) uⁿ⁺¹ = rhs,
		// CAPG remark: this is very very slow.

		// dataOutput.write(timeProvider);

		if(iter == itno - 1){
			std::cout << "Time = " << timeProvider.time() << std::endl;
			InterpolationType::interpolateFunction(initialData, exactSolution);		
			// write_error(l2h1error_ofs);
			write_error(std::cout);
		}
		
		ellipticOp.mass_matrix(U_n, rhs);
		// calculating: rhs = Mⁿ uⁿ  
		// std::cout << "ellipticOp.mass_matrix(U_n, rhs) = " << *rhs.dbegin() << std::endl;

		xi.assign(U_n);

		// debugging
		// std::cout << "U_n = " << *(U_n.dbegin()) << std::endl;

	}
	l2h1error_ofs.close();
	std::cout << std::endl;
}

void BDF::bdf2_heat_algorithm(){
	// get time from parameter file
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
												 "../output/l2h1error"),
			std::ios_base::app};

	Dune::Fem::GridTimeProvider<GridType> timeProvider(t_0, grid);
	timeProvider.init(dT);     // Do your first action before you enter the for loop.
	deformation.set_time_provider(timeProvider);

	DiscreteFunctionType U_n {"U_n",dfSpace};
	DiscreteFunctionType rhs {"rhs", dfSpace};
	DiscreteFunctionType load_vector {"load_vector", dfSpace};
	DiscreteFunctionType exactSolution {"exactSolution",dfSpace};
	DiscreteFunctionType U_nm1 {"U_nm1", dfSpace};
	DiscreteFunctionType M_U_nm1 {"M_U_nm1", dfSpace};
	DiscreteFunctionType U_nm2 {"U_nm2", dfSpace};
	DiscreteFunctionType M_U_nm2 {"M_U_nm2", dfSpace};
	DiscreteFunctionType xi {"xi", dfSpace};
	
	InitialDataType initialData {timeProvider};
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
		os << std::defaultfloat << timeProvider.deltaT() << ' ' 
		<< std::scientific 
		<< l2norm.distance(exactSolution, U_n) << ' '
		<< h1norm.distance(exactSolution, U_n) << std::endl;
	};
	// write_error(l2h1error_ofs);
	
	NonlinearModel model {timeProvider};
	NonlinearOperator ellipticOp {xi, model, 1.5};
	const double solverEps =
		Dune::Fem::Parameter::getValue<double>( "heat.solvereps", 1e-8 );
	LinearInverseOperatorType solver(ellipticOp, solverEps, solverEps);

 	// dataOutput.write(timeProvider);

	// first step
	int iter = 0;

	// debugging
	// std::cout << "\n==========\n"
	// 		  << std::defaultfloat << "time = " << timeProvider.time() 
	// 		  << std::scientific << std::endl;

	InterpolationType::interpolateFunction(initialData, U_n);
	// std::cout << "U_n = " << *(U_n.dbegin()) << std::endl;

	U_nm2.assign(U_n);
	// std::cout << "U_nm2.assign(U_n) = " << *U_nm2.dbegin() << std::endl;
	ellipticOp.mass_matrix(U_nm2, M_U_nm2);
	// std::cout << "ellipticOp.mass_matrix(U_nm2, M_U_nm2) = " 
	// 		  << *M_U_nm2.dbegin() << std::endl;

	// if(iter == itno){
	// 	std::cout << "Time = " << timeProvider.time() << std::endl;
	// 	InterpolationType::interpolateFunction(initialData, exactSolution);		
	// 	write_error(l2h1error_ofs);
	// }

	// second step
	++iter;		// iter == 1
	timeProvider.next(dT);

	std::cout << "----------\n"
			  << std::defaultfloat << "time = " << timeProvider.time()
			  << std::scientific << std::endl;

	InterpolationType::interpolateFunction(initialData, U_n);
	// std::cout << "U_n = " << *(U_n.dbegin()) << std::endl;
	InterpolationType::interpolateFunction(initialData, exactSolution);	

	U_nm1.assign(U_n);
	// std::cout << "U_nm1.assign(U_n) = " << *U_nm1.dbegin() << std::endl;
	ellipticOp.mass_matrix(U_nm1, M_U_nm1);
	// std::cout << "ellipticOp.mass_matrix(U_nm1, M_U_nm1) = " << *M_U_nm1.dbegin()
	// 		  << std::endl;

	// if(iter == itno){
	// 	std::cout << "Time = " << timeProvider.time() << std::endl;
	// 	InterpolationType::interpolateFunction(initialData, exactSolution);	
	// 	write_error(l2h1error_ofs);
	// }

	++iter;
	for(timeProvider.next(dT); iter <= itno; timeProvider.next(dT), ++iter)
		// put your loop-action here, but not the last action
	{
		std::cout << "----------\n"
				  << std::defaultfloat<< "time = " << timeProvider.time()
				  << std::scientific << std::endl;

		// prepare 'rhs'
		rhs.clear();	// rhs.clear() == rhs = 0
		// std::cout << "rhs.clear() = " << *rhs.dbegin() << std::endl;
		rhs.axpy((-1) * .5, M_U_nm2);		// a.axpy(α,v) == a += α·v
		// std::cout << "rhs.axpy((-1) * .5, M_U_nm2) = " << *rhs.dbegin() << std::endl;
		rhs.axpy((-1) * (-2.), M_U_nm1);
		// std::cout << "rhs.axpy((-1) * (-2.), M_U_nm1) = " 
		// 		  << *rhs.dbegin() << std::endl;
		assembleRHS(f, load_vector);
		// std::cout << "assembleRHS(f, load_vector) = " << *load_vector.dbegin() 
		// 		  << std::endl;
		rhs.axpy(timeProvider.deltaT(), load_vector);
		std::cout << "rhs.axpy(timeProvider.deltaT(), load_vector) = " 
				  << *rhs.dbegin() << std::endl;

		// prepare 'solver'
		// InterpolationType::interpolateFunction(initialData, exactSolution);
		// ellipticOp.get_xi(exactSolution);
		// std::cout << "ellipticOp.get_xi(exactSolution) = " << *exactSolution.dbegin()
		// 		  << std::endl;
		xi.clear();
		xi.axpy(-1.,U_nm2);
		xi.axpy(2.,U_nm1);
		std::cout << "xi = " << *xi.dbegin() << std::endl;

		solver(rhs, U_n);	// Solve:  (Mⁿ⁺¹ + Δt A(Uⁿ)) uⁿ⁺¹ = rhs,
		std::cout << "solver(rhs, U_n) = " << *U_n.dbegin() << std::endl;

		// write output
		// dataOutput.write(timeProvider);
		if(iter == itno){
			std::cout << std::defaultfloat
					  << "Time = " << timeProvider.time() 
					  << std::scientific << std::endl;
			InterpolationType::interpolateFunction(initialData, exactSolution);		
			// write_error(l2h1error_ofs);
			write_error(std::cout);
		}
		
		// cycle bdf values
		U_nm2.assign(U_nm1);
		// std::cout << "U_nm2.assign(U_nm1) = " << *U_nm2.dbegin() << std::endl;
		U_nm1.assign(U_n);
		// std::cout << "U_nm1.assign(U_n) = " << *U_nm1.dbegin() << std::endl;
		M_U_nm2.assign(M_U_nm1);
		// std::cout << "M_U_nm2.assign(M_U_nm1) = " << *M_U_nm2.dbegin() << std::endl;
		ellipticOp.mass_matrix(U_n, M_U_nm1);
		// std::cout << "ellipticOp.mass_matrix(U_n, M_U_nm1) = " 
		// 		  << *M_U_nm1.dbegin() << std::endl;

		// debugging
		// std::cout << "U_n = " << *(U_n.dbegin()) << std::endl;

	}
	l2h1error_ofs.close();
	std::cout << std::endl;
}


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
												 "../output/l2h1error"),
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

	auto write_error = [&](std::ostream& os)
	// helper function for 'l2h1error_ofs'
		{
			os << std::defaultfloat << timeProvider.deltaT() << ' ' 
			   << std::scientific 
			   << l2norm.distance(exact_solution, U_np1) << ' '
			   << h1norm.distance(exact_solution, U_np1) << std::endl;
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
