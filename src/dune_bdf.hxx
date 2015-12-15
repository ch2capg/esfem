/*********************************************************************
 *  dune_bdf.hpp, dune_bdf.hxx                                       *
 *                                                                   *
 *  Defines auxiliar classes like a string vector with the correct   *
 *  file names.                                                      *
 *                                                                   *
 *  Revision history:                                                *
 *  none (this is experimental and error phrone)                     *
 *                                                                   *
 *                                                                   *
 *  Created by Christian Power on 19.06.14.                          *
 *  Copyright (c) 2014 Christian Power. All rights reserved.         *
 *                                                                   *
 *********************************************************************/

#ifndef DUNE_BDF_HXX
#define DUNE_BDF_HXX

// all implementation details

// helper function

// Usually in only template definitions are put into a *.hxx file
// (CAPG convention), but the dune makefiles are very complicated,
// so I put also the function definition in there.

// starting the interface for EvoMapType::evolve()
struct SolverType : 
	NUMERIK::NewtonLinearSolverType< Eigen::Vector3d, Eigen::Matrix3d, SolverType>{
	Eigen::Vector3d operator() (const Eigen::Matrix3d& A, const Eigen::Vector3d& b) const {
		return A.fullPivLu().solve( (-1) * b);
	}
};

struct CAPGNormType : NUMERIK::NormFunctorType<double, Eigen::Vector3d, CAPGNormType>{
	double operator() (Eigen::Vector3d v) const{
		return v.norm();
	}
};

struct SurfaceVectorfieldType : 
	NUMERIK::C1SpaceTimeFunctionType
	<Eigen::Vector3d, Eigen::Matrix3d, SurfaceVectorfieldType>{
	typedef C1SpaceTimeFunctionType
	<Eigen::Vector3d, Eigen::Matrix3d, SurfaceVectorfieldType> BaseType;
	using BaseType::get_time;
	using BaseType::get_position;
	using BaseType::get_dT;

	SurfaceVectorfieldType(const double starting_time, const Eigen::Vector3d& starting_pos) : 
		BaseType(starting_time, starting_pos) {}
	Eigen::Vector3d evaluate(Eigen::Vector3d vec) const{
		double x = vec(0), y = vec(1), z = vec(2), t = get_time();
		Eigen::Vector3d f_t_vec( 
			16.*M_PI*pow(x,3)*cos(2*M_PI*t)/( (pow(y,2) + pow(z,2) + 16.*pow(x,2)/pow(sin(2.*M_PI*t) + 4.,2) )*pow(sin(2.*M_PI*t) + 4.,3) ),
			4.*M_PI*pow(x,2)*y*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16.*pow(x,2)/pow(sin(2.*M_PI*t) + 4.,2) )*pow(sin(2.*M_PI*t) + 4.,2)),
			4.*M_PI*pow(x,2)*z*cos(2.*M_PI*t)/((pow(y,2) + pow(z,2) + 16.*pow(x,2)/pow(sin(2*M_PI*t) + 4.,2) )*pow(sin(2.*M_PI*t) + 4,2) ) );
		return  vec - get_dT() * f_t_vec - get_position();
		// x_next - dT * f(t_next, x_next) - x_n
	}
	Eigen::Matrix3d jacobian(Eigen::Vector3d vec) const{
		double x = vec(0), y = vec(1), z = vec(2), t = get_time();
		Eigen::Matrix3d jac_f_t_vec;
		jac_f_t_vec << 
			// first row
			48*M_PI*pow(x,2)*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,3)) - 512*M_PI*pow(x,4)*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,5)), 
			-32*M_PI*pow(x,3)*y*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,3)), 
			-32*M_PI*pow(x,3)*z*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,3)),
			// second row
			8*M_PI*x*y*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)) - 128*M_PI*pow(x,3)*y*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,4)), 
			-8*M_PI*pow(x,2)*pow(y,2)*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)) + 4*M_PI*pow(x,2)*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)), 
			-8*M_PI*pow(x,2)*y*z*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)), 
			// third row
			8*M_PI*x*z*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)) - 128*M_PI*pow(x,3)*z*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,4)), 
			-8*M_PI*pow(x,2)*y*z*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)), 
			-8*M_PI*pow(x,2)*pow(z,2)*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)) + 4*M_PI*pow(x,2)*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2));
		return Eigen::Matrix3d::Identity() - get_dT() * jac_f_t_vec;
		// id_3x3 - dT * jac(f)(t_next, x_next)
	}
};

struct ALESurfaceVectorfieldType : 
	NUMERIK::C1SpaceTimeFunctionType
	<Eigen::Vector3d, Eigen::Matrix3d, ALESurfaceVectorfieldType>{
	typedef C1SpaceTimeFunctionType
	<Eigen::Vector3d, Eigen::Matrix3d, ALESurfaceVectorfieldType> BaseType;
	using BaseType::get_time;
	using BaseType::get_position;
	using BaseType::get_dT;

	ALESurfaceVectorfieldType(const double starting_time, 
							  const Eigen::Vector3d& starting_pos) : 
		BaseType(starting_time, starting_pos) {}

	Eigen::Vector3d evaluate(Eigen::Vector3d vec) const{
		double x = vec(0), y = vec(1), z = vec(2), t = get_time();
		Eigen::Vector3d f_t_vec( 
			.5 * M_PI * cos(2. * M_PI * t) / ( 2. * sqrt( 1. + .25 * sin(2. * M_PI * t) ) ) *x, 0, 0 );
		return  vec - get_dT() * f_t_vec - get_position();
		// x_next - dT * f(t_next, x_next) - x_n
	}
	Eigen::Matrix3d jacobian(Eigen::Vector3d vec) const{
		double x = vec(0), y = vec(1), z = vec(2), t = get_time();
		Eigen::Matrix3d jac_f_t_vec;
		jac_f_t_vec << 
			// first row
			.5 * M_PI * cos(2. * M_PI * t) / ( 2. * sqrt( 1. + .25 * sin(2. * M_PI * t) ) ) , 0. , 0.,
			// second row
			0. , 0. , 0.,
			// third row
			0. , 0. , 0.;
		return Eigen::Matrix3d::Identity() - get_dT() * jac_f_t_vec;
		// id_3x3 - dT * jac(f)(t_next, x_next)
	}
};

struct StationaryType : 
	NUMERIK::C1SpaceTimeFunctionType <Eigen::Vector3d, Eigen::Matrix3d, StationaryType>{
	
	typedef C1SpaceTimeFunctionType
	<Eigen::Vector3d, Eigen::Matrix3d, StationaryType> BaseType;
	using BaseType::get_time;
	using BaseType::get_position;
	using BaseType::get_dT;
	
	StationaryType(const double starting_time, 
				   const Eigen::Vector3d& starting_pos) : 
		BaseType(starting_time, starting_pos) {}

	Eigen::Vector3d evaluate(Eigen::Vector3d vec) const{
		return  vec - get_position();
	}
	Eigen::Matrix3d jacobian(Eigen::Vector3d vec) const{
		return Eigen::Matrix3d::Identity();
	}
};

struct RotatingType : 
	NUMERIK::C1SpaceTimeFunctionType <Eigen::Vector3d, Eigen::Matrix3d, RotatingType>{ 
	
	typedef C1SpaceTimeFunctionType
	<Eigen::Vector3d, Eigen::Matrix3d, RotatingType> BaseType;
	using BaseType::get_position;
	
	RotatingType(const double starting_time, 
				   const Eigen::Vector3d& starting_pos) : 
		BaseType(starting_time, starting_pos) {}

	Eigen::Vector3d evaluate(Eigen::Vector3d vec) const{
		return  vec - get_position();
	}
	Eigen::Matrix3d jacobian(Eigen::Vector3d vec) const{
		return Eigen::Matrix3d::Identity();
	}
};
#endif // #ifndef DUNE_BDF_HXX
