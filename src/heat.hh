#ifndef HEAT_HH
#define HEAT_HH

#include <cmath>

#include <dune/fem/function/common/function.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/quadrature/quadrature.hh>

// EvolutionFunction
// -----------------

template <class FunctionSpace, class Impl>
struct EvolutionFunction : Dune::Fem::Function<FunctionSpace, Impl>{
	EvolutionFunction(const Dune::Fem::TimeProviderBase& time_provider) :
		tp {time_provider} {}
	double get_time() const{ return tp.time(); }
	// double get_time() const{ return t; }
	// void set_time(const double time) { t = time; }
private:
	const Dune::Fem::TimeProviderBase& tp;
	// double t;
};

namespace Problem{
				  	
	// RHSFunction
    // -----------
	template< class FunctionSpace >
	struct RHSFunction : EvolutionFunction<FunctionSpace, 
										   RHSFunction<FunctionSpace> >{
		typedef FunctionSpace FunctionSpaceType;
		typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
		typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
		typedef typename FunctionSpaceType::DomainType DomainType;
		typedef typename FunctionSpaceType::RangeType RangeType;
		static const int dimDomain = DomainType::dimension;
		static const int dimRange = RangeType::dimension;

		typedef EvolutionFunction<FunctionSpace, 
								  RHSFunction<FunctionSpace> >  Base;
		using Base::get_time;
		
		RHSFunction(const Dune::Fem::TimeProviderBase& time_provider) :
			Base {time_provider} {}
		void evaluate(const DomainType &xx, RangeType &phi) const{
			
			const DomainFieldType x = xx[ 0 ], y = xx[ 1 ], z = xx[ 2 ];
			const double t = get_time();

			// ALE paper 2. review
			phi = 4*(4*M_PI*pow(x,2)*pow(y,2)*pow(z,2)*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),3)*pow(sin(2*M_PI*t) + 4,2)) + (2*M_PI*pow(x,2)*pow(y,2)*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)) - M_PI*pow(x,2)*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)))*(pow(y,2)/(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2)) - 1) + (2*M_PI*pow(x,2)*pow(z,2)*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)) - M_PI*pow(x,2)*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)))*(pow(z,2)/(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2)) - 1) - 4*(3*M_PI*pow(x,2)*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,3)) - 32*M_PI*pow(x,4)*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,5)))*(16*pow(x,2)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)) - 1) - 8*x*y*(M_PI*x*y*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)) - 16*M_PI*pow(x,3)*y*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,4)))/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*(sin(2*M_PI*t) + 4)) - 8*x*z*(M_PI*x*z*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)) - 16*M_PI*pow(x,3)*z*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,4)))/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*(sin(2*M_PI*t) + 4)) + 32*M_PI*pow(x,4)*pow(y,2)*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),3)*pow(sin(2*M_PI*t) + 4,4)) + 32*M_PI*pow(x,4)*pow(z,2)*cos(2*M_PI*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),3)*pow(sin(2*M_PI*t) + 4,4)))*x*y*exp(-6*t) + 4*M_PI*pow(x,3)*y*cos(2*M_PI*t)*exp(-6*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)) - 6*x*y*exp(-6*t) + 16*M_PI*pow(x,3)*y*cos(2*M_PI*t)*exp(-6*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,3)) + 2*(x*(pow(y,3)/pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2) - y/(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2)))*exp(-6*t) + 4*x*pow(y,3)*exp(-6*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*(sin(2*M_PI*t) + 4)) - 4*x*y*exp(-6*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*(sin(2*M_PI*t) + 4)))*(pow(y,2)/(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2)) - 1) + (2*x*pow(y,2)*z*exp(-6*t)/pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2) - x*z*exp(-6*t)/(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2)) + 8*x*pow(y,2)*z*exp(-6*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*(sin(2*M_PI*t) + 4)) - 4*x*z*exp(-6*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*(sin(2*M_PI*t) + 4)))*y*z/(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2)) + 2*(x*pow(y,2)*z*exp(-6*t)/pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2) + 4*x*pow(y,2)*z*exp(-6*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*(sin(2*M_PI*t) + 4)))*y*z/(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2)) + (2*x*y*pow(z,2)*exp(-6*t)/pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2) - x*y*exp(-6*t)/(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2)) + 8*x*y*pow(z,2)*exp(-6*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*(sin(2*M_PI*t) + 4)) - 4*x*y*exp(-6*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*(sin(2*M_PI*t) + 4)))*(pow(z,2)/(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2)) - 1) - 8*(4*y*(x/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)) - 16*pow(x,3)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,4)))*exp(-6*t) + x*y*exp(-6*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*(sin(2*M_PI*t) + 4)) - 16*pow(x,3)*y*exp(-6*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,3)))*(16*pow(x,2)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)) - 1) - 4*((pow(y,2)/(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2)) - 1)*exp(-6*t) - 32*pow(x,2)*pow(y,2)*exp(-6*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)) + 4*pow(y,2)*exp(-6*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*(sin(2*M_PI*t) + 4)) - 128*pow(x,2)*pow(y,2)*exp(-6*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,3)))*x*y/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*(sin(2*M_PI*t) + 4)) - 4*((16*pow(x,2)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2)) - 1)*exp(-6*t) - 8*pow(x,2)*pow(y,2)*exp(-6*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*(sin(2*M_PI*t) + 4)) - 32*pow(x,2)*pow(y,2)*exp(-6*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)) + 4*pow(x,2)*exp(-6*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*(sin(2*M_PI*t) + 4)))*x*y/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*(sin(2*M_PI*t) + 4)) - 4*(y*z*exp(-6*t)/(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2)) - 32*pow(x,2)*y*z*exp(-6*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)) + 4*y*z*exp(-6*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*(sin(2*M_PI*t) + 4)) - 128*pow(x,2)*y*z*exp(-6*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,3)))*x*z/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*(sin(2*M_PI*t) + 4)) + 32*(pow(x,2)*y*z*exp(-6*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*(sin(2*M_PI*t) + 4)) + 4*pow(x,2)*y*z*exp(-6*t)/(pow(pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2),2)*pow(sin(2*M_PI*t) + 4,2)))*x*z/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*(sin(2*M_PI*t) + 4));
		}
	};


	// InitialData
	// -----------
	template <class FunctionSpace>
	struct InitialData : EvolutionFunction<FunctionSpace, 
										   InitialData<FunctionSpace> >{
		typedef FunctionSpace FunctionSpaceType;
		typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
		typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
		typedef typename FunctionSpaceType::DomainType DomainType;
		typedef typename FunctionSpaceType::RangeType RangeType;
		static const int dimDomain = DomainType::dimension;
		static const int dimRange  = RangeType::dimension;

		typedef EvolutionFunction<FunctionSpace, 
								  InitialData<FunctionSpace> >  Base;
		using Base::get_time;
		
		InitialData(const Dune::Fem::TimeProviderBase& time_provider) :
			Base {time_provider} {}
		void evaluate(const DomainType& xx, RangeType& phi) const{
			// CAPG 11.06.2014
			const double x = xx[0], y = xx[1]; // z = xx[2]; unused
			const double t = get_time();

			phi = x*y * exp(-6.*t);	// heat problem

			// Other Functions
			// This is the constant one function.
			// phi = RangeType {1};
			// phi = RangeType {0};
			// phi = x * y;	// poisson problem
		}
  };

	// HeatModel
	// ---------
	template <class FunctionSpace>
	struct HeatModel{
		typedef FunctionSpace FunctionSpaceType;
		typedef typename FunctionSpaceType::DomainType DomainType;
		typedef typename FunctionSpaceType::RangeType RangeType;
		typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
		typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
		typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

		HeatModel(const Dune::Fem::TimeProviderBase& input_timeProvider)
			: timeProvider (input_timeProvider) {}

		void set_somefactor(const double d) {somefactor = d;}
		double get_somefactor() const {return somefactor;}


		template< class Entity, class Point >
		void massFlux ( const Entity& entity, const Point& xx,
						const RangeType& value, RangeType& flux ) const;
		// 'massFlux' is the semilinearity in the elliptic part of our PDE.
		// massFlux(entity, xx, value, flux) calculates: flux = m(xx, u(xx) );
		//
		// Legend: value == u( quadratur[pt] ), xx == quadrature[pt].  
		// 'entity' is needed, if you want the global coordinates from xx, c.f. code
		// below.
		

		template< class Entity, class Point >
		void diffusiveFlux ( const Entity &entity, const Point &x, 
							 const JacobianRangeType &gradient, 
							 JacobianRangeType &flux ) const;
		// 'diffusiveFlux' is the multiplication with the diffusive tensor 
		// in the elliptic part of our PDE.
		// diffuxivFlux(enityt, x, gradient, flux) calculates:
		// flux = D(x, u(x)) ∇u(x);
		//
		// c.f. massFlux.
		
		template< class Entity, class Point >
		void aleFlux ( const Entity& entity, const Point& xx, 
					   const RangeType& value, 
					   JacobianRangeType& flux ) const;

		template< class Entity, class Point >
		void massFlux_2 ( const Entity& entity, const Point& xx,
						  const RangeType& value, RangeType& flux ) const;

	private:
		const Dune::Fem::TimeProviderBase& timeProvider;
		double somefactor {1.};
	};

	// NonlinearHeatModel
	// ---------
	template <typename FunctionSpace>
	struct NonlinearHeatModel{
		typedef FunctionSpace FunctionSpaceType;
		typedef typename FunctionSpaceType::DomainType DomainType;
		typedef typename FunctionSpaceType::RangeType RangeType;
		typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
		typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
		typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
		static const int dim_domain = DomainType::dimension;
		static const int dim_range = RangeType::dimension;

		explicit NonlinearHeatModel(const Dune::Fem::TimeProviderBase& time_provider) :
			tp {time_provider} {}

		template <typename Entity, typename Point>
		void jac_ran_to_jac_ran(const Entity& entity, const Point& xx, 
							    const JacobianRangeType& du_x, const RangeType& xi_x, 
								JacobianRangeType& a_xi_du_x) const;
		template <typename Entity, typename Point>
		void ran_to_jac_ran(const Entity& entity, const Point& xx, const RangeType& u_x, 
							JacobianRangeType& vec_u_x) const;
		template <typename Entity, typename Point>
		void jac_ran_to_ran(const Entity& entity, const Point& xx, 
							const JacobianRangeType& du_x, RangeType& du_x_vec) const;
		template <typename Entity, typename Point>
		void ran_to_ran(const Entity& entity, const Point& xx, const RangeType& u_x, 
						RangeType& u_x_chi) const;
	private:
		const Dune::Fem::TimeProviderBase& tp;
	};



	// All implementation details
	
	// Implementation details for 'HeatModel'

    // Assumption: massFlux is (Mf) divergence of the velocity field
    // Modification: added time
	template<class FunctionSpace> template< class Entity, class Point >
    inline void HeatModel<FunctionSpace>
	::massFlux( const Entity &entity, const Point &xx,
				const RangeType &value, RangeType &flux ) const
    {
		// ------------------------------------------------------------
		// Legend for the argument:
		// value = u( quadrature[ pt ] )
		// Point = quadrature[ pt ]
		// ------------------------------------------------------------
		
		flux = value;
		// if(explicitOp == true)
		// 	flux = value;  // comment this out!
		// else
		// 	flux = value;
		
		// CAPG 11.06.2014
		// need some testing if it actually works;
		// don't forget to comment out the time variable
		// const RangeType x = xx[ 0 ], y = xx[ 1 ], z = xx[ 2 ];
		// const double t = timeVariable_;
    }

	template<class FunctionSpace> template<class Entity, class Point>
	inline void HeatModel<FunctionSpace>
	::diffusiveFlux( const Entity &entity, const Point &x, 
					 const JacobianRangeType &gradient, JacobianRangeType &flux ) const
  	{
		flux = gradient;				// Laplace-Beltrami example
		flux *= timeProvider.deltaT();	// moving surfaces example
		/*
		if(explicitOp == true){
			flux *=  0;
			// flux[0][0] *= 0;
			// flux[0][1] *= 0;
			// flux[0][2] *= 0;
			// flux[0][3] *= 0; // danger!!!!! this doesn't throw an error
		}
		else{
			flux *=timeProvider_.deltaT();
			// CAPG:You have multiplied // TODO: he whole equation with \Delta t.
			//  The somefactor = theta. Thus I've commented it out!!
			// flux *= timeProvider_.deltaT();
		}
		// We use the identiy matrix, hence we comment this line out.  
		// DomainType xGlobal = entity.geometry().global ( Dune::coordinate(x) );
		// flux[0][0] *= (xGlobal[1]>0)?0.1:2.1;   //Warning flux is not a Matrix!!
		// flux[0][1] *= (xGlobal[1]>0)?0.1:2.1;
		*/
	}

	template<class FunctionSpace> template< class Entity, class Point >
	inline void HeatModel<FunctionSpace>
	::aleFlux ( const Entity& entity, const Point& xx, 
				const RangeType& value, 
				JacobianRangeType& flux ) const {
		// We assemble the weak equation: ∫ Uʰ (Wʰ - Vʰ) · ∇ϕʰ  ∀ϕʰ FE-functions

		//flux *= 0;
		const DomainType xGlobal = entity.geometry().global ( Dune::coordinate(xx) );
		const double x = xGlobal[0], y = xGlobal[1], z = xGlobal[2], 
			t = timeProvider.time(), dT = timeProvider.deltaT();
		// flux is a R¹ˣ³ Matrix, so flux has two indices.
		// flux = W_h
		flux[0][0] = .5 * M_PI * cos(2. * M_PI * t) / ( 2. * sqrt( 1. + .25 * sin(2. * M_PI * t) ) ) *x;
		flux[0][1] = 0.;
		flux[0][2] = 0.;
		// flux -= V_h
		flux[0][0] -= 16.*M_PI*pow(x,3)*cos(2*M_PI*t)/( (pow(y,2) + pow(z,2) + 16.*pow(x,2)/pow(sin(2.*M_PI*t) + 4.,2) )*pow(sin(2.*M_PI*t) + 4.,3) );
		flux[0][1] -= 4.*M_PI*pow(x,2)*y*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16.*pow(x,2)/pow(sin(2.*M_PI*t) + 4.,2) )*pow(sin(2.*M_PI*t) + 4.,2));
		flux[0][2] -= 4.*M_PI*pow(x,2)*z*cos(2.*M_PI*t)/((pow(y,2) + pow(z,2) + 16.*pow(x,2)/pow(sin(2*M_PI*t) + 4.,2) )*pow(sin(2.*M_PI*t) + 4,2) ); 
		flux *= value * dT;
		// flux[0][0] *= value * dT;
		// flux[0][1] *= value * dT;
		// flux[0][2] *= value * dT;
	}

	template<class FunctionSpace> template< class Entity, class Point >
    inline void HeatModel<FunctionSpace>
	::massFlux_2( const Entity &entity, const Point &xx,
				const RangeType &value, RangeType &flux ) const
    {
		// ------------------------------------------------------------
		// Legend for the argument:
		// value = u( quadrature[ pt ] )
		// Point = quadrature[ pt ]
		// ------------------------------------------------------------
		
		flux = value * timeProvider.deltaT();
		// if(explicitOp == true)
		// 	flux = value;  // comment this out!
		// else
		// 	flux = value;
		
		// CAPG 11.06.2014
		// need some testing if it actually works;
		// don't forget to comment out the time variable
		// const RangeType x = xx[ 0 ], y = xx[ 1 ], z = xx[ 2 ];
		// const double t = timeVariable_;
    }

	// Implementation details for 'NonlinearHeatModel'
	template <typename FunctionSpace> 
	template <typename Entity, typename Point>
	inline void NonlinearHeatModel<FunctionSpace>::
	jac_ran_to_jac_ran(const Entity& entity, const Point& xx, 
					   const JacobianRangeType& du_x, const RangeType& xi_x, 
					   JacobianRangeType& a_xi_du_x) const{
		// DomainType xGlobal = entity.geometry().global(Dune::coordinate(xx));
		// xGlobal.size() == dim_domain == 3;
		// DomainFieldType x = xGlobal[0], y = xGlobal[1], z = xGlobal[2];
		// double t = tp.time(), dT = tp.deltaT();
	    // 'a_xi_du_x' is a R¹ˣ³ Matrix, so it has two indices.		

		// double a_tensor =  1. - ( exp(- (xi_x * xi_x)/4. ) )/2.;
		// double a_tensor =  1.;
		a_xi_du_x = du_x;
		// a_xi_du_x *= tp.deltaT() * a_tensor;
		a_xi_du_x *= tp.deltaT();
	}
	template <typename FunctionSpace> 
	template <typename Entity, typename Point>
	inline void NonlinearHeatModel<FunctionSpace>::
	ran_to_jac_ran(const Entity& entity, const Point& xx, const RangeType& u_x, 
				   JacobianRangeType& vec_u_x) const{
		DomainType xGlobal = entity.geometry().global(Dune::coordinate(xx));
		// xGlobal.size() == dim_domain == 3;
		DomainFieldType x = xGlobal[0], y = xGlobal[1], z = xGlobal[2];
		double t = tp.time(), dT = tp.deltaT();
		
		vec_u_x[0][0] = 16*M_PI*pow(x,3)*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,3));
		vec_u_x[0][1] = 4*M_PI*pow(x,2)*y*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2));
		vec_u_x[0][2] = 4*M_PI*pow(x,2)*z*cos(2*M_PI*t)/((pow(y,2) + pow(z,2) + 16*pow(x,2)/pow(sin(2*M_PI*t) + 4,2))*pow(sin(2*M_PI*t) + 4,2));
		vec_u_x *= (-1);						// vec_u_x = - I_h v_h == - V_h
		vec_u_x[0][0] += x * M_PI * cos(2*M_PI*t) 
			/ (4. + sin(2*M_PI*t));				// vec_u_x = W_h - V_h
		vec_u_x *= u_x * dT;					// vec_u_x = dT * U_h * (W_h - V_h)
	}
	template <typename FunctionSpace> 	
	template <typename Entity, typename Point>
	inline void NonlinearHeatModel<FunctionSpace>::
	jac_ran_to_ran(const Entity& entity, const Point& xx, 
				   const JacobianRangeType& du_x, RangeType& du_x_vec) const{
		// DomainType xGlobal = entity.geometry().global(Dune::coordinate(xx));
		// xGlobal.size() == dim_domain == 3;
		// DomainFieldType x = xGlobal[0], y = xGlobal[1], z = xGlobal[2];
		// double t = tp.time(), dT = tp.deltaT();

		du_x_vec = RangeType {0};
	}
	template <typename FunctionSpace> 
	template <typename Entity, typename Point>
	inline void NonlinearHeatModel<FunctionSpace>::
	ran_to_ran(const Entity& entity, const Point& xx, const RangeType& u_x, 
			   RangeType& u_x_chi) const{
		// DomainType xGlobal = entity.geometry().global(Dune::coordinate(x));
		// xGlobal.size() == dim_domain == 3;
		// DomainFieldType x = xGlobal[0], y = xGlobal[1], z = xGlobal[2];
		// double t = tp.time(), dT = tp.deltaT();		
		u_x_chi = u_x;
	}
}
#endif // #ifndef HEAT_HH
