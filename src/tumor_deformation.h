#ifndef DEFORMATION_HH
#define DEFORMATION_HH

#include <dune/common/exceptions.hh>
#include <dune/grid/geometrygrid/coordfunction.hh>
#include <dune/fem/space/common/functionspace.hh>

#include<cassert>
#include<dune_bdf.hpp>

#include <Eigen/Dense>

// DeformationCoordFunction
// ------------------------

class DeformationCoordFunction
  : public Dune::AnalyticalCoordFunction<double, 3, 3, DeformationCoordFunction>{
  typedef Dune::AnalyticalCoordFunction<double, 3, 3, DeformationCoordFunction> 
  BaseType;
public:
  typedef Dune::Fem::FunctionSpace<double, double, 3, 3> FunctionSpaceType;

  typedef BaseType::DomainVector DomainVector;
  typedef BaseType::RangeVector RangeVector;

  DeformationCoordFunction() {};
  DeformationCoordFunction (BDF::EvoMapType& evoMap)
    : evoMapPt {&evoMap} {}
  void set_time_provider(Dune::Fem::TimeProviderBase& time_provider){
    p_tp = &time_provider;
  }
  void evaluate (const DomainVector& x, RangeVector& y) const{
    // double t = p_tp->time();
    
    // '2015linfty' example
    // y[0] = sqrt(1. + sin(2. * M_PI * t )/4.) * x[0];
    // y[1] = x[1];
    // y[2] = x[2];		
    
    // ODE solver code or nonlinear evolution code
    // std:: cout <<  "Hello\n"
    // 	       << x[0] << ' ' << x[1] << ' ' << x[2] << std::endl;

    std::array<double, 3> domain_node { {x[0], x[1], x[2]} };

    const std::array<double, 3> range_node = (*evoMapPt)[domain_node];
    for(std::size_t s = 0; s < range_node.size(); ++s)
      y[s] = range_node[s];

    // 'baseball bat' example
    //const double newtime
    // = (time < 0.5 ? 0.0 : (time < 0.75 ? 4*(time - 0.5) : 1.0));
    //const double newtime = std::min( time, 1.0 );
    //const double r1 = std::abs( x[ 0 ] );
    //const double target 
    // = (1.0 - (r1*r1))*((r1*r1) + 0.05) + (r1*r1)*sqrt(1.0 - (r1*r1));
			
    //const double r2 = std::sqrt( x[1]*x[1] + x[2]*x[2] );
    //const double factor 
    // = std::exp( -2*newtime )*r2 + (1.0 - std::exp( -2*newtime ))*target;

    //y[ 0 ] = 2 * x[ 0 ] + newtime*(x[ 0 ] > 0 ? 2.0 : -1.0 )*x[ 0 ];
    //y[ 1 ] = factor * x[ 1 ] / (r2 + 0.000001);
    //y[ 2 ] = factor * x[ 2 ] / (r2 + 0.000001);
    
    //y[ 0 ] = x[ 0 ];
    //y[ 1 ] = x[ 1 ]; //(1+ .05 * std::sin(4*M_PI*time));
    //y[ 2 ] = x[ 2 ]; //* (1+ .5 * std::cos(2*M_PI*time));	  
  }
private:
  BDF::EvoMapType* evoMapPt = nullptr;	
	// use this if your surface evolves via an ODE
  Dune::Fem::TimeProviderBase* p_tp = nullptr;
};

// ------------------------------------------------------------
// CAPG: For what is this needed?????
// ------------------------------------------------------------

//! deformation depending on a discrete function 
template <class DiscreteFunctionType>
class DeformationDiscreteFunction
  : public Dune::DiscreteCoordFunction<double, 3, 
				       DeformationDiscreteFunction<DiscreteFunctionType> 
				       >{
  typedef Dune::DiscreteCoordFunction<double, 3, 
				      DeformationDiscreteFunction<DiscreteFunctionType>
				      > BaseType;
  typedef typename DiscreteFunctionType::GridType GridType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename DiscreteFunctionType::RangeType  RangeType;
public:
  DeformationDiscreteFunction (const DiscreteFunctionType& vertices )
    : vertices_( vertices )
  {}

  template< class HostEntity , class RangeVector >
  void evaluate ( const HostEntity& hostEntity, unsigned int corner,
		  RangeVector& y ) const{
    DUNE_THROW(Dune::NotImplemented,"evaluate not implemented for codim > 0");
  }

  template <class RangeVector> 
  void evaluate ( const typename GridType:: template Codim<0>::Entity& entity, 
		  unsigned int corner, RangeVector& y ) const{
    typedef typename GridType::ctype  ctype;
    enum { dim = GridType::dimension };
    
    const Dune::GenericReferenceElement<ctype, dim>& refElement
      = Dune::GenericReferenceElements<ctype, dim>::general( entity.type() );
    
    LocalFunctionType localVertices = vertices_.localFunction( entity );	
    
    localVertices.evaluate( refElement.position( corner, dim ), y );
  }
  // void setTime ( const double time ){}
protected:
  const DiscreteFunctionType& vertices_;
};

#endif // #ifndef DEFORMATION_HH
