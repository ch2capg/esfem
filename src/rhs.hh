#ifndef RHS_HH
#define RHS_HH

#include <dune/fem/quadrature/cachingquadrature.hh>

// assembleRHS
// -----------

template< class Function, class DiscreteFunction >
void assembleRHS ( const Function &function, DiscreteFunction &rhs )
{
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunction::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  rhs.clear();

  const DiscreteFunctionSpaceType &dfSpace = rhs.space();

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();

    LocalFunctionType rhsLocal = rhs.localFunction( entity );
    
    QuadratureType quadrature( entity, 2*dfSpace.order()+1 );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      // obtain quadrature point
      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );

      // evaluate f
      typename Function::RangeType f;
      function.evaluate( geometry.global( x ), f );

      // multiply by quadrature weight
      f *= quadrature.weight( pt ) * geometry.integrationElement( x );

      // add f * phi_i to rhsLocal[ i ]
      rhsLocal.axpy( quadrature[ pt ], f );
    }
  }
  rhs.communicate();
}

#endif // #ifndef RHS_HH
