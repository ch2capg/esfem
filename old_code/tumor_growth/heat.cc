#include <config.h>

// iostream includes
#include <iostream>

// include GeometryGrid
#include <dune/grid/geometrygrid.hh>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// include discrete function space
//#include <dune/fem/space/lagrangespace.hh> veraltet
#include <dune/fem/space/lagrange.hh>

// include discrete function
#include <dune/fem/function/adaptivefunction.hh>

// include solvers
//#include <dune/fem/solver/inverseoperators.hh> veraltet
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/oemsolver.hh>

// include Lagrange interpolation
#include <dune/fem/operator/lagrangeinterpolation.hh>

// include norms
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// include parameter handling
#include <dune/fem/io/parameter.hh>

// include output
#include <dune/fem/io/file/dataoutput.hh>

// local includes
#include "deformation.hh"
#include "elliptic.hh"
#include "heat.hh"
#include "rhs.hh"
// local CAPG includes
#include "dune_typedef_heat.hpp"
// this header only contains lots of typedef's
#include "dune_bdf.hpp"
#include "dune_heat_algorithm.hpp"

// main
// ----

int main ( int argc, char **argv )
try
{
  // initialize MPI, if necessary
  Dune::Fem::MPIManager::initialize(argc, argv);

  // append all given parameters
  Dune::Fem::Parameter::append(argc, argv);
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append(argv[i]);
  Dune::Fem::Parameter::append("../data/parameter");

  // BDF::ie_heat_algorithm();
  // BDF::bdf2_heat_algorithm();  
  // BDF::nonlinear_algorithm();
  // nonlinear_evolution();
  tumor_growth_code();

  return 0;
}
catch( const Dune::Exception& e )
{
  std::cerr << "Error: " << e << std::endl;
  return 1;
}
