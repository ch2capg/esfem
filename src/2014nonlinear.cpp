#include <config.h>

#include <iostream>

// dune headers
#include <dune/grid/geometrygrid.hh>	// geometry grid
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>	// grid part
#include <dune/fem/space/lagrange.hh>	// discrete function space
#include <dune/fem/function/adaptivefunction.hh>	// discrete function
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/oemsolver.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/dataoutput.hh>

// two headers from config.h in fem-school
// #include <dune/grid/alugrid.hh>
// #include <dune/grid/io/file/dgfparser/dgfalu.hh>

// two headers from config.h in dune-alugrid
// #include <dune/alugrid/grid.hh>
// #include <dune/alugrid/dgf.hh>	// this is not working!

// own developped moving surface code 
#include "deformation.hh"
#include "elliptic.hh"
#include "heat.hh"
#include "rhs.hh"
#include "dune_typedef_heat.hpp"
#include "dune_bdf.hpp"
#include "dune_heat_algorithm.hpp"

int main(int argc, char** argv) try {
  Dune::Fem::MPIManager::initialize(argc, argv);

  Dune::Fem::Parameter::append(argc, argv);
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append(argv[i]);

  const std::string parameter_dir {"data/2014nonlinear_parameter.txt"};
  Dune::Fem::Parameter::append(parameter_dir.c_str());
  BDF::nonlinear_algorithm();
}
catch(const Dune::Exception& e){
  std::cerr << "Error: " << e << std::endl;
  return 1;
}
