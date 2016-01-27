#include "2015linfty.h"

int main(int argc, char** argv) try{
  Dune::Fem::MPIManager::initialize(argc, argv);

  Dune::Fem::Parameter::append(argc, argv);
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append(argv[i]);
  Dune::Fem::Parameter::append("../data/parameter");

  BDF::nonlinear_algorithm();

  return 0;
}
catch( const Dune::Exception& e )
{
  std::cerr << "Error: " << e << std::endl;
  return 1;
}
