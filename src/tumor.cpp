#include <config.h>
#include <iostream>

#include <dune/fem/io/parameter.hh>
// #include <dune/fem/io/file/dataoutput.hh>

#include "tumor_code.h"

int main(int argc, char** argv) try{
  // initialize MPI, if necessary
  Dune::Fem::MPIManager::initialize(argc, argv);

  // append all given parameters
  Dune::Fem::Parameter::append(argc, argv);
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append(argv[i]);
  const std::string p_file 
  {"/Users/christianpower/cpp/DISS_surfaces/data/tumor_parameter.txt"};
  Dune::Fem::Parameter::append(p_file.c_str() );

  tumor_growth_code();
}
catch(const Dune::Exception& e){
  std::cerr << "Error: " << e << std::endl;
  return 1;
}
catch(const std::exception& e){
  std::cerr << e.what() << std::endl;
  return 2;
}
