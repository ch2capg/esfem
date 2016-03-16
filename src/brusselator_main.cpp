#include <iostream>
#include <config.h>
#include <dune/fem/io/parameter.hh>
#include <dune/common/exceptions.hh>
// #include "linHeat_algo.h"
#include "brusselator_algo.h"

void print_errMsg(const std::exception&);

int main(int argc, char** argv) try{
  // linHeat_algo(argc, argv);
  Esfem::brusselator_algo(argc, argv);
} 
catch(const Dune::Exception& e){
  std::cerr << "Dune error: " << e << std::endl;
  return 1;
}
catch(const std::exception& e){
  print_errMsg(e);
  return 2;
}

// ----------------------------------------------------------------------
// Internal Implementation

void print_errMsg(const std::exception& e){
  std::cerr << e.what() << '\n';
  try{
    std::rethrow_if_nested(e);
  }
  catch(const std::exception& nested){
    print_errMsg(nested);
  }
}
