/*! \file io_dgf.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Testing the components of io_dgf.h.

     Created by Christian Power on 08.03.2016

     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include "esfem.h"

void print_errMsg(const std::exception&);
template<typename FEFun>
std::string name(const FEFun& fef){
  constexpr auto base_dir = "/Users/christianpower/cpp/DISS_surfaces/test/";
  return base_dir + fef.name() + ".dgf";
}

using namespace std;

int main(int argc, char** argv) try{
  class {
    int i {0};
  public:
    void operator()(const string& m = "Hello"){
      clog << m << ' ' << i++ << endl;
    }
  } hello;
  
  const string this_dir {"/Users/christianpower/cpp/DISS_surfaces/test/"};
  
  Dune::Fem::MPIManager::initialize(argc, argv);
  
  const auto parameter_file = this_dir + "test_parameter.txt";
  Esfem::Io::Parameter data {argc, argv, parameter_file};
  // Two additional parameter
  const auto dgf_inputFile = this_dir + "input_test00.dgf";
  const auto dgf_outputFile = "/Users/christianpower/cpp/DISS_surfaces/test/"
    "output_test00.dgf";
  
  const Esfem::Io::Dgf::Handler dgf_interpreter {dgf_inputFile};
  Esfem::Grid::Grid_and_time grid {data};
  Esfem::Grid::Scal_FEfun fef {"fef", grid};
  Esfem::Grid::Vec_FEfun vfef {"vfef", grid};

  dgf_interpreter.read(dgf_inputFile, vfef);
  vfef *= -1.;
  dgf_interpreter.write(name(vfef), vfef);

  fef += 1.;
  dgf_interpreter.write(name(fef), fef);

  try{
    hello(name(fef));
    dgf_interpreter.read(name(fef), vfef); }
  catch(const exception& e){
    clog << "Caught error 0 as expected." << endl;
    print_errMsg(e);
  }  
  try{
    hello(name(vfef));
    dgf_interpreter.read(name(vfef), fef); }
  catch(const exception& e){
    clog << "Caught error 1 as expected." << endl;
    print_errMsg(e);
  }
 }
 catch(const exception& e){
   print_errMsg(e);
   return 1;
 }
 catch(const Dune::Exception& e){
   cerr << "Dune error: " << e << endl;
   return 2;
 }
 catch(...){
   cerr << "Unkown error in main().\n";
   return 3;
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

/*! Log:
 */
