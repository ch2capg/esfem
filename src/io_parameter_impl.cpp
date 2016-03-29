/*! \file io_parameter_impl.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Implementation details for io_parameter_impl.h

     Created by Christian Power on 15.03.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
*/

#include <fstream>
#include <dassert.h>
#include <config.h>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include "io_parameter_impl.h"
#include "grid.h"
#include "esfem_error.h"

using namespace std;

const std::string& Esfem::Impl::project_dir(){
  // static const std::string project_dir {"/Users/christianpower/cpp/DISS_surfaces/"}; 
  static const std::string project_dir {"/home/power/cpp/DISS_surfaces/"}; 
  // os depending
  return project_dir;
}
void Esfem::Impl::dune_fem_parameter_append(int argc, char** argv, const std::string& file){
  using Dune::Fem::Parameter;
  Parameter::append(argc, argv);
  for( int i = 1; i < argc; ++i )
    Parameter::append(argv[i]);
#ifndef NDEBUG
  std::clog << "Using parameter file: " << file << std::endl;
#endif 
  Parameter::append(file);
}
std::string Esfem::Impl::get_gridKey(){
  const auto grid_key = Dune::Fem::IOInterface::
    defaultGridKey(Esfem::Grid::Grid_and_time::Grid::dimension);
#ifndef NDEBUG
  std::clog << "grid_key: " << grid_key << std::endl;
#endif
  return grid_key;
}
std::string Esfem::Impl::get_macroGrid(){
  const auto grid_key = get_gridKey();
  const auto macro_grid = Dune::Fem::Parameter::getValue<std::string>(grid_key);
  if(Dune::Fem::MPIManager::rank() == 0)
    std::clog << "Loading macro grid: " << macro_grid << std::endl;
  return macro_grid;
}
std::string Esfem::Impl::doubleVector_to_string(const std::vector<double>& vd){
  constexpr size_t max_vec_size = 100;
  Assert::dynamic<Assert::level(1), Esfem::Parameter_error>
    (vd.size() <= max_vec_size, __FILE__, __LINE__, 
     "Vector is to big for doubleVector_to_string()");
  std::ostringstream oss;
  oss << '{';
  for(std::size_t it = 0; it < vd.size() - 1; ++it)
    oss << vd[it] << ", ";
  oss <<  vd.back() << '}';
  return oss.str();
}
void Esfem::Impl::file_check(const std::vector<std::string>& file_list){  
  for(const auto& file_name : file_list){
    fstream fs {file_name};
    Assert::dynamic<Assert::level(1), Esfem::Parameter_error>
      (fs.is_open(), __FILE__, __LINE__, "Could not open file " + file_name);
  }
}
