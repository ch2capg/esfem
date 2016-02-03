/*! \file parameter.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     Implementation details for esfem.h
     Created by Christian Power on 20.01.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */
#include "parameter.h"
#include <string>
#include <sstream>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include <bdf.h>
#include "fe_grid.h"


void dune_fem_parameter_append(int argc, char** argv, const char* file);
std::string get_macroGrid();
std::string doubleVector_to_string(const std::vector<double>&);

const std::string project_dir {"/Users/christianpower/cpp/DISS_surfaces/"};

// ----------------------------------------------------------------------
// Implementation of namespace Parameter

struct Parameter::PDE_data::Data{
  double t_0 {0.0};
  double dT {0.1};
  double t_end {0.6};
  std::string grid_dgf {project_dir + "data/sphere/Sphere_1.dgf"};
  std::string error_log {project_dir + "output/norm_errors.log"};
  std::string paraview {project_dir + "output/video_"};
  double eps {1e-8};
  std::vector<double> bdf_alphas {-1., 1.};
  std::vector<double> bdf_gammas {1.};
  // bdf coefficients ordered w.r.t. polynomial order (from X^0 to X^n)
};
Parameter::PDE_data::PDE_data(int argc, char** argv,
			      const char* const parameter_file_name){
  dune_fem_parameter_append(argc, argv, parameter_file_name);
  auto bdf_no = Dune::Fem::Parameter::getValue<int>("heat.bdf", 1);
  
  d_ptr = new Data {
    Dune::Fem::Parameter::getValue<double>("heat.starttime",0.0),
    Dune::Fem::Parameter::getValue<double>("heat.timestep",0.1),
    Dune::Fem::Parameter::getValue<double>("heat.endtime",0.6),
    get_macroGrid(),
    Dune::Fem::Parameter::
    getValue<std::string>("fem.io.errorFile",
			  "/Users/christianpower/cpp/DISS_surfaces/"
			  "output/l2h1error"),
    Dune::Fem::Parameter::
    getValue<std::string> ("fem.io.outputName",
			   "/Users/christianpower/cpp/DISS_surfaces/"
			   "output/ALE_LiteratureExample-"),
    Dune::Fem::Parameter::getValue<double>("heat.solvereps", 1e-8),
    NUMERIK::bdf_alpha_coeff(bdf_no),
    NUMERIK::bdf_gamma_coeff(bdf_no)
  };
}
Parameter::PDE_data::~PDE_data(){
  delete d_ptr;
#ifdef DEBUG 
  std::cerr << "~PDE_data(): deleted d_ptr\n";
#endif
}
const std::string& Parameter::PDE_data::grid() const noexcept{
  return d_ptr -> grid_dgf;
}
const std::string& Parameter::PDE_data::error_log() const noexcept{
  return d_ptr -> error_log;
}
double Parameter::PDE_data::start_time() const noexcept{
  return d_ptr -> t_0;
}
double Parameter::PDE_data::global_timeStep() const{
  return d_ptr -> dT;
}
long Parameter::PDE_data::max_timeSteps() const{
  const auto t_n = d_ptr->t_end;
  const auto t_0 = d_ptr->t_0;
  if(t_0 > t_n){
    std::ostringstream err_msg;
    err_msg << "Error in max_timeSteps().  t_start > t_end, "
	    << t_0 << " > " << t_n;
    throw std::runtime_error {err_msg.str()};
  }
  return (t_n - t_0)/d_ptr->dT + .1;
}
double Parameter::PDE_data::eps() const noexcept{
  return d_ptr -> eps;
}
const std::vector<double>& Parameter::PDE_data::bdf_alphas() const noexcept{
  return d_ptr -> bdf_alphas;
}
const std::vector<double>& Parameter::PDE_data::bdf_gammas() const noexcept{
  return d_ptr -> bdf_gammas;
}
std::ostream& Parameter::operator<<(std::ostream& os, const PDE_data& d){
  const auto p = d.d_ptr;
  os << "t_0: " << p->t_0 << '\n'
     << "dT: " << p->dT << '\n'
     << "t_end: " << p->t_end << '\n'
     << "grid_dgf: " << p->grid_dgf << '\n'
     << "error_log: " << p->error_log << '\n'
     << "paraview: " << p->paraview << '\n'
     << "eps: " << p->eps << '\n'
     << "bdf_alphas: " << doubleVector_to_string(p->bdf_alphas) << '\n'
     << "bdf_gammas: " << doubleVector_to_string(p->bdf_gammas);
  return os;
}

Parameter::Error_stream::Error_stream(const PDE_data& d)
  : ofs {d.error_log(), std::ios_base::app}
    {}
Parameter::Error_stream& Parameter::Error_stream::operator<<(StdManip manip){
#if DEBUG > 8
  std::clog << "Using operator<<(StdManip)." << std::endl;
#endif
  manip(ofs);
  return *this;
}
Parameter::Error_stream& Parameter::Error_stream::
operator<<(const std::vector<double>& vd){
  for(std::size_t it = 0; it < vd.size() -1; ++it)
    ofs << vd[it] << ' ';
  ofs << vd.back();
  return *this;
}
Parameter::Error_stream& Parameter::Error_stream::close(){
  ofs.close();
  return *this;
}

// ----------------------------------------------------------------------
// Internal implementation

void dune_fem_parameter_append(int argc, char** argv, const char* file){
  using Dune::Fem::Parameter;
  Parameter::append(argc, argv);
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append(argv[i]);
  std::clog << "Using parameter file: " << file << std::endl;
  Dune::Fem::Parameter::append(file);
}
std::string get_macroGrid(){
  const auto grid_key =
    Dune::Fem::IOInterface::defaultGridKey(FE_grid::Grid::dimension);
#ifdef DEBUG
  std::clog << "grid_key: " << grid_key << std::endl;
#endif
  
  const auto macro_grid = Dune::Fem::Parameter::getValue<std::string>(grid_key);
  if(Dune::Fem::MPIManager::rank() == 0)
    std::clog << "Loading macro grid: " << macro_grid << std::endl;

  return macro_grid;
}
std::string doubleVector_to_string(const std::vector<double>& vd){
  std::ostringstream oss;
  oss << '{';
  for(std::size_t it = 0; it < vd.size() - 1; ++it)
    oss << vd[it] << ", ";
  oss <<  vd.back() << '}';
  return oss.str();
}
/*! Log:
 */
