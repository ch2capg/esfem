/*! \file io_parameter.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 27. Januar 2016

     Implementation details for io_parameter.h
     Created by Christian Power on 27.01.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <config.h>
#include <sstream>
#include <stdexcept>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include <bdf.h>
#include "io_parameter.h"
#include "grid.h"

void dune_fem_parameter_append(int argc, char** argv, const std::string& file);
std::string get_macroGrid();
std::string doubleVector_to_string(const std::vector<double>&);

static const std::string project_dir {"/Users/christianpower/cpp/DISS_surfaces/"};

// ----------------------------------------------------------------------
// Implementation of Esfem::Io::Parameter

struct Esfem::Io::Parameter::Data{
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
  double tg_a {.1};
  double tg_b {.9};
  double tg_Dc {10.};
  double tg_gamma {30.};
};
Esfem::Io::Parameter::Parameter(int argc, char** argv,
				const std::string& parameter_file_name){
  dune_fem_parameter_append(argc, argv, parameter_file_name);
  auto bdf_no = Dune::Fem::Parameter::getValue<int>("heat.bdf", 1);
  if(bdf_no > 6 || bdf_no < 1)
    throw std::runtime_error{"Bad parameter 'heat.bdf'.  Enter a value in [1,6]."};
  
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
    NUMERIK::bdf_gamma_coeff(bdf_no),
    Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.a", .1),
    Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.b", .9),
    Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.Dc", 10.),
    Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.gamma", 30.)
  };
}
Esfem::Io::Parameter::~Parameter(){
  delete d_ptr;
  d_ptr = nullptr;
#ifdef DEBUG 
  std::cerr << "~Parameter(): deleted d_ptr\n";
#endif
}
const std::string& Esfem::Io::Parameter::grid() const noexcept{
  return d_ptr -> grid_dgf;
}
const std::string& Esfem::Io::Parameter::error_log() const noexcept{
  return d_ptr -> error_log;
}
const std::string& Esfem::Io::Parameter::paraview() const noexcept{
  return d_ptr -> paraview;
}
double Esfem::Io::Parameter::start_time() const noexcept{
  return d_ptr -> t_0;
}
double Esfem::Io::Parameter::global_timeStep() const{
  return d_ptr -> dT;
}
long Esfem::Io::Parameter::max_timeSteps() const{
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
double Esfem::Io::Parameter::eps() const noexcept{
  return d_ptr -> eps;
}
const std::vector<double>& Esfem::Io::Parameter::bdf_alphas() const noexcept{
  return d_ptr -> bdf_alphas;
}
const std::vector<double>& Esfem::Io::Parameter::bdf_gammas() const noexcept{
  return d_ptr -> bdf_gammas;
}
double Esfem::Io::Parameter::tg_a() const noexcept{
  return d_ptr -> tg_a;
}
double Esfem::Io::Parameter::tg_b() const noexcept{
  return d_ptr -> tg_b;
}
double Esfem::Io::Parameter::tg_Dc() const noexcept{
  return d_ptr -> tg_Dc;
}
double Esfem::Io::Parameter::tg_gamma() const noexcept{
  return d_ptr -> tg_gamma;
}
std::ostream& Esfem::Io::operator<<(std::ostream& os, const Parameter& d){
  const auto p = d.d_ptr;
  os << "t_0: " << p->t_0 << '\n'
     << "dT: " << p->dT << '\n'
     << "t_end: " << p->t_end << '\n'
     << "grid_dgf: " << p->grid_dgf << '\n'
     << "error_log: " << p->error_log << '\n'
     << "paraview: " << p->paraview << '\n'
     << "eps: " << p->eps << '\n'
     << "bdf_alphas: " << doubleVector_to_string(p->bdf_alphas) << '\n'
     << "bdf_gammas: " << doubleVector_to_string(p->bdf_gammas) << '\n'
     << "tg_a: " << p->tg_a << '\n'
     << "tg_b: " << p->tg_b << '\n'
     << "tg_Dc: " << p->tg_Dc << '\n'
     << "tg_gamma: " << p->tg_gamma;
  return os;
}

// ----------------------------------------------------------------------
// Internal Implementation

void dune_fem_parameter_append(int argc, char** argv, const std::string& file){
  using Dune::Fem::Parameter;
  Parameter::append(argc, argv);
  for( int i = 1; i < argc; ++i )
    Parameter::append(argv[i]);
#ifdef DEBUG
  std::clog << "Using parameter file: " << file << std::endl;
#endif 
  Parameter::append(file);
}
std::string get_macroGrid(){
  const auto grid_key = Dune::Fem::IOInterface::
    defaultGridKey(Esfem::Grid::Grid_and_time::Grid::dimension);
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
