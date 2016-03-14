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

namespace Esfem{
  void dune_fem_parameter_append(int argc, char** argv, const std::string& file);
  std::string get_gridKey();
  std::string get_macroGrid();
  std::string doubleVector_to_string(const std::vector<double>&);
  void check_bdfNo(const int);
}

static const std::string project_dir {"/Users/christianpower/cpp/DISS_surfaces/"};

// ----------------------------------------------------------------------
// Implementation of Esfem::Io::Parameter

struct Esfem::Io::Parameter::Data{
  const double t_0 {0.0};
  const double dT {0.1};
  const double t_end {0.6};
  const std::string grid_dgf {project_dir + "data/sphere/Sphere_1.dgf"};
  const std::string error_log {project_dir + "output/norm_errors.log"};
  const std::string paraview {project_dir + "output/video_"};
  const double eps {1e-8};
  const std::vector<double> bdf_alphas {-1., 1.};
  const std::vector<double> bdf_gammas {1.};
  // bdf coefficients ordered w.r.t. polynomial order (from X^0 to X^n)
  const double tg_a {.1};
  const double tg_b {.9};
  const double tg_Dc {10.};
  const double tg_gamma {30.};
  const double u_hom_value {1.};
  const double w_hom_value {.9};
  const double u_pertubation {.01};
  const double w_pertubation {.01};
  const std::string u_init_dof {project_dir + "output/u_dof.log"};
  const std::string w_init_dof {project_dir + "output/w_dof.log"};
};

Esfem::Io::Parameter::Parameter(int argc, char** argv,
				const std::string& parameter_file_name) try{
  using Dune::Fem::Parameter;
    
  Esfem::dune_fem_parameter_append(argc, argv, parameter_file_name);
  
  const auto bdf_no = Dune::Fem::Parameter::getValue<int>("heat.bdf", 1);
  check_bdfNo(bdf_no);

  const auto a = Parameter::getValue<double>("tumor_growth.heat.a", .1);
  const auto b = Parameter::getValue<double>("tumor_growth.heat.b", .9);
  const auto u_hom = Parameter::getValue<double>("tumor_growth.heat.u_hom",
						 a + b);
  const auto w_hom = Parameter::getValue<double>("tumor_growth.heat.w_hom",
						 b  / ( u_hom * u_hom ) );
  
  d_ptr = new Data {
    Parameter::getValue<double>("heat.starttime",0.0),
    Parameter::getValue<double>("heat.timestep",0.1),
    Parameter::getValue<double>("heat.endtime",0.6),
    Esfem::get_macroGrid(),
    Parameter::getValue<std::string>("fem.io.errorFile",
				      project_dir + "output/l2h1error"),
    Parameter::getValue<std::string> ("fem.io.outputName",
				       project_dir + "output/ALE_LiteratureExample-"),
    Dune::Fem::Parameter::getValue<double>("heat.solvereps", 1e-8),
    NUMERIK::bdf_alpha_coeff(bdf_no),
    NUMERIK::bdf_gamma_coeff(bdf_no),
    a,
    b,
    Parameter::getValue<double>("tumor_growth.heat.Dc", 10.),
    Parameter::getValue<double>("tumor_growth.heat.gamma", 30.),
    u_hom,
    w_hom,
    Parameter::getValue<double>("tumor_growth.heat.u_pertubation", .01),
    Parameter::getValue<double>("tumor_growth.heat.w_pertubation", .01),
    Parameter::getValue<std::string>("tumor_growth.io.u_init_dof",
				     project_dir + "output/u_dof.log"),
    Parameter::getValue<std::string>("tumor_growth.io.w_init_dof",
				     project_dir + "output/w_dof.log")
  };
 }
 catch(const std::exception&){
   std::throw_with_nested(std::runtime_error
			  {"Error in constructor of Parameter."});
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
double Esfem::Io::Parameter::u_hom_value() const noexcept{
  return d_ptr -> u_hom_value;
}
double Esfem::Io::Parameter::w_hom_value() const noexcept{
  return d_ptr -> w_hom_value;
}
double Esfem::Io::Parameter::u_pertubation() const noexcept{
  return d_ptr -> u_pertubation;
}
double Esfem::Io::Parameter::w_pertubation() const noexcept{
  return d_ptr -> w_pertubation;
}
const std::string& Esfem::Io::Parameter::u_init_dof() const noexcept{
  return d_ptr -> u_init_dof;
}
const std::string& Esfem::Io::Parameter::w_init_dof() const noexcept{
  return d_ptr -> w_init_dof;
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
     << "bdf_alphas: " << Esfem::doubleVector_to_string(p->bdf_alphas) << '\n'
     << "bdf_gammas: " << Esfem::doubleVector_to_string(p->bdf_gammas) << '\n'
     << "tg_a: " << p->tg_a << '\n'
     << "tg_b: " << p->tg_b << '\n'
     << "tg_Dc: " << p->tg_Dc << '\n'
     << "tg_gamma: " << p->tg_gamma << '\n'
     << "u_hom_value: " << p -> u_hom_value << '\n'
     << "w_hom_value: " << p -> w_hom_value << '\n'
     << "u_pertubation: " << p -> u_pertubation << '\n'
     << "w_pertubation: " << p -> w_pertubation << '\n'
     << "u_init_dof: " << p -> u_init_dof << '\n'
     << "w_init_dof: " << p -> w_init_dof
    ;
  return os;
}

// ----------------------------------------------------------------------
// Internal Implementation

void Esfem::dune_fem_parameter_append(int argc, char** argv, const std::string& file){
  using Dune::Fem::Parameter;
  Parameter::append(argc, argv);
  for( int i = 1; i < argc; ++i )
    Parameter::append(argv[i]);
#ifdef DEBUG
  std::clog << "Using parameter file: " << file << std::endl;
#endif 
  Parameter::append(file);
}
std::string Esfem::get_gridKey(){
  const auto grid_key = Dune::Fem::IOInterface::
    defaultGridKey(Esfem::Grid::Grid_and_time::Grid::dimension);
#ifdef DEBUG
  std::clog << "grid_key: " << grid_key << std::endl;
#endif
  return grid_key;
}
std::string Esfem::get_macroGrid(){
  const auto grid_key = get_gridKey();
  const auto macro_grid = Dune::Fem::Parameter::getValue<std::string>(grid_key);
  if(Dune::Fem::MPIManager::rank() == 0)
    std::clog << "Loading macro grid: " << macro_grid << std::endl;
  return macro_grid;
}
std::string Esfem::doubleVector_to_string(const std::vector<double>& vd){
  std::ostringstream oss;
  oss << '{';
  for(std::size_t it = 0; it < vd.size() - 1; ++it)
    oss << vd[it] << ", ";
  oss <<  vd.back() << '}';
  return oss.str();
}
void Esfem::check_bdfNo(const int i) try{
  if(i > 6 || i < 1)
    throw std::runtime_error{"Bad parameter \"heat.bdf\".  Enter a value in [1,6]."};
 }
 catch(const std::exception&){
   std::throw_with_nested(std::runtime_error {"Error in check_bdfNo()."});
 }
/*! Log:
 */
