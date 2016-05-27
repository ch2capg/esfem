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
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include <bdf.h>
#include <dassert.h>
#include "esfem_error.h"
#include "io_parameter.h"
#include "io_parameter_impl.h"

using namespace std;
// using Dune::Fem::Parameter;
using Esfem::Io::Parameter;
// ----------------------------------------------------------------------
// Implementation of Esfem::Io::Parameter

struct Parameter::Data{
  const double t_0;
  const double t_end;
  const double t_end_pattern;
  const double dT;
  const std::string grid_dgf;
  const std::string error_log;
  const std::string paraview;
  const double eps;
  const std::vector<double> bdf_alphas;
  const std::vector<double> bdf_gammas;
  // bdf coefficients ordered w.r.t. the polynomial order: from X^0 to X^n
  const double tg_a;
  const double tg_b;
  const double tg_Dc;
  const double tg_gamma;
  const double velocity_regularization;
  const double surface_growthFactor;
  const double mcf_regularization;
  const double u_hom_value;
  const double w_hom_value;
  const double u_pertubation;
  const double w_pertubation;
  const std::string u_init_dof;
  const std::string w_init_dof;
  Data();
};

Esfem::Io::Parameter::Data::Data()
  : t_0 {Dune::Fem::Parameter::getValue<double>("heat.starttime", 0.0)},
  t_end {Dune::Fem::Parameter::getValue<double>("heat.endtime", 0.6)},
  t_end_pattern {Dune::Fem::Parameter::getValue<double>("heat.pattern.endtime", t_end)},
  dT {Dune::Fem::Parameter::getValue<double>("heat.timestep",0.1)},
  grid_dgf {Esfem::Impl::get_macroGrid()},
  error_log {Dune::Fem::Parameter::getValue<std::string>
      ("fem.io.errorFile", Esfem::Impl::project_dir() + "output/l2h1error")},
  paraview {Dune::Fem::Parameter::getValue<std::string> 
      ("fem.io.outputName", Esfem::Impl::project_dir() + 
       "output/ALE_LiteratureExample-")},
  eps {Dune::Fem::Parameter::getValue<double>("heat.solvereps", 1e-8)},
  bdf_alphas {NUMERIK::bdf_alpha_coeff
      (Dune::Fem::Parameter::getValue<int>("heat.bdf", 1))},
  bdf_gammas {NUMERIK::bdf_gamma_coeff
      (Dune::Fem::Parameter::getValue<int>("heat.bdf", 1))},
  tg_a {Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.a", .1)},
  tg_b {Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.b", .9)},
  tg_Dc {Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.Dc", 10.)},
  tg_gamma {Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.gamma", 30.)},
  velocity_regularization
  {Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.alpha", 1e-3)},
  surface_growthFactor
  {Dune::Fem::Parameter::getValue<double>("tumor_growth.surface.delta", .4)},
  mcf_regularization
  {Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.epsilon", .01)},
  u_hom_value {Dune::Fem::Parameter::getValue<double>
      ("tumor_growth.heat.u_hom", tg_a + tg_b)},
  w_hom_value {Dune::Fem::Parameter::getValue<double>
      ("tumor_growth.heat.w_hom", tg_b  / ( u_hom_value * u_hom_value ) )},
  u_pertubation {Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.u_pertubation", .01)},
  w_pertubation {Dune::Fem::Parameter::getValue<double>("tumor_growth.heat.w_pertubation", .01)},
  u_init_dof {Dune::Fem::Parameter::getValue<std::string>
      ("tumor_growth.io.u_init_dof", Esfem::Impl::project_dir() + "output/u_dof.log")},
  w_init_dof {Dune::Fem::Parameter::getValue<std::string>
      ("tumor_growth.io.w_init_dof", Esfem::Impl::project_dir() + "output/w_dof.log")}
{
  Assert::dynamic<Assert::level(1), Esfem::Parameter_error>
    (eps > 0, Assert::compose(__FILE__, __LINE__, "Non positive tolerance."));
  Assert::dynamic<Assert::level(1), Esfem::Parameter_error>
    ( dT > eps, Assert::compose(__FILE__, __LINE__, "Time step too small."));
  Assert::dynamic<Assert::level(1), Esfem::Parameter_error>
    ((t_0 <= t_end_pattern) && (t_end_pattern <= t_end),
     Assert::compose
     (__FILE__, __LINE__, "Something is wrong with heat.starttime, "
      "heat.endtime or heat.pattern.starttime."));
  Esfem::Impl::file_check({grid_dgf, error_log, u_init_dof, w_init_dof});
}

Esfem::Io::Parameter::Parameter(int argc, char** argv,
				const std::string& parameter_file_name){
  Esfem::Impl::dune_fem_parameter_append(argc, argv, parameter_file_name);
  d_ptr =  make_unique<Data>();

  // const auto bdf_no = Dune::Fem::Parameter::getValue<int>("heat.bdf", 1);
  // check_bdfNo(bdf_no);

  // const auto t0 = Parameter::getValue<double>("heat.starttime", 0.0);
  // const auto t_end = Parameter::getValue<double>("heat.endtime", 0.6);
  // const auto t0_pat = Parameter::getValue<double>("heat.pattern.starttime", t_end);
  // check t0, t0_pattern, t_end

  // const auto a = Parameter::getValue<double>("tumor_growth.heat.a", .1);
  // const auto b = Parameter::getValue<double>("tumor_growth.heat.b", .9);
  // const auto u_hom = Parameter::getValue<double>("tumor_growth.heat.u_hom",
  // 						 a + b);
  // const auto w_hom = Parameter::getValue<double>("tumor_growth.heat.w_hom",
  // 						 b  / ( u_hom * u_hom ) );
  
  // d_ptr =  make_unique<Data>();
    // (Data{
    //   t0,
    // 	Parameter::getValue<double>("heat.timestep",0.1),
    // 	t0_pat,
    // 	t_end,
    // 	Esfem::Impl::get_macroGrid(),
    // 	Parameter::getValue<std::string>
    // 	("fem.io.errorFile", Esfem::Impl::project_dir() + "output/l2h1error"),
    // 	Parameter::getValue<std::string> 
    // 	("fem.io.outputName", Esfem::Impl::project_dir() + 
    // 	 "output/ALE_LiteratureExample-"),
    // 	Dune::Fem::Parameter::getValue<double>("heat.solvereps", 1e-8),
    // 	NUMERIK::bdf_alpha_coeff(bdf_no),
    // 	NUMERIK::bdf_gamma_coeff(bdf_no),
    // 	a,
    // 	b,
    // 	Parameter::getValue<double>("tumor_growth.heat.Dc", 10.),
    // 	Parameter::getValue<double>("tumor_growth.heat.gamma", 30.),
    // 	u_hom,
    // 	w_hom,
    // 	Parameter::getValue<double>("tumor_growth.heat.u_pertubation", .01),
    // 	Parameter::getValue<double>("tumor_growth.heat.w_pertubation", .01),
    // 	Parameter::getValue<std::string>("tumor_growth.io.u_init_dof",
    // 					 Esfem::Impl::project_dir() + "output/u_dof.log"),
    // 	Parameter::getValue<std::string>("tumor_growth.io.w_init_dof",
    // 					 Esfem::Impl::project_dir() + "output/w_dof.log")
    // 	}
    //   );
}
Esfem::Io::Parameter::~Parameter() = default;
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
double Esfem::Io::Parameter::global_timeStep() const noexcept{
  return d_ptr -> dT;
}
long Esfem::Io::Parameter::max_timeSteps() const{
  const auto t_n = d_ptr->t_end;
  const auto t_0 = d_ptr->t_0;
  return (t_n - t_0)/d_ptr->dT + .1;
}
long Esfem::Io::Parameter::prePattern_timeSteps() const{
  const auto t_n = d_ptr->t_end_pattern;
  const auto t_0 = d_ptr->t_0;
  return (t_n - t_0)/d_ptr->dT + .1;
}
long Esfem::Io::Parameter::pattern_timeSteps() const{
  const auto t_n = d_ptr->t_end;
  const auto t_0 = d_ptr->t_end_pattern;
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
double Parameter::velocity_regularization() const noexcept{
  return d_ptr -> velocity_regularization;
}
double Parameter::surface_growthFactor() const noexcept{
  return d_ptr -> surface_growthFactor;
}
double Parameter::mcf_regularization() const noexcept{
  return d_ptr -> mcf_regularization;
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
  const auto& p = d.d_ptr;
  os << "t_0: " << p->t_0 << '\n'
     << "dT: " << p->dT << '\n'
     << "t_end_pattern: " << p-> t_end_pattern << '\n'
     << "t_end: " << p->t_end << '\n'
     << "grid_dgf: " << p->grid_dgf << '\n'
     << "error_log: " << p->error_log << '\n'
     << "paraview: " << p->paraview << '\n'
     << "eps: " << p->eps << '\n'
     << "bdf_alphas: " << Esfem::Impl::doubleVector_to_string(p->bdf_alphas) << '\n'
     << "bdf_gammas: " << Esfem::Impl::doubleVector_to_string(p->bdf_gammas) << '\n'
     << "tg_a: " << p->tg_a << '\n'
     << "tg_b: " << p->tg_b << '\n'
     << "tg_Dc: " << p->tg_Dc << '\n'
     << "tg_gamma: " << p->tg_gamma << '\n'
     << "velocity_regularization: " << p->velocity_regularization << '\n'
     << "surface_growthFactor: " << p->surface_growthFactor << '\n'
     << "mcf_regularization: " << p->mcf_regularization << '\n'
     << "u_hom_value: " << p -> u_hom_value << '\n'
     << "w_hom_value: " << p -> w_hom_value << '\n'
     << "u_pertubation: " << p -> u_pertubation << '\n'
     << "w_pertubation: " << p -> w_pertubation << '\n'
     << "u_init_dof: " << p -> u_init_dof << '\n'
     << "w_init_dof: " << p -> w_init_dof
    ;
  return os;
}
