/*! \file brusselator_algo_impl.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Implementation details for brusselator_algo_impl.h

     Created by Christian Power on 16.03.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
*/

#include "brusselator_algo_impl.h"

using namespace Esfem;
using Scal_FEfun_set = FEfun_set<Esfem::Grid::Scal_FEfun>;
using Vec_FEfun_set = FEfun_set<Esfem::Grid::Vec_FEfun>;
using Identity = SecOrd_op::Identity;

// ----------------------------------------------------------------------
// Implementation of structs 

Rhs::Rhs(const Io::Parameter& p, const Grid::Grid_and_time& gt)
  : u {p, gt}, w {p, gt}
{}

// Init_data::Init_data(const Grid::Grid_and_time& gt)
//   : u {gt}, w {gt}
// {}
Init_data::Init_data(const Esfem::Io::Parameter& p)
  : u {p, Growth::promoting},
  w {p, Growth::inhibiting}
{}

Solver::Solver(const Esfem::Io::Parameter& p, const Esfem::Grid::Grid_and_time& g,
	       const Scal_FEfun_set& u_set,
	       const Scal_FEfun_set& w_set)
  : u {p, g, Growth::promoting, u_set.app, w_set.app},
  w {p, g, Growth::inhibiting, u_set.app, u_set.app}
{}

Err_cal::Err_cal(const Esfem::Grid::Grid_and_time& g,
		 const Scal_FEfun_set& u_set,
		 const Scal_FEfun_set& w_set)
  : u {g, u_set.exact, u_set.fun},
  w {g, w_set.exact, w_set.fun}
{}

Err_stream::Err_stream(const Esfem::Io::Parameter& p)
  : u {"_u", p}, w {"_w", p}
{}

// ----------------------------------------------------------------------
// helper functions

void Esfem::first_interpolate(const Identity& identity,
			      const Init_data& id,
			      Scal_FEfun_set& u,
			      Scal_FEfun_set& w,
			      Vec_FEfun_set& surface){
  id.u.interpolate(u.fun);
  id.u.interpolate(u.app);
  id.u.interpolate(u.exact);
  id.w.interpolate(w.fun);
  id.w.interpolate(w.app);
  id.w.interpolate(w.exact);
  identity.interpolate(surface.fun);
}
void Esfem::update_exactSolution(const Init_data& id, Scal_FEfun_set& u, Scal_FEfun_set& w){
  id.u.interpolate(u.exact);
  id.w.interpolate(w.exact);
}
void Esfem::massMatrix_rhsLes(const Solver& s, Scal_FEfun_set& u, Scal_FEfun_set& w){
  s.u.mass_matrix(u.fun, u.rhs_les);
  s.w.mass_matrix(w.fun, w.rhs_les);  
}
void Esfem::massMatrixConstOne_rhsLes(const Solver& s, Scal_FEfun_set& u, Scal_FEfun_set& w){
  s.u.add_massMatrixConstOne_to(u.rhs_les);
  s.w.add_massMatrixConstOne_to(w.rhs_les);
}
void Esfem::assemble_and_addScaled_rhsLes(const Rhs& rhs, Scal_FEfun_set& u, Scal_FEfun_set& w){
  rhs.u.assemble_and_addScaled_to(u.rhs_les);
  rhs.w.assemble_and_addScaled_to(w.rhs_les);  
}
void Esfem::solve_pde(const Solver& s, Scal_FEfun_set& u, Scal_FEfun_set& w){
  s.u.solve(u.rhs_les, u.fun);
  u.app = u.fun;    
  s.w.solve(w.rhs_les, w.fun);
  w.app = w.fun;  
}

// ------------------------------------------------------------
// io helper functions

void Esfem::generate_header_line(Err_stream& es){
  generate_header_line(es.u);
  generate_header_line(es.w);
}
void Esfem::generate_header_line(Esfem::Io::Error_stream& es){
  es << "timestep" << "\t"
     << "L2err" << "\t\t"
     << "H1err" << std::endl;
  es << std::scientific;
}
void Esfem::write_error_line(Err_stream& es,
			     const Esfem::Grid::Grid_and_time& gt,
			     const Err_cal& ec){
  write_error_line(es.u, gt, ec.u);
  write_error_line(es.w, gt, ec.w);
}
void Esfem::write_error_line(Io::Error_stream& es,
			     const Esfem::Grid::Grid_and_time& g,
			     const Io::L2H1_calculator& cal){
  const auto& tp = g.time_provider();
  es << tp.deltaT() << '\t'
     << cal.l2_err() << '\t'
     << cal.h1_err() << std::endl; 
}
std::string Esfem::compose_dgfName(const std::string& fun_name,
				   const std::string& dir){
  constexpr auto suffix= ".dgf";
  return dir + fun_name + suffix;
}
void Esfem::clog_uw(const Scal_FEfun_set& u, const Scal_FEfun_set& w){
  const Grid::Scal_FEfun::Dune_FEfun& u_fun = u.fun;
  const Grid::Scal_FEfun::Dune_FEfun& u_rhsLes = u.rhs_les;
  const Grid::Scal_FEfun::Dune_FEfun& w_fun = w.fun;
  const Grid::Scal_FEfun::Dune_FEfun& w_rhsLes = w.rhs_les;
  std::clog << *u_fun.dbegin() << '\t'
	    << *u_rhsLes.dbegin() << '\n'    
	    << *w_fun.dbegin() << '\t'
	    << *w_rhsLes.dbegin() << std::endl;  
}

