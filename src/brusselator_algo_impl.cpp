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

using Esfem::Rhs;
using Esfem::Init_data;
using Esfem::Scalar_solver;
using Esfem::Err_cal;
using Esfem::Err_stream;
using Esfem::PreLoop_helper;
using Esfem::RhsAndSolve_helper;
using Esfem::Helper_surface;
using Esfem::SecOrd_op::Identity;
using Scal_FEfun_set = Esfem::FEfun_set<Esfem::Grid::Scal_FEfun>;
using Vec_FEfun_set = Esfem::FEfun_set<Esfem::Grid::Vec_FEfun>;

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

Scalar_solver::Scalar_solver
(const Esfem::Io::Parameter& p,
 const Esfem::Grid::Grid_and_time& g,
 const Scal_FEfun_set& u_set,
 const Scal_FEfun_set& w_set)
  : u {p, g, Growth::promoting, u_set.app, w_set.app},
  w {p, g, Growth::inhibiting, u_set.app, u_set.app}
{}

Scalar_solver::Scalar_solver
(const Io::Parameter& p,
 const Grid::Grid_and_time& g,
 const Grid::Scal_tiny_FEfun_set& u_set,
 const Grid::Scal_tiny_FEfun_set& w_set)
  : u {p, g, Growth::promoting, u_set.fun, w_set.fun},
  w {p, g, Growth::inhibiting, u_set.fun, u_set.fun}
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
// Implementation PreLoop_helper

PreLoop_helper::PreLoop_helper(Brusselator_scheme& bs_ref)
  : bs {bs_ref},
  io {bs.io},
  fef {bs.fef},
  init_data {bs.data},
  err_cal {bs.fix_grid, fef.u, fef.w},
  paraview {bs.data, bs.fix_grid, fef.u.fun, fef.w.fun},
  solver {bs.data, bs.fix_grid, fef.u, fef.w}
{}

void PreLoop_helper::first_interpolate(){
  interpolate(init_data.u, fef.u);
  interpolate(init_data.w, fef.w);
  bs.identity.interpolate(fef.surface.fun);
}
void PreLoop_helper::headLine_in_errFile(){
  head_line(io.u);
  head_line(io.w);
}
void PreLoop_helper::plot_errors_in_errFile(){
  const auto& g = bs.fix_grid;
  write_error_line(io.u, g.time_provider(), err_cal.u);
  write_error_line(io.w, g.time_provider(), err_cal.w);
}
void PreLoop_helper::plot_paraview(){
  paraview.write();
}
void PreLoop_helper::prepare_rhs(){
  auto& u = fef.u;
  auto& w = fef.w;
  solver.u.mass_matrix(u.fun, u.rhs_les);
  solver.w.mass_matrix(w.fun, w.rhs_les);
}

// ----------------------------------------------------------------------
// Implementation PrePattern_helper 

PrePattern_helper::PrePattern_helper(Brusselator_scheme& bs_input)
  : bs {bs_input},
  io {bs.io},
  fef {bs.fef},
  err_cal {bs.fix_grid, fef.u, fef.w},
  paraview {bs.data, bs.fix_grid, fef.u.fun, fef.w.fun},
  solver {bs.data, bs.fix_grid, fef.u, fef.w}
{}

void PrePattern_helper::finalize_rhs(){
  auto& u_rhs = fef.u.rhs_les;
  auto& w_rhs = fef.w.rhs_les;
  solver.u.add_massMatrixConstOne_to(u_rhs);
  solver.u.add_massMatrixConstOne_to(w_rhs);
}
void PrePattern_helper::solve_pde(){
  auto& u = fef.u;
  auto& w = fef.w;
  solver.u.solve(u.rhs_les, u.fun);
  u.app = u.fun;    
  solver.w.solve(w.rhs_les, w.fun);
  w.app = w.fun;    
}
void PrePattern_helper::prepare_rhs(){
  auto& u = fef.u;
  auto& w = fef.w;
  solver.u.mass_matrix(u.fun, u.rhs_les);
  solver.w.mass_matrix(w.fun, w.rhs_les);
}
void PrePattern_helper::plot_errors_in_errFile(){
  const auto& g = bs.fix_grid;
  write_error_line(io.u, g.time_provider(), err_cal.u);
  write_error_line(io.w, g.time_provider(), err_cal.w);
}
void PrePattern_helper::plot_paraview(){
  paraview.write();
}

// ----------------------------------------------------------------------
// Implementation RhsAndSolve_helper

RhsAndSolve_helper::RhsAndSolve_helper(Brusselator_scheme& bs)
  : fef {bs.fef},
  grid {bs.data, compose_dgfName(fef.surface.fun.name()), 
      bs.fix_grid.time_provider().time()},
  u {fef.u, grid},
  w {fef.w, grid},
  X {fef.surface, grid},
  ss  {data, grid, u, w},
  vs {data, grid, u.fun}
{}
RhsAndSolve_helper::scalar_massMatrix(){
  ss.u.mass_matrix(u.rhs_les);
  ss.w.mass_matrix(w.rhs_les);
  fef.u.rhs_les = u.rhs_les;
  fef.w.rhs_les = w.rhs_les;
}
RhsAndSolve_helper::solve_surface_and_save(){
  vs.rhs(X.fun, X.rhs_les);
  vs.solve(X.rhs_les, X.fun);
  fef.surface.fun = X.fun;
  fef.surface.write(io.dgf_handler, "./");     
}

// ----------------------------------------------------------------------
// Implementation Helper_surface

Helper_surface::
Helper_surface(const Scal_FEfun_set& u_input, Vec_FEfun_set& surface_input,
	       const Io::Parameter& p, const Grid::Grid_and_time& g)
  : u {u_input},
  surface {surface_input},
  gt {p, compose_dgfName(surface_input.fun.name()), g.time_provider().time()},
  surface_loc {"surface_loc", g},
  solver {p, gt, u}
{}
void Helper_surface::solve_for_surface(){
}

// ----------------------------------------------------------------------
// Implementation Helper_uw

Helper_uw::Helper_uw(Scal_FEfun_set& u_input, Scal_FEfun_set& w_input,
		     Vec_FEfun_set& surface, const Io::Parameter& p,
		     const Grid::Grid_and_time& g)
  : u {u_input},
  w {w_input},
  gt {p, compose_dgfName(surface.fun.name()), g.time_provider().time()},
  u_loc {"u_loc", g},
  w_loc {"w_loc", g},
  errCal {p, u_loc, w_loc},
  solver {p, gt, u_loc, w_loc},
  paraview {p, g, u_loc.fun, w_loc.fun}
{}
void Helper_uw::solve_pde(){
}
void Helper_uw::write_error_line(){
}
void Helper_uw::paraview_plot(){
}

// ----------------------------------------------------------------------
// helper functions

void Esfem::interpolate(const SecOrd_op::Init_data& id, Scalar_FEfun_set& f){
  id.interpolate(f.fun);
  f.app = f.fun;
  f.exact = f.fun;
}
void Esfem::head_line(Io::Error_stream& file){
  es << "timestep" << "\t"
     << "L2err" << "\t\t"
     << "H1err" << std::endl;
  es << std::scientific;
}
void Esfem::write_error_line(Io::Error_stream& file,
			     const Dune::Fem::TimeProviderBase& tp,
			     const Io::L2H1_calculator& cal){
  file << tp.deltaT() << '\t'
       << cal.l2_err() << '\t'
       << cal.h1_err() << std::endl; 
}

// ----------------------------------------------------------------------
// old code

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
