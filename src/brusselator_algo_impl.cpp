/*! \file brusselator_algo_impl.cpp
    \brief Implementation details for brusselator_algo_impl.h

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Idea
     --------------------------------------------------

     Providing many struct, such that functions like initial value,
     right-hand side function or PDE solver do not have to be specified for 
     \f$ u\f$ and \f$ w\f$ separately.  

    \author Christian Power
    \date 23. April 2016
    \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
*/

#include "brusselator_algo_impl.h"
#include "esfem_error.h"

/*! \name Convenience typedefs */
//@{
using Bruss_error = Esfem::BrusselatorScheme_error;
using Esfem::SecOrd_op::Identity;
using Scal_FEfun_set = Esfem::Grid::Scal_FEfun_set;
using Vec_FEfun_set = Esfem::Grid::Vec_FEfun_set;
//@}

/*! \name Implemented in this file */
//@{
using Esfem::Rhs;
using Esfem::Init_data;
using Esfem::Scalar_solver;
using Esfem::Err_cal;
using Esfem::PrePattern_helper;
using Esfem::PreLoop_helper;
using Esfem::RhsAndSolve_helper;
using Esfem::Pattern_helper;
//@}

// ----------------------------------------------------------------------
// Implementation of structs 

Rhs::Rhs(const Grid::Grid_and_time& gt)
  : u {gt, Growth::promoting}, w {gt, Growth::inhibiting}
{}

// Init_data::Init_data(const Grid::Grid_and_time& gt)
//   : u {gt}, w {gt}
// {}
Init_data::Init_data(const Esfem::Io::Parameter& p)
  : u {p, Growth::promoting}, w {p, Growth::inhibiting}
{}

Scalar_solver::Scalar_solver
(const Esfem::Io::Parameter& p,
 const Esfem::Grid::Grid_and_time& g,
 const Esfem::Grid::Scal_FEfun_set& u_set,
 const Esfem::Grid::Scal_FEfun_set& w_set)
  :u {p, g, Growth::promoting, u_set.app, w_set.app},
  w {p, g, Growth::inhibiting, u_set.app, u_set.app}
{}

Scalar_solver::Scalar_solver
(const Esfem::Io::Parameter& p,
 const Esfem::Grid::Grid_and_time& g,
 const Esfem::Grid::Scal_tiny_FEfun_set& u_set,
 const Esfem::Grid::Scal_tiny_FEfun_set& w_set)
  : u {p, g, Growth::promoting, u_set.fun, w_set.fun},
  w {p, g, Growth::inhibiting, u_set.fun, u_set.fun}
{}

Err_cal::Err_cal(const Esfem::Grid::Grid_and_time& g,
		 const Scal_FEfun_set& u_set,
		 const Scal_FEfun_set& w_set)
  : u {g, u_set.exact, u_set.fun},
  w {g, w_set.exact, w_set.fun}
{}

// Err_stream::Err_stream(const Esfem::Io::Parameter& p)
//   : u {"_u", p}, w {"_w", p}
// {}

// ----------------------------------------------------------------------
// Implementation PreLoop_helper

PreLoop_helper::PreLoop_helper(Brusselator_scheme& bs_input)
  : bs {bs_input},
  init_data {bs.data},
  err_cal {bs.fix_grid, bs.fef.u, bs.fef.w},
  paraview {bs.data, bs.fix_grid, bs.fef.u.fun, bs.fef.w.fun},
  solver {bs.data, bs.fix_grid, bs.fef.u, bs.fef.w}
{}

void PreLoop_helper::first_interpolate(){
  interpolate(init_data.u, bs.fef.u);
  interpolate(init_data.w, bs.fef.w);
  bs.io.identity.interpolate(bs.fef.surface.fun);
}
void PreLoop_helper::headLine_in_errFile(){
  head_line(bs.io.u);
  head_line(bs.io.w);
}
void PreLoop_helper::plot_errors_in_errFile(){
  const auto& tp = bs.fix_grid.time_provider();
  write_error_line(bs.io.u, tp, err_cal.u);
  write_error_line(bs.io.w, tp, err_cal.w);
}
void PreLoop_helper::prepare_rhs(){
  auto& u = bs.fef.u;
  auto& w = bs.fef.w;
  solver.u.mass_matrix(u.fun, u.rhs_les);
  solver.w.mass_matrix(w.fun, w.rhs_les);
}

// ----------------------------------------------------------------------
// Implementation PrePattern_helper 

PrePattern_helper::PrePattern_helper(Brusselator_scheme& bs)
  : io {bs.io},
  u {bs.fef.u},
  w {bs.fef.w},
  tp {bs.fix_grid.time_provider()},
  err_cal {bs.fix_grid, u, w},
  paraview {bs.data, bs.fix_grid, u.fun, w.fun},
  solver {bs.data, bs.fix_grid, u, w}
{}

// void PrePattern_helper::finalize_rhs(){
void PrePattern_helper::rhs(){
  solver.u.mass_matrix(u.fun, u.rhs_les);
  solver.w.mass_matrix(w.fun, w.rhs_les);
  solver.u.add_massMatrixConstOne_to(u.rhs_les);
  solver.w.add_massMatrixConstOne_to(w.rhs_les);
}
void PrePattern_helper::solve_pde(){
  solver.u.solve(u.rhs_les, u.fun);
  u.app = u.fun;    
  solver.w.solve(w.rhs_les, w.fun);
  w.app = w.fun;    
}
void PrePattern_helper::prepare_rhs(){
  solver.u.mass_matrix(u.fun, u.rhs_les);
  solver.w.mass_matrix(w.fun, w.rhs_les);
}
void PrePattern_helper::plot_errors_in_errFile(){
  write_error_line(io.u, tp, err_cal.u);
  write_error_line(io.w, tp, err_cal.w);
}

// ----------------------------------------------------------------------
// Implementation RhsAndSolve_helper

RhsAndSolve_helper::RhsAndSolve_helper(Brusselator_scheme& bs_input)
try :bs {bs_input},
  fef {bs.fef},
  grid {bs.data,
      Grid::compose_dgfName(fef.surface.fun.name(), fef.tmpFile_path), 
      bs.fix_grid.time_provider().time()},
  u {fef.u, grid},
  w {fef.w, grid},
  X {fef.surface, grid},
  ss {bs.data, grid, u, w},
  vs {bs.data, grid, u.fun},
  v_rhs {bs.fix_grid}
{}
catch(std::exception&){
  std::throw_with_nested(Bruss_error {"RhsAndSolve_helper()"});
 }
 catch(...){
   throw Bruss_error {"RhsAndSolve_helper(), unknown error"};
 }
void RhsAndSolve_helper::scalar_massMatrix(){
  ss.u.mass_matrix(u.fun, u.rhs_les);
  ss.w.mass_matrix(w.fun, w.rhs_les);
  fef.u.rhs_les = u.rhs_les;
  fef.w.rhs_les = w.rhs_les;
}
void RhsAndSolve_helper::brusselator_rhs(){
  vs.brusselator_rhs(X.fun, X.rhs_les);	
}
void RhsAndSolve_helper::addScaled_surfaceLoadVector(){
  v_rhs.assemble_and_addScaled_to(X.rhs_les);
}
void RhsAndSolve_helper::solve_surface_and_save(){
  vs.solve(X.rhs_les, X.fun);
  fef.surface.fun = X.fun;
  fef.surface.write(bs.io.dgf_handler, fef.tmpFile_path);
}

// ----------------------------------------------------------------------
// Implementation Pattern_helper

Pattern_helper::Pattern_helper(Brusselator_scheme& bs_input)
  :bs {bs_input},
   grid
  {bs.data,
      Grid::compose_dgfName(bs.fef.surface.fun.name(), bs.fef.tmpFile_path ), 
      bs.fix_grid.time_provider().time()},
   u {bs.fef.u, grid},
   w {bs.fef.w, grid},
   err_cal {grid, u, w},
   paraview {bs.data, grid, u.fun, w.fun},
   solver {bs.data, grid, u, w}
{}
void Pattern_helper::finalize_scalarPDE_rhs(){
  solver.u.add_massMatrixConstOne_to(u.rhs_les);
  // rhs.u.add_scaled(u.rhs_les);
  solver.w.add_massMatrixConstOne_to(w.rhs_les);
  // rhs.w.add_scaled(w.rhs_les);
}
void Pattern_helper::solve_scalarPDE(){
  solver.u.solve(u.rhs_les, u.fun);
  u.app = u.fun;    
  solver.w.solve(w.rhs_les, w.fun);
  w.app = w.fun;
  bs.fef.u = u;
  bs.fef.w = w;
}
void Pattern_helper::plot_errors_in_errFile(){
  const auto& tp = grid.time_provider();
  write_error_line(bs.io.u, tp, err_cal.u);
  write_error_line(bs.io.w, tp, err_cal.w);  
}

// ----------------------------------------------------------------------
// helper functions

void Esfem::interpolate(const SecOrd_op::Init_data& id, Scal_FEfun_set& f){
  id.interpolate(f.fun);
  f.app = f.fun;
  f.exact = f.fun;
}
void Esfem::head_line(Io::Error_stream& file){
  file << "timestep" << "\t"
       << "L2err" << "\t\t"
       << "H1err" << std::endl;
  file << std::scientific;
}
void Esfem::write_error_line(Io::Error_stream& file,
			     const Dune::Fem::TimeProviderBase& tp,
			     const Io::L2H1_calculator& cal){
  file << tp.deltaT() << '\t'
       << cal.l2_err() << '\t'
       << cal.h1_err() << std::endl; 
}
