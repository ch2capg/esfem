/*! \file brusselator_algo.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Februar 2016

     Implementation details for brusselator_algo.h
     Created by Christian Power on 04.02.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <config.h>
#include "brusselator_algo.h"

using namespace Esfem;

using Scal_FEfun_set = FEfun_set<Esfem::Grid::Scal_FEfun>;
using Vec_FEfun_set = FEfun_set<Esfem::Grid::Vec_FEfun>;

void clog_uw(const Scal_FEfun_set& u, const Scal_FEfun_set& w);

void brusselator_algo(int argc, char** argv){
  Dune::Fem::MPIManager::initialize(argc, argv);

  const auto parameter_file =
    "/home/power/cpp/DISS_surfaces/data/tumor_parameter.txt";  

  FEM_data fem {argc, argv, parameter_file};

  pre_loop_action(fem);
  fem.next_timeStep();
  for(long it = 0; it < fem.prePattern_timeSteps(); ++it){
    pre_pattern_action(fem);
    fem.next_timeStep();
  }
  intermediate_action(fem);
  // fem.next_timeStep(); // Do we need this?
  for(long it = 0; it < fem.pattern_timeSteps(); ++it){
    pattern_action(fem);
    fem.next_timeStep();
  }
  final_action(fem);

  Io::Parameter data {argc, argv, parameter_file};
  const Io::Dgf::Handler dgf_interpreter {data.grid()};
  const Init_data init_data {data};
  Err_stream err_log {data};
  Grid::Grid_and_time grid {data};
  Scal_FEfun_set u {"u", grid};
  Scal_FEfun_set w {"w", grid};
  const Err_cal err_cal {grid, u, w};
  Solver solver {data, grid, u, w};
  Io::Paraview paraview_plot {data, grid, u.fun, w.fun};

  pre_loop_action(dgf_interpreter, init_data, err_log,
		  grid, u, w, err_cal, solver, paraview_plot);

  grid.next_timeStep(data.global_timeStep());
  for(long it =0; it < data.max_timeSteps(); ++it){
    std::cout << "step number " << it << std::endl;
    pre_pattern_action(data, dgf_interpreter, err_log, grid,
		       u, w, err_cal, solver, paraview_plot);
    grid.next_timeStep(data.global_timeStep());
  }

  u >> dgf_interpreter;
  // Io::dof_to_file(u.fun.cbegin(), u.fun.cend(), data.u_init_dof());
}

// ----------------------------------------------------------------------
// Implementation structs 
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

FEM_data::FEM_data(int argc, char** argv, const std::string& parameter_fname)
try : data {argc, argv, parameter_fname},
		 dgf_handler {data.grid()},
		 estream {data},
		 fix_grid {data},
		 u {"u", fix_grid},
		 w {"w", fix_grid}
{}
catch(const std::exception&){
  throw_with_nested(std::runtime_error {"Error in constructor of FEM_data."});
 }
 catch(...){
   throw std::runtime_error {"Unkown error in constructor of FEM_data."};
 }
void FEM_data::next_timeStep(){
  fix_grid.next_timeStep(data.global_timeStep());
}
long FEM_data::prePattern_timeSteps() const{
  return data.prePattern_timeSteps();
}
long FEM_data::pattern_timeSteps() const{
  return data.pattern_timeSteps();
}

// ----------------------------------------------------------------------
// helper functions

void first_interpolate(const Init_data& id,
		       Scal_FEfun_set& u,
		       Scal_FEfun_set& w){
  id.u.interpolate(u.fun);
  id.u.interpolate(u.app);
  id.u.interpolate(u.exact);
  id.w.interpolate(w.fun);
  id.w.interpolate(w.app);
  id.w.interpolate(w.exact);
}
void update_exactSolution(const Init_data& id, Scal_FEfun_set& u, Scal_FEfun_set& w){
  id.u.interpolate(u.exact);
  id.w.interpolate(w.exact);
}
void massMatrix_rhsLes(const Solver& s, Scal_FEfun_set& u, Scal_FEfun_set& w){
  s.u.mass_matrix(u.fun, u.rhs_les);
  s.w.mass_matrix(w.fun, w.rhs_les);  
}
void massMatrixConstOne_rhsLes(const Solver& s, Scal_FEfun_set& u, Scal_FEfun_set& w){
  s.u.add_massMatrixConstOne_to(u.rhs_les);
  s.w.add_massMatrixConstOne_to(w.rhs_les);
}
void assemble_and_addScaled_rhsLes(const Rhs& rhs, Scal_FEfun_set& u, Scal_FEfun_set& w){
  rhs.u.assemble_and_addScaled_to(u.rhs_les);
  rhs.w.assemble_and_addScaled_to(w.rhs_les);  
}
void solve_pde(const Solver& s, Scal_FEfun_set& u, Scal_FEfun_set& w){
  s.u.solve(u.rhs_les, u.fun);
  u.app = u.fun;    
  s.w.solve(w.rhs_les, w.fun);
  w.app = w.fun;  
}

// ----------------------------------------------------------------------
// loop action

void pre_loop_action(FEM_data&){
}
void pre_pattern_action(FEM_data&){
}
void intermediate_action(FEM_data&){
}
void pattern_action(FEM_data&){
}
void final_action(FEM_data&){
}

void pre_loop_action(const Esfem::Io::Dgf::Handler& h,
		     const Init_data& id,
		     Err_stream& es,
		     const Esfem::Grid::Grid_and_time& gt,
		     FEfun_set<Esfem::Grid::Scal_FEfun>& u,
		     FEfun_set<Esfem::Grid::Scal_FEfun>& w,
		     const Err_cal& ec, 
		     const Solver& s, 
		     Esfem::Io::Paraview& p
		     ){
  first_interpolate(id, u, w);

  generate_header_line(es);
  write_error_line(es, gt, ec);
  p.write();

  massMatrix_rhsLes(s, u, w);
  u >> h;
  w >> h;
  u.write(h, "./");
  w.write(h, "./");
}
void pre_pattern_action(const Esfem::Io::Parameter& data,
			const Esfem::Io::Dgf::Handler& h,
			Err_stream& es,
			const Esfem::Grid::Grid_and_time& gt,
			FEfun_set<Esfem::Grid::Scal_FEfun>& u,
			FEfun_set<Esfem::Grid::Scal_FEfun>& w,
			const Err_cal& ec,
			const Solver& s,
			Esfem::Io::Paraview& p
			){
  Grid::Grid_and_time grid_tmp {data};
  Scal_FEfun_set u_tmp {u, grid_tmp};
  u_tmp = u;
  Scal_FEfun_set w_tmp {w, grid_tmp};
  w_tmp = w;
  const Err_cal errCal_tmp {grid_tmp, u_tmp, w_tmp};
  Solver solver {data, grid_tmp, u_tmp, w_tmp};
  Io::Paraview paraview_plot {data, grid_tmp, u_tmp.fun, w_tmp.fun};

  massMatrixConstOne_rhsLes(s, u, w);
  // assemble_and_addScaled_rhsLes(rhs, u, w);	// testing rhs
  solve_pde(s, u, w);
  massMatrix_rhsLes(s, u, w);
  // update_exactSolution(init_data, u, w);	// testing rhs

  write_error_line(es, gt, ec);
  p.write();  
}

// ------------------------------------------------------------
// io helper functions

void generate_header_line(Err_stream& es){
  generate_header_line(es.u);
  generate_header_line(es.w);
}
void generate_header_line(Esfem::Io::Error_stream& es){
  es << "timestep" << "\t"
     << "L2err" << "\t\t"
     << "H1err" << std::endl;
  es << std::scientific;
}
void write_error_line(Err_stream& es,
		      const Esfem::Grid::Grid_and_time& gt,
		      const Err_cal& ec){
  write_error_line(es.u, gt, ec.u);
  write_error_line(es.w, gt, ec.w);
}
void write_error_line(Io::Error_stream& es,
		      const Esfem::Grid::Grid_and_time& g,
		      const Io::L2H1_calculator& cal){
  const auto& tp = g.time_provider();
  es << tp.deltaT() << '\t'
     << cal.l2_err() << '\t'
     << cal.h1_err() << std::endl; 
}

// ----------------------------------------------------------------------
// Internal implementation

void clog_uw(const Scal_FEfun_set& u, const Scal_FEfun_set& w){
  const Grid::Scal_FEfun::Dune_FEfun& u_fun = u.fun;
  const Grid::Scal_FEfun::Dune_FEfun& u_rhsLes = u.rhs_les;
  const Grid::Scal_FEfun::Dune_FEfun& w_fun = w.fun;
  const Grid::Scal_FEfun::Dune_FEfun& w_rhsLes = w.rhs_les;
  std::clog << *u_fun.dbegin() << '\t'
	    << *u_rhsLes.dbegin() << '\n'    
	    << *w_fun.dbegin() << '\t'
	    << *w_rhsLes.dbegin() << std::endl;  
}

/*! Log:
 */
