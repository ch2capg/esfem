/*! \file brusselator_algo.cpp
    \brief Implementation for `brusselator_algo.h`

     Revision history
     --------------------------------------------------

          Revised by Christian Power June 2016
          Revised by Christian Power May 2016
          Revised by Christian Power April 2016
          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Idea
     --------------------------------------------------
     
     Implementation of Esfem::brusselator_algo() and
     the class `Esfem::Brusselator_scheme`

     \author Christian Power 
     \date 16. June 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
*/

#include <config.h>
#include "brusselator_algo.h"
#include "brusselator_algo_impl.h"
#include "secOrd_op_solutionDriven.h"
#include "esfem_error.h"

#ifndef PFILE
#error Give a complete path for parameter file in macro variable PFILE.
#endif

using Bruss_error = Esfem::BrusselatorScheme_error;
/*!< \brief Shorter name for convenience */
using Esfem::Brusselator_scheme;
using Scal_FEfun_set = Esfem::Grid::Scal_FEfun_set;
/*!< \brief Four functions of type \f$ f\colon \R^3 \to \R \f$ */
using Vec_FEfun_set = Esfem::Grid::Vec_FEfun_set;
/*!< \brief Four functions of type \f$ f\colon \R^3 \to \R^3 \f$ */
using Vector_solver = Esfem::SecOrd_op::Solution_driven;
/*!< \brief Solver for the `X` respectively `surface` variable */

void Esfem::brusselator_algo(int argc, char** argv){
  Dune::Fem::MPIManager::initialize(argc, argv);
  constexpr auto parameter_file = PFILE;
  Brusselator_scheme fem {argc, argv, parameter_file};
  // fem.prePattern_loop();
  // fem.intermediate_action(); 
  // fem.pattern_loop();
  // fem.final_action();
  // fem.standard_esfem(); // c++ code works flawless
  // fem.eoc_logisticSphere(); // does not work
  // fem.eoc_mcf(); // code works
  // fem.eoc_sls(); // code works
  // fem.sd(); // code works
  fem.eoc_sdp(); // not finished
}

// ----------------------------------------------------------------------
// Brusselator_scheme implementation

// ------------------------------------------------------------
// Brusselator_scheme public 

/*! \todo Test if folder for tmp file exists. */
Brusselator_scheme::
Brusselator_scheme(int argc, char** argv,
		   const std::string& parameter_fname)
try :data {argc, argv, parameter_fname},
  io {data},
  fix_grid {data},
  norm {fix_grid},  
  fef {fix_grid}, 
  // exact {fix_grid, data} // constructor for random initial data
  exact {fix_grid} // constructor for analytic initial data
{
  pre_loop_action(); // initialize member fef
}
catch(...){
  throw_with_nested(Bruss_error {"Constructor."});    
 }

void Brusselator_scheme::standard_esfem(){
  Rhs load_vector {fix_grid};
  Scalar_solver solver {data, fix_grid, fef.u, fef.w};

  error_on_intSurface(); // error on t_0
  for(long it = 0; it < time_steps(); ++it){
    // rhs = M^n u^n
    solver.u.mass_matrix(fef.u.fun, fef.u.rhs_les);
    solver.w.mass_matrix(fef.w.fun, fef.w.rhs_les);
    next_timeStep(); // next surface

    // rhs += M^{n+1} 1
    solver.u.add_massMatrixConstOne_to(fef.u.rhs_les);
    solver.w.add_massMatrixConstOne_to(fef.w.rhs_les);
    // rhs += load_vector
    load_vector.u.assemble_and_addScaled_to(fef.u.rhs_les);
    load_vector.w.assemble_and_addScaled_to(fef.w.rhs_les);

    solver.u.solve(fef.u.rhs_les, fef.u.fun);
    fef.u.app = fef.u.fun;    
    solver.w.solve(fef.w.rhs_les, fef.w.fun);
    fef.w.app = fef.w.fun;

    exact.u.interpolate(fef.u.exact);
    exact.w.interpolate(fef.w.exact);
    error_on_intSurface(); // error on t_n+1
  }
}

void Brusselator_scheme::eoc_logisticSphere(){  
  io.identity.interpolate(fef.surface.fun);

  for(long it = 0; it < pattern_timeSteps(); ++it){
    Grid::Grid_and_time grid {data,
	Grid::compose_dgfName(fef.surface.fun.name(), fef.tmpFile_path), 
	fix_grid.time_provider().time()};
    Grid::Scal_FEfun_set u {fef.u, grid};
    Grid::Vec_FEfun_set X {fef.surface, grid};
    SecOrd_op::Init_data u_init {grid, Growth::promoting}; // perhaps not efficient
    SecOrd_op::Vec_rhs X_loadVector {grid};
    SecOrd_op::Solution_driven X_solver {data, grid, u.app};
    
    // calcul1ate X.fun
    u_init.interpolate(u.app);
    X.app = fef.surface.fun;
    X_solver.brusselator_rhs(X.app, X.rhs_les); 
    // X_loadVector.assemble_and_addScaled_to(X.rhs_les);

    X_solver.solve(X.rhs_les, X.fun);

    next_timeStep();
    
    // save surface
    fef.surface.fun = X.fun; // swap would be more efficient
    fef.surface.write(io.dgf_handler, fef.tmpFile_path);    

    // calculate error 
    io.identity.interpolate(fef.surface.exact);    
    io.surface << fix_grid.time_provider().deltaT() << ' '
	       << norm.l2_err(fef.surface.fun, fef.surface.exact) << ' '
	       << norm.h1_err(fef.surface.fun, fef.surface.exact) 
	       << std::endl;
    // I've changed order to get a stationary surface
    // io.identity.interpolate(fef.surface.fun);
    // fef.surface.write(io.dgf_handler, fef.tmpFile_path);


    
    // update_surface(); // calculate exact surface
    // update_scalar_solution(); // on the exact surface 
    // rhs_and_solve_SPDE(); // also scalar rhs
    // update_velocity();    // exact and approximation 
    // error_on_intSurface(); // Error on surface(t_n)
    // next_timeStep(); // next surface
    // Pattern_helper helper {*this};
    // helper.finalize_scalarPDE_rhs();
    // helper.solve_scalarPDE();
    // helper.errors_on_numSurface();
    // helper.plot_paraview();
  }
}
/*! \pre delta == 0, epsilon == 1, alpha == 0. */
void Brusselator_scheme::eoc_mcf(){
  std::unique_ptr<SecOrd_op::vIdata> ex_ptr 
  {SecOrd_op::vIdata::new_sms(fix_grid)};
  SecOrd_op::Solution_driven X_solver {data, fix_grid, fef.u.app};
  ex_ptr->interpolate(fef.surface.fun);
  fef.surface.exact = fef.surface.fun;
  for(long it = 0; it < pattern_timeSteps(); ++it){
    X_solver.brusselator_rhs(fef.surface.fun, fef.surface.rhs_les);
    // les_rhs = M_3 u^n
    X_solver.solve(fef.surface.rhs_les, fef.surface.fun);
    // (M_3 + tau A_3) u^n+1 = les_rhs
    next_timeStep();
    fix_grid.new_nodes(fef.surface.exact);
    exact.X_ptr->interpolate(fef.surface.exact);
    io.surface << fix_grid.time_provider().deltaT() << ' '
	       << norm.l2_err(fef.surface.fun, fef.surface.exact) << ' '
	       << norm.h1_err(fef.surface.fun, fef.surface.exact) 
	       << std::endl;
    fix_grid.new_nodes(fef.surface.fun);
  }
}

void Brusselator_scheme::eoc_sls(){
  using namespace SecOrd_op;
  using std::unique_ptr;
  unique_ptr<vRhs> vRhs_ptr {vRhs::new_sls(fix_grid)};
  unique_ptr<vIdata> ex_ptr {vIdata::new_sls(fix_grid)};
  unique_ptr<sIdata> u_ptr {sIdata::new_1ssef(fix_grid)};
  Solution_driven X_solver {data, fix_grid, fef.u.app};
  ex_ptr->interpolate(fef.surface.fun);
  u_ptr->interpolate(fef.u.app);
  fef.surface.exact = fef.surface.fun;
  for(long it = 0; it < pattern_timeSteps(); ++it){
    X_solver.brusselator_rhs(fef.surface.fun, fef.surface.rhs_les);
    // les_rhs = (M_3 + alpha A_3) X^n
    vRhs_ptr->addScaled_to(fef.surface.rhs_les);
    // les_rhs += tau load_vector^n
    X_solver.solve(fef.surface.rhs_les, fef.surface.fun);
    // (M_3 + (alpha + tau * epsilon) A_3) X^n+1 = les_rhs
    next_timeStep();
    fix_grid.new_nodes(fef.surface.exact);
    ex_ptr->interpolate(fef.surface.exact);
    u_ptr->interpolate(fef.u.app);
    io.surface << fix_grid.time_provider().deltaT() << ' '
	       << norm.l2_err(fef.surface.fun, fef.surface.exact) << ' '
	       << norm.h1_err(fef.surface.fun, fef.surface.exact) 
	       << std::endl;
    fix_grid.new_nodes(fef.surface.fun);
  }
}

void Brusselator_scheme::sd(){
  using namespace SecOrd_op;
  std::unique_ptr<vRhs> vRhs_ptr {vRhs::new_sd(fix_grid)};
  std::unique_ptr<vIdata> ex_ptr {vIdata::new_sd(fix_grid)};
  Solution_driven X_solver {data, fix_grid, fef.u.app};
  ex_ptr->interpolate(fef.surface.fun);
  fef.surface.exact = fef.surface.fun;
  for(long it = 0; it < pattern_timeSteps(); ++it){
    X_solver.brusselator_rhs(fef.surface.fun, fef.surface.rhs_les);
    vRhs_ptr->addScaled_to(fef.surface.rhs_les);
    X_solver.solve(fef.surface.rhs_les, fef.surface.fun);
    next_timeStep();
    fix_grid.new_nodes(fef.surface.exact);
    ex_ptr->interpolate(fef.surface.exact); 
    io.surface << fix_grid.time_provider().deltaT() << ' '
	       << norm.l2_err(fef.surface.fun, fef.surface.exact) << ' '
	       << norm.h1_err(fef.surface.fun, fef.surface.exact) 
	       << std::endl;
    fix_grid.new_nodes(fef.surface.fun);
  }
}

void Brusselator_scheme::eoc_sdp(){
  using namespace SecOrd_op;
  using std::unique_ptr;
  unique_ptr<vRhs> g_load {vRhs::new_sls(fix_grid)};
  unique_ptr<sRhs> f_load {sRhs::new_sdp_u(fix_grid)};
  unique_ptr<vIdata> X_ex {vIdata::new_sls(fix_grid)};
  unique_ptr<vIdata> V_ex {vIdata::new_v_sls(fix_grid)};
  unique_ptr<sIdata> u_ex {sIdata::new_1ssef(fix_grid)};
  Scalar_solver solver {data, fix_grid, fef.u, fef.w}; // I just need `solver.u`
  Solution_driven X_solver {data, fix_grid, fef.u.app};
  X_ex->interpolate(fef.surface.fun);
  u_ex->interpolate(fef.u.app);
  fef.velocity.rhs_les = fef.surface.exact = fef.surface.fun;
  fef.u.fun = fef.u.app;
  for(long it = 0; it < pattern_timeSteps(); ++it){
    X_solver.brusselator_rhs(fef.surface.fun, fef.surface.rhs_les);
    g_load->addScaled_to(fef.surface.rhs_les);
    X_solver.solve(fef.surface.rhs_les, fef.surface.fun);
    solver.u.mass_matrix(fef.u.fun, fef.u.rhs_les);
    next_timeStep(); // next surface
    fix_grid.new_nodes(fef.surface.fun);    
    f_load->addScaled_to(fef.u.rhs_les);
    solver.u.solve(fef.u.rhs_les, fef.u.fun);
    calculate_velocity();
    fef.u.app = fef.u.fun;

    X_ex->interpolate(fef.surface.exact);
    fix_grid.new_nodes(fef.surface.exact); 
    u_ex->interpolate(fef.u.exact);
    V_ex->interpolate(fef.velocity.exact);
    print(io.u, fef.u);
    print(io.surface, fef.surface);
    print(io.velocity, fef.velocity);
    fix_grid.new_nodes(fef.surface.fun);
  }
}

void Brusselator_scheme::calculate_velocity(){
  using dunefef = Esfem::Grid::Vec_FEfun::Dune_FEfun;
  dunefef& v = fef.velocity.fun, x_old = fef.velocity.rhs_les;
  const dunefef& x_new = fef.surface.fun;
  const double dT = fix_grid.time_provider().deltaT();
  v.clear();
  v.axpy(1./dT, x_new);
  v.axpy(1./dT, x_old);
  x_old.assign(x_new);
}

// --------------------------------------------------
// Brusselator_scheme loop action

void Brusselator_scheme::prePattern_loop(){
  PrePattern_helper helper {*this};
  for(long it = 0; it < prePattern_timeSteps(); ++it, next_timeStep()){
    helper.rhs();
    // helper.solve_pde(); // Does this make any sense?
    helper.plot_paraview();
  }
}
void Brusselator_scheme::intermediate_action(){
  const auto uw_path = fef.tmpFile_path + "intermediate_";
  switch(prePattern_timeSteps()){
  case 0: // heat.starttime == heat.pattern.endtime
    fef.u.read(io.dgf_handler, uw_path);
    fef.w.read(io.dgf_handler, uw_path);
    break;
  default:
    fef.u.write(io.dgf_handler, uw_path);
    fef.w.write(io.dgf_handler, uw_path);
    // fef.surface.write(io.dgf_handler, fef.tmpFile_path);
    break;
  };
}
void Brusselator_scheme::pattern_loop(){
  for(long it = 0; it < pattern_timeSteps(); ++it){
    update_surface();
    update_scalar_solution();
    rhs_and_solve_SPDE();
    update_velocity();    
    error_on_intSurface(); // Error on surface(t_n)
    next_timeStep();
    Pattern_helper helper {*this};
    helper.finalize_scalarPDE_rhs();
    helper.solve_scalarPDE();
    // helper.errors_on_numSurface();
    // helper.plot_paraview();
  }
}
void Brusselator_scheme::final_action(){
  fef.u.write(io.dgf_handler, "./final_");
  fef.w.write(io.dgf_handler, "./final_");
  fef.surface.write(io.dgf_handler, "./final_"); 
}

// ------------------------------------------------------------
// Brusselator_scheme private

void Brusselator_scheme::update_surface(){
  io.identity.interpolate(fef.surface.exact);

  // Save old surface for the velocity
  fef.surface.app = fef.surface.fun; 
}
void Brusselator_scheme::update_velocity(){
  exact.v.interpolate(fef.velocity.exact);

  // update fef.velocity.fun
  const auto dT = fix_grid.time_provider().deltaT();
  auto xNew_ptr = fef.surface.fun.cbegin();
  auto xOld_ptr = fef.surface.app.cbegin();
  for(auto v_ptr = fef.velocity.fun.begin();
      v_ptr != fef.velocity.fun.end();
      ++v_ptr, ++xNew_ptr, ++xOld_ptr)
    *v_ptr = (*xNew_ptr - *xOld_ptr)/ dT;
}
void Brusselator_scheme::update_scalar_solution(){
  exact.u.interpolate(fef.u.exact);
  exact.w.interpolate(fef.w.exact);
}

void Brusselator_scheme::error_on_intSurface(){
  const auto dT = fix_grid.time_provider().deltaT();
  // scalar_error
  io.u << dT << '\t'
       << norm.l2_err(fef.u.exact, fef.u.fun) << '\t'
       << norm.h1_err(fef.u.exact, fef.u.fun) << std::endl;
  io.w << dT << '\t'
       << norm.l2_err(fef.w.exact, fef.w.fun) << '\t'
       << norm.h1_err(fef.w.exact, fef.w.fun) << std::endl;
  // surface_error
  io.surface << dT << '\t'
	     << norm.l2_err(fef.surface.exact, fef.surface.app)
	     << '\t'
	     << norm.h1_err(fef.surface.exact, fef.surface.app)
	     << std::endl;
  // velocity error
  io.velocity << dT << '\t'
	      << norm.l2_err(fef.velocity.fun, fef.velocity.exact)
	      << '\t'
	      << norm.h1_err(fef.velocity.fun, fef.velocity.exact)
	      << std::endl;
}

void Brusselator_scheme::pre_loop_action(){
  PreLoop_helper helper {*this};
  // helper.random_initialValues();
  helper.analytic_initialValues();
  exact.X_ptr->interpolate(fef.surface.fun);
  helper.headLine_in_errFile();
  helper.save_surface();
  // helper.plot_errors_in_errFile();
  // helper.plot_paraview();
  // next_timeStep();
}
void Brusselator_scheme::rhs_and_solve_SPDE(){
  RhsAndSolve_helper helper {*this};
  helper.scalar_massMatrix();
  helper.brusselator_rhs();
  helper.addScaled_surfaceLoadVector();
  helper.solve_surface_and_save();
}
// ----------------------------------------------------------------------
// Implementation of structs
// Brusselator_scheme::Fef and Brusselator_scheme::Io

Brusselator_scheme::Fef::Fef(const Esfem::Grid::Grid_and_time& gt)
  :u {"u", gt}, w {"w", gt},
   surface {"surface", gt}, velocity {"velocity", gt}
{}

Brusselator_scheme::Io::Io(const Esfem::Io::Parameter& p)
  :dgf_handler {p.grid()}, u {"_u", p}, w {"_w", p},
   surface {"_X", p}, velocity {"_v",p}
{}

Brusselator_scheme::Init_data::Init_data(const Esfem::Grid::Grid_and_time& gt)
  :u {gt, Growth::promoting}, w {gt, Growth::inhibiting}, v {gt},
   X_ptr {// SecOrd_op::vIdata::new_ssef(gt)
     SecOrd_op::vIdata::new_sms(gt)}
{}
Brusselator_scheme::Init_data::Init_data(const Esfem::Grid::Grid_and_time& gt, const Esfem::Io::Parameter& p)
  :u {p, Growth::promoting}, w {p, Growth::inhibiting}, v {gt},
   X_ptr {SecOrd_op::vIdata::new_ssef(gt)}
{}
