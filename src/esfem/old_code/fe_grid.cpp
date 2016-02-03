/*! \file fe_grid.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 21. Januar 2016

     Implementation details for esfem.h
     Created by Christian Power on 21.01.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */
#include "fe_grid.h"
#include <sstream>
#include <stdexcept>
#include "check_pseudo_externals.h"

#ifdef DEBUG
#include <iostream>
#endif

constexpr int no_externalVars = 4;
static Check_pseudo_externals::Control_table<no_externalVars> control_vars {};

void FE_grid::init_grid_and_timeProvider(const Parameter::PDE_data& d) try{
  auto& deformation = Grid_evolution::deformation(d);
  const auto& grid_name = d.grid(); 
  auto& hGrid_ptr = FE_grid::host_grid(&grid_name);
  hGrid_ptr -> loadBalance();
  auto& grid = FE_grid::grid(&hGrid_ptr, &deformation);
  FE_grid::grid_part(&grid);
  const auto t_0 = d.start_time();
  auto& tp = FE_grid::time_provider(&t_0, &grid);
  deformation.set_timeProvider(tp);
  tp.init(d.global_timeStep());
}
catch(const std::exception& e){
  std::string err_msg {"Error in init_grid().\n"};
  err_msg += e.what();
  throw std::runtime_error {err_msg};
}
FE_grid::HGrid_ptr& FE_grid::host_grid(const std::string* const gridfile_ptr){
#ifdef DEBUG
  std::clog << "host_grid(\"" << *gridfile_ptr << "\")" << std::endl;
#endif
  constexpr int var_no = 0;
  if(!gridfile_ptr && control_vars.not_init<var_no>())
    throw std::logic_error
    {"Error in host_grid().  This pseudo exteral variable has not been initialized, "
	"but function is invoked with default arguments."};
  static Dune::GridPtr<FE_grid::Host_grid> host_grid {*gridfile_ptr};
  control_vars.just_init<var_no>();
  return host_grid;
}
FE_grid::Grid& FE_grid::grid(FE_grid::HGrid_ptr* const hg_ptr_ptr,
			     Grid_evolution::DeformationCoordFunction* const d_ptr){
#ifdef DEBUG
  std::clog << "grid(" << hg_ptr_ptr << ", " << d_ptr << ')' << std::endl;
#endif
  constexpr int var_no = 1;
  if((!hg_ptr_ptr || !d_ptr) && control_vars.not_init<var_no>())
    throw std::logic_error
    {"Error in grid().  This pseudo exteral variable has not been initialized, "
	"but function is invoked with default arguments."};
  static FE_grid::Grid grid(**hg_ptr_ptr, *d_ptr);
  control_vars.just_init<var_no>();
  return grid;
}
FE_grid::Grid_part& FE_grid::grid_part(FE_grid::Grid* const g_ptr){
#ifdef DEBUG
  std::clog << "grid_part(" << g_ptr << ")" << std::endl;
#endif
  constexpr int var_no = 2;
  if(!g_ptr && control_vars.not_init<var_no>())
    throw std::logic_error
    {"Error in grid_part().  This pseudo exteral variable has not been initialized, "
	"but function is invoked with default arguments."};
  static FE_grid::Grid_part grid_part(*g_ptr);
  control_vars.just_init<var_no>();
  return grid_part;
}
FE_grid::Time_provider& FE_grid::time_provider(const double* const t0_ptr ,
					       FE_grid::Grid* const g_ptr){
#ifdef DEBUG
  std::clog << "time_provider(" << t0_ptr << ", " << g_ptr << ")" << std::endl;
#endif
  constexpr int var_no = 3;
  if((!t0_ptr || !g_ptr) && control_vars.not_init<var_no>())
    throw std::logic_error
    {"Error in grid_part().  This pseudo exteral variable has not been initialized, "
	"but function is invoked with default arguments."};
  static FE_grid::Time_provider tp {*t0_ptr, *g_ptr};
  control_vars.just_init<var_no>();
  return tp;
}
/*! Log:
 */
