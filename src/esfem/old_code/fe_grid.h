/*! \file fe_grid.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 25. Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 25.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef FE_GRID_H
#define FE_GRID_H 

#include <config.h>
#include <dune/grid/geometrygrid.hh>	// geometry grid
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>	// grid part
#include <dune/fem/solver/timeprovider.hh>
#include "parameter.h"
#include "grid_evolution.h"

namespace FE_grid{
  using Host_grid = Dune::GridSelector::GridType;
  using HGrid_ptr = Dune::GridPtr<Host_grid>;
  using Grid = Dune::GeometryGrid<Host_grid, Grid_evolution::DeformationCoordFunction>;
  using Grid_part
  = Dune::Fem::AdaptiveLeafGridPart<Grid, Dune::InteriorBorder_Partition>;
  using Time_provider = Dune::Fem::GridTimeProvider<Grid>;
  
  void init_grid_and_timeProvider(const Parameter::PDE_data&);

  // pseudo external variables
  HGrid_ptr& host_grid(const std::string* const gridfile_ptr = nullptr);
  Grid& grid(HGrid_ptr* const = nullptr,
	     Grid_evolution::DeformationCoordFunction* const = nullptr);
  Grid_part& grid_part(Grid* const = nullptr);
  Time_provider& time_provider(const double* const t0_ptr = nullptr,
			       Grid* const = nullptr);
}

#endif // FE_GRID_H

/*! Log:
 */
