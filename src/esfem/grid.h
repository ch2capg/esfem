/*! \file grid.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 27. Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 27.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef GRID_H
#define GRID_H 

#include <config.h>
#include <dune/grid/geometrygrid.hh>	// geometry grid
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>	// grid part
#include <dune/grid/geometrygrid/coordfunction.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>	// discrete function space
#include <dune/fem/function/adaptivefunction.hh>	// discrete function
#include "esfem_fwd.h"

namespace Esfem{
  namespace Grid{
    class Deformation
      : public Dune::AnalyticalCoordFunction<double, 3, 3, Deformation>
    {
    public:
      using Base = Dune::AnalyticalCoordFunction<double, 3, 3, Deformation>;
      using Domain = Base::DomainVector;
      using Range = Base::RangeVector;

      Deformation();
      explicit Deformation(const Io::Parameter&) = delete;
      ~Deformation();
      Deformation(const Deformation&) = delete;
      Deformation& operator=(const Deformation&) = delete;

      void evaluate(const Domain&, Range&) const;
      void set_timeProvider(const Dune::Fem::TimeProviderBase&);
    private:
      struct Data;
      Data* d_ptr {nullptr};
    };
    class Grid_and_time{
    public:
      using Host_grid = Dune::GridSelector::GridType;
      using Grid = Dune::GeometryGrid<Host_grid, Deformation>;
      using Grid_part
      = Dune::Fem::AdaptiveLeafGridPart<Grid, Dune::InteriorBorder_Partition>;
      using Function_space = Dune::Fem::
	FunctionSpace<double, double, Grid::dimensionworld, 1>;
      using FE_space = Dune::Fem::
	LagrangeDiscreteFunctionSpace<Function_space, Grid_part, POLORDER>;
      
      explicit Grid_and_time(const Io::Parameter&);
      ~Grid_and_time();

      void next_timeStep(const double);
      Dune::Fem::TimeProviderBase& time_provider();
      const Dune::Fem::TimeProviderBase& time_provider() const;
      Grid& grid() const;
      Grid_part& grid_part() const;
      FE_space& fe_space() const;
    private:
      struct Data;
      Data* d_ptr {nullptr};
    };
    class Scal_FEfun{
    public:
      using Dune_FEfun = Dune::Fem::AdaptiveDiscreteFunction<Grid_and_time::FE_space>;

      explicit Scal_FEfun(const std::string& fun_name, const Grid_and_time& gt)
	: fun {fun_name, gt.fe_space()}
      { fun *= 0; }
      Scal_FEfun& operator=(const Scal_FEfun& other){
	fun.assign(other.fun);
      }
      Scal_FEfun(const Scal_FEfun&) = delete;
      
      operator Dune_FEfun&(){ return fun; }
      operator const Dune_FEfun&() const { return fun; }
    private:
      Dune_FEfun fun;
    };
  }
}

#endif // GRID_H

/*! Log:
 */
