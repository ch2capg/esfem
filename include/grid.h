/*! \file grid.h
    \author Christian Power
    \date 22. March 2016

    \brief Fundamental header providing grid, time provider and finite element functions

    Revision history
    --------------------------------------------------

          Revised by Christian Power March 2016
	  Revised by Christian Power February 2016
          Originally written by Christian Power
               (power22c@gmail.com) January 2016

     Idea
     --------------------------------------------------

     Wrapper classes for the dune implementations.  This is very usefull
     since the construction of those classes is quite non trivial.

         Created by Christian Power on 17.03.2016
         Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef GRID_H
#define GRID_H 

#include <memory>
#include <config.h>
#include <dune/grid/geometrygrid.hh>	// geometry grid
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>	// grid part
#include <dune/grid/geometrygrid/coordfunction.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>	// discrete function space
#include <dune/fem/function/adaptivefunction.hh>	// discrete function
#include <dune/fem/function/blockvectorfunction.hh>	// 2. discrete function

#include "esfem_fwd.h"

namespace Esfem{
  namespace Grid{
    inline constexpr int grid_dim(){
      return Dune::GridSelector::dimgrid;
    }
    inline constexpr int world_dim(){
      return Dune::GridSelector::dimworld;
    }

    static_assert( grid_dim() == 2, "Bad grid dimension.");
    static_assert( world_dim() == 3, "Bad world dimension.");
    
    class Deformation
      : public Dune::
               AnalyticalCoordFunction<double, world_dim(), world_dim(), Deformation>
    {
    public:
      using Base = Dune::AnalyticalCoordFunction
	<double, world_dim(), world_dim(), Deformation>;
      using Domain = Base::DomainVector;
      using Range = Base::RangeVector;

      Deformation();
      explicit Deformation(const Impl::Evolving_grid&);
      ~Deformation();

      void evaluate(const Domain&, Range&) const;
      void set_timeProvider(const Dune::Fem::TimeProviderBase&);
    private:
      struct Data;
      std::unique_ptr<Data> d_ptr {};
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
      using Vec_Function_space = Dune::Fem::
	FunctionSpace<double, double, Grid::dimensionworld, Grid::dimensionworld>;
      using Vec_FE_space = Dune::Fem::
	LagrangeDiscreteFunctionSpace<Vec_Function_space, Grid_part, POLORDER>;
      
      explicit Grid_and_time(const Io::Parameter&);
      explicit Grid_and_time(const Io::Parameter&, const std::string& dgf_file,
			     const double t0);
      ~Grid_and_time();

      void next_timeStep(const double);
      void new_nodes(const Vec_FEfun&) = delete;
      
      Dune::Fem::TimeProviderBase& time_provider();
      const Dune::Fem::TimeProviderBase& time_provider() const;
      
      Grid& grid() const;
      Grid_part& grid_part() const;
      FE_space& fe_space() const;
      Vec_FE_space& vec_fe_space() const;
    private:
      struct Data;
      std::unique_ptr<Data> d_ptr {};
    };
    class Scal_FEfun{
    public:
      // using Dune_FEfun = Dune::Fem::AdaptiveDiscreteFunction<Grid_and_time::FE_space>;
      using Dune_FEfun 
      = Dune::Fem::ISTLBlockVectorDiscreteFunction<Grid_and_time::FE_space>;
      
      explicit Scal_FEfun(const std::string& name, const Grid_and_time& gt);
      explicit Scal_FEfun(const Scal_FEfun&, const Grid_and_time&);
      Scal_FEfun& operator=(const Scal_FEfun&);
      
      operator Dune_FEfun&(){ return fun; }
      operator const Dune_FEfun&() const { return fun; }

      auto size() const { return fun.size(); }
      auto begin() { return fun.dbegin(); }
      auto end() { return fun.dend(); }
      auto begin() const { return fun.dbegin(); }
      auto end() const { return fun.dend(); }
      auto cbegin() const { return fun.dbegin(); }
      auto cend() const { return fun.dend(); }
      auto name() const { return fun.name(); }

      Scal_FEfun& operator+=(const double d);
      Scal_FEfun& operator*=(const double d);
    private:
      Dune_FEfun fun;
    };
    class Vec_FEfun{
    public:
      using Dune_FEfun
      = Dune::Fem::AdaptiveDiscreteFunction<Grid_and_time::Vec_FE_space>;
      // using Dune_FEfun 
      // = Dune::Fem::ISTLBlockVectorDiscreteFunction<Grid_and_time::Vec_E_space>;
      
      explicit Vec_FEfun(const std::string& name, const Grid_and_time&);
      explicit Vec_FEfun(const Vec_FEfun&, const Grid_and_time&);
      Vec_FEfun& operator=(const Vec_FEfun&);

      operator Dune_FEfun&(){ return fun; }
      operator const Dune_FEfun&() const { return fun; }

      auto size() const { return fun.size(); }
      auto begin() { return fun.dbegin(); }
      auto end() { return fun.dend(); }
      auto begin() const { return fun.dbegin(); }
      auto end() const { return fun.dend(); }
      auto cbegin() const { return fun.dbegin(); }
      auto cend() const { return fun.dend(); }
      auto name() const { return fun.name(); }

      Vec_FEfun& operator+=(const double);
      Vec_FEfun& operator*=(const double);
    private:
      Dune_FEfun fun;
    };
  }	// namespace Grid
}	// namespace Esfem

#endif // GRID_H
