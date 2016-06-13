/*! \file grid.h
    \brief Fundamental header providing grid, time provider and finite element functions

    Revision history
    --------------------------------------------------

          Revised by Christian Power June 2016
          Revised by Christian Power Mai 2016
          Revised by Christian Power March 2016
	  Revised by Christian Power February 2016
          Originally written by Christian Power
               (power22c@gmail.com) January 2016

    Idea
    --------------------------------------------------

    Wrapper classes for the dune implementations.  This is very usefull
    since the construction of those classes is quite non trivial.

    \author Christian Power
    \date 7. June 2016
    \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
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

      //! Analytic formula
      Deformation();
      //! Hash map for grid evolution
      explicit Deformation(const std::string& fname);
      //! Pointer to implementation requires this.
      ~Deformation();

      //! Explicit mapping
      void evaluate(const Domain&, Range&) const;
      //! Finish cyclic dependency
      void set_timeProvider(const Dune::Fem::TimeProviderBase&);
      //! Update hash map
      Deformation& operator=(const Vec_FEfun&);
    private:
      struct Data;
      //! Pointer to implementation
      std::unique_ptr<Data> d_ptr;
    };
    //! FE grid, FE space, function space
    class Grid_and_time{
    public:
      //! Important for dune 
      using Host_grid = Dune::GridSelector::GridType;
      //! \copybrief Host_grid
      using Grid = Dune::GeometryGrid<Host_grid, Deformation>;
      //! ALU-Grid
      using Grid_part
      = Dune::Fem::AdaptiveLeafGridPart<Grid, Dune::InteriorBorder_Partition>;
      //! \f$f\colon \R^3\to \R\f$
      using Function_space = Dune::Fem::
	FunctionSpace<double, double, Grid::dimensionworld, 1>;
      //! Scalar valued finite element space
      using FE_space = Dune::Fem::
	LagrangeDiscreteFunctionSpace<Function_space, Grid_part, POLORDER>;
      //! \f$ f\colon \R^3\to \R^3\f$
      using Vec_Function_space = Dune::Fem::
	FunctionSpace<double, double, Grid::dimensionworld, Grid::dimensionworld>;
      //! Vector valued finite element space
      using Vec_FE_space = Dune::Fem::
	LagrangeDiscreteFunctionSpace<Vec_Function_space, Grid_part, POLORDER>;

      //! Construct first grid
      explicit Grid_and_time(const Io::Parameter&);
      //! Read from an dgf
      explicit Grid_and_time(const Io::Parameter&, const std::string& dgf_file,
			     const double t0);
      //! Pointer to implementation requires this.
      ~Grid_and_time();

      //! ++time
      void next_timeStep(const double);
      //! Update hash grid
      void new_nodes(const Vec_FEfun&);

      //! Get time provider
      Dune::Fem::TimeProviderBase& time_provider();
      //! Get time provider
      const Dune::Fem::TimeProviderBase& time_provider() const;

      //! Grid for finite element functions
      Grid& grid() const;
      //! I believe for the norm
      Grid_part& grid_part() const;
      //! Finite element space
      /*! Used for dune operator */
      FE_space& fe_space() const;
      //! Vector valued finite element space
      /*! Used for dune operator */
      Vec_FE_space& vec_fe_space() const;
    private:
      struct Data;
      //! Pointer to implementation
      std::unique_ptr<Data> d_ptr;
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
      // = Dune::Fem::ISTLBlockVectorDiscreteFunction<Grid_and_time::Vec_FE_space>;
      
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
