/*! \file grid_evolution.h

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

#ifndef GRID_EVOLUTION_H
#define GRID_EVOLUTION_H 

#include "parameter.h"
#include <dune/grid/geometrygrid/coordfunction.hh>
#include <dune/fem/solver/timeprovider.hh>
// #include <dune/fem/function/common/function.hh>

namespace Grid_evolution{
  class DeformationCoordFunction
    : public Dune::AnalyticalCoordFunction<double, 3, 3, DeformationCoordFunction>
  {
  public:
    using Base = Dune::AnalyticalCoordFunction<double, 3, 3, DeformationCoordFunction>;
    using Domain = Base::DomainVector;
    using Range = Base::RangeVector;
  
    DeformationCoordFunction() = default;
    explicit DeformationCoordFunction(const Parameter::PDE_data&);
    DeformationCoordFunction(const DeformationCoordFunction&) = delete;
    DeformationCoordFunction(DeformationCoordFunction&&) = delete;
    DeformationCoordFunction& operator=(const DeformationCoordFunction&) = delete;
    DeformationCoordFunction& operator=(DeformationCoordFunction&&) = delete;
    ~DeformationCoordFunction() = default;
  
    void set_timeProvider(Dune::Fem::TimeProviderBase&);
    void evaluate(const Domain&, Range& y) const;
  private:
    class Grid_values;
    Grid_values* gv_ptr {nullptr}; 
    Dune::Fem::TimeProviderBase* tp_ptr {nullptr}; 
  };

  // pseudo external variables
  DeformationCoordFunction& deformation(const Parameter::PDE_data&);
}

#endif // GRID_EVOLUTION_H

/*! Log:
 */
