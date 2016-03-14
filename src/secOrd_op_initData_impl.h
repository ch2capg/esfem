/*! \file secOrd_op_initData_impl.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Februar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 21.02.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef SECORD_OP_INITDATA_IMPL_H
#define SECORD_OP_INITDATA_IMPL_H 

#include <random>
#include <functional>
#include <string>
#include <config.h>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include "grid.h"

namespace Esfem{
  class Explicit_initial_data
    : public Dune::Fem::Function<Esfem::Grid::Grid_and_time::Function_space,
				 Explicit_initial_data>
  {
  public:
    using FE_space = Esfem::Grid::Grid_and_time::Function_space;
    using Domain = FE_space::DomainType;
    using Range = FE_space::RangeType;

    explicit Explicit_initial_data(const Esfem::Grid::Grid_and_time&);
    Explicit_initial_data(const Explicit_initial_data&) = delete;
    Explicit_initial_data& operator=(const Explicit_initial_data&) = delete;

    void evaluate(const Domain&, Range&) const;
  private:
    const Dune::Fem::TimeProviderBase& tp;
  };
  class Random_initial_data
    : public Dune::Fem::Function<Esfem::Grid::Grid_and_time::Function_space,
				 Random_initial_data>
  {
  public:
    using FE_space = Esfem::Grid::Grid_and_time::Function_space;
    using Domain = typename FE_space::DomainType;
    using Range = typename FE_space::RangeType;

    explicit Random_initial_data(const Esfem::Io::Parameter&,
				  const Esfem::Growth);
    void evaluate(const Domain&, Range&) const;
  
  private:
    using Random_dist = std::uniform_real_distribution<>;
    using Random_engine = std::default_random_engine;
    std::function<double()> random_fun;

    explicit Random_initial_data(const double hom_value, const double pertubation);
  };

  // ----------------------------------------------------------------------
  // helper functions 

  // ------------------------------------------------------------
  // Random_init_data helper
  double hom_value(const Esfem::Io::Parameter&, const Esfem::Growth);
  double pertubation(const Esfem::Io::Parameter&, const Esfem::Growth);
  std::string print_configuration(const Esfem::Io::Parameter&, const Esfem::Growth);

  // ------------------------------------------------------------
  // Init_data::Data helper
  std::string dof_filename(const Io::Parameter&, const Growth);
}


#endif // SECORD_OP_INITDATA_IMPL_H

/*! Log:
 */
