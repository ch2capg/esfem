/*! \file secOrd_op_solutionDriven.cpp
    \author Christian Power
    \date 17. March 2016

    \brief Implementing `Esfem::SecOrd_op::Solution_driven`

     Revision history
     --------------------------------------------------

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Idea
     --------------------------------------------------

     Implementation of `secOrd_op_solutionDriven.h`

         Created by Christian Power on 17.03.2016
         Copyright (c) 2016 Christian Power.  All rights reserved.

*/

#include "secOrd_op_solutionDriven.h"
#include "io_parameter.h"
#include "grid.h"

using Esfem::SecOrd_op::Solution_driven;
using FEfun = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using Vec_FEfun = Esfem::Grid::Vec_FEfun::Dune_FEfun;

// ----------------------------------------------------------------------
// Implementation of Solution_driven::Data

struct Solution_driven::Data{
  const double alpha;
  const double delta;
  const double epsilon;
  const Dune::Fem::TimeProviderBase& tp;
  const FEfun& u;
  Data(const Io::Parameter&, const Grid::Grid_and_time&
       const Grid::Scal_FEfun& u_wrapper);
};

Solution_driven::Data::Data(const Io::Parameter& p, const Grid::Grid_and_time& g,
			    const Grid::Scal_FEfun& u_input)
  : alpha {p.velocity_regularization()},
  delta {p.surface_growthFactor()},
  epsilon {p.mcf_regularization()},
  tp {g.time_provider()},
  u {static_cast<FEfun>(u_input)}
{}

// ----------------------------------------------------------------------
// Implementation of Solution_driven

Solution_driven::Solution_driven(const Io::Parameter& p, const Grid::Grid_and_time& g,
				 const Grid::Scal_FEfun& u)
  : d_ptr {std::make_unique<Data>(p, g, u)}
{}

Solution_driven::~Solution_driven() = default;

void Solution_driven::solve(const Grid::Vec_FEfun& rhs, Grid::Vec_FEfun& lhs) const{
  const FEfun& fef1 = rhs;
  FEfun& fef2 = lhs;
  d_ptr -> mcf_solver(fef1, fef2);
}

void Solution_driven::rhs(Grid::Vec_FEfun&) const{
}
