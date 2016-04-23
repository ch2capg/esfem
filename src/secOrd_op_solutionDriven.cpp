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

#include <config.h>
#include <dune/fem/solver/cginverseoperator.hh>
#include "secOrd_op_solutionDriven.h"
#include "secOrd_op_solutionDriven_impl.h"
#include "io_parameter.h"
#include "grid.h"


using Esfem::SecOrd_op::Solution_driven;
using Esfem::Impl::MCF_op;
using FEfun = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using Vec_FEfun = Esfem::Grid::Vec_FEfun::Dune_FEfun;
using Vec_cg_solver = Dune::Fem::CGInverseOperator<Vec_FEfun>;

// ----------------------------------------------------------------------
// Implementation of Solution_driven::Data

struct Solution_driven::Data{
  MCF_op mcf_op;
  Vec_cg_solver cg_solver;
  Data(const Io::Parameter&, const Grid::Grid_and_time&,
       const Grid::Scal_FEfun& u_wrapper);
};

Solution_driven::Data::Data(const Io::Parameter& p, const Grid::Grid_and_time& g,
			    const Grid::Scal_FEfun& u_wrapper)
  :mcf_op {p, g, u_wrapper},
   cg_solver {mcf_op, p.eps(), p.eps()}
{}

// ----------------------------------------------------------------------
// Implementation of Solution_driven

Solution_driven::Solution_driven(const Io::Parameter& p,
				 const Grid::Grid_and_time& g,
				 const Grid::Scal_FEfun& u)
  :d_ptr {std::make_unique<Data>(p, g, u)}
{}

Solution_driven::~Solution_driven() = default;

void Solution_driven::
solve(const Grid::Vec_FEfun& rhs, Grid::Vec_FEfun& lhs) const{
  const Vec_FEfun& vfef1 = rhs;
  Vec_FEfun& vfef2 = lhs;
  d_ptr -> cg_solver(vfef1, vfef2);
}

void Solution_driven::rhs(const Grid::Vec_FEfun& rhs, Grid::Vec_FEfun& lhs) const{
  const Vec_FEfun& vfef1 = rhs;
  Vec_FEfun& vfef2 = lhs;
  d_ptr -> mcf_op.rhs(vfef1, vfef2);
}
