/*! \file io_l2h1Calculator.cpp
    \brief Implementation of io_l2h1Calculator.h

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     Idea
     --------------------------------------------------

     No special idea.

     \author Christian Power
     \date 27. April 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#include <config.h>
#include <stdexcept>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include "io_l2h1Calculator.h"
#include "grid.h"

using namespace std;
using Esfem::Io::L2H1_calculator;

//! \f$f\colon \R^3 \to \R\f$
using FEfun = Esfem::Grid::Scal_FEfun::Dune_FEfun;
//! \f$f\colon \R^3 \to \R^3\f$
using Vec_FEfun = Esfem::Grid::Vec_FEfun::Dune_FEfun;
//! \f$L^2\f$-norm
using L2_norm = Dune::Fem::L2Norm<Esfem::Grid::Grid_and_time::Grid_part>;
//! \f$H^1\f$-norm
using H1_norm = Dune::Fem::H1Norm<Esfem::Grid::Grid_and_time::Grid_part>;

// ----------------------------------------------------------------------
// Helper template function

template<typename FEF, typename Norm>
double norm_err_helper(const Norm& n, const FEF& u1, const FEF& u2){
  using Dfef = typename FEF::Dune_FEfun;
  const Dfef& du1 = u1;
  const Dfef& du2 = u2;
  return n.distance(du1, du2);
}

// ----------------------------------------------------------------------
// Actual implementations 

//! %Data members of L2H1_calculator
struct L2H1_calculator::Data{
  //! Dune \f$L^2\f$-norm functor
  L2_norm l2;
  //! Dune \f$H^1\f$-norm functor
  H1_norm h1;
  //! Get grid
  /*! \post Grid_and_time must outlive this object. */
  Data(const Grid::Grid_and_time& gt) :l2 {gt.grid_part()}, h1{gt.grid_part()} {}
};

L2H1_calculator::L2H1_calculator(const Grid::Grid_and_time& gt)
  :d_ptr {make_unique<Data>(gt)} {}
Esfem::Io::L2H1_calculator::~L2H1_calculator() = default;

double L2H1_calculator::l2_err(const Grid::Scal_FEfun& u, const Grid::Scal_FEfun& uN) const{
  return norm_err_helper(d_ptr->l2, u, uN);
}
double L2H1_calculator::l2_err(const Grid::Vec_FEfun& u, const Grid::Vec_FEfun& uN) const{
  return norm_err_helper(d_ptr->l2, u, uN);
}
double L2H1_calculator::h1_err(const Grid::Scal_FEfun& u, const Grid::Scal_FEfun& uN) const{
  return norm_err_helper(d_ptr->h1, u, uN);
}
double L2H1_calculator::h1_err(const Grid::Vec_FEfun& u, const Grid::Vec_FEfun& uN) const{
  return norm_err_helper(d_ptr->h1, u, uN);
}
