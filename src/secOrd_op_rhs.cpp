/*! \file secOrd_op_rhs.cpp
    \brief Implementing `secOrd_op_rhs.h`

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) January 2016

     Idea
     --------------------------------------------------

     Implementation of classes `Rhs` and `Vec_rhs`.

     \author Christian Power
     \date 16. April 2016
     \copyright Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <cmath>
#include <config.h>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include "secOrd_op_rhs.h"
#include "secOrd_op_rhs_impl.h"

//! Implementing this
using Esfem::SecOrd_op::Rhs;
//! Implementing this
using Esfem::SecOrd_op::Vec_rhs;
//! Dune scalar valued finite element function
using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
//! Dune vector valued finite element function
using Vec_FE_function = Esfem::Grid::Vec_FEfun::Dune_FEfun;

// ----------------------------------------------------------------------
// Implementation of Rhs

Rhs::Rhs(const Grid::Grid_and_time& gt)
  :d_ptr {std::make_unique<Data>(gt)}
{}
Rhs::~Rhs() = default;
void Rhs::assemble_and_addScaled_to(Grid::Scal_FEfun& fef){
  static FE_function load_vector {"local_tmp", d_ptr -> fe_space};
  assemble_RHS(d_ptr -> rhs, load_vector);
  FE_function& dune_fef = fef;
  dune_fef.axpy(d_ptr -> tp.deltaT(), load_vector); 
}

// ----------------------------------------------------------------------
// Implementation of Vec_rhs

Vec_rhs::Vec_rhs(const Grid::Grid_and_time& gt)
  :d_ptr {std::make_unique<Data>(gt)}
{}
Vec_rhs::~Vec_rhs = default;
void Vec_rhs::assemble_and_addScaled_to(Grid::Vec_FEfun& vfef){
  assemble_vecRhs(d_ptr -> vec_rhs, vec_load_vector);
  Vec_FE_function& dune_vfef = vfef;
  dune_vfef.axpy(d_ptr -> tp.deltaT(), vfef_rhs);
}

// ----------------------------------------------------------------------
// Implementation of Rhs::Data and Vec_rhs::Data

Rhs::Data::Data(const Grid::Grid_and_time& gt)
  :rhs {gt.time_provider()}, tp {gt.time_provider()},
   fe_space {gt.fe_space()}
{}

Vec_rhs::Data::Data(const Grid::Grid_and_time& gt)
  :rhs {gt.time_provider()}, tp {gt.time_provider()},
   fe_space {gt.fe_space()}
{}
