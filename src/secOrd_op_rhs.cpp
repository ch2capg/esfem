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
     \date 23. April 2016
     \copyright Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <config.h>
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

Rhs::Rhs(const Grid::Grid_and_time& gt, const Growth type)
  :d_ptr {std::make_unique<Data>(gt, type)}
{}
Rhs::~Rhs() = default;
void Rhs::assemble_and_addScaled_to(Grid::Scal_FEfun& fef){
  assemble_RHS(d_ptr -> rhs, d_ptr -> load_vector);
  FE_function& dune_fef = fef;
  dune_fef.axpy(d_ptr -> tp.deltaT(), d_ptr -> load_vector); 
}

// ----------------------------------------------------------------------
// Implementation of Vec_rhs

Vec_rhs::Vec_rhs(const Grid::Grid_and_time& gt)
  :d_ptr {std::make_unique<Data>(gt)}
{}
Vec_rhs::~Vec_rhs() = default;
void Vec_rhs::assemble_and_addScaled_to(Grid::Vec_FEfun& vfef){
  assemble_RHS(d_ptr -> rhs, d_ptr -> load_vector);
  Vec_FE_function& dune_vfef = vfef;
  dune_vfef.axpy(d_ptr -> tp.deltaT(), d_ptr -> load_vector);
}

// ----------------------------------------------------------------------
// Implementation of Rhs::Data and Vec_rhs::Data

Rhs::Data::Data(const Grid::Grid_and_time& gt, const Growth type)
  :tp {gt.time_provider()}, rhs {tp, type},
   load_vector {"load_vector", gt.fe_space()}
{}

Vec_rhs::Data::Data(const Grid::Grid_and_time& gt)
  :tp {gt.time_provider()}, rhs {tp},
   load_vector {"vec_load_vector", gt.vec_fe_space()}
{}
