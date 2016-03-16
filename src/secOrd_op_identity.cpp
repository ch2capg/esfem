/*! \file secOrd_op_identity.cpp
    \author Christian Power
    \date 9. March 2016

    \brief Implementing secOrd_op_identity.h. 

     Revision history:
     --------------------------------------------------

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Idea
     --------------------------------------------------

     The actual dune stuff is in secOrd_op_identity_impl.{h,cpp}.

         Created by Christian Power on 16.03.2016
         Copyright (c) 2016 Christian Power.  All rights reserved.
     
*/

#include "secOrd_op_identity.h"
#include "secOrd_op_identity_impl.h"

// ----------------------------------------------------------------------
// Implementation of Identity::Data

using namespace std;
using Esfem::SecOrd_op::Identity;
using Esfem::Impl::Identity_impl;

struct Identity::Data{
  Identity_impl {};
};

// ----------------------------------------------------------------------
// Implementation of Identity

Identity::Identity()
  : d_ptr {make_unique<Data>()}
{}

Identity::~Identity() = default;

void Identity::operator()(Grid::Vec_FEfun& vfef) const{
  using FE_function = Grid::Vec_FEfun::Dune_FEfun;
  Dune::LagrangeInterpolation<FE_function>::
    interpolateFunction(d_ptr -> Identity_impl, vfef);
}
