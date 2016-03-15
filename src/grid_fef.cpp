/*! \file grid.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 27. Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 27.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#include "grid.h"

using Scal_FEfun = Esfem::Grid::Scal_FEfun;
using Vec_FEfun = Esfem::Grid::Vec_FEfun;

// ----------------------------------------------------------------------
// Implementation Scal_FEfun

Scal_FEfun::Scal_FEfun(const std::string& fun_name, const Grid_and_time& gt)
  : fun {fun_name, gt.fe_space()}
 { 
   fun.clear(); 
 }
Scal_FEfun& Scal_FEfun::operator=(const Scal_FEfun& other){
  fun.assign(other.fun);
  return *this;
}
Scal_FEfun& Scal_FEfun::operator+=(const double d){
  for(auto it = fun.dbegin(); it != fun.dend(); ++it) *it += d;
  return *this;
}
Scal_FEfun& Scal_FEfun::operator*=(const double d){
  for(auto it = fun.dbegin(); it != fun.dend(); ++it) *it *= d;
  return *this;
}
Scal_FEfun::Scal_FEfun(const Scal_FEfun& other, const Grid_and_time& gt)
  : fun {other.name() + "+", gt.fe_space()}
{ 
  fun.assign(other.fun); 
}

// ----------------------------------------------------------------------
// Implemenation Vec_FEfun

Vec_FEfun::Vec_FEfun(const std::string& fun_name, const Grid_and_time& gt)
  : fun {fun_name, gt.vec_fe_space()}
{ 
  fun.clear(); 
}
Vec_FEfun::Vec_FEfun(const Vec_FEfun& other, const Grid_and_time& gt)
  : fun {other.name() + "+", gt.vec_fe_space()}
{ 
  fun.assign(other.fun); 
}
Vec_FEfun& Vec_FEfun::operator=(const Vec_FEfun& other){
  fun.assign(other.fun);
  return *this;
}

Vec_FEfun& Vec_FEfun::operator+=(const double d){
  for(auto it = fun.dbegin(); it != fun.dend(); ++it) *it += d;
  return *this;
}
Vec_FEfun& Vec_FEfun::operator*=(const double d){
  for(auto it = fun.dbegin(); it != fun.dend(); ++it) *it *= d;
  return *this;
}
