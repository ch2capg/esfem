/*! \file 2015linfty.h

    \brief Numerical experiments for the 2015linfty paper

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 17. Dezember 2015

     In the paper we considered space discretization of the linear 
     heat equation on evolving surfaces

     ∂•u + div(v) u - Δu = f

     with linear finite elements.  The purpose is to measure the error in 
     maximum norm.  Although the paper does not cover full discretization,  
     we discretized in time with BDF 4.  

     The header w1l.norm.hh implements
     
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 17.12.2015
     Copyright (c) 2015 Christian Power.  All rights reserved.
 */

#ifndef 2015LINFTY_H
#define 2015LINFTY_H 


#include "dune_headers.h"

#include "deformation.hh"
#include "elliptic.hh"
#include "heat.hh"
#include "rhs.hh"
// local CAPG includes
#include "dune_typedef_heat.hpp"
// this header only contains lots of typedef's
#include "dune_bdf.hpp"
#include "dune_heat_algorithm.hpp"


#endif // 2015LINFTY_H

/*! Log:
 */
