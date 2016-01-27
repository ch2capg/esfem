/*! \file dune_headers.h

    \brief Include all standard dune headers

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 17. Dezember 2015

     This header includes config.h and all standard dune header used by my FEM code.  

     Created by Christian Power on 17.12.2015
     Copyright (c) 2015 Christian Power.  All rights reserved.
 */

#ifndef DUNE_HEADERS_H
#define DUNE_HEADERS_H 



#include <config.h>

// iostream includes
#include <iostream>

// include GeometryGrid
#include <dune/grid/geometrygrid.hh>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// include discrete function space
//#include <dune/fem/space/lagrangespace.hh> veraltet
#include <dune/fem/space/lagrange.hh>

// include discrete function
#include <dune/fem/function/adaptivefunction.hh>

// include solvers
//#include <dune/fem/solver/inverseoperators.hh> veraltet
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/oemsolver.hh>

// include Lagrange interpolation
#include <dune/fem/operator/lagrangeinterpolation.hh>

// include norms
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// include parameter handling
#include <dune/fem/io/parameter.hh>

// include output
#include <dune/fem/io/file/dataoutput.hh>

#endif // DUNE_HEADERS_H

/*! Log:
 */
