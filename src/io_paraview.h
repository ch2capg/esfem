/*! \file io_paraview.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 31. Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 31.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef IO_PARAVIEW_H
#define IO_PARAVIEW_H 

#include <memory>
#include "esfem_fwd.h"

namespace Esfem{
  namespace Io{
    class Paraview{
    public:
      explicit Paraview(const Parameter&, const Grid::Grid_and_time&,
			Grid::Scal_FEfun&, Grid::Scal_FEfun&,
			const int refinement_label = 0);
      ~Paraview();
      Paraview(const Paraview&) = delete;
      Paraview& operator=(const Paraview&) = delete;

      void write();
    private:
      struct Data;
      std::unique_ptr<Data> d_ptr;
    };
  } 
}

#endif // IO_PARAVIEW_H

/*! Log:
 */
