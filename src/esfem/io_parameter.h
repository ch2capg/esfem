/*! \file io_parameter.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 27.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef IO_PARAMETER_H
#define IO_PARAMETER_H 

#include <string>
#include <vector>
#include "esfem_fwd.h"

namespace Esfem{
  namespace Io{
    class Parameter{
    public:
      explicit Parameter(int argc, char** argv, const std::string& parameter_file_name);
      ~Parameter();
      Parameter(const Parameter&) = delete;
      Parameter& operator=(const Parameter&) = delete;

      // ----------------------------------------------------------------------
      const std::string& grid() const noexcept;
      const std::string& error_log() const noexcept;
      const std::string& paraview() const noexcept;
      
      double start_time() const noexcept;
      double global_timeStep() const;
      long max_timeSteps() const;
    
      double eps() const noexcept;
    
      const std::vector<double>& bdf_alphas() const noexcept;
      const std::vector<double>& bdf_gammas() const noexcept;

      double tg_a() const noexcept;
      double tg_b() const noexcept;
      double tg_Dc() const noexcept;
      double tg_gamma() const noexcept;
      
      // ----------------------------------------------------------------------
      friend std::ostream& operator<<(std::ostream&, const Parameter&);
    private:
      struct Data;
      Data* d_ptr {nullptr};
    };
  }
}

#endif // IO_PARAMETER_H

/*! Log:
 */
