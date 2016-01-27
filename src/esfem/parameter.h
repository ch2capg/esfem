/*! \file parameter.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 25.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef PARAMETER_H
#define PARAMETER_H 

#include <config.h>
#include <iostream>
#include <vector>

namespace Parameter{
  class PDE_data{
  public:
    PDE_data() = delete;
    explicit PDE_data(int argc, char** argv, const char* const parameter_file_name);
    PDE_data(const PDE_data&) = delete;
    PDE_data(PDE_data&&) = delete;
    PDE_data& operator=(const PDE_data&) = delete;
    PDE_data& operator=(PDE_data&&) = delete;
    ~PDE_data();

    const std::string& grid() const noexcept;
    const std::string& error_log() const noexcept;
    double start_time() const noexcept;
    double global_timeStep() const;
    long max_timeSteps() const;
    double eps() const noexcept;
    const std::vector<double>& bdf_alphas() const noexcept;
    const std::vector<double>& bdf_gammas() const noexcept;
    
    friend std::ostream& operator<<(std::ostream&, const PDE_data&);
  private:
    struct Data;
    Data* d_ptr;
  };
  std::ofstream& err_stream(const PDE_data&);
}

#endif // PARAMETER_H

/*! Log:
 */
