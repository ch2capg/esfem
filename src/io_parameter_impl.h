/*! \file io_parameter_impl.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Implementation details for io_parameter.h

     Created by Christian Power on 15.03.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
*/

#ifndef IO_PARAMETER_IMPL_H
#define IO_PARAMETER_IMPL_H

#include <string>
#include <vector>

namespace Esfem{
  namespace Impl{
    const std::string& project_dir();
    void dune_fem_parameter_append(int argc, char** argv, const std::string& file);
    std::string get_gridKey();
    std::string get_macroGrid();
    std::string doubleVector_to_string(const std::vector<double>&);
    void file_check(const std::vector<std::string>& file_list);
  }
}

#endif // IO_PARAMETER_IMPL_H
