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
#include <dassert.h>

namespace Esfem{
  namespace Impl{
    struct Parameter_error : std::runtime_error{
      explicit Parameter_error(const std::string& msg) 
	: std::runtime_error {"Error in Esfem::Io::Parameter.\n" + msg}
      {}
    };
    template<int level>
    void parameter_assert(const bool assertion, const char* fname, const int line, 
    			  const std::string& msg){
      Assert::dynamic<level, Parameter_error>(assertion, fname, line, msg);
    }

    const std::string& project_dir();
    void dune_fem_parameter_append(int argc, char** argv, const std::string& file);
    std::string get_gridKey();
    std::string get_macroGrid();
    std::string doubleVector_to_string(const std::vector<double>&);
    void file_check(const std::vector<std::string>& file_list);
  }
}

template<>
inline void Assert::dynamic<false, Esfem::Impl::Parameter_error>(const bool, const std::string&){};

#endif // IO_PARAMETER_IMPL_H
