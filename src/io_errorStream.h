/*! \file io_errorStream.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power February 2016
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     Provides a wrapper class for fstream, that is compatible with the 
     Esfem::Io::Parameter class.  

     Created by Christian Power on 19.02.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef IO_ERRORSTREAM_H
#define IO_ERRORSTREAM_H 

#include <fstream>
#include "esfem_fwd.h"

namespace Esfem{
  namespace Io{
    class Error_stream{
    public:
      explicit Error_stream(const Parameter&);
      // filename is parameter.error_log()
      explicit Error_stream(const std::string& suffix, const Parameter&);
      // filename is suffix + parameter.error_log()
      Error_stream(const Error_stream&) = delete;
      Error_stream& operator=(const Error_stream&) = delete;

      // ----------------------------------------------------------------------
      // Providing 'endl', 'scientific' and 'defaultfloat'
      using Basic_ostream = std::basic_ostream<char, std::char_traits<char> >;
      using StdManip = Basic_ostream& (*)(Basic_ostream&);
      Error_stream& operator<<(StdManip);
    
      template<typename T>
      Error_stream& operator<<(const T&);
      Error_stream& close();
    private:
      std::ofstream ofs;

      explicit Error_stream(const std::string& fname);
      // filename is fname
    };

    // ----------------------------------------------------------------------
    // Template implementation
    template<typename T>
    Error_stream& Error_stream::operator<<(const T& t){
      ofs << t;
      return *this;
    }
  } 
}

#endif // IO_ERRORSTREAM_H

/*! Log:
 */
