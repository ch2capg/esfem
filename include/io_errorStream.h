/*! \file io_errorStream.h
    \brief Providing a costume stream 

    Revision history:
    --------------------------------------------------

         Revised by Christian Power April 2016
         Revised by Christian Power February 2016
         Originally written by Christian Power
              (power22c@gmail.com) Januar 2016

    Idea
    --------------------------------------------------
    Provides a wrapper class for fstream, that is compatible with the 
    Esfem::Io::Parameter class.  

    \author Christian Power
    \date 25. April 2016
    \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef IO_ERRORSTREAM_H
#define IO_ERRORSTREAM_H 

#include <fstream>
#include "esfem_fwd.h"

namespace Esfem{
  namespace Io{
    //! Wrapper class for a std::stream
    class Error_stream{
    public:
      //! `filename` will be `Io::Parameter::error_log` 
      explicit Error_stream(const Parameter&);
      //! `filename` will be `Io::Parameter::error_log + suffix` 
      explicit Error_stream(const std::string& suffix, const Parameter&);

      // ----------------------------------------------------------------------
      // Providing 'endl', 'scientific' and 'defaultfloat'
      //! Auxiliary typedef
      using Basic_ostream = std::basic_ostream<char, std::char_traits<char> >;
      //! The type of std::endl, std::float, etc.
      using StdManip = Basic_ostream& (*)(Basic_ostream&);
      //! Now you can use std::endl, etc.
      Error_stream& operator<<(StdManip);

      //! Stream type `T` out
      template<typename T>
      Error_stream& operator<<(const T&);
      //! Close the stream
      Error_stream& close();
    private:
      //! The actual stream implementation
      std::ofstream ofs;

      //! Constructor helper function
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
