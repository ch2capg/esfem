/*! \file io_errorStream.cpp
    \brief Implementation of io_errorStream.h 

    Revision history
    --------------------------------------------------

         Revised by Christian Power April 2016
         Revised by Christian Power February 2016
         Originally written by Christian Power
              (power22c@gmail.com) Januar 2016

    Idea
    --------------------------------------------------

    No special idea

    \author Christian Power
    \date 25. April 2016
    \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#include <sstream>
#include "io_errorStream.h"
#include "io_parameter.h"

#ifdef DEBUG
#include <iostream>
#endif 

void cannot_open_file(const std::string& fname);

Esfem::Io::Error_stream::Error_stream(const std::string& fname)
  :ofs {fname, std::ios_base::app}
{
#ifdef DEBUG
  std::clog << "Opening " << fname << " for Error_stream" << std::endl;
#endif
  if(!ofs) cannot_open_file(fname);
}
Esfem::Io::Error_stream::Error_stream(const Parameter& d)
  :Error_stream {d.error_log()}
{}
Esfem::Io::Error_stream::Error_stream(const std::string& suffix, const Parameter& d)
  :Error_stream {d.error_log() + suffix}
{}
Esfem::Io::Error_stream& Esfem::Io::Error_stream::operator<<(StdManip manip){
#if DEBUG > 8
  std::clog << "Using operator<<(StdManip)." << std::endl;
#endif
  manip(ofs);
  return *this;
}
Esfem::Io::Error_stream& Esfem::Io::Error_stream::close(){
  ofs.close();
  return *this;
}

// ----------------------------------------------------------------------
// Internal implementation

void cannot_open_file(const std::string& fname){
  std::ostringstream err_msg;
  err_msg << "Error in constructor of Error_stream. "
    "Could not open file: " << fname;
  throw std::runtime_error {err_msg.str()};
}

/*! Log:
    19.02.16: Fixed bug in with cannot_open_file().  runtime_error was not
     throwing.  
 */
