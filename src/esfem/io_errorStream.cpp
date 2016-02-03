/*! \file io_errorStream.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 27. Januar 2016

     Implementation details for io_errorStream.h
     Created by Christian Power on 27.01.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <sstream>
#include "io_errorStream.h"
#include "io_parameter.h"

#ifdef DEBUG
#include <iostream>
#endif 

Esfem::Io::Error_stream::Error_stream(const Parameter& d)
  : ofs {d.error_log(), std::ios_base::app}
    {
#ifdef DEBUG
      std::clog << "Opening " << d.error_log() << " for Error_stream" << std::endl;
#endif
      if(!ofs){
	std::ostringstream err_msg;
	err_msg << "Error in constructor of Error_stream. "
	  "Could not open file: " << d.error_log();
	std::runtime_error {err_msg.str()};
      }
    }
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

/*! Log:
 */
