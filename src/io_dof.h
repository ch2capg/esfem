/*! \file io_dof.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 21. Februar 2016

     Routine to read and write dof of finite element functions to files.

     Created by Christian Power on 21.02.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef IO_DOF_H
#define IO_DOF_H 

#include <fstream>
#include "esfem_fwd.h"

namespace Esfem{
  namespace Io{
    template<typename It>
    void file_to_dof(It first, It last,
		     const std::string& filename);
    // Bad style above.  Output algorithm only have on iterator not two.
    template<typename It>
    void file_to_dof(It first, const std::string& filename, std::size_t dof_no);
    template<typename It>
    void dof_to_file(It first, It last,
		     const std::string& filename);

    // ----------------------------------------------------------------------
    // Template implementation

    // ------------------------------------------------------------
    // Error helper functions
    
    inline void err_open_file(const std::string& fname){
      throw std::logic_error {"Could not open file " + fname};
    }
    inline void err_file_to_small(const std::string& fname){
      throw std::runtime_error {"File \"" + fname + "\" not big enough."};
    }
    
    template<typename It>
    void file_to_dof(It first, It last,
		     const std::string& filename) try{      
      std::ifstream ifs {filename};
      if(!ifs) err_open_file(filename);
      while(first != last){
	double value {0.};
	if(ifs >> value) *first = value;
	else err_file_to_small(filename);
	++first;
      }
    }
    catch(const std::exception&){
      std::throw_with_nested(std::runtime_error
			     {"Error in file_to_dof()."});
    }
    template<typename It>
    void file_to_dof(It first, const std::string& filename,
		     std::size_t dof_no) try{      
      std::ifstream ifs {filename};
      if(!ifs) err_open_file(filename);
      for( ; dof_no > 0; --dof_no){
	double value {0.};
	if(ifs >> value) *first = value;
	else err_file_to_small(filename);
	++first;
      }
      if(!(ifs >> std::ws).eof()) throw std::runtime_error
	{"Stuff left in file \"" + filename + "\"."};
    }
    catch(const std::exception&){
      std::throw_with_nested(std::runtime_error
			     {"Error in file_to_dof()."});
    }
    template<typename It>
    void dof_to_file(It first, It last,
		     const std::string& filename) try{
      std::ofstream ofs {filename};
      if(!ofs) err_open_file(filename);
      while(first != last){
	if( !(ofs << *first << std::endl) )
	  err_file_to_small(filename);
	++first;
      }
    }
    catch(const std::exception&){
      std::throw_with_nested(std::runtime_error
			     {"Error in dof_to_file()."});
    }
  }
}


#endif // IO_DOF_H

/*! Log:
 */
