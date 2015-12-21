/*! \file iodof.h

    \brief Input and output for finite element function

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 17. Dezember 2015

     Giving interface for input and/or output of nodal values for a dune 
     finite elemente function. 

     Created by Christian Power on 17.12.2015
     Copyright (c) 2015 Christian Power.  All rights reserved.
 */

#ifndef IODOF_H
#define IODOF_H 

#include <string>
#include <fstream>

enum class IO_direction {write, read};

class IO_dune_fem{
public:
  IO_dune_fem(const std::string& filename_series_prefix,
	      const IO_direction write_or_read,
	      const size_t starting_no = 0);
  template<typename DiscreteFunction>
  void operator()(DiscreteFunction& fe);	
private:
  std::string filename;
  IO_direction io_dir;
  size_t f_no;
};

//======================================================================
// All implementation details

//----------------------------------------------------------------------
// Implementation details for IO_dune_fem
IO_dune_fem::IO_dune_fem(const std::string& filename_series_prefix,
			 const IO_direction write_or_read,
			 const size_t starting_no)
  : filename {filename_series_prefix}, io_dir {write_or_read}, f_no {starting_no}
    {
      switch(io_dir){
      case IO_direction::write: case IO_direction::read:
	break;
      default:
	throw std::runtime_error 
	{"Error in the constructor of IO_dune_fem.  "
	    "Encountered unknown IO_direction type."};
      }
    }

template<typename DiscreteFunction>
void IO_dune_fem::operator()(DiscreteFunction& fe){
  auto err_handling = [](const std::string& msg){
    std::string full_msg {"Error in IO_dune_fem::operator().\n"};
    full_msg += msg;
    throw std::runtime_error {full_msg};
  };
  
  const std::string curr_file {filename + std::to_string(f_no)};
  
  switch(io_dir){
  case IO_direction::write: {
    std::ofstream ofs {curr_file};
    if(!ofs) err_handling("Could not open file: " + curr_file);
    for(auto it = fe.dbegin(); it != fe.dend(); ++it)
      ofs << *it << std::endl;
  }
    break;
  case IO_direction::read: {
    std::ifstream ifs {curr_file};
    if(!ifs) err_handling("Could not open file: " + curr_file);
    for(auto it = fe.dbegin(); it != fe.dend(); ++it){
      double value;
      ifs >> value;
      *it = value;
    }      
  }
    break;
  default:
    err_handling("Impossible error in the switch statement.");
    break;
  }
  ++f_no;
}

#endif // IODOF_H

/*! Log:
 */
