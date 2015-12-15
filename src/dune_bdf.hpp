#ifndef DUNE_BDF_HPP
#define DUNE_BDF_HPP

/*********************************************************************
 *  dune_bdf.hpp, dune_bdf.hxx                                       *
 *                                                                   *
 *  Defines auxiliar classes like a string vector with the correct   *
 *  file names.                                                      *
 *                                                                   *
 *  Revision history:                                                *
 *  none (this is experimental and error phrone)                     *
 *                                                                   *
 *                                                                   *
 *  Created by Christian Power on 09.03.15.                          *
 *  Copyright (c) 2015 Christian Power. All rights reserved.         *
 *                                                                   *
 *********************************************************************/

#include <iostream>		// debugging
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <exception>
#include <unordered_map>
//#include <map>
#include <array>
#include <algorithm>
#include <Eigen/Dense>

namespace BDF
{
  class EvoMapType{
  public:
	  EvoMapType(const std::string& input_filename);

	  template<typename In, typename Integrator>
	   // requires Input_iterator<In>()
	  void evolve(In p_t_0, In p_t_end, Integrator integrator);
	  
	  const std::array<double,3>& operator[] (const std::array<double,3>& vec) const;
  private:
	  std::unordered_map<std::string, std::array<double, 3> > evoMap;
	   // evoMap: Γʰ(0) → Γʰ(n), evoMap[node⁰] = nodeⁿ;
	  std::string filename;
	  size_t precision;
  };

	void ie_heat_algorithm();		// BDF1 algorithm
	void bdf2_heat_algorithm();		// debugging
	void nonlinear_algorithm();		// the actuall BDF algorithm;
	// void algorithm();
	// void poisson_problem();	// for testing reasons

	// All implementations details

	// Implementation details for 'EvoMapType's helper functions

	std::vector<std::string> getVertexList(const std::string& filename)
	 // opens the file 'filename' and return a vector of lines, where each line
	 // represents a vertex of the grid.  The lines are whitespace cleaned.
	{
		// local helper functions
		auto error = [](const std::string& msg){ 
			throw 
			std::runtime_error("ERROR: BDF::getVertexList()\n"
							   + msg ); 
		};
		auto clean_whitespace = [](const std::string& line)
		// A 'line' may begin with whitespace. Clean it!
			{
				std::string clean_line;
				std::istringstream iss {line};
				for(std::string s; iss >> s; )
					clean_line += s + ' ';
				clean_line.pop_back();	// delete last whitespace
				return clean_line;
			};
		// end local helper functions

		std::ifstream ifs {filename};
		if(!ifs) error("File \"" + filename + "\" does not exist!");

		std::vector<std::string> lines;
		for(std::string line; getline(ifs, line); )
			lines.push_back(clean_whitespace(line));
		ifs.close();

		// get list of vertices in terms of iterators
		std::vector<std::string>::iterator first =
			std::find(lines.begin(), lines.end(), "VERTEX");
		if( first == lines.end() ) error("Missing keyword \"VERTEX\" in file \""
										 + filename + "\"!");
		++first;	// "VERTEX" is not part of the list
		std::vector<std::string>::iterator last =
			std::find(first, lines.end(), "#");
		if( last == lines.end() ) error("Missing ending symbol \"#\" in file \""
										+ filename + "\"!");
	
		return std::vector<std::string> {first, last};
	
	}

	size_t get_precision(const std::vector<std::string>& vertex_list){
		if(vertex_list.size() == 0) return 0;

		std::vector<size_t> precision_vec;
		for(const std::string& line : vertex_list){
			std::istringstream interpreter {line};
			for(std::string s; interpreter >> s; ){
				precision_vec.push_back(s.size());
			}
		}

		std::vector<size_t>::iterator  p_precision = 
			std::max_element(precision_vec.begin(), precision_vec.end());
		return *p_precision;
	}

	template<typename T, size_t N >
	std::string array_to_string(const std::array<T, N> ad,
								const unsigned int precision = 8){
		std::ostringstream oss;
		oss << std::scientific;
		 // setprecision and scientific is essential! For setprecision it is mandatory
		 // that the number  be greater then precision in the dgf file.
		for(std::size_t s = 0; s < ad.size() -1 ; ++s)
			oss << std::setprecision(precision) << ad[s] << ' ';
		oss << ad.back();	// the last one without ' '
		return oss.str();
	}

	template<typename T, size_t N >
	std::array<T, N> string_to_array(const std::string s){
		auto error = [] (const std::string& msg) { 
			throw std::runtime_error("ERROR: BDF::string_to_array()\n" + msg); 
		};
		std::istringstream interpreter {s};
		std::array<T, N> return_array;
		for(auto& it : return_array)
			if(!(interpreter >> it)) error("operator>>() failed!");
		if(!(interpreter >> std::ws).eof())	// stuff left in stream?
			error("whitespace left!");
		return return_array;
	}

	// Implementation details for 'EvoMapType'

	BDF::EvoMapType::EvoMapType( const std::string& input_filename ) :
		filename {input_filename}
	{
		std::vector<std::string> lines = getVertexList(filename);
		precision = get_precision(lines);
		//std::cout << "precision = " << precision << std::endl; // debugging
		for(auto& line : lines)
			evoMap[ array_to_string( string_to_array<double, 3>(line), precision) ] = 
				string_to_array<double, 3>(line);
	}

	const std::array<double, 3>&
	BDF::EvoMapType::operator[] (const std::array<double, 3>& node) const { 
		return evoMap.at(array_to_string(node, precision));
	}

	template<typename In, typename Integrator>
	 // requires Input_iterator<In>() && 
	 // Vector3d Integrator.integrate(In first, In last, Vector3d x_0)
	 // Vector3d from the 'Eigen' library
	void BDF::EvoMapType::evolve(In p_t_0, In p_t_end, Integrator integrator){
		for(auto& it : evoMap){
			Eigen::Map<Eigen::Vector3d> x_0 {it.second.data()};	// get vector
			x_0 = integrator.integrate(p_t_0, p_t_end, x_0);	// perform action
			std::array<double, 3> it_new;						// translate it back
			std::copy(x_0.data(), x_0.data() + x_0.size(), it_new.begin());
			it.second = it_new;

			// nice testing idea by Balázs
			// std::cout << "surface error: " 
			// 		  << (x_0(0)*x_0(0) / (1. + (sin(2. * M_PI * time_vec.back() ) / 4.)
			//				) ) + x_0(1) * x_0(1) + x_0(2) * x_0(2) - 1
			// 		  << std::endl;
		}
	}

}




// #include "dune_bdf.hxx"
// CAPG remark: (CAPG convention) usually we include in *.hxx files only
// template related implementation, but because the dune make files are
// complicated, we put all the *.cpp stuff there

#endif // #ifndef DUNE_BDF_HPP
