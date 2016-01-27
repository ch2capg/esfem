/*
  For a given level set function this program tests the vertices of a dgf file
  to see if they are precise enough.

  Created by Christian Power on 11. January 2016.
  Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <valarray>

// cpp 2
// #include <algorithm>

// cpp 3
// #include <numeric>
// #include <cmath>

// header 1
void check_usage(int argc, char** argv);

// header 2
using Vertex = std::valarray<double>;
using Vertices = std::vector<Vertex>;
const Vertices get_vertices(const std::string& filename);

// cpp 2 -> header 2.1
// const Vertex convert_to_vertex(const std::string& str_vertex);

// header 3
double level_set_function(const Vertex& vertex);

int main(int argc, char** argv) try{
  check_usage(argc, argv);

  Vertices v_list = get_vertices(argv[1]);

  std::valarray<double> lsf_values (v_list.size());
  for(size_t it = 0; it < v_list.size(); ++it)
    lsf_values[it] = level_set_function(v_list[it]);

  // std::cout << "All values" << std::endl;
  // for(const auto value : lsf_values)
  //   std::cout << value << std::endl;
  
  std::cout << std::scientific;
  std::cout << "Mean value: " << lsf_values.sum() / v_list.size() << std::endl;

  std::cout << "Median: ";
  std::sort(std::begin(lsf_values), std::end(lsf_values));
  if(v_list.size() % 2)
    std::cout << lsf_values[ v_list.size()/2 ];
  else
    std::cout << (lsf_values[ v_list.size()/2 - 1 ] + lsf_values[ v_list.size()/2])/2;
  std::cout << std::endl;
}
 catch(std::exception& e){
   std::cerr << e.what() << '\n';
   return 1;
 }


// cpp 1
void check_usage(int argc, char** argv){
  if(argc != 2){
    std::ostringstream err_msg;
    err_msg << "Usage: " << argv[0] << " filename " << '\n'
	    << "Program terminating due to invalid input.";
    throw std::runtime_error{err_msg.str()};
  }
}

// cpp 2
#include <algorithm>

const Vertex convert_to_vertex(const std::string& str_vertex);

const Vertices get_vertices(const std::string& filename){
  std::ifstream dgf_file {filename};
  if(!dgf_file) throw std::runtime_error {"File \"" + filename + "\" does not exist."};

  std::vector<std::string> lines;
  for(std::string line; std::getline(dgf_file,line); )
    lines.push_back(line);

  const std::string begin_vertices = "VERTEX";
  auto first_it = std::find(std::cbegin(lines), std::cend(lines), begin_vertices);
  if(first_it == std::cend(lines))
    throw std::runtime_error
    {"Could not find " + begin_vertices + " in file \"" + filename + "\"."};
  ++first_it;	// skip "VERTEX"
  
  const std::string end_vertices = "#";
  const auto second_it = std::find(first_it, std::cend(lines), end_vertices);
  if(second_it == first_it)
    throw std::runtime_error
    {"Empty vertex list detected in file \"" + filename + "\"."};

  Vertices vertex_list;
  int line_no = 1;
  for( ; first_it != second_it; ++first_it){
    try{
      const auto numerical_vertex = convert_to_vertex(*first_it);
      vertex_list.push_back(numerical_vertex);
    }
    catch(const std::exception& e){
      std::cerr << e.what() << '\n';
      const std::string err_msg {"line" + std::to_string(line_no)};
      throw std::runtime_error {err_msg};
    }
    ++line_no;
  }

  return vertex_list;
}

const Vertex convert_to_vertex(const std::string& str_vertex){
  std::istringstream iss {str_vertex};
  double d1, d2, d3;
  if(!(iss >> d1 >> d2 >> d3) || !(iss >> std::ws).eof() )
    throw std::runtime_error {"convert_to_vertex() failed."};
  return {d1, d2, d3};
}

// cpp 3
#include <numeric>
#include <cmath>

double level_set_function(const Vertex& vertex){
  // x**2 + y**2 + z**2 - 1;
  using std::begin;
  using std::end;
  return std::inner_product(begin(vertex), end(vertex), begin(vertex), 0.);
}
