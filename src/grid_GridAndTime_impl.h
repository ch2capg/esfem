/*! \file grid_GridAndTime_impl.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 25. Februar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 25.02.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef GRID_GRIDANDTIME_IMPL_H
#define GRID_GRIDANDTIME_IMPL_H 

#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include "grid.h"

namespace Esfem{
  namespace Impl{
    class Evolving_grid{
    public:
      using Node = Grid::Deformation::Domain;

      Evolving_grid(const std::string& filename);
    
      const Node& operator[](const Node&) const;
      Evolving_grid& operator=(const Grid::Vec_FEfun& new_nodes);

      using Nodes_key = std::vector<std::string>;
      using Map = std::unordered_map<std::string, Node>;
    private:    
      const Nodes_key original_vertices {};
      std::size_t digit_precision {8};
      Map vertices_map {};

      template<typename It, typename Ot>
      void update_vertices(It first, It last, Ot out, const int dim);
      // Assumption: distance(first, last) == dim * vertices_map.size();
    };

    // ----------------------------------------------------------------------
    // Constructor helper functions

    Evolving_grid::Nodes_key get_vertexList(const std::string& filename);
    // Extract from a dgf all vertices.  The resulting vector is then compatible
    // with a Dune discrete finite element function.
    std::size_t get_relevant_precision(const Evolving_grid::Nodes_key& vertex_list);
    // Calculate minimum precision needed for all Nodes_key
    Evolving_grid::Map
    identity_map(const Evolving_grid::Nodes_key&, const std::size_t precision);
    
    std::size_t get_relevant_precision(const std::string& number);
    // Assuming scientific format, e.g., 0.9298e-10.
    // For this example the function would return 4.
    std::vector<std::size_t>
    create_precision_vec(const std::vector<std::string>& vertex_list);
    // Applies the function above to a whole vector, where each line contains
    // three numbers.
    
    std::size_t max(const std::vector<std::size_t>&);
    
    std::string node_to_string(const Evolving_grid::Node&, const int precision);
    Evolving_grid::Node string_to_node(const std::string&);
    std::string normalize_node(const std::string&, const int precision);
    
    std::string clean_whitespace(const std::string&);
    
    std::vector<std::string> dgfFile_to_vec(const std::string& filename);
    
    template<typename It>
    It dgf_find_vertex(It first, It last, const std::string& filename);
    template<typename It>
    It dgf_close_list(It first, It last, const std::string& filename);
    
    // ----------------------------------------------------------------------
    // template implementation
    
    template<typename It, typename Ot>
    void Evolving_grid::update_vertices(It first, It last, Ot out, const int dim){
      while(first != last){
	Node new_vertex {};
	for(int i = 0; i < dim; ++i){
	  new_vertex[i] = *first;
	  ++first;
	}
	try{
	  vertices_map.at(*out) = new_vertex;
	}
	catch(...){
	  const auto err_msg = "Bad node: " + *out;
	  throw std::logic_error {err_msg};
	}
	++out;
      }
    }

    template<typename It>
    It dgf_find_vertex(It first, It last, const std::string& filename){
      std::string starting_point {"VERTEX"};
      It rv = std::find(first, last, "VERTEX");
      if(rv == last)
	throw std::runtime_error
	{"Missing keyword \"" + starting_point
	    + "\" in file " + filename + "."};
      ++rv; // "VERTEX" is not part of the list
      return rv;
    }
    template<typename It>
    It dgf_close_list(It first, It last, const std::string& filename){
      std::string ending_point {"#"};
      It rv = std::find(first, last, ending_point);
      if(rv == last)
	throw std::runtime_error
	{"Missing ending symbol \"" + ending_point
	    + "\" in file " + filename + "."};
      return rv;
    }
  }
}

#endif // GRID_GRIDANDTIME_IMPL_H

/*! Log:
 */
