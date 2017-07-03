/*! \file grid_GridAndTime_impl.h
    \brief Hash maps for evolving grids

     Revision history
     --------------------------------------------------

          Revised by Christian Power June 2016
          Originally written by Christian Power
               (power22c@gmail.com) 25. Februar 2016

    \author Christian Power
    \date 7. June 2016
    \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef GRID_GRIDANDTIME_IMPL_H
#define GRID_GRIDANDTIME_IMPL_H 

#include <vector>
#include <array>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <boost/functional/hash.hpp>
#include "grid.h"

//! Template specializations
namespace std{
  //! Hash map for Domain points
  /*! \remark Unfortunately does dune change the last precision digits. */
  template<> struct std::hash<Esfem::Grid::Deformation::Domain>{
    //! Produce hash value via boost helper
    /*! \pre Dimension is 3. */
    std::size_t operator()(const Esfem::Grid::Deformation::Domain& n) const{
      using std::hash;
      // return hash<double>{}(n[0])
      //        ^ (hash<double>{}(n[1])<<1)>>1
      //        ^ (hash<double>{}(n[2])<<1)>>1;
      using boost::hash_value;
      using boost::hash_combine;
      size_t seed {0};
      hash_combine(seed, hash_value(n[0]));
      hash_combine(seed, hash_value(n[1]));
      hash_combine(seed, hash_value(n[2]));
      return seed;
    }
  };
  //! Hash map for an array of integers 
  template<> struct std::hash<array<int, Esfem::Grid::world_dim()> >{
    //! Convenience
    using Base = array<int, Esfem::Grid::world_dim()>;
    //! Produce hash value via boost helper
    /*! \pre Dimension is 3. */
    std::size_t operator()(const Base& ai) const{
      // using std::hash;
      // return hash<double>{}(n[0])
      //        ^ (hash<double>{}(n[1])<<1)>>1
      //        ^ (hash<double>{}(n[2])<<1)>>1;
      using boost::hash_value;
      using boost::hash_combine;
      size_t seed {0};
      for(int i = 0; i < Esfem::Grid::world_dim(); ++i)
	hash_combine(seed, hash_value(ai[i]));
      return seed;
    }
  };
  //! operator==() for Domain points
  template<> struct std::equal_to<Esfem::Grid::Deformation::Domain>{
    //! Component-wise comparison
    /*! \pre Dimension is 3. */
    bool operator()(const Esfem::Grid::Deformation::Domain& lhs,
		    const Esfem::Grid::Deformation::Domain& rhs) const{
      bool rv = true;
      for(int i = 0; i < Esfem::Grid::Deformation::Domain::dimension; ++i)
	rv = rv && lhs[i] == rhs[i];
      return rv;
    }
  };
  //! operator==() for Domain points
  template<> struct std::equal_to<array<int, Esfem::Grid::world_dim()> >{
    //! Short alias
    using Base = array<int, Esfem::Grid::world_dim()>;
    //! Component-wise comparison
    /*! \pre Dimension is 3. */
    bool operator()(const Base& lhs, const Base& rhs) const{
      bool rv = true;
      for(int i = 0; i < Esfem::Grid::world_dim(); ++i)
	rv = rv && lhs[i] == rhs[i];
      return rv;
    }
  };  
}

namespace Esfem{
  namespace Impl{    
    //! Evolution via a hash map
    /*! I use std::string as key. */
    class Evolving_grid{
    public:
      //! What is a node for us
      using Node = Grid::Deformation::Domain;
      //! Empty map
      Evolving_grid() = default;
      //! Construct grid from dgf
      explicit Evolving_grid(const std::string& filename);

      //! Access node
      const Node& operator[](const Node&) const;
      //! Get new nodes
      Evolving_grid& operator=(const Grid::Vec_FEfun& new_nodes);

      //! List of keys for hash map
      using Nodes_key = std::vector<std::string>;
      //! The actual map
      using Map = std::unordered_map<std::string, Node>;
    private:
      //! The list
      const Nodes_key original_vertices;
      //! Precision of the strings
      std::size_t digit_precision {8};
      //! Hash map
      Map vertices_map;

      //! Helper for operator=()
      /*! \pre distance(first, last) == dim * vertices_map.size().*/
      template<typename It, typename Ot>
      void update_vertices(It first, It last, Ot out, const int dim);
    };

    //! hash grid plus helper functions
    namespace hash{
      //! Basic dune entity
      using range = Grid::Deformation::Domain;
      //! Int because of rounding errors
      using key = std::array<int, range::dimension>;
      
      //! Evolution via a hash map
      /*! I use boost functions to combine hash values of doubles. */
      class grid{
	//! Short alias for dimension
	static constexpr auto dim = range::dimension;
	//! Hash map type
	using map = std::unordered_map<key, range>;
	//! Sequence
	using seq = std::vector<key>;
	//! Map 
	map m;
	//! Original list
	seq ol;
      public:
	//! Restricted hash map iterator
	class iterator{
	public:
	  //! Actual iterator
	  using iterator_impl = map::iterator;
	  //! Pseudo copy iterator
	  explicit iterator(const iterator_impl& it_in) :it {it_in} {}
	  //! Get range value
	  range& operator*() const{
	    return (*it).second;
	  }
	  //! Increment
	  iterator& operator++(){
	    ++it;
	    return *this;
	  }
	  //! Is equal?
	  bool operator==(const iterator& other){
	    return this->it == other.it;
	  }
	  //! Is different?
	  bool operator!=(const iterator& other){
	    return !(*this == other);
	  }
	private:
	  //! Instance of the actual iterator
	  iterator_impl it;
	};
	//! Range for-loop
	iterator begin(){
	  return iterator {m.begin()};
	}
	//! Range for-loop
	iterator end(){
	  return iterator {m.end()};
	}	
	//! Reports errors
	struct bad : std::runtime_error{
	  //! Forward error message
	  explicit bad(const std::string& msg) :std::runtime_error {msg} {}
	};
	//! Empty map
	grid() = default;
	//! Read first list from file name
	explicit grid(const std::string& fname);
	//! Get first list
	/*! \remark Use Lagrange interpolation of the identity. */
	explicit grid(const Grid::Vec_FEfun& init_keys);
	//! Get new list
	grid& operator=(const Grid::Vec_FEfun& value_list);
	//! Checked access
	/*! \throws bad as a nested object if the point does not exists*/
	const range& operator[](const range&) const;
      };
      //! Make a string out of a range
      std::string to_string(const range&);
      //! Multiply each entry with 1e5 and then truncate
      key to_key(const range&);
    }
    
    // ----------------------------------------------------------------------
    // Evolving grid constructor helper functions

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
