/*! \file io_dgf.h
    \author Christian Power
    \date 22. March 2016

    \brief Provides class Io::Dgf::Handler to communicate
           between dgf files and finite element functions

    Revision history
    --------------------------------------------------

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Idea
     --------------------------------------------------

     Providing routines to input resp. output nodal values of a scalar valued
     finite element function to files according to the dune grid format.

         Created by Christian Power on 17.03.2016
         Copyright (c) 2016 Christian Power.  All rights reserved.
*/

#ifndef IO_DGF_H
#define IO_DGF_H 

#include <sstream>
#include <iterator>
#include <iomanip>
#include "esfem_fwd.h"

#include <iostream>

namespace Esfem{
  namespace Io{
    namespace Dgf{
      class Handler{
      public:
	Handler(const std::string& dgf_filename);
	~Handler();
      
	void write(const std::string& out_filename, const Grid::Scal_FEfun&) const;
	void write(const std::string& out_filename, const Grid::Vec_FEfun&) const;

	void read(const std::string& in_filename, Grid::Vec_FEfun&) const;
	void read(const std::string& in_filename, Grid::Scal_FEfun&) const;
	std::size_t no_of_nodes() const noexcept;
	
	using Vertex = std::vector<double>;
	using Vertices = std::vector<Vertex>;
	using Simplex = std::vector<short>;
	using Simplices = std::vector<Simplex>;
	using Lines = std::vector<std::string>;
      private:
	struct Data;
	std::unique_ptr<Data> d_ptr;
      };

      // ----------------------------------------------------------------------
      // Constructor helper functions

      Handler::Simplices create_simplices(const std::string& filename);
      Handler::Lines create_stringSimplices(const std::string& filename);

      std::size_t get_precision(const std::string& filename);
      std::size_t calculate_dof(const std::string& dgf_file);
      
      Handler::Vertices create_vertices(const std::string& filename);
      Handler::Lines create_stringVertices(const std::string& filename);
      
      // ----------------------------------------------------------------------
      // Helper functions for dgf files

      // ------------------------------------------------------------
      // Navigation 

      template<typename It>
      It find(It first, It last, const std::string& keyword);
      template<typename It>
      It close_list(It first, It last, const std::string& keyword);

      // ------------------------------------------------------------
      // Read and write

      template<typename FEfun>
      FEfun& operator<<(FEfun&, const Handler::Vertices&);
      
      template<typename CC>	// We assume that CC::value_type is a container.
      void to_dgfFile(const CC&, std::ostream&, const std::size_t precision);
      
      template<typename C>
      C to(const std::string&);
      template<typename C>
      std::string from(const C&, const std::size_t precision);

      template<typename C>
      Handler::Vertices to_vertices(const C&, const int dim = 3);
      template<typename It>
      Handler::Vertices
      to_vertices(It first, It last, const int dim, const std::size_t size);

      // ----------------------------------------------------------------------
      // template implementation

      template<typename It>
      It find(It first, It last, const std::string& keyword){
	It rv = std::find(first, last, keyword);
	if(rv == last) throw std::runtime_error
	  {"Could not find \"" + keyword + "\"."};
	++rv;
	return rv;
      }    
      template<typename It>
      It close_list(It first, It last, const std::string& keyword){
	It rv = std::find(first, last, keyword);
	if(rv == last) throw std::runtime_error
	  {"Could not find closing \"" + keyword + "\"."};
	return rv;
      }
      template<typename CC>
      void to_dgfFile(const CC& ccontainer, std::ostream& os,
		      const std::size_t precision){
	for(const auto& container : ccontainer)
	  os << from(container, precision) << std::endl;
      }
      template<typename C>
      C to(const std::string& s){
	using Number = typename C::value_type;
	std::istringstream iss {s};
	C rv {};
	for(Number n {0}; iss >> n; ) rv.push_back(n);
	return rv;
      }
      template<typename C>
      std::string from(const C& container, const std::size_t precision){
	std::ostringstream oss {};
	oss << std::scientific;
	for(const auto it : container)
	  if( !(oss << std::setprecision(precision) << it << ' ') )
	    throw std::runtime_error {"Esfem::Io::Dgf::from()."};
	return oss.str();
      }
      template<typename C>
      Handler::Vertices to_vertices(const C& container, const int dim){
	using std::begin;
	using std::end;
	const auto size = container.size();
	if(size % dim) throw std::logic_error
	  {"Bad degree of freedom dimension."};
	return to_vertices(begin(container), end(container), dim, size);
      }
      template<typename It>
      Handler::Vertices to_vertices(It first, It last,
				    const int dim, const std::size_t size) try{
	using Vertex = Handler::Vertices::value_type;
	Handler::Vertices rv {};
	rv.reserve(size);
	while(first != last){
	  Vertex v(dim);
	  for(auto& it : v){
	    it = *first;
	    ++first;
	  }
	  rv.emplace_back(std::move(v));
	}
	return rv;
      }
      catch(const std::exception&){
	std::throw_with_nested(std::runtime_error {"Esfem::Io::Dgf::to_vertices()."});
      }

      template<typename FEfun>
      FEfun& operator<<(FEfun& fef, const Handler::Vertices& v){
	using std::begin;
	using std::end;
	auto o_it = begin(fef);
	auto first = cbegin(v);
	auto last = cend(v);
	while(first != last){
	  const auto vertex = *first;
	  for(const auto d : vertex){
	    *o_it = d;
	    ++o_it;
	  }
	  ++first;
	}
	if(o_it != end(fef)) throw std::logic_error
	  {"Esfem::Io::Dgf::operator<<().  Dof to small"};
	return fef;
      }
    }	// namespace Dgf
  }	// namespace Io
}	// namespace Esfem

#endif // IO_DGF_H
