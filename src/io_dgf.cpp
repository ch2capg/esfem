/*! \file io_dgf.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Implementation details for io_dgf.h
     Created by Christian Power on 08.03.2016
     
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <fstream>
#include <sstream>
#include <vector>
#include "io_dgf.h"
#include "grid_GridAndTime_impl.h"

using namespace std;

using Esfem::Io::Dgf::Handler;
using Simplex = Handler::Simplex;
using Simplices = Handler::Simplices;
using Vertex = Handler::Vertex;
using Vertices = Handler::Vertices;
using Lines = Handler::Lines;
using Vec_FEfun = Esfem::Grid::Vec_FEfun;
using Scal_FEfun = Esfem::Grid::Scal_FEfun;

static constexpr int dim = 3;

// ----------------------------------------------------------------------
// Handler implementation

struct Handler::Data{
  size_t digit_precision {8};
  const size_t scal_dim_dof;
  const size_t vec_dim_dof;
  const Simplices simplices;
  Data(const string&);

  // --------------------------------------------------
  // Helper function
  // void to_dgfFile(const Vertices&, std::ostream&);
  // string from(const Vertex&);
};

// ------------------------------------------------------------
// Handler::Data Implementation

Handler::Data::Data(const string& fname)
  : digit_precision {get_precision(fname)},
    scal_dim_dof {calculate_dof(fname)},
    vec_dim_dof {dim * scal_dim_dof},
    simplices {create_simplices(fname)}
{}
// void Handler::Data::to_dgfFile(const Vertices& v, std::ostream& os){
//   for(const auto& vertex : v)
//     os << from(vertex, digit_precision) << endl;
// }
// string Handler::Data::from(const Vertex& v){
//   std::ostringstream oss {};
//   oss << scientific;
//   for(const auto& d : v)
//     if( !(oss << setprecision(digit_precision) << d << ' ') )
//       throw runtime_error {"Esfem::Io::Dgf::Handler::Data::from()."};
//   return oss.str();
// }

// ------------------------------------------------------------
// Handler member functions

Handler::Handler(const std::string& dgf_filename)
try : d_ptr {make_unique<Data>(dgf_filename)}
{}
 catch(const exception&){
   throw_with_nested
     (runtime_error {"Error in constructor of Esfem::Io::Handler."});
 }
 catch(...){
   throw runtime_error {"Unknown error in constructor of Esfem::Io::Handler."};
 }

Handler::~Handler() = default;

void Handler::write(const std::string& out_filename, const Scal_FEfun& fef) const try{
  if(fef.size() != no_of_nodes()) throw logic_error
    {"Finit element function has a bad degree of freedom dimension."};
  
  ofstream dgf_file {out_filename};
  if(!dgf_file) throw runtime_error {"Could not open file " + out_filename};

  dgf_file << "DGF" << std::endl;
  const auto vertices = Esfem::Io::Dgf::to_vertices(fef, 1); 
  dgf_file << "VERTEX" << std::endl;
  to_dgfFile(vertices, dgf_file, d_ptr -> digit_precision);
  dgf_file << "#" << '\n'
	   << "SIMPLEX" << std::endl;
  Esfem::Io::Dgf::to_dgfFile(d_ptr -> simplices, dgf_file, d_ptr -> digit_precision);
  dgf_file << "#" << std::endl;
  dgf_file << '#' << std::endl;
 }
 catch(const exception&){
   throw_with_nested
     (runtime_error {"Error in Handler::write()."});
 }
 catch(...){
   throw runtime_error {"Unknown error in Handler::write()."};
 }

void Handler::write(const std::string& out_filename, const Vec_FEfun& vfef) const try{
  if(vfef.size() != d_ptr -> vec_dim_dof) throw logic_error
    {"Finit element function has a bad degree of freedom dimension."};
  
  ofstream dgf_file {out_filename};
  if(!dgf_file) throw runtime_error {"Could not open file " + out_filename};

  dgf_file << "DGF" << std::endl;
  const auto vertices = Esfem::Io::Dgf::to_vertices(vfef); 
  dgf_file << "VERTEX" << std::endl;
  to_dgfFile(vertices, dgf_file, d_ptr -> digit_precision);
  dgf_file << "#" << '\n'
	   << "SIMPLEX" << std::endl;
  to_dgfFile(d_ptr -> simplices, dgf_file, d_ptr -> digit_precision);
  dgf_file << "#" << std::endl;
  dgf_file << '#' << std::endl;
 }
 catch(const exception&){
   throw_with_nested
     (runtime_error {"Error in Handler::write()."});
 }
 catch(...){
   throw runtime_error {"Unknown error in Handler::write()."};
 }

void Handler::read(const string& in_filename, Vec_FEfun& vfef) const try{
  const auto vertices = create_vertices(in_filename);

  // We assume that vertices[i] always has the same size.
  const auto dim_of_vertices = vertices.size() * vertices.front().size();
  const auto dim_of_fef = vfef.size();
  if(dim_of_fef != dim_of_vertices) throw logic_error
    {"Input file is incompatible with finite element function."};
  
  vfef << vertices;  
 }
 catch(const exception&){
   throw_with_nested
     (runtime_error {"Error in handler::read()."});
 }
 catch(...){
   throw runtime_error {"Unknown error in Handler::write()."};
 }

void Handler::read(const string& in_filename, Scal_FEfun& fef) const try{
  const auto vertices = create_vertices(in_filename);

  // We assume that vertices[i] always has the same size.
  const auto dim_of_vertices = vertices.size() * vertices.front().size();
  const auto dim_of_fef = fef.size();
  if(dim_of_fef != dim_of_vertices) throw logic_error
    {"Input file is incompatible with finite element function."};
  
  fef << vertices;  
 }
 catch(const exception&){
   throw_with_nested
     (runtime_error {"Error in handler::read()."});
 }
 catch(...){
   throw runtime_error {"Unknown error in Handler::write()."};
 }


size_t Handler::no_of_nodes() const noexcept{
  return d_ptr -> scal_dim_dof;
}

// ----------------------------------------------------------------------
// Constructor helper functions

Simplices Esfem::Io::Dgf::create_simplices(const std::string& filename) try{
  Simplices rv {};
  const auto simplices_list = create_stringSimplices(filename);
  rv.reserve(simplices_list.size());
  for(const auto& line : simplices_list)
    rv.emplace_back(to<Simplex>(line));
  return rv;
 }
 catch(const exception&){
   throw_with_nested
     (runtime_error {"Error in create_simplices()."});
 }
Lines Esfem::Io::Dgf::create_stringSimplices(const string& filename) try{
  const auto lines = Esfem::Impl::dgfFile_to_vec(filename);
  auto first = Esfem::Io::Dgf::find(cbegin(lines), cend(lines), "SIMPLEX");
  auto last = Esfem::Io::Dgf::close_list(first, cend(lines), "#");
  return {first, last};
 }
 catch(const exception&){
   throw_with_nested
     (runtime_error {"Error in create_stringSimplices()."});
 }
size_t Esfem::Io::Dgf::get_precision(const string& dgf_file){
  const auto vertices_list = create_stringVertices(dgf_file);
  return Esfem::Impl::get_relevant_precision(vertices_list);
}
size_t Esfem::Io::Dgf::calculate_dof(const string& dgf_file){
  const auto lines = Esfem::Impl::dgfFile_to_vec(dgf_file);
  auto first = Esfem::Io::Dgf::find(cbegin(lines), cend(lines), "VERTEX");
  auto last = Esfem::Io::Dgf::close_list(first, cend(lines), "#");
  return distance(first, last);
}
Vertices Esfem::Io::Dgf::create_vertices(const string& filename) try{
  Vertices rv {};
  const auto vertices_list = create_stringVertices(filename);
  rv.reserve(vertices_list.size());
  for(const auto& line : vertices_list)
    rv.emplace_back(to<Vertex>(line));
  return rv;  
 }
 catch(const exception&){
   throw_with_nested
     (runtime_error {"Error in Esfem::Io::Dgf::create_vertices()."});
 }
 catch(...){
   throw runtime_error {"Unknown error in Esfem::Io::Dgf::create_vertices()."};
 }
Lines Esfem::Io::Dgf::create_stringVertices(const string& filename){
  const auto lines = Esfem::Impl::dgfFile_to_vec(filename);
  auto first = Esfem::Io::Dgf::find(cbegin(lines), cend(lines), "VERTEX");
  auto last = Esfem::Io::Dgf::close_list(first, cend(lines), "#");
  return {first, last};
}

/*! Log:
 */
