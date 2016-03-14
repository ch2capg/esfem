/*! \file grid_GridAndTime_impl.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 25. Februar 2016

     Implementation details for grid_GridAndTime_impl.h
     Created by Christian Power on 25.02.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <iterator>
#include "grid_GridAndTime_impl.h"

using namespace std;
using Esfem::Impl::Evolving_grid;
using Node = Evolving_grid::Node;
using Nodes_key = Evolving_grid::Nodes_key;
using Map = Evolving_grid::Map;

Esfem::Impl::Evolving_grid::Evolving_grid(const std::string& filename)
try : original_vertices {get_vertexList(filename)},
  digit_precision {get_relevant_precision(original_vertices)},
  vertices_map {identity_map(original_vertices, digit_precision)}
{}
catch(const std::exception&){
  std::throw_with_nested
    (logic_error {"Error in constructor of Evolving_grid."});
 }
 catch(...){
   throw logic_error {"Unknown error in constructor of Evolving_grid."};
 }

const Node& Esfem::Impl::Evolving_grid::operator[](const Node& node) const try{
  const auto node_as_string = node_to_string(node, digit_precision);
  const Node* n_ptr {nullptr};
  try{ n_ptr = &vertices_map.at(node_as_string); }
  catch(...){
    throw runtime_error {"Bad input: " + node_as_string};
  }
  return *n_ptr;
 }
 catch(const std::exception&){
   std::throw_with_nested
     (runtime_error {"Error in Evolving_grid::operator[]."});
 }
 catch(...){
   throw runtime_error {"Unkown error in Evolving_grid::operator[]."};
 }

Esfem::Impl::Evolving_grid& Esfem::Impl::Evolving_grid::
operator=(const Grid::Vec_FEfun& new_nodes) try{
  using std::cbegin;
  using std::cend;
    
  constexpr std::size_t dim = 3;
  static_assert(dim == Node::dimension, "Bad dimension.");
  const auto no_vertices = original_vertices.size();  
  const auto dof_times_dim = new_nodes.size();
  if(no_vertices * dim != dof_times_dim) throw logic_error {"Bad dimension."};
  
  update_vertices(cbegin(new_nodes), cend(new_nodes), cbegin(original_vertices), dim);
 }
 catch(const std::exception&){
   std::throw_with_nested
     (runtime_error {"Error in Evolving_grid::operator=()."});
 }
 catch(...){
   throw runtime_error {"Unkown error in Evolving_grid::operator=()."};
 }

// ----------------------------------------------------------------------
// helper functions

Nodes_key Esfem::Impl::
get_vertexList(const std::string& filename) try{
  using std::cbegin;
  using std::cend;
  const auto whole_file = dgfFile_to_vec(filename);

  auto first = dgf_find_vertex(cbegin(whole_file), cend(whole_file), filename);
  auto last = dgf_close_list(first, cend(whole_file), filename);  
  return {first, last};
 }
 catch(const std::exception&){
   std::throw_with_nested
     (logic_error {"Error in get_vertexList()."});
 }
 catch(...){
   throw runtime_error {"Unknown error in get_vertexList()."};
 }
size_t Esfem::Impl::
get_relevant_precision(const vector<string>& vertex_list) try{
  if(vertex_list.empty()) throw logic_error {"Empty input"};
  const auto precision_vec = create_precision_vec(vertex_list);
  return Esfem::Impl::max(precision_vec);
 }
 catch(const std::exception&){
   std::throw_with_nested
     (runtime_error {"Error in get_relevant_precision(const vector<string>&)"});
 }
 catch(...){
   throw runtime_error
   {"Unknown error in get_relevant_precision(const vector<string>&)"};
 }

Map Esfem::Impl::
identity_map(const Nodes_key& nodes_lines, const std::size_t precision){
  Map rv {};
  rv.reserve(nodes_lines.size());
  for(const auto& line : nodes_lines)    
    rv[normalize_node(line, precision)] = string_to_node(line);
  return rv;
}

std::size_t Esfem::Impl::get_relevant_precision(const std::string& number) try{
  const auto dot_iter = number.find(".");
  const auto exp_iter = number.find("e");
  const auto dist = exp_iter - dot_iter ;
  if(dist < 1) throw logic_error {"Input has not scientific format."};
  return dist;
 }
 catch(const std::exception&){
   std::throw_with_nested(logic_error
			  {"Error in get_relevant_precision(const string&)"});
 }
 catch(...){
   throw runtime_error {"Unknown error in get_relevant_precision(const string&)"};
 }
vector<size_t> Esfem::Impl::create_precision_vec(const vector<string>& vertex_list){
  vector<size_t> rv {};
  
  constexpr int world_dim {3};
  rv.reserve(vertex_list.size() * world_dim);
  // We assume that each line of vertex_list contains three numbers.
  
  for(const auto& line : vertex_list){
    istringstream interpreter {line};
    for(string s; interpreter >> s; )
      rv.push_back(get_relevant_precision(s));
  }
  return rv;
}

std::size_t Esfem::Impl::max(const std::vector<std::size_t>& v){
  const auto rv_ptr = max_element(cbegin(v), cend(v));
  if(rv_ptr == cend(v)) throw logic_error {"Error in max()"};
  return *rv_ptr;
}

std::string Esfem::Impl::node_to_string(const Node& node, const int precision) try{
  // assuming that the dgf file uses std::scientific as it should be.
  std::ostringstream oss;
  oss << std::scientific;
  for(std::size_t i = 0; i < Node::dimension - 1; ++i)
    if( !(oss << std::setprecision(precision) << node[i] << ' ') )
      throw runtime_error {"Node is too big."};
  oss << std::setprecision(precision) << node[Node::dimension - 1];

  // std::clog << "node_to_string()" << '\n'
  // 	    << "precision: " << precision << '\n'
  // 	    << "input: " << node[0] << ' ' << node[1] <<  ' ' << node[2] << std::endl;
  // std::clog <<  "output: " << oss.str() << std::endl;
  
  return oss.str();
 }
 catch(const std::exception&){
   std::throw_with_nested(runtime_error {"Error in node_to_string()."});
 }
 catch(...){
   throw runtime_error {"Unknown error in node_to_string()."};
 }

Node Esfem::Impl::string_to_node(const std::string& str_node) try{
  std::istringstream interpreter {str_node};
  Node rv;
  double value {0.};
  for(std::size_t it = 0; it < Node::dimension; ++it){
    if(interpreter >> value) rv[it] = value;
    else throw runtime_error("operator>>() failed!");
  }
  if( !(interpreter >> std::ws).eof()) // stuff left in stream
    throw runtime_error("Stuff left in node.");
  return rv;
 }
 catch(const std::exception&){
   std::throw_with_nested
     (runtime_error {"Error in string_to_node()."});
 }
 catch(...){
   throw runtime_error {"Unknown error in string_to_node()."};
 }

std::string Esfem::Impl::
normalize_node(const std::string& str_node, const int precision){
  return node_to_string(string_to_node(str_node), precision);
}

std::string Esfem::Impl::clean_whitespace(const std::string& line){
  std::string clean_line {};
  std::istringstream iss {line};
  for(std::string s; iss >> s; )
    clean_line += s + ' ';
  clean_line.pop_back();	// delete last whitespace
  return clean_line;
}

std::vector<std::string> Esfem::Impl::
dgfFile_to_vec(const std::string& filename) try{
  std::ifstream dgf_file {filename};
  if(!dgf_file) throw logic_error {"Could not open \"" + filename + "\"."};
  
  std::vector<std::string> lines;
  for(std::string line; std::getline(dgf_file, line); )
    lines.push_back(clean_whitespace(line));
  dgf_file.close();
  return lines;
 }
 catch(const std::exception&){
   std::throw_with_nested
     (runtime_error {"Error in dgfFile_to_vec()."});
 }
 catch(...){
   throw runtime_error {"Unknown error in dgfFile_to_vec()."};
 }

/*! Log:
 */
