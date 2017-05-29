/*! \file io_paraview.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 31. Januar 2016

     Implementation details for io_paraview.h
     Created by Christian Power on 31.01.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <sstream>
#include <config.h>
#include <dune/fem/io/file/dataoutput.hh>
#include "io_paraview.h"
#include "io_parameter.h"
#include "grid.h"

using namespace std;

using Grid = Esfem::Grid::Grid_and_time::Grid;
using FEfun = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using Dune_tuple = Dune::tuple<FEfun*, FEfun*>;
using DataOutput = Dune::Fem::DataOutput<Grid, Dune_tuple>;

class DataOutputParameters 
  : public Dune::Fem::LocalParameter<Dune::Fem::DataOutputParameters,
				     DataOutputParameters> {
public:
  explicit DataOutputParameters(const Esfem::Io::Parameter&, const int = 0);
  DataOutputParameters(const DataOutputParameters&);
  std::string prefix() const;
private:
  std::string filename;
  int step;
};

// ----------------------------------------------------------------------
// Implementation of io_paraview.h

struct Esfem::Io::Paraview::Data{
  const Dune::Fem::TimeProviderBase& tp;
  Dune_tuple tuple;
  DataOutputParameters dop;
  DataOutput plotter;
  Data(const Parameter& p, const Grid::Grid_and_time& gt,
       Grid::Scal_FEfun& fef1, Grid::Scal_FEfun& fef2,
       const int refinement_label)
    : tp {gt.time_provider()},
      tuple {fef1.data_ptr(), fef2.data_ptr()}, // &static_cast<FEfun&>(fef2)
      dop {p, refinement_label}, plotter {gt.grid(), tuple, dop}
  {}
};

Esfem::Io::Paraview::Paraview(const Parameter& p , const Grid::Grid_and_time& gt,
			      Grid::Scal_FEfun& fef1, Grid::Scal_FEfun& fef2,
			      const int refinement_label)
try : d_ptr {make_unique<Data>(p, gt, fef1, fef2, refinement_label)}
{}
catch(const std::exception&){
  std::throw_with_nested(std::logic_error
			 {"Error in constructor of Paraview."});
 }
 catch(...){
   throw std::logic_error{"Unkown error in constructor of Paraview."};
 }
Esfem::Io::Paraview::~Paraview() = default;
// {
//   delete d_ptr;
//   d_ptr = nullptr;
// #ifdef DEBUG
//   std::cerr << "~Paraview(): delete d_ptr.\n";
// #endif
// }
void Esfem::Io::Paraview::write(){
  d_ptr -> plotter.write( d_ptr -> tp);
}
// ----------------------------------------------------------------------
// Internal implementation

DataOutputParameters::DataOutputParameters(const Esfem::Io::Parameter& p, const int i)
  : filename {p.paraview()}, step {i}
{}
DataOutputParameters::DataOutputParameters(const DataOutputParameters& other)
  : step {other.step}
{}
inline std::string DataOutputParameters::prefix() const {
  std::ostringstream s;
  s << filename << step << "-";
  return s.str();
}

/*! Log:
 */
