/*! \file esfem.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     This programm implements a basic expression calculator.
     Input from cin; output to cout.  The grammar for input is: etc.

     Created by Christian Power on 19.01.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef ESFEM_H
#define ESFEM_H 

#include "parameter.h"
#include "grid_evolution.h"
#include "fe_grid.h"
#include "discrete_function.h"
#include "operator.h"

/*
class DataOutputParameters 
  : public Dune::Fem::LocalParameter<Dune::Fem::DataOutputParameters,
				     DataOutputParameters> {
public:
  DataOutputParameters(const std::string name, const int step);
  DataOutputParameters(const DataOutputParameters&);
  std::string prefix() const;
private:
  std::string name_;
  int step_;
};
*/

// ----------------------------------------------------------------------
// Implementations

/*
Dune::GridPtr<Host_grid> init_dgf_file(){
  const std::string gridkey =
    Dune::Fem::IOInterface::defaultGridKey( Grid::dimension );
  std::cout << "gridkey: " << gridkey << std::endl;
  const std::string gridfile =
    Dune::Fem::Parameter::getValue< std::string >( gridkey );
  if( Dune::Fem::MPIManager::rank() == 0 )
    std::cout << "Loading macro grid: " << gridfile << std::endl;
  Dune::GridPtr< Host_grid > hostGrid {gridfile};
}

DataOutputParameters::DataOutputParameters(const std::string name, const int step)
  : name_(name), step_( step )
{}
DataOutputParameters::DataOutputParameters(const DataOutputParameters& other)
  : step_( other.step_ )
{}
inline std::string DataOutputParameters::prefix() const {
  std::stringstream s;
  s << name_ << step_ << "-";
  return s.str();
}
*/

#endif // ESFEM_H

/*! Log:
 */
