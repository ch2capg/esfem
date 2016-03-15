/*! \file io_l2h1Calculator.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     Implementation details for io_l2h1Calculator.h
     Created by Christian Power on 27.01.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <config.h>
#include <stdexcept>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include "io_l2h1Calculator.h"
#include "grid.h"

using namespace std;
using FEfun = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using L2_norm = Dune::Fem::L2Norm<Esfem::Grid::Grid_and_time::Grid_part>;
using H1_norm = Dune::Fem::H1Norm<Esfem::Grid::Grid_and_time::Grid_part>;

struct Esfem::Io::L2H1_calculator::Data{
  const FEfun& u;
  const FEfun& uN;
  L2_norm l2;
  H1_norm h1;
};

Esfem::Io::L2H1_calculator::
L2H1_calculator(const Grid::Grid_and_time& gt,
		const Grid::Scal_FEfun& exact_solution,
		const Grid::Scal_FEfun& numerical_solution)
  : d_ptr {make_unique<Data>
    (exact_solution, numerical_solution,
     L2_norm {gt.grid_part()}, 
     H1_norm {gt.grid_part()})}
{}
catch(const std::exception&){
  std::throw_with_nested(std::logic_error
			 {"Error in constructor of L2H1_calculator."});
 }
 catch(...){
   throw std::logic_error{"Unknown error in constructor of L2H1_calculator."};
 }
}

Esfem::Io::L2H1_calculator::
~L2H1_calculator() = default;
// {
//   delete d_ptr;
//   d_ptr = nullptr;
// #ifdef DEBUG
//   std::cerr << "~L2H1_calculator(): delete d_ptr\n";
// #endif
// }
double Esfem::Io::L2H1_calculator::l2_err() const{
  return d_ptr -> l2.distance(d_ptr -> u, d_ptr -> uN);
}
double Esfem::Io::L2H1_calculator::h1_err() const{
  return d_ptr -> h1.distance(d_ptr -> u, d_ptr -> uN);
}

/*! Log:
 */
