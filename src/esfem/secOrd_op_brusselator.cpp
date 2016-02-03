/*! \file secOrd_op_brusselator.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 01. Februar 2016

     Implementation details for secOrd_op_brusselator.h
     Created by Christian Power on 01.02.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <stdexcept>
#include <config.h>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/oemsolver.hh>
#include "secOrd_op_brusselator.h"
#include "secOrd_op_brusselator_impl.h"
#include "io_parameter.h"
#include "grid.h"

#ifdef DEBUG
#include <iostream>
#endif

using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using Solver = Dune::Fem::CGInverseOperator<FE_function>;

struct Esfem::SecOrd_op::Brusselator::Data{
  bool owns;
  FE_function* fef1_ptr {nullptr};
  FE_function* fef2_ptr {nullptr};
  const FE_function& fef1_ref;
  const FE_function& fef2_ref;
  Brusselator_op bruss_op;
  Solver bruss_solver;
  explicit Data(const Io::Parameter&, const Grid::Grid_and_time&,
		const Growth);
  explicit Data(const Io::Parameter&, const Grid::Grid_and_time&,
		const Growth, const FE_function&, const FE_function&);
  ~Data();
  Data(const Data&) = delete;
  Data& operator=(const Data&) = delete;
};

Esfem::SecOrd_op::Brusselator::
Brusselator(const Io::Parameter& p, const Grid::Grid_and_time& gt, const Growth type){
  try{
    d_ptr = new Data {p, gt, type};
  }
  catch(const std::exception&){
    std::throw_with_nested(std::logic_error
			   {"Error in constructor of Brusselator."});
  }
  catch(...){
    throw std::logic_error {"Unkown error in constructor of Brusselator."};
  }
}
Esfem::SecOrd_op::Brusselator::
Brusselator(const Io::Parameter& p, const Grid::Grid_and_time& gt, const Growth type,
	    const Grid::Scal_FEfun& fef1, const Grid::Scal_FEfun& fef2){
  try{
    const FE_function& fef1_ref = fef1;
    const FE_function& fef2_ref = fef2;
    d_ptr = new Data {p, gt, type, fef1_ref, fef2_ref};
  }
  catch(const std::exception&){
    std::throw_with_nested(std::logic_error
			   {"Error in constructor of Brusselator."});
  }
  catch(...){
    throw std::logic_error {"Unkown error in constructor of Brusselator."};
  }
}
Esfem::SecOrd_op::Brusselator::~Brusselator(){
  delete d_ptr;
  d_ptr = nullptr;
#ifdef DEBUG
  std::cerr << "~Brusselator(): delete d_ptr.\n";
#endif
}
void Esfem::SecOrd_op::Brusselator::
solve(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const{
  const FE_function& rhs_ref = rhs;
  FE_function& lhs_ref = lhs;
  d_ptr -> bruss_solver(rhs_ref, lhs_ref);
}
void Esfem::SecOrd_op::Brusselator::
mass_matrix(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const{
  const FE_function& rhs_ref = rhs;
  FE_function& lhs_ref = lhs;
  d_ptr -> bruss_op.mass_matrix(rhs_ref, lhs_ref);
}
void Esfem::SecOrd_op::Brusselator::massMatrix_constOne(Grid::Scal_FEfun& fef) const{
  FE_function& fef_ref = fef;
  d_ptr -> bruss_op.massMatrix_constOne(fef_ref);
}
void Esfem::SecOrd_op::Brusselator::
assign_firstArg_quadMassMatrix(const Grid::Scal_FEfun& fef){
  if(d_ptr -> owns){
    const FE_function& fef_ref = fef;
    d_ptr -> fef1_ptr -> assign(fef_ref);
  }
  else
    throw std::logic_error
    {"Error in Brusselator::assign_firstArg_quadMassMatrix().  "
	"*this has used a wrong constructor."};
}
void Esfem::SecOrd_op::Brusselator::
assign_secondArg_quadMassMatrix(const Grid::Scal_FEfun& fef){
  if(d_ptr -> owns){
    const FE_function& fef_ref = fef;
    d_ptr -> fef2_ptr -> assign(fef_ref);
  }
  else
    throw std::logic_error
    {"Error in Brusselator::assign_secondArg_quadMassMatrix().  "
	"*this has used a wrong constructor."};
}

// ----------------------------------------------------------------------
// Internal implementation

Esfem::SecOrd_op::Brusselator::Data::
Data(const Io::Parameter& p, const Grid::Grid_and_time& gt, const Growth type)
  : owns {true},
    fef1_ptr {new FE_function {"fef1", gt.fe_space()}},
    fef2_ptr {new FE_function {"fef2", gt.fe_space()}},
    fef1_ref {*fef1_ptr},
    fef2_ref {*fef2_ptr},
    bruss_op {p, gt, type, fef1_ref, fef2_ref},
    bruss_solver {bruss_op, p.eps(), p.eps()}
  {}
Esfem::SecOrd_op::Brusselator::Data::
Data(const Io::Parameter& p, const Grid::Grid_and_time& gt,
     const Growth type, const FE_function& fef1,
     const FE_function& fef2)
  : owns {false}, fef1_ref {fef1}, fef2_ref {fef2},
    bruss_op {p, gt, type, fef1_ref, fef2_ref},
    bruss_solver {bruss_op, p.eps(), p.eps()}
  {}
Esfem::SecOrd_op::Brusselator::Data::
~Data(){
  if(owns){
    delete fef1_ptr;
    fef1_ptr = nullptr;
    delete fef2_ptr;
    fef2_ptr = nullptr;
#ifdef DEBUG
    std::cerr << "~Brusselator::Data(): delete fef1_ptr, delete fef2_ptr.\n";
#endif
  }
}

/*! Log:
 */
