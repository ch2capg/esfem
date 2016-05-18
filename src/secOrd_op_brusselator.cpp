/*! \file secOrd_op_brusselator.cpp
    \brief Implementation of secOrd_op_brusselator.h

     Revision history
     --------------------------------------------------

          Revised by Christian Power May 2016
          Originally written by Christian Power
               (power22c@gmail.com) Februar 2016

     \author Christian Power 
     \date 18. May 2016
     \copyright Copyright (c) 2016 Christian Power. All rights reserved.
 */

#include <stdexcept>
#include <config.h>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/oemsolver.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlsolver.hh>
#include "secOrd_op_brusselator.h"
#include "secOrd_op_brusselator_impl.h"
#include "io_parameter.h"
#include "grid.h"

using namespace std;
using FE_function = Esfem::Grid::Scal_FEfun::Dune_FEfun;
using CG_method = Dune::Fem::CGInverseOperator<FE_function>;

using Linear_operator = Dune::Fem::ISTLLinearOperator<FE_function, FE_function>;
using GMRes
= Dune::Fem::ISTLGMResOp<FE_function, Linear_operator>;
// = Dune::Fem::ISTLCGOp<FE_function, LinearOperatorType>;

struct Esfem::SecOrd_op::Brusselator::Data{
  bool owns;
  FE_function* fef1_ptr {nullptr};
  FE_function* fef2_ptr {nullptr};
  const FE_function& fef1_ref;
  const FE_function& fef2_ref;
  FE_function tmp_var;
  Brusselator_op bruss_op;
  Linear_operator bruss_matrix;
  CG_method bruss_cg;
  GMRes bruss_gmres;
  explicit Data(const Io::Parameter&, const Grid::Grid_and_time&,
		const Growth);
  explicit Data(const Io::Parameter&, const Grid::Grid_and_time&,
		const Growth, const FE_function&, const FE_function&);
  ~Data();
  Data(const Data&) = delete;
  Data& operator=(const Data&) = delete;
};

Esfem::SecOrd_op::Brusselator::
Brusselator(const Io::Parameter& p, const Grid::Grid_and_time& gt, const Growth type)
try : d_ptr {make_unique<Data>(p, gt, type)}
{}
catch(const std::exception&){
  std::throw_with_nested(std::logic_error
			 {"Error in constructor of Brusselator."});
 }
 catch(...){
   throw std::logic_error {"Unkown error in constructor of Brusselator."};
 }

Esfem::SecOrd_op::Brusselator::
Brusselator(const Io::Parameter& p, const Grid::Grid_and_time& gt, const Growth type,
	    const Grid::Scal_FEfun& fef1, const Grid::Scal_FEfun& fef2)
try{
  const FE_function& fef1_ref = fef1;
  const FE_function& fef2_ref = fef2;
  d_ptr = make_unique<Data>
    (p, gt, type, fef1_ref, fef2_ref);
 }
 catch(const std::exception&){
   std::throw_with_nested(std::logic_error
			  {"Error in constructor of Brusselator."});
  }
 catch(...){
   throw std::logic_error {"Unkown error in constructor of Brusselator."};
 }

Esfem::SecOrd_op::Brusselator::~Brusselator() 
= default;
// {
//   delete d_ptr;
//   d_ptr = nullptr;
// #ifdef DEBUG
//   std::cerr << "~Brusselator(): delete d_ptr.\n";
// #endif
// }
void Esfem::SecOrd_op::Brusselator::
solve(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const{
  const FE_function& rhs_ref = rhs;
  FE_function& lhs_ref = lhs;

  d_ptr -> bruss_cg(rhs_ref, lhs_ref);
  // d_ptr -> bruss_op.jacobian(rhs_ref, d_ptr -> bruss_matrix);
  // d_ptr -> bruss_gmres(rhs_ref, lhs_ref);
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
add_massMatrixConstOne_to(Grid::Scal_FEfun& fef) const{
  FE_function& fef_ref = fef;
  d_ptr -> bruss_op.massMatrix_constOne(d_ptr -> tmp_var);
  fef_ref.axpy(1., d_ptr -> tmp_var);
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

// testing
void Esfem::SecOrd_op::Brusselator::
operator()(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const{
  const FE_function& rhs_ref = rhs;
  FE_function& lhs_ref = lhs;
  d_ptr -> bruss_op(rhs,lhs);
  // d_ptr -> bruss_op.jacobian(rhs_ref, d_ptr -> bruss_matrix);
  // d_ptr -> bruss_matrix(rhs_ref, lhs_ref);
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
    tmp_var {"tmp_var", gt.fe_space()},
    bruss_op {p, gt, type, fef1_ref, fef2_ref},
    bruss_matrix {"assempled elliptic operator", gt.fe_space(), gt.fe_space()},
    bruss_cg {bruss_op, p.eps(), p.eps()},
    bruss_gmres {bruss_matrix, p.eps(), p.eps()}
  {}
Esfem::SecOrd_op::Brusselator::Data::
Data(const Io::Parameter& p, const Grid::Grid_and_time& gt,
     const Growth type, const FE_function& fef1,
     const FE_function& fef2)
  : owns {false}, fef1_ref {fef1}, fef2_ref {fef2},
    tmp_var {"tmp_var", gt.fe_space()},
    bruss_op {p, gt, type, fef1_ref, fef2_ref},
    bruss_matrix {"assempled elliptic operator", gt.fe_space(), gt.fe_space()},
    bruss_cg {bruss_op, p.eps(), p.eps()},
    bruss_gmres {bruss_matrix, p.eps(), p.eps()}
  {}
Esfem::SecOrd_op::Brusselator::Data::
~Data(){
  if(owns){
    delete fef1_ptr;
    delete fef2_ptr;
    fef1_ptr = fef2_ptr = nullptr;
  }
}

/*! Log:
 */
