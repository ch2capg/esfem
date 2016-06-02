/*! \file io_l2h1Calculator.cpp
    \brief Implementation of io_l2h1Calculator.h

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     \author Christian Power
     \date 27. April 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#include <config.h>
#include <stdexcept>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dassert.h>
#include "io_l2h1Calculator.h"
#include "grid.h"


using namespace std;
using Esfem::Io::L2H1_calculator;

//! \f$f\colon \R^3 \to \R\f$
using FEfun = Esfem::Grid::Scal_FEfun::Dune_FEfun;
//! \f$f\colon \R^3 \to \R^3\f$
using Vec_FEfun = Esfem::Grid::Vec_FEfun::Dune_FEfun;
//! \f$L^2\f$-norm
using L2_norm = Dune::Fem::L2Norm<Esfem::Grid::Grid_and_time::Grid_part>;
//! \f$H^1\f$-norm
using H1_norm = Dune::Fem::H1Norm<Esfem::Grid::Grid_and_time::Grid_part>;

// ----------------------------------------------------------------------
// Helper template function

template<typename FEF, typename Norm>
double norm_err_helper(const Norm& n, const FEF& u1, const FEF& u2){
  using Dfef = typename FEF::Dune_FEfun;
  const Dfef& du1 = u1;
  const Dfef& du2 = u2;
  return n.distance(du1, du2);
}

// ----------------------------------------------------------------------
// Actual implementations 

//! %Data members of L2H1_calculator
struct L2H1_calculator::Data{
  //! Dune \f$L^2\f$-norm functor
  L2_norm l2;
  //! Dune \f$H^1\f$-norm functor
  H1_norm h1;
  //! First vector valued finite element function
  /*! __Invariant:__ `size() == worlddim`*/
  vector<FEfun> vec_fefun1;
  //! Second vector valued finite element function
  /*! \copydetails vec_fefun1 */
  vector<FEfun> vec_fefun2;
  //! Get grid
  /*! \post Grid_and_time must outlive this object. 
   \note I could add try and catch block, to avoid cryptic 
   error messages. */
  Data(const Grid::Grid_and_time& gt);
  //! Change value of vec_fefun1
  /*! \copydetails assign() */
  void assign_v1(const Vec_FEfun& v){ assign(v, vec_fefun1); }
  //! Change_value of vec_fefun2
  /*! \copydetails assign() */
  void assign_v2(const Vec_FEfun& v){ assign(v, vec_fefun2); }
private:
  //! Constructor helper
  void init_vec(vector<FEfun>& v, const Grid::Grid_and_time::FE_space& fes);
  //! Actual implementation of assign_v1() and assign_v2()
  /*! \post Dimension of Range should be `world_dim()`
    and the underlying grid should be the same. */
  void assign(const Vec_FEfun& v, vector<FEfun>& vec);
};

L2H1_calculator::Data::Data(const Grid::Grid_and_time& gt) 
:l2 {gt.grid_part()}, h1{gt.grid_part()} {
  init_vec(vec_fefun1, gt.fe_space());
  init_vec(vec_fefun2, gt.fe_space());
}
void L2H1_calculator::Data::init_vec
(vector<FEfun>& v, const Grid::Grid_and_time::FE_space& fes){
  constexpr auto d = Grid::world_dim();
  v.reserve(d);
  for(int it = 0; it < d; ++it){
    v.emplace_back("vec_fefun", fes);
    v.back().clear();
  }
}
void L2H1_calculator::Data::assign
(const Vec_FEfun& v, vector<FEfun>& vec){
  constexpr auto d = Grid::world_dim();
  size_t no_nodes = v.size() / d; // assuming the post condition
  Assert::dynamic(v.size() == no_nodes * d, 
		  Assert::compose(__FILE__, __LINE__, "assign()"));
  for(int jt = 0; jt < d; ++jt){
    auto vfef_ptr = v.dbegin();
    auto fef_ptr = vec[jt].dbegin();
    vfef_ptr += jt;
    for(size_t it = 0; it < no_nodes; ++it, vfef_ptr += d, ++fef_ptr)
      *fef_ptr = *vfef_ptr;
  }
}


L2H1_calculator::L2H1_calculator(const Grid::Grid_and_time& gt)
  :d_ptr {make_unique<Data>(gt)} {}
Esfem::Io::L2H1_calculator::~L2H1_calculator() = default;

double L2H1_calculator::l2_err(const Grid::Scal_FEfun& u, const Grid::Scal_FEfun& uN) const{
  const FEfun& u1 = u;
  const FEfun& u2 = uN;
  return d_ptr->l2.distance(u1, u2);
}
double L2H1_calculator::l2_err(const Grid::Vec_FEfun& u, const Grid::Vec_FEfun& uN) const{
  d_ptr->assign_v1(u);
  d_ptr->assign_v2(uN);
  auto rv = 0.;
  for(size_t it = 0; it < d_ptr->vec_fefun1.size(); ++it){
    rv += d_ptr->l2.distance(d_ptr->vec_fefun1[it], d_ptr->vec_fefun2[it]);
  }
  return rv;
  // return norm_err_helper(d_ptr->l2, u, uN);
}
double L2H1_calculator::h1_err(const Grid::Scal_FEfun& u, const Grid::Scal_FEfun& uN) const{
  const FEfun& u1 = u;
  const FEfun& u2 = uN;
  return d_ptr->h1.distance(u1, u2);
}
double L2H1_calculator::h1_err(const Grid::Vec_FEfun& u, const Grid::Vec_FEfun& uN) const{
  d_ptr->assign_v1(u);
  d_ptr->assign_v2(uN);
  auto rv = 0.;
  for(size_t it = 0; it < d_ptr->vec_fefun1.size(); ++it){
    rv += d_ptr->h1.distance(d_ptr->vec_fefun1[it], d_ptr->vec_fefun2[it]);
  }
  return rv;
  // return norm_err_helper(d_ptr->h1, u, uN);
}
