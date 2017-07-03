/*! \file grid_ode.cpp
    \brief Implementation for grid_ode.h

     Revision history
     --------------------------------------------------

          Revised by Christian Power Juli 2017
          Originally written by Christian Power
               (power22c@gmail.com) Juli 2017

    \author Christian Power 
    \date 02. Juli 2017
    \copyright Copyright (c) 2017 Christian Power.  All rights reserved.
 */

#include "grid_ode.h"
#include <boost/numeric/odeint.hpp>
#include <Eigen/Dense>
#include <implicit_euler.h>
#include <bdf.h>
#include "alePaper_impl.h"

using namespace std;
namespace ode = Esfem::ODE;
namespace ss = Esfem::singleStep_integrator;
using namespace boost::numeric::odeint;
using namespace Eigen;
using namespace NUMERIK;
using Esfem::Grid::Deformation;
//! \f$\R^3\f$
using Domain = Esfem::Grid::Deformation::Domain;
//! \f$\R^3\f$
using Range = Esfem::Grid::Deformation::Range;

static_assert(Esfem::Grid::Deformation::Domain::dimension == 3,
	      "Bad domain dimension.");
static_assert(Esfem::Grid::Deformation::Range::dimension == 3,
	      "Bad range dimension.");

//! Normal and ALE movement 
/*! \sa Brusselator_scheme::ale_normalMovement(), 
  Brusselator_scheme::ale_aleMovement()
 */
namespace alePaper{
  // ODE part
  // --------------------------------------------------
  struct LinSolv{
    Vector3d linear_solve(const Matrix3d& A, const Vector3d& b) const {
      return A.fullPivLu().solve(b);  
    }
  };
  struct MyNorm{
    double norm(const Vector3d& v) const {
      return v.norm();
    }
  };
  struct MyFun{
    double t;
    void set_time(const double t_in){
      t = t_in;
    }
    Vector3d evaluate(Vector3d vec) const{
      vector_type vec_in(3);
      for(int i{}; i < dimWorld; ++i) vec_in(i) = vec(i);
      vector_type vec_out(3);
      levelset_grad lg;
      lg(vec_in, vec_out, t);
      return Vector3d(vec_out(0), vec_out(1), vec_out(2));
    }
    Matrix3d jacobian(Vector3d vec) const{
      vector_type vec_in(3);
      for(int i{}; i < dimWorld; ++i) vec_in(i) = vec(i);
      matrix_type mat_out(3,3);
      jac_levelset_grad jlg;
      jlg(vec_in, mat_out, t);
      Matrix3d rv;
      for(int i{}; i < dimWorld; ++i)
	for(int j{}; j < dimWorld; ++j)
	  rv(i,j) = mat_out(i,j);
      return rv;
    }
  };  
  struct MyIdentityMat{
    Matrix3d identity_matrix() const { return Matrix3d::Identity(); }
  };

  //! Home made implicit Euler routine
  struct capg_implicitEuler : ss::ss_int{
    //! Use my integrator 
    void integrate(const Domain& x, Range& y, 
		   const double t, const double dT) const override{
      Vector3d x_0(x[0], x[1], x[2]);
      array<double, 2> time_vec {t-dT, t};
      ImplicitEuler<LinSolv, MyNorm, MyFun, MyIdentityMat, Matrix3d, Vector3d>
	ie {NewtonParameters<> {1e-6, // tolerance
	    50, // max Newton steps
	    10, // max damping steps
	    .5, // damping factor for steps
	    .5 // Guaranteed decrease factor
	    }
      };
      Vector3d ans = ie.integrate(time_vec.begin(), time_vec.end(), x_0);
      for(int i{0}; i < dimWorld; ++i) y[i] = ans(i);      
    }
  };
  
  // ToDo: Can be used for another routine
  // Arbitrary-Lagrangian-Eulerian movement of a problem
  // void normalMovement(const Domain&, Range& y, 
  // 		      const double t, const double dT) noexcept{
  //   // vector_type start_val(dimWorld);
  //   // for(int i{0}; i < dimWorld; ++i) start_val[i] = y[i];
  //   // integrate_n_steps(implicit_euler<double> {}, 
  //   // 		      make_pair(levelset_grad{}, jac_levelset_grad {}),
  //   // 		      start_val, t, dT, 1);
  //   // for(int i{0}; i < dimWorld; ++i) y[i] = start_val[i];
  // }
} // namespace alePaper

unique_ptr<ode::ss_int> ss::capg_ie(){
  return make_unique<alePaper::capg_implicitEuler>();
}
