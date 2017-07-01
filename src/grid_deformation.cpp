/*! \file grid_deformation.cpp
    \brief Implementation of Deformation class

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     Idea
     --------------------------------------------------

     Explicit flow functions are coded as inline functions.

    \author Christian Power
    \date 23. April 2016
    \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#include <cmath>
#include <numeric>
#include "grid.h"
#include "grid_GridAndTime_impl.h"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace boost::numeric::odeint;
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;
template<class T>
using matrix_column = boost::numeric::ublas::matrix_column<T>;
using Esfem::Grid::Deformation;
//! \f$\R^3\f$
using Domain = Esfem::Grid::Deformation::Domain;
//! \f$\R^3\f$
using Range = Esfem::Grid::Deformation::Range;

static_assert(Esfem::Grid::Deformation::Domain::dimension == 3,
	      "Bad domain dimension.");
static_assert(Esfem::Grid::Deformation::Range::dimension == 3,
	      "Bad range dimension.");

namespace{
  //! x-axis is a bouncing sinus
  void bouncing_ellipsoid(const Domain& x, Range& y, double t) noexcept{
    const double at = 1. + .25 * sin(M_PI * 2.0 * t);
    y[0] = at * x[0];
    y[1] = x[1];
    y[2] = x[2];
  }
  //! final picture looks like a baseball bat
  void baseball_bat(const Domain& x, Range& y, double t) noexcept{
    // const double newtime = 
    //   (t < 0.5 ? 0.0 : (t < 0.75 ? 4*(t - 0.5) : 1.0));
    const double newtime = std::min( t, 1.0 );

    const double r1 = std::abs( x[ 0 ] );
    const double target = (1.0 - (r1*r1))*((r1*r1) + 0.05) 
      + (r1*r1)*sqrt(1.0 - (r1*r1)); 

    const double r2 = std::sqrt( x[1]*x[1] + x[2]*x[2] );
    const double factor = std::exp( -2*newtime )*r2 
      + (1.0 - std::exp( -2*newtime ))*target;

    y[ 0 ] = 
      1.0 * x[ 0 ] + newtime*(x[ 0 ] > 0 ? 1.0 : -.5 )*x[ 0 ];
    y[ 1 ] = factor * x[ 1 ] / (r2 + 0.000001);
    y[ 2 ] = factor * x[ 2 ] / (r2 + 0.000001);
  }
  //! Ball is shrinking
  void linear_decrease(const Domain& x, Range& y, double t) noexcept{
    y = x;
    y *= 1. - t/2.;
  }
}
//! Normal and ALE movement 
/*! \sa Brusselator_scheme::ale_normalMovement(), 
  Brusselator_scheme::ale_aleMovement()
 */
namespace alePaper{
  //! Arbitrary-Lagrangian-Eulerian movement of a problem
  void aleMovement(const Domain& x, Range& y, double t) noexcept{
    // auto G = [](double s){ return 200. * s * (s - 199./200.); };
    auto L = [](double t){ return 1. + .2 * sin(4.*M_PI*t); };
    auto K = [](double t){ return .1 + .05 * sin(2.*M_PI*t); };
    y[0]= x[0] * K(t)/K(0);
    y[1]= x[1] * K(t)/K(0);
    y[2]= x[2] * L(t)/L(0);
  }
  // ODE part
  // --------------------------------------------------
  //! Dimension of the world 
  constexpr int dimWorld {3};
  //! Step size for finite difference
  constexpr double fd_step {1e-6};
  
  //! Zero set is the manifold
  /*! \f$ d(x,t) := x_1^2 + x_2^2 + K(t)^2 G\Big( \frac{x_3^2}{L(t)^2} \Big)-K(t)^2. \f$
   */
  struct levelset_function{
    //! Apply it
    /*! \pre x.size() == 3*/
    double operator()(const vector_type& x, const double t) const noexcept{
      auto G = [](double s){ return 200. * s * (s - 199./200.); };
      auto L = [](double t){ return 1. + .2 * sin(4.*M_PI*t); };
      auto K = [](double t){ return .1 + .05 * sin(2.*M_PI*t); };
      return x[0]*x[0] + x[1]*x[1] + K(t)*K(t)*
	G( x[2]*x[2] / (L(t) * L(t))) - K(t)*K(t);
    }
  };
  //! This is the vector field on the right-hand side
  /*! It reads as
      \f{equation*}
      \dot{a}_j = V_j \nu_j, \quad \textnormal{where} 
      \quad V_j=\frac{-\partial_t d(a_j,t)}{|\nabla d(a_j,t)|}, 
      \quad \nu_j=\frac{\nabla d(a_j,t)}{|\nabla d(a_j,t)|}.
      \f}
   */
  struct levelset_grad{
    //! Use finite differences to calculate the gradient
    void update_spatial_grad(const vector_type& x,
			     const double t, vector_type& dx) const noexcept{
      levelset_function lf;
      vector_type xh {x};
      for(int i {0}; i < dimWorld; ++i){
	xh[i] += fd_step;
	dx[i] = (lf(xh, t) - lf(x, t))/fd_step;
	xh[i] = x[i];
      }
      cout << "levelset_grad::update_spatial_grad: " << dx << endl;
    }
    //! Use a finite difference to calculate the time derivative 
    double time_grad(const vector_type& x, const double t) const noexcept{
      cout << "levelset_grad::time_grad()" << endl;
      levelset_function lf;
      return (lf(x, t + fd_step) - lf(x,t))/fd_step;
    }
    //! Calculate the vector field
    void operator()(const vector_type& x, vector_type& dxdt, 
		    const double t) const noexcept{
      vector_type grad_d(dimWorld);      
      update_spatial_grad(x, t, grad_d);
      auto ddist_dt = time_grad(x, t);
      double norm_d = norm_2(grad_d);
      dxdt = - ddist_dt / norm_d * grad_d / norm_d;
      cout << "levelset_grad: " << dxdt << endl;
    }
  };
  //! Calculate Jacobian of the vector field above
  struct jac_levelset_grad{
    //! Calculate the Jacobian
    void operator()(const vector_type& x, matrix_type& dxdt, const double t) 
      const noexcept{
      cout << "jac_levelset_grad() begin" << endl;
      static vector_type fxh(dimWorld);
      static vector_type fx(dimWorld);
      vector_type xh {x};
      levelset_grad lgrad;
      for(int j{0}; j < dimWorld; ++j){
	xh(j) += fd_step;
	lgrad(x, fx, t);
	lgrad(xh, fxh, t);
	matrix_column<matrix_type> dxdt_j(dxdt, j);
	dxdt_j = (fxh - fx)/fd_step;
	xh(j) = x(j);
      }
      cout << dxdt << endl;
    }
  };
  //! Arbitrary-Lagrangian-Eulerian movement of a problem
  void normalMovement(const Domain&, Range& y, 
		      const double t, const double dT) noexcept{
    vector_type start_val(dimWorld);
    for(int i{0}; i < dimWorld; ++i) start_val[i] = y[i];
    integrate_n_steps(implicit_euler<double> {}, 
		      make_pair(levelset_grad{}, jac_levelset_grad {}),
		      start_val, t, dT, 1);
    cout << "Integrate once: " << start_val << endl;
    for(int i{0}; i < dimWorld; ++i) y[i] = start_val[i];
  }
} // namespace alePaper


//! \f$id\colon \R^3 \to \R^3\f$
static inline void identity(const Domain& x, Range& y) noexcept{  
  y[0] = x[0]; 
  y[1] = x[1]; 
  y[2] = x[2]; 
}


//! \f$r(t) = \frac{r_{end} r_0}{r_{end} e^{-kt} + r_0 (1 - e^{-kt})}\f$
/*! \param t Current time
  \param x Point from the initial surface
  \retval y \f$ y = r(t) x\f$
  \pre Initial surface is a sphere.
 */
static inline void logistic_growth(const double t, const Domain& x, Range& y) noexcept{
  const double r_end = 2., r0 = 1., k = .5; // logistic function parameter
  const double r = r_end * r0 / (r_end*exp(-k*t) + r0*(1-exp(-k*t)));
  y = x;
  y *= r;
}

//! Dalquist test equation with \f$\lambda=1\f$
static inline void dalquist(const double t, const Domain& x, Range& y){
  const double factor = exp(-t);
  y = x;
  y *= factor;
}

//! \f$ R(t) = \sqrt{ R_0^2 - 2nt}\f$
static inline void mcf_sphere(const double t, const Domain& x, Range& y){
  const double norm_square
    // = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    // = 1.;
    = std::inner_product(&x[0], &x[0] + Domain::dimension, &x[0], 0.);
  const double factor = sqrt( norm_square - 2 * Esfem::Grid::grid_dim() * t);
  y = x;
  y *= factor;
}

// ----------------------------------------------------------------------
// Implementaion of Deformation

struct Esfem::Grid::Deformation::Data{
  // Old hash map
  // Impl::Evolving_grid eg_ptr;
  //! New hash map
  Impl::hash::grid hg;
  //! I do not assume ownership
  const Dune::Fem::TimeProviderBase* tp_ptr {nullptr};
  //! All pointer to zero
  Data() = default;
  //! Use hash map
  Data(const std::string& fname) :hg {fname} {}
};

Esfem::Grid::Deformation::Deformation() :d_ptr {std::make_unique<Data>()} {}
Esfem::Grid::Deformation::Deformation(const std::string& fname)
  :d_ptr {std::make_unique<Data>(fname)} {}
Esfem::Grid::Deformation::~Deformation() = default;

void Esfem::Grid::Deformation::
set_timeProvider(const Dune::Fem::TimeProviderBase& tp){
  d_ptr -> tp_ptr = &tp;
}
void Esfem::Grid::Deformation::evaluate(const Domain& x, Range& y) const{
  const double t = d_ptr -> tp_ptr->time();

  // mcf_sphere(t, x, y);

  // dalquist(t, x, y);

  // logistic_growth(t, x, y);

  // identity(x, y);
  // bouncing_ellipsoid(x, y, t);
  // baseball_bat(x, y, t);
  // linear_decrease(x, y, t);

  // alePaper::aleMovement(x, y, t);
  y = d_ptr->hg[x];
  alePaper::normalMovement(x, y, t, d_ptr->tp_ptr->deltaT());
  
  // y = d_ptr->hg[x];

  // const auto eg = *(d_ptr -> eg_ptr);
  // y = eg[x];
}
Deformation& Deformation::operator=(const Vec_FEfun& rhs){
  d_ptr->hg = rhs;
  return *this;
}
