/*! \file alePaper_impl.h
    \brief Levelset function and some derivatives

     Revision history
     --------------------------------------------------

          Revised by Christian Power July 2017
          Originally written by Christian Power
               (power22c@gmail.com) July 2017

    \author Christian Power 
    \date 03. Juli 2017
    \copyright Copyright (c) 2017 Christian Power.  All rights reserved.
 */

#ifndef ALEPAPER_IMPL_GUARD
#define ALEPAPER_IMPL_GUARD 1

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace alePaper{
  typedef boost::numeric::ublas::vector< double > vector_type;
  typedef boost::numeric::ublas::matrix< double > matrix_type;
  template<class T>
  using matrix_column = boost::numeric::ublas::matrix_column<T>;

  //! Dimension of the world 
  constexpr int dimWorld {3};
  //! Step size for finite difference
  constexpr double fd_step {1e-6};

  //! Part of the levelset function 
  auto G = [](double s){ return 200. * s * (s - 199./200.); };
  //! Part of the levelset function 
  auto L = [](double t){ return 1. + .2 * sin(4.*M_PI*t); };
  //! Part of the levelset function 
  auto K = [](double t){ return .1 + .05 * sin(2.*M_PI*t); };

  //! Zero set is the manifold
  /*! \f$ d(x,t) := x_1^2 + x_2^2 + K(t)^2 G\Big( \frac{x_3^2}{L(t)^2} \Big)-K(t)^2.\f$ 
   */
  struct levelset_function{
    //! Apply it
    /*! \pre x.size() == 3*/
    double operator()(const vector_type& x, const double t) const noexcept{
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
    }
    //! Use a finite difference to calculate the time derivative 
    double time_grad(const vector_type& x, const double t) const noexcept{
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
    }
  };
  //! Calculate Jacobian of the vector field above
  struct jac_levelset_grad{
    //! Calculate the Jacobian
    void operator()(const vector_type& x, matrix_type& dxdt, const double t) 
      const noexcept{
      static vector_type fxh(dimWorld);
      static vector_type fx(dimWorld);
      vector_type xh {x};
      levelset_grad lgrad;
      for(int j{0}; j < dimWorld; ++j){
	xh(j) += fd_step;
	lgrad(x, fx, t);
	lgrad(xh, fxh, t);
	matrix_column<matrix_type> 
	  dxdt_j(dxdt, j);
	dxdt_j = (fxh - fx)/fd_step;
	xh(j) = x(j);
      }
    }
  };  
} // namespace alePaper

#endif // ALEPAPER_IMPL_GUARD
