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
#include "grid.h"
#include "grid_GridAndTime_impl.h"

//! \f$\R^3\f$
using Domain = Esfem::Grid::Deformation::Domain;
//! \f$\R^3\f$
using Range = Esfem::Grid::Deformation::Range;

static_assert(Esfem::Grid::Deformation::Domain::dimension == 3,
	      "Bad domain dimension.");
static_assert(Esfem::Grid::Deformation::Range::dimension == 3,
	      "Bad range dimension.");

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
  y[0] = r * x[0]; 
  y[1] = r * x[1]; 
  y[2] = r * x[2]; 
}

//! Dalquist test equation with \f$\lambda=1\f$
static inline void dalquist(const double t, const Domain& x, Range& y){
  const double factor = exp(t);
  y[0] = factor * x[0]; 
  y[1] = factor * x[1];
  y[2] = factor * x[2];
}

// ----------------------------------------------------------------------
// Implementaion of Deformation

struct Esfem::Grid::Deformation::Data{
  const Impl::Evolving_grid* eg_ptr {nullptr};
  const Dune::Fem::TimeProviderBase* tp_ptr {nullptr};
  Data() = default;
  Data(const Impl::Evolving_grid& eg) :eg_ptr {&eg} {}
};

Esfem::Grid::Deformation::Deformation()
  :d_ptr {std::make_unique<Data>()}
{}

Esfem::Grid::Deformation::Deformation(const Impl::Evolving_grid& eg)
  :d_ptr {std::make_unique<Data>(eg)}
{}

Esfem::Grid::Deformation::~Deformation() = default;

void Esfem::Grid::Deformation::
set_timeProvider(const Dune::Fem::TimeProviderBase& tp){
  d_ptr -> tp_ptr = &tp;
}
void Esfem::Grid::Deformation::evaluate(const Domain& x, Range& y) const{
  const double t = d_ptr -> tp_ptr->time();
  dalquist(t, x, y);

  // logistic_growth(t, x, y);

  // identity(x, y);

  // const auto eg = *(d_ptr -> eg_ptr);
  // y = eg[x];
}
