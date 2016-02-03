/*! \file operator_initData.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     Implementation details for operator_initData.h
     Created by Christian Power on 25.01.2016
     Copyright (c) 2016 Christian Power. All rights reserved.
 */
#include "operator.h"
#include <cmath>
#include <dune/fem/operator/lagrangeinterpolation.hh>

class Initial_data
  : public Dune::Fem::Function<Discrete_function::Function_space, Initial_data>
{
public:
  using Base = Discrete_function::Function_space;
  using Domain = Base::DomainType;
  using Range = Base::RangeType;
  
  Initial_data() = delete;
  explicit Initial_data(const Dune::Fem::TimeProviderBase&);
  Initial_data(const Initial_data&) = delete;
  Initial_data(Initial_data&&) = delete;
  Initial_data& operator=(const Initial_data&) = delete;
  Initial_data& operator=(Initial_data&&) = delete;
  ~Initial_data() = default;

  void evaluate(const Domain&, Range&) const;
private:
  const Dune::Fem::TimeProviderBase& tp;
};

// ----------------------------------------------------------------------
// Implementation esfem.h

void Operator::interpolate_initData(const Dune::Fem::TimeProviderBase& tp,
				    Discrete_function::FE_function& fef){
  static const Initial_data initD {tp};
  Dune::LagrangeInterpolation<Discrete_function::FE_function>::
    interpolateFunction(initD, fef);
}

// ----------------------------------------------------------------------
// Implementation Initial_data

Initial_data::Initial_data(const Dune::Fem::TimeProviderBase& tpb)
  : tp {tpb}
{}
void Initial_data::evaluate(const Domain& d, Range& r) const{
  static_assert(Domain::dimension == 3, "Bad domain dimension.");
  static_assert(Range::dimension == 1, "Bad range dimension.");
  const double x = d[0];
  const double y = d[1];
  const double t = tp.time();
  r = std::exp(-6.*t)*x*y;
}

/*! Log:
 */



