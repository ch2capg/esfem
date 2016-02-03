#include "operator.h"
#include <cmath>
// #include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

using Discrete_function::FE_function;
using Geometry
= FE_function::DiscreteFunctionSpaceType::IteratorType::Entity::Geometry;
using Grid_part = FE_function::GridPartType;
using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;

class RHS_data
  : public Dune::Fem::Function<Discrete_function::Function_space, RHS_data>
{
public:
  using Base = Discrete_function::Function_space;
  using Domain = Base::DomainType;
  using Range = Base::RangeType;
  
  RHS_data() = delete;
  explicit RHS_data(const Dune::Fem::TimeProviderBase&);
  RHS_data(const RHS_data&) = delete;
  RHS_data(RHS_data&&) = delete;
  RHS_data& operator=(const RHS_data&) = delete;
  RHS_data& operator=(RHS_data&&) = delete;
  ~RHS_data() = default;

  void evaluate(const Domain&, Range&) const;
  Range operator()(const Domain&) const;
private:
  const Dune::Fem::TimeProviderBase& tp;
};

void matrixFree_assembly(const Dune::Fem::TimeProviderBase&, const Geometry&,
			 const Quadrature&, FE_function::LocalFunctionType&);

// ----------------------------------------------------------------------
// Implementation esfem.h

void Operator::assemble_RHS(const Dune::Fem::TimeProviderBase& tp,
			    Discrete_function::FE_function& fef){  
  fef.clear();
  const auto& df_space = fef.space();
  
  for(const auto& entity : df_space){
    const auto& geometry = entity.geometry();
    const Quadrature quad {entity, 4 * df_space.order() + 1};
    auto fef_local = fef.localFunction(entity);
    matrixFree_assembly(tp, geometry, quad, fef_local);
  }
  fef.communicate();
}

// ----------------------------------------------------------------------
// Implementation this file

RHS_data::RHS_data(const Dune::Fem::TimeProviderBase& tpb)
  : tp {tpb}
{}
void RHS_data::evaluate(const Domain& d, Range& r) const{
  static_assert(Domain::dimension == 3, "Bad domain dimension.");
  static_assert(Range::dimension == 1, "Bad range dimension.");
  const double x = d[0];
  const double y = d[1];
  const double z = d[2];
  const double t = tp.time();
  r =
#include "u_rhs.txt"
    ;
}
RHS_data::Range RHS_data::operator()(const Domain& d) const{
  Range r {0};
  evaluate(d,r);
  return r;
}

void matrixFree_assembly(const Dune::Fem::TimeProviderBase& tp, const Geometry& g,
			 const Quadrature& q, FE_function::LocalFunctionType& f_loc){
  static const RHS_data rhs_fun {tp};
  for(std::size_t pt = 0; pt < q.nop(); ++pt){
    const auto& x = q.point(pt);
    RHS_data::Range fx {rhs_fun(g.global(x))};
    fx *= q.weight(pt) * g.integrationElement(x);
    f_loc.axpy(q[pt], fx);
  }
}
