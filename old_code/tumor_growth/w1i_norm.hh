#ifndef W1I_NORM_HH
#define W1I_NORM_HH

#include <dune/fem/quadrature/cachingquadrature.hh>

template<typename In>
double linfty_norm(In first, In last) noexcept {
  In max = first;
  while(first != last){
    if(std::abs(*max) < std::abs(*first)) 
      max = first;
    ++first;
  }
  if(max == last) return 0.;
  return std::abs(*max);
}

template<typename DiscreteFunction>
double linfty_norm(DiscreteFunction df){
  return linfty_norm(df.dbegin(), df.dend());
}

template<typename T>
class TD;

template<class DiscreteFunction>
double w1infty_norm(const DiscreteFunction& df){
  // calculates the W^{1,\infty} seminorm of 'df'
  using DiscreteFunction_space
    = typename DiscreteFunction::DiscreteFunctionSpaceType;
  using LocalFunction
    = typename DiscreteFunction::LocalFunctionType;

  using Iterator
    = typename DiscreteFunction_space::IteratorType;
  using Entity
    = typename Iterator::Entity;
  using Geometry
    = typename Entity::Geometry;

  using GridPart
    = typename DiscreteFunction_space::GridPartType;
  using Quadrature
    = Dune::Fem::CachingQuadrature<GridPartType, 0>;
  // ----------------------------------------------------------------------

  auto max_value = 0.;

  const DiscreteFunction_space& df_space = df.space();
  for(const auto& entity : df_space){
    const Geometry& geo = entity.geometry();
    const LocalFunction df_local = df.localFunction(entity);
    // Debug: The following code works only for df_local.order() == 1
    // if(df_local.order() != 1)
    //   throw std::runtime_error {"Polynomial order is to not 1."};
    Quadrature quadrature(entity, 1);
    const auto x = quadrature[0];
    
    typename LocalFunction::JacobianRangeType jac_df;
    df_local.jacobian(x,jac_df);

    auto w1infty_of_jac_df 
      = std::sqrt( jac_df[0][0] * jac_df[0][0] + jac_df[0][1] * jac_df[0][1]
		   + jac_df[0][2] * jac_df[0][2]);

    if(w1infty_of_jac_df > max_value) max_value = w1infty_of_jac_df;
  }
  return max_value;
}

#endif // #ifndef W1I_NORM_HH
