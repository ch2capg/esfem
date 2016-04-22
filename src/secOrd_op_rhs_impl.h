/*! \file secOrd_op_rhs_impl.h
    \brief Helper classes for `secOrd_op_rhs.cpp`

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) April 2016

     Idea
     --------------------------------------------------

     Implement `Vec_rhs::Data`

     \author Christian Power 
     \date 16. April 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef SECORD_OP_RHS_IMPL_H
#define SECORD_OP_RHS_IMPL_H 

#include <functional>
#include <config.h>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include "grid.h"
#include "secOrd_op_rhs.h"


namespace Esfem{
  //! Implementation details
  namespace Impl{
    //! Scalar valued right-hand side for the surface PDE
    class Rhs_fun
      :public Dune::Fem::Function
    <Esfem::Grid::Grid_and_time::Function_space, Rhs_fun>
    {
    public:
      //! Template argument
      using Base = Esfem::Grid::Grid_and_time::Function_space;
      //! \f$ \R^3\f$
      using Domain = Base::DomainType;
      //! \f$ \R^1\f$
      using Range = Base::RangeType;

      static_assert(Domain::dimension == 3, "Bad Domain dimension");
      static_assert(Range::dimension == 1, "Bad Range dimension");
      
      //! Get time, time step and indicate right-hand side for `u` or `w`
      /*! \warning If time provider outlives `this`, then you are in problem. */
      Rhs_fun(const Dune::Fem::TimeProviderBase&,
	      const Growth);
      //! No copy constructor
      Rhs_fun(const Rhs_fun&) = delete;
      //! No copy assignment 
      Rhs_fun& operator=(const Rhs_fun&) = delete;

      //! y = f(x)
      void evaluate(const Domain&, Range&) const;
      //! \copybrief evaluate()
      Range operator()(const Domain&) const;
    private:
      //! Current time and time step
      const Dune::Fem::TimeProviderBase& tp;
      //! Right-hand side for \f$ u\f$ or \f$ w\f$
      std::function<void(const Domain&,Range&)> fun_impl;
    };

    //! Vector valued right-hand side for the mean curvature equation
    class Vec_rhs_fun
      :public Dune::Fem::Function
    <Esfem::Grid::Grid_and_time::Vec_Function_space, Vec_rhs_fun>
    {
    public:
      //! Template argument
      using Base = Esfem::Grid::Grid_and_time::Vec_Function_space;
      //! \f$ \R^3 \f$
      using Domain = Base::DomainType;
      //! \f$ \R^3 \f$
      using Range = Base::RangeType;

      static_assert(Domain::dimension == 3, "Bad Domain dimension");
      static_assert(Range::dimension == 3, "Bad Range dimension");

      //! Get time and time step
      /*! \copydetails Rhs_fun::Rhs_fun() */
      explicit Vec_rhs_fun(const Dune::Fem::TimeProviderBase&);
      //! No copy construct
      Vec_rhs_fun(const Vec_rhs_fun&) = delete;
      //! No copy assignment
      Vec_rhs_fun& operator=(const Vec_rhs_fun&) = delete;

      //! \copybrief Rhs_fun::evaluate()
      void evaluate(const Domain&, Range&) const;
      //! \copybrief evaluate()
      Range operator()(const Domain&) const;
    private:
      //! Time and time step
      const Dune::Fem::TimeProviderBase& tp;
    };    

    //! Assemble load vector
    /*! \tparam Rhs Insert Rhs_fun or Vec_rhs_fun
      \tparam Fef Insert 
      Scal_FEfun::Dune_FEfun or Vec_FEfun::Dune_FEfun
     */
    template<typename Rhs, typename Fef>
    void assemble_RHS(const Rhs& rhs, Fef& fef);
  } // namespace Impl

  //! %Data details of Rhs
  struct SecOrd_op::Rhs::Data{
    //! Time and time step
    const Dune::Fem::TimeProviderBase& tp;
    //! Scalar valued right-hand side function 
    Impl::Rhs_fun rhs;
    //! Finite element load vector
    Esfem::Grid::Scal_FEfun::Dune_FEfun load_vector;

    //! Get space-time manifold and finite element space
    Data(const Grid::Grid_and_time&, const Growth);
  };

  //! %Data details of Vec_rhs
  struct SecOrd_op::Vec_rhs::Data{
    //! \copybrief Rhs::Data::tp
    const Dune::Fem::TimeProviderBase& tp;
    //! Vector valued right-hand side function
    Impl::Vec_rhs_fun rhs;
    //! \copybrief Rhs::Data::load_vector
    Esfem::Grid::Vec_FEfun::Dune_FEfun load_vector;
    
    //! \copybrief Rhs::Data::Data()
    Data(const Grid::Grid_and_time&);
  };
} // namespace Esfem

// ----------------------------------------------------------------------
// Template implementation

template<typename Rhs, typename Fef>
void Esfem::Impl::assemble_RHS(const Rhs& rhs, Fef& fef){
  using Range = typename Rhs::Range;
  using Grid_part = typename Fef::GridPartType;
  using Quadrature = Dune::Fem::CachingQuadrature<Grid_part, 0>;
  
  fef.clear();
  const auto& df_space = fef.space();
  for(const auto& entity : df_space){
    const auto& geometry = entity.geometry();
    const Quadrature quad {entity, 2 * df_space.order() + 1};
    auto fef_local = fef.localFunction(entity);
    for(std::size_t pt = 0; pt < quad.nop(); ++pt){
      const auto& x = quad.point(pt);
      Range fx {rhs(geometry.global(x))};
      // rhs.evaluate(geometry.global(x), fx);
      fx *= quad.weight(pt) * geometry.integrationElement(x);
      fef_local.axpy(quad[pt], fx);
    }  
  }
  fef.communicate();    
}

// ----------------------------------------------------------------------
// Inline implementation

inline void Esfem::Impl::Rhs_fun::evaluate(const Domain& d, Range& r) const{
  fun_impl(d,r);
  // const double x = d[0];
  // const double y = d[1];
  // // const double z = d[2];
  // const double t = tp.time();
  // r = std::exp(-6.*t)*x*y;
}
inline Esfem::Impl::Rhs_fun::Range
Esfem::Impl::Rhs_fun::operator()(const Domain& d) const{
  Range r {0};
  evaluate(d,r);
  return r;
}

#endif // SECORD_OP_RHS_IMPL_H
