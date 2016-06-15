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
      //! Dynamic assert
      void dassert(const bool assertion, const std::string& msg);
    };

    //! Vector valued right-hand side for the mean curvature equation
    /*! This is my right-hand function:
      \f[
      \biggl[ k \Bigl( 1 - \frac{r(t)}{r_{end}}\Bigr) \lvert X \rvert
      + 2 \Bigl(\alpha k \Bigl( 1 - \frac{r(t)}{r_{end}}\Bigr) 
      + \varepsilon \Bigr) \frac{1}{\lvert X \rvert} 
      - \delta xy e^{-6t}\biggr]
      \f]
      Currently I am testing.  So, I use an easier right-hand side.
      I am doing classic mean curvature flow on the sphere.  
      Right-hand side must be zero for this example.  Also I expect
      unit sphere as starting value, alpha = 0 and epsilon = 1.

      Old test, which I want to check again later: 
      - \f$ \surface_0 = S^2\f$ 
      - Exact flow: \f$ \Phi(x,t) = e^{-t} x\f$
      - Velocity: \f$ v(x,t) = -x\f$ 
      - Normal: \f$ n = x/ |x|\f$
      - Mean curvature: \f$ H = 2/ |x|\f$
      - \f$\delta = \varepsilon = 0\f$, hence the operator is
      \f[
      v - \alpha \Delta v = g(x,t) = -(1 + \alpha 2/ |x|^2) x.
      \f]
      \todo Finish test and clean up comments.
     */
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
      //! \f$\alpha \Delta v\f$
      double alpha;
      //! \f$ \varepsilon \Delta X\f$
      double epsilon;
      //! Initial size of sphere
      double r_start;
      //! Final size of sphere
      double r_end;
      //! Steepness of logistic growh
      double k;
      //! \f$\delta u\f$
      double delta;
      //! Helper variable, which save computations
      mutable double cache[4];
      //! Conditional update member `cache`
      /*! `cache[0]` holds the last time.  If the current time does not equals `cache[0]`,
	then update the cache.  The content of cache is 
	- `cache[0] = tp.time()`,
	- \f$ k( 1 - \frac{r(t)}{r_{end}})\f$, 
	- \f$2 ( \alpha k( 1 - \frac{r(t)}{r_{end}}) + \varepsilon)\f$, 
	- \f$\delta e^{-6t}\f$,

	where \f$r(t) = \frac{r_{end}r_0}{r_{end} e^{-kt} + r_0 (1 - e^{-kt})}\f$
      */      
      void update_cache() const;
    };    

    //! Right-hand side for surface logistic sphere experiment
    /*! \pre Use evaluate() on the exact surface.  
      I assume that it is a sphere, such that \f$ r(t) = |x|\f$.
      \sa Esfem::Brusselator_scheme::eoc_sls() */
    struct sls_rhs 
      : Dune::Fem::Function
         <Esfem::Grid::Grid_and_time::Vec_Function_space, sls_rhs>,
	SecOrd_op::vRhs{
      //! Dune function 
      using dBase = Esfem::Grid::Grid_and_time::Vec_Function_space;
      //! \f$ \R^3\f$
      using dom = dBase::DomainType;
      //! \f$ \R^3\f$
      using ran = dBase::RangeType;
      //! For my generic algorithm
      using Range = ran;

      //! Get time and finit element space
      /*! \post Grid and time outlive this object. */
      sls_rhs(const Grid::Grid_and_time&);
      sls_rhs* clone() override{ return new sls_rhs {*this}; }
      void addScaled_to(Grid::Vec_FEfun& rhs) override;
      //! Needed for interpolation 
      ran operator()(const dom&) const;
    private:
      //! Manifold dimension
      static constexpr int dim {2};
      //! Time step
      const Dune::Fem::TimeProviderBase& tp;      
      //! Load vector
      Esfem::Grid::Vec_FEfun::Dune_FEfun lvec;
      //! Carrying capacity
      double r_end;
      //! \f$\alpha\f$
      double a;
      //! \f$\varepsilon\f$
      double e;
      //! Growth rate
      double k;
    };

    //! Assemble load vector
    /*! \tparam Rhs Deduce Rhs_fun or Vec_rhs_fun
      \tparam Fef Deduce 
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
}
inline Esfem::Impl::Rhs_fun::Range
Esfem::Impl::Rhs_fun::operator()(const Domain& d) const{
  Range r {0};
  evaluate(d,r);
  return r;
}

#endif // SECORD_OP_RHS_IMPL_H
