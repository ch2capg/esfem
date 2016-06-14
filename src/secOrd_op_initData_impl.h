/*! \file secOrd_op_initData_impl.h
    \brief Implementation details for `secOrd_op_initData.cpp`

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Idea
     --------------------------------------------------

     Actual implementation of Initial data.

     \author Christian Power 
     \date 26. April 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef SECORD_OP_INITDATA_IMPL_H
#define SECORD_OP_INITDATA_IMPL_H 

#include <random>
#include <functional>
#include <string>
#include <config.h>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include "secOrd_op_initData.h"
#include "grid.h"

namespace Esfem{
  namespace Impl{
    //! Initial data is given via an analytic formula
    class Explicit_initial_data
      : public Dune::Fem::Function<Esfem::Grid::Grid_and_time::Function_space,
				   Explicit_initial_data>
    {
    public:
      //! \f$ f\colon\R^3\to \R\f$
      using Fun_space = Esfem::Grid::Grid_and_time::Function_space;
      //! \f$ \R^3\f$
      using Domain = Fun_space::DomainType;
      //! \f$ \R\f$
      using Range = Fun_space::RangeType;

      static_assert(Domain::dimension == 3, "Bad Domain dimension");
      static_assert(Range::dimension == 1, "Bad Range dimension");
    
      //! Time provider and determine if we want initial data for \f$ u\f$ or \f$ w\f$
      /*! \post Input Grid_and_time must outlive this object!
	\todo We have hard coded the parameter into the function.  Change
	the constructor so that it reads the parameter.
      */
      Explicit_initial_data(const Esfem::Grid::Grid_and_time&,
			    const Esfem::Growth);
      //! No copy constructor
      Explicit_initial_data(const Explicit_initial_data&) = delete;
      //! No copy assignment 
      Explicit_initial_data& operator=(const Explicit_initial_data&) = delete;

      //! Required for the dune Lagrange interpolation class
      void evaluate(const Domain&, Range&) const;
    private:
      //! Current time 
      const Dune::Fem::TimeProviderBase& tp;
      //! Initial data function for \f$ u\f$ or \f$ w\f$
      std::function<void(const Domain&,Range&)> fun_impl;
    };

    //! Random perturbation around an equilibrium point
    class Random_initial_data
      : public Dune::Fem::Function<Esfem::Grid::Grid_and_time::Function_space,
				   Random_initial_data>
    {
    public:
      //! \copybrief Explicit_initial_data::Fun_space
      using Fun_space = Esfem::Grid::Grid_and_time::Function_space;
      //! \copybrief Explicit_initial_data::Domain
      using Domain = typename Fun_space::DomainType;
      //! \copybrief Explicit_initial_data::Range
      using Range = typename Fun_space::RangeType;

      static_assert(Domain::dimension == 3, "Bad Domain dimension");
      static_assert(Range::dimension == 1, "Bad Range dimension");
    
      //! \copybrief Explicit_initial_data()
      Random_initial_data(const Esfem::Io::Parameter&,
			  const Esfem::Growth);
    
      //! \copybrief Explicit_initial_data::evaluate()
      void evaluate(const Domain&, Range&) const;
    private:
      //! Random distribution
      using Random_dist = std::uniform_real_distribution<>;
      //! Seed 
      using Random_engine = std::default_random_engine;
    
      //! Random function for evaluate()
      std::function<double()> random_fun;

      //! Constructs `random_fun`
      Random_initial_data(const double hom_value, const double pertubation);
    };

    //! First eigen function of the sphere
    /*! \f$xye^{-6t}\f$*/
    struct sphere_1EF
      : Dune::Fem::Function<Esfem::Grid::Grid_and_time::Function_space,
			    sphere_1EF>,
	SecOrd_op::sIdata{
      //! \f$f\colon \R^3\to\R^1\f$
      using fspace = Esfem::Grid::Grid_and_time::Function_space;
      //! \f$ \R^3\f$
      using domain = typename fspace::DomainType;
      //! \f$\R^1\f$
      using range = typename fspace::RangeType;
      static_assert(domain::dimension == 3, "Bad Domain dimension");
      static_assert(range::dimension == 1, "Bad Range dimension");
      //! Get time provider
      sphere_1EF(const Grid::Grid_and_time&);
      sphere_1EF* clone() override;
      void interpolate(Grid::Scal_FEfun&) const override; 
      //! \copybrief Explicit_initial_data::evaluate()
      void evaluate(const domain&, range&) const;
    private:
      //! Current time
      const Dune::Fem::TimeProviderBase& tp;      
    };
    //! Second eigen function of the sphere
    /*! \f$yze^{-6t}\f$*/
    struct sphere_2EF
      : Dune::Fem::Function<Esfem::Grid::Grid_and_time::Function_space,
			    sphere_2EF>,
	SecOrd_op::sIdata{
      //! \f$f\colon \R^3\to\R^1\f$
      using fspace = Esfem::Grid::Grid_and_time::Function_space;
      //! \f$ \R^3\f$
      using domain = typename fspace::DomainType;
      //! \f$\R^1\f$
      using range = typename fspace::RangeType;
      static_assert(domain::dimension == 3, "Bad Domain dimension");
      static_assert(range::dimension == 1, "Bad Range dimension");
      //! Get time provider
      sphere_2EF(const Grid::Grid_and_time&);
      sphere_2EF* clone() override;
      void interpolate(Grid::Scal_FEfun&) const override; 
      //! \copybrief Explicit_initial_data::evaluate()
      void evaluate(const domain&, range&) const;
    private:
      //! Current time
      const Dune::Fem::TimeProviderBase& tp;      
    };
    //! Thrid eigen function of the sphere
    /*! \f$xze^{-6t}\f$*/
    struct sphere_3EF
      : Dune::Fem::Function<Esfem::Grid::Grid_and_time::Function_space,
			    sphere_3EF>,
	SecOrd_op::sIdata{
      //! \f$f\colon \R^3\to\R^1\f$
      using fspace = Esfem::Grid::Grid_and_time::Function_space;
      //! \f$ \R^3\f$
      using domain = typename fspace::DomainType;
      //! \f$\R^1\f$
      using range = typename fspace::RangeType;
      static_assert(domain::dimension == 3, "Bad Domain dimension");
      static_assert(range::dimension == 1, "Bad Range dimension");
      //! Get time provider
      sphere_3EF(const Grid::Grid_and_time&);
      sphere_3EF* clone() override;
      void interpolate(Grid::Scal_FEfun&) const override; 
      //! \copybrief Explicit_initial_data::evaluate()
      void evaluate(const domain&, range&) const;
    private:
      //! Current time
      const Dune::Fem::TimeProviderBase& tp;      
    };


    //! Implementation of SecOrd_op::vIdata::ssef()
    struct sphere_eigenFun
      : Dune::Fem::Function<Esfem::Grid::Grid_and_time::Vec_Function_space,
			    sphere_eigenFun>,
        SecOrd_op::vIdata{
      //! \f$f\colon \R^3\to\R^3\f$
      using Fun_space = Esfem::Grid::Grid_and_time::Vec_Function_space;
      //! \f$\R^3\f$
      using Domain = typename Fun_space::DomainType;
      //! \f$\R^3\f$
      using Range = typename Fun_space::RangeType;

      static_assert(Domain::dimension == 3, "Bad Domain dimension");
      static_assert(Range::dimension == 3, "Bad Range dimension");
      //! Get time provider
      sphere_eigenFun(const Grid::Grid_and_time&);
      sphere_eigenFun* clone() override;
      void interpolate(Grid::Vec_FEfun&) const override; 
      //! \copybrief Explicit_initial_data::evaluate()
      void evaluate(const Domain&, Range&) const;
    private:
      //! Current time
      const Dune::Fem::TimeProviderBase& tp;
    };

    //! Implementation of SecOrd_op::vIdata::ssef()
    /*! \pre The initial sphere is the unit sphere \f$S^2\f$.  */
    struct sphere_mcf_sol
      : Dune::Fem::Function<Esfem::Grid::Grid_and_time::Vec_Function_space,
			    sphere_mcf_sol>,
        SecOrd_op::vIdata{
      //! \f$f\colon \R^3\to\R^3\f$
      using Fun_space = Esfem::Grid::Grid_and_time::Vec_Function_space;
      //! \f$\R^3\f$
      using Domain = typename Fun_space::DomainType;
      //! \f$\R^3\f$
      using Range = typename Fun_space::RangeType;

      static_assert(Domain::dimension == 3, "Bad Domain dimension");
      static_assert(Range::dimension == 3, "Bad Range dimension");
      //! Get time provider
      sphere_mcf_sol(const Grid::Grid_and_time&);
      sphere_mcf_sol* clone() override;
      void interpolate(Grid::Vec_FEfun&) const override; 
      //! \copybrief Explicit_initial_data::evaluate()
      void evaluate(const Domain&, Range&) const;
    private:
      //! Current time
      const Dune::Fem::TimeProviderBase& tp;
    };    

    //! Sphere growing via na logistic growth function
    /*! \pre I assume that I am on the exact surface, which is a sphere.
      \sa Esfem::Brusselator_scheme::eoc_sls() */
    struct sls_iData
      : Dune::Fem::Function<Esfem::Grid::Grid_and_time::Vec_Function_space,
			    sls_iData>,
        SecOrd_op::vIdata{
      //! \f$f\colon \R^3\to\R^3\f$
      using Fun_space = Esfem::Grid::Grid_and_time::Vec_Function_space;
      //! \f$\R^3\f$
      using Domain = typename Fun_space::DomainType;
      //! \f$\R^3\f$
      using Range = typename Fun_space::RangeType;

      static_assert(Domain::dimension == 3, "Bad Domain dimension");
      static_assert(Range::dimension == 3, "Bad Range dimension");
      //! Get time provider and logistic growth parameter
      sls_iData(const Grid::Grid_and_time&);
      sls_iData* clone() override{ return new sls_iData {*this}; }
      void interpolate(Grid::Vec_FEfun&) const override; 
      //! \copybrief Explicit_initial_data::evaluate()
      void evaluate(const Domain&, Range&) const;
    private:
      //! Current time
      const Dune::Fem::TimeProviderBase& tp;
      //! Initial radius
      double rA;
      //! End radius or carrying capcity
      double rE;
      //! Growth rate
      double k;
    };
    
    //! The actual implementation of the velocity
    class Analytic_velocity
      : public Dune::Fem::Function<Esfem::Grid::Grid_and_time::Vec_Function_space,
				   Analytic_velocity>
    {
    public:
      //! \f$f\colon \R^3\to\R^3\f$
      using Fun_space = Esfem::Grid::Grid_and_time::Vec_Function_space;
      //! \f$\R^3\f$
      using Domain = typename Fun_space::DomainType;
      //! \f$\R^3\f$
      using Range = typename Fun_space::RangeType;

      static_assert(Domain::dimension == 3, "Bad Domain dimension");
      static_assert(Range::dimension == 3, "Bad Range dimension");

      //! Get time provider
      Analytic_velocity(const Esfem::Grid::Grid_and_time&);
      //! No copy constructor
      Analytic_velocity(const Analytic_velocity&) = delete;
      //! No copy assignment 
      Analytic_velocity& operator=(const Analytic_velocity&) = delete;

      //! \copybrief Explicit_initial_data::evaluate()
      void evaluate(const Domain&, Range&) const;
    private:
      //! Current time
      const Dune::Fem::TimeProviderBase& tp;
    };

    // ----------------------------------------------------------------------
    // helper functions 

    // ------------------------------------------------------------
    // Random_init_data helper
    double hom_value(const Esfem::Io::Parameter&, const Esfem::Growth);
    double pertubation(const Esfem::Io::Parameter&, const Esfem::Growth);
    std::string print_configuration(const Esfem::Io::Parameter&, const Esfem::Growth);

    // ------------------------------------------------------------
    // Init_data::Data helper
    std::string dof_filename(const Io::Parameter&, const Growth);
  } // namespace Impl

  //! %Data members of Init_data
  /*! \note If we want to get rid of the unique pointer, we may use inheritance. */
  struct SecOrd_op::Init_data::Data{
    //! Name for a log file 
    const std::string dof_io_filename {};
    //! Analytically given function.
    std::unique_ptr<Impl::Explicit_initial_data> eid_ptr;
    //! Function with random distributed coefficients.
    std::unique_ptr<Impl::Random_initial_data> rid_ptr;
    
    //! `eid_ptr` constructor
    /*! \post Grid_and_time must outlive this object. */
    Data(const Grid::Grid_and_time&, const Growth);
    //! `rid_ptr` constructor
    Data(const Io::Parameter&, const Growth);
  };

  //! %Data members of Exact_velocity
  struct SecOrd_op::Exact_velocity::Data{
    //! Functor of the analytic function
    Impl::Analytic_velocity v_fun;
    
    //! Get time provider
    Data(const Grid::Grid_and_time&);
  };
} // namespace Esfem

// ----------------------------------------------------------------------
// Inline implementation

inline void Esfem::Impl::Explicit_initial_data::
evaluate(const Domain& d, Range& r) const{
  fun_impl(d,r);
}    
inline void Esfem::Impl::Random_initial_data::
evaluate(const Domain&, Range& q) const{
  q = random_fun(); 
}  

#endif // SECORD_OP_INITDATA_IMPL_H
