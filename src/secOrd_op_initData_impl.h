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
     \date 22. April 2016
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

      //! Required for the interpolation class
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

/*! Log:
 */
