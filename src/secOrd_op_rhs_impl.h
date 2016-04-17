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

#include "secOrd_op_rhs.h"
#include "grid.h"

namespace Esfem{
  namespace Impl{
    class RHS_data
      :public Dune::Fem::Function
    <Esfem::Grid::Grid_and_time::Function_space, RHS_data>
    {
    public:
      using Base = Esfem::Grid::Grid_and_time::Function_space;
      using Domain = Base::DomainType;
      using Range = Base::RangeType;
  
      explicit RHS_data(const Dune::Fem::TimeProviderBase&);
      RHS_data(const RHS_data&) = delete;
      RHS_data& operator=(const RHS_data&) = delete;

      void evaluate(const Domain&, Range&) const;
      Range operator()(const Domain&) const;
    private:
      const Dune::Fem::TimeProviderBase& tp;
    };

    class Vec_rhs_data
      :public Dune::Fem::Function
    <Esfem::Grid::Grid_and_time::Function_space, Vec_rhs_data>
    {
    public:
      using Base = Esfem::Grid::Grid_and_time::Function_space;
      using Domain = Base::DomainType;
      using Range = Base::RangeType;
  
      explicit Vec_rhs_data(const Dune::Fem::TimeProviderBase&);
      Vec_rhs_data(const Vec_rhs_data&) = delete;
      Vec_rhs_data& operator=(const Vec_rhs_data&) = delete;

      void evaluate(const Domain&, Range&) const;
      Range operator()(const Domain&) const;
    private:
      const Dune::Fem::TimeProviderBase& tp;
    };    
    
    void assemble_RHS(const RHS_data&, Esfem::Grid::Scal_FEfun::Dune_FEfun&);
    void massMatrix_for_entity
    (const Esfem::Grid::Scal_FEfun::Dune_FEfun::
     DiscreteFunctionSpaceType::IteratorType::Entity::Geometry&,
     const Dune::Fem::CachingQuadrature<Grid_part, 0>&,
     const RHS_data&,
     Esfem::Grid::Scal_FEfun::Dune_FEfun::LocalFunctionType&);
  } // namespace Impl

  struct SecOrd_op::Rhs::Data{
    RHS_data rhs;
    const Dune::Fem::TimeProviderBase& tp;
    const Grid::Grid_and_time::FE_space& fe_space;
    Data(const Grid::Grid_and_time& gt);
  };

  struct SecOrd_op::Vec_rhs::Data{
    Vec_rhs_fun rhs;
    const Dune::Fem::TimeProviderBase& tp;
    const Grid::Grid_and_time::FE_space& fe_space;
    Data(const Grid::Grid_and_time& gt)
  };  
} // namespace Esfem

#endif // SECORD_OP_RHS_IMPL_H
