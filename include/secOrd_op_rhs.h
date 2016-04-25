/*! \file secOrd_op_rhs.h
    \brief Load vector for ESFEM discretization of a scalar or vector valued PDE 

     Revision history
     --------------------------------------------------

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     Idea
     --------------------------------------------------

     Wrapper class for the dune implementation.

     \author Christian Power
     \date 23. April 2016
     \copyright Copyright (c) 2016 Christian Power. All rights reserved.
 */

#ifndef SECORD_OP_RHS_H
#define SECORD_OP_RHS_H 

#include <memory>
#include "esfem_fwd.h"

//! The evolving surface finite element method
namespace Esfem{
  //! Parabolic and elliptic second order operators
  namespace SecOrd_op{
    //! Load vector for scalar valued PDE
    /*! 
      Useful formulas
      --------------------------------------------------

      \f{gather*}{
      \diver_{\surface}(X) = \diver_{\R^{n+1}}(\bar{X}) 
      - \surfaceNormal^T (D_{\R^{n+1}}\bar{X}) \surfaceNormal, \\
      \laplaceBeltrami f = \Delta_{\R^{n+1}}\bar{f} - \surfaceNormal^T
      (D_{\R^{n+1}}^2\bar{f}) \surfaceNormal - H \surfaceNormal^T D_{\R^{n+1}}f.
      \f}
     */
    class Rhs{
    public:
      //! Get time, finite element space and an indicator 
      /*! \post Grid_and_time must outlive this object. */
      explicit Rhs(const Grid::Grid_and_time&, const Growth);
      //! Needed for the pointer to data implementation technique.
      ~Rhs();

      //! Assemble the load vector and add with some scale to the input.      
      void assemble_and_addScaled_to(Grid::Scal_FEfun&);
    private:
      struct Data;
      //! Data pointer
      std::unique_ptr<Data> d_ptr;
    };

    //! Load vector for vector valued PDE
    class Vec_rhs{
    public:
      //! Get time and finite element space
      /*! \post Grid_and_time must outlive this object. */
      explicit Vec_rhs(const Grid::Grid_and_time&);
      //! Needed for the pointer to data implementation technique.
      ~Vec_rhs();

      //! Assemble the load vector and add with some scale to the input.      
      void assemble_and_addScaled_to(Grid::Vec_FEfun&);
    private:
      struct Data;
      //! Data pointer
      std::unique_ptr<Data> d_ptr;
    };
  }
}

#endif // SECORD_OP_RHS_H
