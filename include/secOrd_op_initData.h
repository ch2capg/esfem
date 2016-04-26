/*! \file secOrd_op_initData.h
    \brief Providing initial data for the experiment

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     Idea
     --------------------------------------------------

     Wrapper class for the dune function class.

     \author Christian Power 
     \date 22. April 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef SECORD_OP_INITDATA_H
#define SECORD_OP_INITDATA_H 

#include <memory>
#include "esfem_fwd.h"

namespace Esfem{
  namespace SecOrd_op{
    //! Initial data for scalar ESFEM experiments
    class Init_data{
    public:
      //! Constructor for an analytic given initial function
      /*! \post Grid_and_time must outlive this object. */
      Init_data(const Grid::Grid_and_time&, const Growth);
      //! Constructor for an random initial function
      Init_data(const Io::Parameter&, const Growth);
      //! Required for the pointer to implementation technique.
      ~Init_data();

      //! Lagrange interpolation
      void interpolate(Grid::Scal_FEfun&) const;
    private:
      struct Data;
      //! Pointer to data members
      std::unique_ptr<Data> d_ptr;
    };

    //! Analytically given velocity
    class Exact_velocity{
    public:
      //! Fetch time provider
      /*! \post Gird_and_time must outlive this object.
	\todo Change constructor so that it also reads parameter.
       */
      Exact_velocity(const Grid::Grid_and_time&);
      //! Required for the pointer to implementation technique
      ~Exact_velocity();

      //! Lagrange interpolation
      void interpolate(Grid::Vec_FEfun&) const;
    private:
      struct Data;
      //! Pointer to data members
      std::unique_ptr<Data> d_ptr;
    };
  }
}

#endif // SECORD_OP_INITDATA_H
