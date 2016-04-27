/*! \file io_l2h1Calculator.h
    \brief \f$L^2\f$- and \f$H^1\f$-norm

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Originally written by Christian Power
               (power22c@gmail.com) Januar 2016

     Idea
     --------------------------------------------------

     Build a wrapper class for the dune implementation. 

     \author Christian Power
     \date 27. April 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef IO_L2H1CALCULATOR_H
#define IO_L2H1CALCULATOR_H 

#include <memory>
#include "esfem_fwd.h"

namespace Esfem{
  //! Input output routines
  namespace Io{
    //! Wrapper class for dune \f$L^2\f$- and \f$H^1\f$-norms
    class L2H1_calculator{
    public:
      // explicit L2H1_calculator(const Grid::Grid_and_time&,
      // 			       const Grid::Scal_FEfun& exact_solution,
      // 			       const Grid::Scal_FEfun& numerical_solution);
      //! Get grid.
      /*! \post Grid_and_time outlives this object. */
      explicit L2H1_calculator(const Grid::Grid_and_time&);
      //! Required by pointer to implementation technique.
      ~L2H1_calculator();

      //! \f$L^2\f$-distance of two finite element functions
      double l2_err(const Grid::Scal_FEfun&, const Grid::Scal_FEfun&) const;
      //! \copybrief l2_err()
      double l2_err(const Grid::Vec_FEfun&, const Grid::Vec_FEfun&) const;
      //! \f$H^1\f$-distance of two finite element functions
      double h1_err(const Grid::Scal_FEfun&, const Grid::Scal_FEfun&) const;
      //! \copybrief h1_err()
      double h1_err(const Grid::Vec_FEfun&, const Grid::Vec_FEfun&) const;
    private:
      struct Data;
      //! Pointer to data member
      std::unique_ptr<Data> d_ptr;
    };
  }
}

#endif // IO_L2H1CALCULATOR_H
