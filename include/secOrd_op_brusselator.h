/*! \file secOrd_op_brusselator.h
    \author Christian Power
    \date 23. March 2016
  
    \brief Provides the class `Esfem::SecOrd_op::Brusselator`

     Revision history:

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) Februar 2016

     Solves the scalar surface equation using the brusselator model.
     You have to use the enum class `Esfem::Growth` to specify which
     model you want.

     Created by Christian Power on 23.03.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef SECORD_OP_BRUSSELATOR_H
#define SECORD_OP_BRUSSELATOR_H 

#include <memory>
#include "esfem_fwd.h"

namespace Esfem{
  namespace SecOrd_op{
    class Brusselator{
    public:
      Brusselator(const Io::Parameter&, const Grid::Grid_and_time&,
		  const Growth);n
      /*!< \brief This constructor creats two internal FE functions. */
      Brusselator(const Io::Parameter&, const Grid::Grid_and_time&,
		  const Growth, const Grid::Scal_FEfun&,
		  const Grid::Scal_FEfun&);
      /*!< \brief This constructor keeps references of
	          two external FE functions.
      */
      ~Brusselator();

      /*! \name Numerical interface for the finite element code */
      //@{
      void solve(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const;
      void mass_matrix(const Grid::Scal_FEfun& rhs, Grid::Scal_FEfun& lhs) const;
      void massMatrix_constOne(Grid::Scal_FEfun&) const;
      void add_massMatrixConstOne_to(Grid::Scal_FEfun&) const;
      //@}
      
      /*! \name Assignment member functions
	\warning The following two methods throw exceptions if
	         you have constructed
	         this object with the second constructor.  
       */
      //@{
      void assign_firstArg_quadMassMatrix(const Grid::Scal_FEfun&);
      void assign_secondArg_quadMassMatrix(const Grid::Scal_FEfun&);
      //@}

      void operator()(const Grid::Scal_FEfun&, Grid::Scal_FEfun&) const;
      /*! \brief Included for testing reasons. */
    private:
      struct Data;
      std::unique_ptr<Data> d_ptr;
    };
  }	// namespace Grid
}	// namespace Esfem

#endif // SECORD_OP_BRUSSELATOR_H
