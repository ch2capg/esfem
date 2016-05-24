/*! \file secOrd_op_identity.h
    \brief Providing an identity function for vector valued finite element functions

     Revision history
     --------------------------------------------------

          Revised by Christian Power Mai 2016
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Idea
     --------------------------------------------------

     This class provides a way to save the current grid in a vector valued finite 
     element function.  With a different function you can save the nodal values into
     a file with the dune grid format.

     \author Christian Power
     \date 23. Mai 2016
     \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
*/

#ifndef SECORD_OP_IDENTITY_H
#define SECORD_OP_IDENTITY_H

#include <memory>
#include "esfem_fwd.h"

namespace Esfem{
  namespace SecOrd_op{
    //! Interface for the identity function.
    class Identity{
    public:
      //! Set up 
      Identity();
      //! Required by pointer to implementation technique.
      ~Identity();

      void interpolate(Grid::Vec_FEfun&) const;
      /*!< \brief The nodal values of the finite element functions 
	          will be the nodal values itself.
       */
    private:
      struct Data;
      //! Pointer to implementation
      std::unique_ptr<Data> d_ptr;
    };
  }
}

#endif // SECORD_OP_IDENTITY_H
