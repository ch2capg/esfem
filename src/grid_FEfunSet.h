/*! \file grid_FEfunSet.h
    \author Christian Power
    \date 22. March 2016

    \brief Basic struct with input output 
           member functions for a generic numerical experiment

    Revision history
    --------------------------------------------------

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) March 2016

     Idea
     --------------------------------------------------

     This header provides model classes and operator classes to solve 
     a tumor growth model proposed by Elliott and Styles via the ESFEM.  

         Created by Christian Power on 17.03.2016
         Copyright (c) 2016 Christian Power.  All rights reserved.
*/

#ifndef GRID_FEFFUNSET_H
#define GRID_FEFFUNSET_H

#include "esfem_fwd.h"
#include "io_dgf.h"
#include "grid.h"

namespace Esfem{
  namespace Grid{
    template<typename FEfun>
    struct FEfun_set{
      FEfun fun; /*!< Numerical solution */
      FEfun app; /*!< BDF approximation to the numerical solution */
      FEfun exact; /*!< Reference solution */
      FEfun rhs_les; /*!< Right-hand side for the solver */
      FEfun_set(const std::string& name, const Grid::Grid_and_time&);
      /*!< \brief Standard constructor
	\param name Member get concated names like name + "_app" etc.
      */
      FEfun_set(const FEfun_set&, const Grid::Grid_and_time&);
      /*!< \brief Pseudo copy constructor
	
	`Grid::Grid_and_time` is needed to get the correct finite element space.
      */

      void write(const Esfem::Io::Dgf::Handler&,
		 const std::string& dir = "/tmp/");
      /*! \brief Save nodal values in a dgf file. */
      void read(const Esfem::Io::Dgf::Handler&,
		const std::string& dir = "/tmp/");
      /*! \brief Read nodal values from a dgf file. */
    };
    /*!< \brief Minimal set of finite element functions to perform
      higher order BDF-ESFEM with input output helper functions
    */

    inline std::string compose_dgfName(const std::string& fun_name, const std::string& dir = "./"){
      constexpr auto suffix= ".dgf";
      return dir + fun_name + suffix;
    }
    /*!< \brief Returns dune grid format filename
      \param fun_name Expects the result of `Grid::Scal_FEfun::name` or `Grid::Vec_FEfun::name`
      \param dir Directory with trailing backslash
    */  

    // ----------------------------------------------------------------------
    // Template implementation

    template<typename FEfun>
    FEfun_set<FEfun>::
    FEfun_set(const std::string& name,
	      const Esfem::Grid::Grid_and_time& gt)
      : fun {name, gt}, app {name + "_app", gt},
      nexact {name + "_exact", gt}, rhs_les {name + "_rhs_les", gt}
    {}
    template<typename FEfun>
    FEfun_set<FEfun>::
    FEfun_set(const FEfun_set& other, const Esfem::Grid::Grid_and_time& gt)
      : fun {other.fun, gt}, app {other.app, gt}, 
      exact {other.exact, gt}, rhs_les {other.rhs_les, gt}
    {}

    template<typename FEfun>
    void FEfun_set<FEfun>::
    write(const Esfem::Io::Dgf::Handler& h,
	  const std::string& dir){
      h.write(compose_dgfName(fun.name(), dir), fun);
      h.write(compose_dgfName(app.name(), dir), app);
      h.write(compose_dgfName(exact.name(), dir), exact);
      h.write(compose_dgfName(rhs_les.name(), dir), rhs_les);
    }
    template<typename FEfun>
    void FEfun_set<FEfun>::
    read(const Esfem::Io::Dgf::Handler& h, const std::string& dir){
      h.read(compose_dgfName(fun.name(), dir), fun);
      h.read(compose_dgfName(app.name(), dir), app);
      h.read(compose_dgfName(exact.name(), dir), exact);
      h.read(compose_dgfName(rhs_les.name(), dir), rhs_les);
    }
  }	// namespace Grid
}	// namespace Esfem

#endif // GRID_FEFFUNSET_H
