/*! \file brusselator_algo_impl.h

    \brief Implementation details for brusselator_algo.cpp

     Revision history
     --------------------------------------------------

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Idea
     --------------------------------------------------

     To do.

         Created by Christian Power on 17.03.2016
         Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef BRUSSELATOR_ALGO_IMPL_H
#define BRUSSELATOR_ALGO_IMPL_H 

#include <string>
#include "esfem.h"
#include "brusselator_algo.h"

namespace Esfem{
  struct Rhs{
    Esfem::SecOrd_op::Rhs_u u;
    Esfem::SecOrd_op::Rhs_w w;
    explicit Rhs(const Esfem::Io::Parameter&, const Esfem::Grid::Grid_and_time&);
  };
  struct Init_data{
    // Esfem::SecOrd_op::Init_data_u u;
    // Esfem::SecOrd_op::Init_data_w w;
    Esfem::SecOrd_op::Init_data u;
    Esfem::SecOrd_op::Init_data w;
    explicit Init_data(const Esfem::Grid::Grid_and_time&) = delete;
    explicit Init_data(const Esfem::Io::Parameter&);
  };
  
  struct Solver{
    Esfem::SecOrd_op::Brusselator u;
    Esfem::SecOrd_op::Brusselator w;
    explicit Solver(const Esfem::Io::Parameter&, const Esfem::Grid::Grid_and_time&,
		    const FEfun_set<Esfem::Grid::Scal_FEfun>& u_set,
		    const FEfun_set<Esfem::Grid::Scal_FEfun>& w_set);
  };
  struct Err_cal{
    Esfem::Io::L2H1_calculator u;
    Esfem::Io::L2H1_calculator w;
    explicit Err_cal(const Esfem::Grid::Grid_and_time&,
		     const FEfun_set<Esfem::Grid::Scal_FEfun>& u_set,
		     const FEfun_set<Esfem::Grid::Scal_FEfun>& w_set);
  };
  struct Err_stream{
    Esfem::Io::Error_stream u;
    Esfem::Io::Error_stream w;
    explicit Err_stream(const Esfem::Io::Parameter&);
  };

  class PrePattern_helper{
  };
  class Pattern_helper{
  };
  class Helper_surface{
  public:
    Helper_surface(const Scal_FEfun_set& u_input, Vec_FEfun_set& surface,
		   const Io::Parameter&, const Grid::Grid_and_time&);
    Helper_surface(const Helper_surface&) = delete;
    Helper_surface& operator=(const Helper_surface&) = delete;

    void solve_for_surface();
  private:
    /*! \name Reference to container */
    //@{
    const Scal_FEfun_set& u;
    Vec_FEfun_set& surface;
    //@}
    /*! \name Local loop variable */
    //@{
    const Grid::Grid_and_time gt;
    Vec_FEfun_set surface_loc;
    SecOrd_op::Solution_driven solver;
    //@}
  };
  /*!< \brief `Esfem::Brusselator_scheme::pattern_action`
               pattern action implementation details
  */
  class Helper_uw{
  public:
    Helper_uw(Scal_FEfun_set& u_input, Scal_FEfun_set& w_input,
	      Vec_FEfun_set& surface, const Io::Parameter&,
	      const Grid::Grid_and_time&);
    Helper_uw(const Helper_uw&) = delete;
    Helper_uw& operator=(const Helper_uw&) = delete;

    void solve_pde();
    void write_error_line();
    void paraview_plot();
  private:
    /*! \name Reference to container */
    //@{
    Scal_FEfun_set& u; 
    Scal_FEfun_set& w;
    //@}
    /*! \name Local loop variable */
    //@{
    const Grid::Grid_and_time gt; 
    Scal_FEfun_set u_loc;
    Scal_FEfun_set w_loc;
    const Err_cal errCal;
    Solver solver;
    Io::Paraview paraview;
    //@}
  };
  /*!< \brief `Esfem::Brusselator_scheme::pattern_action`
               pattern action implementation details
  */
  
  // ----------------------------------------------------------------------
  // helper functions

  void first_interpolate(const SecOrd_op::Identity&,
			 const Init_data&,
			 FEfun_set<Esfem::Grid::Scal_FEfun>& u,
			 FEfun_set<Esfem::Grid::Scal_FEfun>& w,
			 FEfun_set<Esfem::Grid::Vec_FEfun>& surface);
  void update_exactSolution(const Init_data&,
			    FEfun_set<Esfem::Grid::Scal_FEfun>& u,
			    FEfun_set<Esfem::Grid::Scal_FEfun>& w);
  void massMatrix_rhsLes(const Solver&,
			 FEfun_set<Esfem::Grid::Scal_FEfun>& u,
			 FEfun_set<Esfem::Grid::Scal_FEfun>& w);
  void massMatrixConstOne_rhsLes(const Solver&,
				 FEfun_set<Esfem::Grid::Scal_FEfun>& u,
				 FEfun_set<Esfem::Grid::Scal_FEfun>& w);
  void assemble_and_addScaled_rhsLes(const Rhs&,
				     FEfun_set<Esfem::Grid::Scal_FEfun>& u,
				     FEfun_set<Esfem::Grid::Scal_FEfun>& w);
  void solve_pde(const Solver&,
		 FEfun_set<Esfem::Grid::Scal_FEfun>& u,
		 FEfun_set<Esfem::Grid::Scal_FEfun>& w);


  // ------------------------------------------------------------
  // io helper functions

  void generate_header_line(Err_stream&);
  void generate_header_line(Esfem::Io::Error_stream&);
  void write_error_line(Err_stream&, const Esfem::Grid::Grid_and_time&, const Err_cal&);
  void write_error_line(Esfem::Io::Error_stream&,
			const Esfem::Grid::Grid_and_time&,
			const Esfem::Io::L2H1_calculator&);
  std::string compose_dgfName(const std::string& fun_name, const std::string& dir = "./");
  void clog_uw(const FEfun_set<Esfem::Grid::Scal_FEfun>& u, 
	       const FEfun_set<Esfem::Grid::Scal_FEfun>& w);
} // namespace Esfem

#endif // BRUSSELATOR_ALGO_IMPL_H
