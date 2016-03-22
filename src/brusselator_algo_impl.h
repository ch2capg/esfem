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
    SecOrd_op::Rhs_u u; /*!< Right-hand side for u */
    SecOrd_op::Rhs_w w; /*!< Right-hand side for w */
    explicit Rhs(const Io::Parameter&, const Grid::Grid_and_time&);
  };
  struct Init_data{
    // Esfem::SecOrd_op::Init_data_u u;
    // Esfem::SecOrd_op::Init_data_w w;
    SecOrd_op::Init_data u; /*!< Initial data for u */
    SecOrd_op::Init_data w; /*!< Initial data for w */
    explicit Init_data(const Esfem::Grid::Grid_and_time&) = delete; // to be revised 
    explicit Init_data(const Esfem::Io::Parameter&);
    /*!< \brief Provides uniform distributed random inital values. */
  };
  struct Scalar_solver{
    using Scal_FEfun_set = Grid::FEfun_set<Grid::Scal_FEfun>;
    using Vec_FEfun_set = Grid::FEfun_set<Grid::Vec_FEfun>;

    SecOrd_op::Brusselator u; /*!< PDE with solver for u */
    SecOrd_op::Brusselator w; /*!< PDE with solver for w */
    explicit Solver(const Io::Parameter&,
		    const Grid::Grid_and_time&,
		    const Scal_FEfun_set& u_set,
		    const Scal_FEfun_set& w_set);
  };
  struct Err_cal{    
    using Scal_FEfun_set = Grid::FEfun_set<Grid::Scal_FEfun>;
    using Vec_FEfun_set = Grid::FEfun_set<Grid::Vec_FEfun>;

    Esfem::Io::L2H1_calculator u; /*!< Calculates the error norm for u. */
    Esfem::Io::L2H1_calculator w; /*!< Calculates the error norm for w. */
    Err_cal(const Grid::Grid_and_time&,
	    const Scal_FEfun_set& u_set,
	    const Vec_FEfun_set& w_set);
  };
  class PreLoop_helper{
  public:
    PreLoop_helper(Brusselator_scheme&);
    /*!< \brief Modifies private data members of `Brusselator_scheme` */
    void first_interpolate();
    void headLine_in_errFile();
    void plot_errors_in_errFile();
    void plot_paraview();
    void prepare_rhs();
    /*!< \brief Applies the mass matrix on the
                old surface on `u` and `w`.
	 \warning Overwrites the value of rhs.
    */
  private:
    Brusselator_scheme& bs; /*!< Contains numerical solution. */
    Io& io; /*!< Reference to `bs.io` */
    Fef& fef; /*!< Reference to `bs.fef` */
    const Init_data init_data;
    Err_cal err_cal; /*!< Contains two output files. */
    Io::Paraview paraview;
    Scalar_solver solver; /*!< Brusselator solver */
  };
  class PrePattern_helper{
  public:
    PrePattern_helper(Brusselator_scheme&);
    /*!< \brief Modifies private data members of `Brusselator_scheme` */
    void finalize_rhs();
    /*!< \brief Adds to the rhs, which has been prepared with
                `PreLoop_helper::finalize_rhs` or
                `PrePattern_helper::prepare_rhs`,
		tau * gamma * (a|b) * M * 1.  
     */
    void solve_pde();
    void prepare_rhs();
    void plot_errors_in_errFile();
    void plot_paraview();
  private:
    Brusselator_scheme& bs; /*!< Contains numerical solution. */
    Io& io; /*!< Reference to bs.io */
    Fef& fef; /*!< Reference to bs.fef */
    Err_cal err_cal; /*!< Contains two output files. */
    Io::Paraview paraview; 
    Scalar_solver solver; /*!< Brusselator solver */
  };
  /*!< \brief Implementation details for the pre
              pattern loop in `Brusselator_scheme`
  */
  class Pattern_helper{
  public:
    Pattern_helper(Brusselator_scheme&);
    /*!< \brief Modifies private data members of `Brusselator_scheme` */
  private:
    Grid::Grid_and_time fix_grid;    
  };
  /*!< \brief Implementation details for the 
              pattern loop in `Brusselator_scheme`
  */
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

  void interpolate(const SecOrd_op::Init_data&, Grid::FEfun_set<Grid::Scal_FEfun>&);
  /*!< \brief Gives initial values for the members `fun`, `app` and `exact` */
  void head_line(Esfem::Io::Error_stream&);
  /*!< \brief The head line is `timestep L2err H1err` with tab alignment. */
  void write_error_line(Io::Error_stream& file, const Dune::Fem::TimeProviderBase& ,
		      const Io::L2H1_calculator&);
  /*!< \brief Prints time step, L2 and H1 error to `file` with proper
              tab alignment.
   */


  // ----------------------------------------------------------------------
  // old code
  
  void first_interpolate(const SecOrd_op::Identity&,
			 const Init_data&,
			 FEfun_set<Grid::Scal_FEfun>& u,
			 FEfun_set<Grid::Scal_FEfun>& w,
			 FEfun_set<Grid::Vec_FEfun>& surface);
  void update_exactSolution(const Init_data&,
			    FEfun_set<Grid::Scal_FEfun>& u,
			    FEfun_set<Grid::Scal_FEfun>& w);
  void massMatrix_rhsLes(const Solver&,
			 FEfun_set<Grid::Scal_FEfun>& u,
			 FEfun_set<Grid::Scal_FEfun>& w);
  void massMatrixConstOne_rhsLes(const Solver&,
				 FEfun_set<Grid::Scal_FEfun>& u,
				 FEfun_set<Grid::Scal_FEfun>& w);
  void assemble_and_addScaled_rhsLes(const Rhs&,
				     FEfun_set<Grid::Scal_FEfun>& u,
				     FEfun_set<Grid::Scal_FEfun>& w);
  void solve_pde(const Solver&,
		 FEfun_set<Grid::Scal_FEfun>& u,
		 FEfun_set<Grid::Scal_FEfun>& w);


  // ------------------------------------------------------------
  // io helper functions

  void generate_header_line(Err_stream&);

  void write_error_line(Err_stream&, const Esfem::Grid::Grid_and_time&, const Err_cal&);
  void write_error_line(Io::Error_stream&,
			const Grid::Grid_and_time&,
			const Io::L2H1_calculator&);
  void clog_uw(const FEfun_set<Grid::Scal_FEfun>& u, 
	       const FEfun_set<Grid::Scal_FEfun>& w);
} // namespace Esfem

#endif // BRUSSELATOR_ALGO_IMPL_H
