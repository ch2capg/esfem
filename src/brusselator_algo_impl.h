/*! \file brusselator_algo_impl.h
    \author Christian Power
    \date 29. March 2016
    
    \brief Implementation details for brusselator_algo.cpp

     Revision history
     --------------------------------------------------

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Idea
     --------------------------------------------------

     To do.

         Created by Christian Power on 29.03.2016
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
  /*!< \brief Used to test the other components. */
  struct Init_data{
    // Esfem::SecOrd_op::Init_data_u u;
    // Esfem::SecOrd_op::Init_data_w w;
    SecOrd_op::Init_data u; /*!< Initial data for u */
    SecOrd_op::Init_data w; /*!< Initial data for w */
    explicit Init_data(const Grid::Grid_and_time&) = delete;
    /*!< \todo Revise this constructor.  */
    explicit Init_data(const Io::Parameter&);
    /*!< \brief Provides uniform distributed random inital values. */
  };
  /*!< \brief Initial data for the numerical experiment. */
  struct Scalar_solver{
    SecOrd_op::Brusselator u; /*!< PDE with solver for u */
    SecOrd_op::Brusselator w; /*!< PDE with solver for w */
    Scalar_solver(const Io::Parameter&,
		  const Grid::Grid_and_time&,
		  const Grid::Scal_FEfun_set& u_set,
		  const Grid::Scal_FEfun_set& w_set);
    /*!< \brief Use this constructor if you need to solve a system. */
    Scalar_solver(const Io::Parameter&,
		  const Grid::Grid_and_time&,
		  const Grid::Scal_tiny_FEfun_set& u_set,
		  const Grid::Scal_tiny_FEfun_set& w_set);
    /*!< \brief Use this constructor if you only need a mass matrix.
      Never use this constructor if you need to solve a system.
    */
  };
  struct Err_cal{    
    Esfem::Io::L2H1_calculator u;
    /*!< \brief Calculates the error norm for u. */
    Esfem::Io::L2H1_calculator w;
    /*!< \brief Calculates the error norm for w. */
    Err_cal(const Grid::Grid_and_time&,
	    const Grid::Scal_FEfun_set& u_set,
	    const Grid::Scal_FEfun_set& w_set);
  };
  /*!< \brief Class that calculates errors in the \f$L^2\f$-
    and \f$H^1\f$-norm.
  */
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
    Brusselator_scheme& bs; /*!< \brief Contains numerical solution. */
    const Init_data init_data;
    Err_cal err_cal; /*!< \brief Contains two output files. */
    Io::Paraview paraview;
    /*!< \brief Has reference to `bs.fef.u.fun` and
                `bs.fef.w.fun`
     */
    Scalar_solver solver; /*!< \brief Brusselator solver */
  };
  class PrePattern_helper{
  public:
    PrePattern_helper(Brusselator_scheme&);
    /*!< \brief Modifies private data members of `Brusselator_scheme` */
    void finalize_rhs();
    /*!< \brief Adds to member `rhs_les` from
                `Brusselator_scheme::fef.u` and
		`Brusselator_scheme::fef.w`
      
      We assume that the `rhs_les` has been prepared with
      PreLoop_helper::prepare_rhs() or
      prepare_rhs().  We add
      the following finite element function to the member
      `rhs_les` of `Brusselator_scheme::fef.u` respectively
      `Brusselator_scheme::fef.w`
      \f{gather*}{
        \tau \gamma a M^{n+1} \nodalValue{1}, \\
	\tau \gamma b M^{n+1} \nodalValue{1}.
      \f}
     */
    void solve_pde();
    /*!< \brief New value for member `fun` from
                `Brusselator_scheme::fef.u` and
		`Brusselator_scheme::fef.w`
	 
      We assume that the `rhs_les` has been prepared with
      finalize_rhs().  The vaule of `fun` will be overwritten
      with the solution of the following system for \f$u\f$ and \f$w\f$:
      \f{gather*}{
        (1+ \tau \gamma) (M\nodalValue{u})^{n+1} + \tau (A \nodalValue{u})^{n+1}
        - \tau \gamma M^{n+1}(\nodalValue{u}^n, \nodalValue{w}^n)
        \nodalValue{u}^{n+1} = \nodalValue{y}_{1}, \\
        (M\nodalValue{w})^{n+1} + \tau D_c (A \nodalValue{w})^{n+1} 
        + \tau \gamma M^{n+1}(\nodalValue{u}^{n+1}, \nodalValue{u}^{n+1})
        \nodalValue{w}^{n+1} = \nodalValue{y}_{2},
      \f}
      where \f$ \nodalValue{y}_{1,2}\f$ the right-hand side of the system is.
     */
    void prepare_rhs();
    /*!< \brief New value for member `rhs_les` from
                `Brusselator_scheme::fef.u` and
		`Brusselator_scheme::fef.w`
      
      `rhs_les` will have the following value for \f$u\f$ respectively
      \f$w\f$:
      \f{equation*}{
       (M\nodalValue{u})^n \quad \lor\quad (M \nodalValue{w})^n.
      \f}
     */
    void plot_errors_in_errFile();
    /*!< \brief Prints out time step, \f$L^2\f$-error,
                and \f$H^1\f$-error of (`fun` - `exact`)
		from `Brusselator_scheme::fef.u` and
		`Brusselator_scheme::fef.w`.
     */
    void plot_paraview();
    /*!< \brief Prints out `fun` of `Brusselator_scheme::fef.u` and
		`Brusselator_scheme::fef.w`.
     */
  private:
    Brusselator_scheme& bs; /*!< \copydoc #PreLoop_helper::bs. */
    Err_cal err_cal; /*!< \brief Contains two output files. */
    Io::Paraview paraview;
    /*!< \brief Has reference to member #fun from #bs.fef.u and
                #bs.fef.w
     */
    Scalar_solver solver; /*!< \brief Brusselator solver */
  };
  /*!< \brief Implementation details for
    Brusselator_scheme::prePattern_loop().
  */
  class RhsAndSolve_helper{
  public:
    RhsAndSolve_helper(Brusselator_scheme&);
    void scalar_massMatrix();
    void solve_surface_and_save();
  private:
    using Vector_solver = Esfem::SecOrd_op::Solution_driven;
    /*!< \brief Solver for #X */
    
    Brusselator_scheme& bs; /*!< \copydoc #PreLoop_helper::bs */
    Brusselator_scheme::Fef& fef; /*!< \brief Reference to #bs.fef */
    const Grid::Grid_and_time grid;
    /*!< \brief Temporally grid, hence `const` */
    Grid::Scal_tiny_FEfun_set u;
    /*!< \brief #fef.u.fun and #fef.u.rhs_les on #grid. */
    Grid::Scal_tiny_FEfun_set w;
    /*!< \brief #fef.w.fun and #fef.w.rhs_les on #grid. */
    Grid::Vec_tiny_FEfun_set X;
    /*!< \brief #bs_fef.surface.fun and #bs_fef.surface.rhs_les on #grid. */
    Scalar_solver ss;
    /*!< \brief Mass matrix for #u and #w */
    Vector_solver vs;
    /*!< \brief Solver for #X */
  };
  /*!< \brief Implementation details for 
    Brusselator_scheme::rhs_and_solve_SPDE().
   */
  class Pattern_helper{
  public:
    Pattern_helper(Brusselator_scheme&);
    /*!< \copydoc PrePattern_helper() */
    void finalize_scalarPDE_rhs();
    /*!< \copydoc PrePattern_helper::finalize_rhs() */
    void solve_scalarPDE();
    /*!< \copydoc PrePattern_helper::solve_pde() */
    void plot_errors_in_errFile();
    /*!< \copydoc PrePattern_helper::plot_errors_in_errFile() */
    void plot_paraview();
    /*!< \copydoc PrePattern_helper::plot_paraview() */
  private:
    Brusselator_scheme& bs; /*!< \copydoc #PreLoop_helper::bs */
    const Grid::Grid_and_time grid;
    /*!< \copydoc #RhsAndSolve_helper::grid brief */    
    Grid::Scal_FEfun_set u; /*!< \brief #fef.u on #grid */
    Grid::Scal_FEfun_set w; /*!< \brief #fef.w on #grid */
    Err_cal err_cal; /*!< \copydoc #PrePattern_helper::err_cal */
    Io::Paraview paraview;
    /*!< \brief Has reference to member #u.fun and #w.fun. */
    Scalar_solver solver; /*!< \brief Solver for #u and #w */
  };
  /*!< \brief Implementation details for 
    Brusselator_scheme::pattern_loop().
  */
  // class Helper_surface{
  // public:
  //   Helper_surface(const Grid::Scal_FEfun_set& u_input, Grid::Vec_FEfun_set& surface,
  // 		   const Io::Parameter&, const Grid::Grid_and_time&);
  //   Helper_surface(const Helper_surface&) = delete;
  //   Helper_surface& operator=(const Helper_surface&) = delete;
  // 
  //   void solve_for_surface();
  // private:
  //   /*! \name Reference to container */
  //   //@{
  //   const Grid::Scal_FEfun_set& u;
  //   Grid::Vec_FEfun_set& surface;
  //   //@}
  //   /*! \name Local loop variable */
  //   //@{
  //   const Grid::Grid_and_time gt;
  //   Vec_FEfun_set surface_loc;
  //   SecOrd_op::Solution_driven solver;
  //   //@}
  // };
  // /*!< \brief `Esfem::Brusselator_scheme::pattern_action`
  //              pattern action implementation details
  // */
  // class Helper_uw{
  // public:
  //   Helper_uw(Scal_FEfun_set& u_input, Scal_FEfun_set& w_input,
  // 	      Vec_FEfun_set& surface, const Io::Parameter&,
  // 	      const Grid::Grid_and_time&);
  //   Helper_uw(const Helper_uw&) = delete;
  //   Helper_uw& operator=(const Helper_uw&) = delete;
  // 
  //   void solve_pde();
  //   void write_error_line();
  //   void paraview_plot();
  // private:
  //   /*! \name Reference to container */
  //   //@{
  //   Scal_FEfun_set& u; 
  //   Scal_FEfun_set& w;
  //   //@}
  //   /*! \name Local loop variable */
  //   //@{
  //   const Grid::Grid_and_time gt; 
  //   Scal_FEfun_set u_loc;
  //   Scal_FEfun_set w_loc;
  //   const Err_cal errCal;
  //   Solver solver;
  //   Io::Paraview paraview;
  //   //@}
  // };
  // /*!< \brief `Esfem::Brusselator_scheme::pattern_action`
  //              pattern action implementation details
  // */
  
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


// // ----------------------------------------------------------------------
// // old code
// 
// void first_interpolate(const SecOrd_op::Identity&,
// 			 const Init_data&,
// 			 FEfun_set<Grid::Scal_FEfun>& u,
// 			 FEfun_set<Grid::Scal_FEfun>& w,
// 			 FEfun_set<Grid::Vec_FEfun>& surface);
// void update_exactSolution(const Init_data&,
// 			    FEfun_set<Grid::Scal_FEfun>& u,
// 			    FEfun_set<Grid::Scal_FEfun>& w);
// void massMatrix_rhsLes(const Solver&,
// 			 FEfun_set<Grid::Scal_FEfun>& u,
// 			 FEfun_set<Grid::Scal_FEfun>& w);
// void massMatrixConstOne_rhsLes(const Solver&,
// 				 FEfun_set<Grid::Scal_FEfun>& u,
// 				 FEfun_set<Grid::Scal_FEfun>& w);
// void assemble_and_addScaled_rhsLes(const Rhs&,
// 				     FEfun_set<Grid::Scal_FEfun>& u,
// 				     FEfun_set<Grid::Scal_FEfun>& w);
// void solve_pde(const Solver&,
// 		 FEfun_set<Grid::Scal_FEfun>& u,
// 		 FEfun_set<Grid::Scal_FEfun>& w);


  // ------------------------------------------------------------
  // io helper functions

  // void generate_header_line(Esfem::Io::Err_stream&);
  // void generate_header_line(Esfem::Io::Error_stream&);
  // void write_error_line(Io::Err_stream&, const Esfem::Grid::Grid_and_time&, const Io::Err_cal&);
  // void write_error_line(Io::Error_stream&,
  // 			const Grid::Grid_and_time&,
  // 			const Io::L2H1_calculator&);
  // void clog_uw(const FEfun_set<Grid::Scal_FEfun>& u, 
  // 	       const FEfun_set<Grid::Scal_FEfun>& w);

  // ----------------------------------------------------------------------
  // Inline implementation

  inline void PreLoop_helper::plot_paraview(){
    paraview.write();
  }
  inline void PrePattern_helper::plot_paraview(){
    paraview.write();
  }
  inline void Pattern_helper::plot_paraview(){
    paraview.write();
  }
  
} // namespace Esfem

#endif // BRUSSELATOR_ALGO_IMPL_H
