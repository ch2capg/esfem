/*! \file brusselator_algo_impl.h
    \brief Implementation details for brusselator_algo.cpp

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Idea
     --------------------------------------------------

     Providing many helper classes for the Brusselator_scheme.  
     Since the we solve two scalar equation I have provided several
     struct to avoid word repetition.  

    \author Christian Power
    \date 23. April 2016
    \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef BRUSSELATOR_ALGO_IMPL_H
#define BRUSSELATOR_ALGO_IMPL_H 

#include <string>
#include "esfem.h"
#include "brusselator_algo.h"

namespace Esfem{
  //! For the EOC experiments
  struct Rhs{
    //! Right-hand side for \f$ u \f$
    SecOrd_op::Rhs u;
    //! Right-hand side for \f$ w \f$
    SecOrd_op::Rhs w;
    //! Get time provider
    /*! \todo Right-hand side is currently so coded, that the parameter 
      are not used.  Change this!
    */
    explicit Rhs(const Grid::Grid_and_time&);
    // Rhs(const Io::Parameter&, const Grid::Grid_and_time&);
  };
  
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

  // struct Err_cal{    
  //   Esfem::Io::L2H1_calculator u;
  //   /*!< \brief Calculates the error norm for u. */
  //   Esfem::Io::L2H1_calculator w;
  //   /*!< \brief Calculates the error norm for w. */
  //   Err_cal(const Grid::Grid_and_time&,
  // 	    const Grid::Scal_FEfun_set& u_set,
  // 	    const Grid::Scal_FEfun_set& w_set);
  // };
  // /*!< \brief Class that calculates errors in the \f$L^2\f$-
  //   and \f$H^1\f$-norm.
  // */

  class PreLoop_helper{
  public:
    //! Get access to 
    explicit PreLoop_helper(Brusselator_scheme&);
    /*!< \brief Modifies private data members of `Brusselator_scheme` */
    //! Initial values via an analytic expression
    /*! \sa Rhs */
    void analytic_initialValues();
    //! Initial values via a random distribution
    /*! \sa Rhs */
    void random_initialValues();
    //! Save the initial surface into a temporary file.
    void save_surface();
    //! First line in the error file
    void headLine_in_errFile();
    //! Paraview output 
    void plot_paraview();
    void prepare_rhs();
    /*!< \brief Applies the mass matrix on the
                old surface on `u` and `w`.
	 \warning Overwrites the value of rhs.
    */
  private:
    //! Contains numerical solution
    Brusselator_scheme& bs; 
    //! Reference to `bs.exact`
    const Brusselator_scheme::Init_data& init_data;
    // Err_cal err_cal; /*!< \brief Contains two output files. */
    Io::Paraview paraview;
    /*!< \brief Has reference to `bs.fef.u.fun` and
                `bs.fef.w.fun`
     */
    Scalar_solver solver; /*!< \brief Brusselator solver */
  };
  class PrePattern_helper{
  public:
    explicit PrePattern_helper(Brusselator_scheme&);
    /*!< \brief Modifies private data members of `Brusselator_scheme` */
    // void finalize_rhs();
    void rhs();
    /*!< \brief Adds to member `rhs_les` from
                `Brusselator_scheme::fef.u` and
		`Brusselator_scheme::fef.w`
      
      The new value of 
      `rhs_les` from `Brusselator_scheme::fef.u` respectively
      `Brusselator_scheme::fef.w` will be
      \f{gather*}{
        (M\nodalValue{u})^n + \tau \gamma a M^{n+1} \nodalValue{1}, \\
	(M\nodalValue{w})^n + \tau \gamma b M^{n+1} \nodalValue{1}.
      \f}
      Note that the surface is not changing at this stage.
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
    void plot_paraview();
    /*!< \brief Prints out `fun` of `Brusselator_scheme::fef.u` and
		`Brusselator_scheme::fef.w`.
     */
  private:
    Brusselator_scheme::Io& io;
    /*!< \brief Reference to `Brusselator_scheme::io` */
    Grid::Scal_FEfun_set& u;
    /*!< \brief Reference to `Brusselator_scheme::fef.u` */
    Grid::Scal_FEfun_set& w;
    /*!< \brief Reference to `Brusselator_scheme::fef.w` */
    const Dune::Fem::TimeProviderBase& tp;
    /*!< From `Brusselator_scheme::fix_grid` */
    // Err_cal err_cal; /*!< \brief Contains two output files. */
    Io::Paraview paraview;
    /*!< \brief Has reference to member `fun` from `u` and `w` */
    Scalar_solver solver; /*!< \brief Brusselator solver */
  };
  /*!< \brief Implementation details for
    Brusselator_scheme::prePattern_loop().
  */
  class RhsAndSolve_helper{
  public:
    //! Get private members 
    explicit RhsAndSolve_helper(Brusselator_scheme&);
    //! Assemble \f$ (M\nodalValue{u})^n \f$ and \f$ (M\nodalValue{w})^n \f$
    void scalar_massMatrix();
    /*! \brief \f$ 
          (M_3^n + \alpha A_3^n) \nodalValue{X}^n
	  + \tau \delta M_3^n(\nodalValue{u}^n, 
	  \nodalValue{\surfaceNormal}) \f$
     */
    void brusselator_rhs();
    //! Adds \f$ \tau G^n\f$
    void addScaled_surfaceLoadVector();
    /*! \brief \f$ 
      \parentheses[\big]{M_3^n + (\alpha + \varepsilon\tau) A_3^n} \nodalValue{X}^{n+1}
      = \nodalValue{Y}
      \f$
     */
    void solve_surface_and_save();
  private:
    //! Solver for `X`
    using Vector_solver = Esfem::SecOrd_op::Solution_driven;

    //! \copybrief PreLoop_helper::bs
    Brusselator_scheme& bs;
    //! Reference to `bs.fef`
    Brusselator_scheme::Fef& fef;
    //! Temporally grid, hence `const` 
    const Grid::Grid_and_time grid;
    //! `fef.u.fun` and `fef.u.rhs_les` on `grid`
    Grid::Scal_tiny_FEfun_set u;
    //! `fef.w.fun` and `fef.w.rhs_les` on `grid` 
    Grid::Scal_tiny_FEfun_set w;
    //! `bs_fef.surface.fun` and `bs_fef.surface.rhs_les` on `grid` 
    Grid::Vec_tiny_FEfun_set X;
    //! Mass matrix for `u` and `w`
    Scalar_solver ss;
    //! Solver for `X`
    Vector_solver vs;
    //! Right-hand side for the surface PDE
    SecOrd_op::Vec_rhs v_rhs;
  };
  /*!< \brief Implementation details for 
    Brusselator_scheme::rhs_and_solve_SPDE().
   */
  class Pattern_helper{
  public:
    explicit Pattern_helper(Brusselator_scheme&);
    /*!< \copydoc PrePattern_helper() */
    void finalize_scalarPDE_rhs();
    /*!< \copydoc PrePattern_helper::finalize_rhs() */
    void solve_scalarPDE();
    /*!< \copydoc PrePattern_helper::solve_pde() */
    //! Interpolate exact solutions scalar and vector valued PDE.
    /*! Roughly I interpolate the ESFEM solution on the exact grid and then
      assign nodal values to the current grid. 
     */
    void update_exactSolutions();
    // \f$L^2\f$- and \f$H^1\f$-errors on the numerical surface
    // Not used anymore // void errors_on_numSurface();
    void plot_paraview();
    /*!< \copydoc PrePattern_helper::plot_paraview() */
  private:
    Brusselator_scheme& bs; /*!< \copydoc PreLoop_helper::bs */
    const Grid::Grid_and_time grid;
    /*!< \copydoc RhsAndSolve_helper::grid brief */    
    Grid::Scal_FEfun_set u; /*!< \brief `fef.u` on `grid` */
    Grid::Scal_FEfun_set w; /*!< \brief `fef.w` on `grid` */
    //! Norms on the calculated Grid
    Io::L2H1_calculator norm; 
    Io::Paraview paraview;
    /*!< \brief Has reference to member `u.fun` and `w.fun`. */
    Scalar_solver solver; /*!< \brief Solver for `u` and `w` */
    //! Provides assembly and add_scaled methods.
    Rhs load_vector;
  };
  /*!< \brief Implementation details for 
    Brusselator_scheme::pattern_loop().
  */
  
  // ----------------------------------------------------------------------
  // helper functions

  void interpolate(const SecOrd_op::Init_data&, Grid::FEfun_set<Grid::Scal_FEfun>&);
  /*!< \brief Gives initial values for the members `fun`, `app` and `exact` */
  void head_line(Esfem::Io::Error_stream&);
  /*!< \brief The head line is `timestep L2err H1err` with tab alignment. */
  // void write_error_line(Io::Error_stream& file, const Dune::Fem::TimeProviderBase&,
  // 			const Io::L2H1_calculator&);
  /*!< \brief Prints time step, L2 and H1 error to `file` with proper
              tab alignment.
   */

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
