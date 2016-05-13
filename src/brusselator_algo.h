/*! \file brusselator_algo.h
    \brief Numerical experiment for the solution driven paper

     Revision history
     --------------------------------------------------

          Revised by Christian Power April 2016
          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Idea
     --------------------------------------------------

     This header provides model classes and operator classes to solve 
     a tumor growth model proposed by Elliott and Styles via the ESFEM.
     Enter path to writable directory in the macro variable FEF_PATH.


     Partial differential equation
     ==================================================

     Parameter
     --------------------------------------------------

     - \f$a,b, D_{c} \in \R\f$ for the scalar equation with  
       \f$\gamma \sim \vol(\surface)\f$.
     - \f$\varepsilon, \delta, \alpha \in \R\f$ for the surface equation, where 
       \f$\varepsilon\f$ and \f$\alpha\f$ are small regularization parameter.

     Smooth problem
     --------------------------------------------------

     Find \f$u\colon \surface \to \R\f$ (growth-promoting), 
     \f$w\colon \surface \to \R\f$ (growth-inhibiting) and
     \f$X\colon \surface_{0} \times [0,T] \to \R^{m+1}\f$ such that
     \f{gather*}{
       \matd u + u \diver(v) - \laplaceBeltrami u = f_1(u,w) + f_{(1)}, \\
       \matd w + w \diver(v) - D_c \laplaceBeltrami w = f_2(u,w) + f_{(2)},
     \f}
     where for \f$f_1,\, f_2\f$ we use the Brusselator model
     \f{equation*}{
       f_1(u,w) = \gamma (a - u + u^2 w) \quad \land \quad f_2(u,w)
       = \gamma (b - u^2 w),
     \f}
     and \f$ f_{(1)}\f$ and \f$ f_{(2)}\f$ are some forcing terms,
     and for the surface 
     \f{align*}{
       v - \alpha \laplaceBeltrami v = {} &  
       \parentheses[\big]{\varepsilon (-\meanCurvature) + \delta u} \surfaceNormal
       = \varepsilon \laplaceBeltrami X + \delta u \surfaceNormal + g, \\
       \dell_{t} X = {} & v(X),
     \f}
     where \f$ g\f$ is a forcing term.

     Finite element discretization
     --------------------------------------------------

     Find \f$\nodalValue{u}\colon I \to \R^{N}\f$ (growth-promoting nodal values), 
     \f$\nodalValue{w}\colon I \to \R^{N}\f$ (growth-inhibiting 
     nodal values) and 
     \f$\nodalValue{X}\colon I \to \R^{3N}\f$ (surface nodal values) such that
     \f{gather*}{
       \parentheses[\big]{M(X) + \alpha A(X)} \dell_t X = 
       \varepsilon A(X)X + \delta M(u,\nodalValue{\surfaceNormal})
       + G, \\
       \dell_t \parentheses[\big]{M(X) \nodalValue{u} } + A(X) \nodalValue{u}
       = \gamma \parentheses[\big]{a M(X) \nodalValue{1} 
       - M(X)\nodalValue{u}
       + M(X; \nodalValue{u},\nodalValue{w}) \nodalValue{u}}
       + F_{(1)} \\
       \dell_t \parentheses[\big]{M(X) \nodalValue{w}} 
       + D_c A(X) \nodalValue{w}
       = \gamma \parentheses[\big]{b M(X) \nodalValue{1} 
       - M(X; \nodalValue{u},\nodalValue{u}) \nodalValue{w}}
       + F_{(2)}.
     \f}
     We note that instead of the \f$ L^2\f$-projection we use the interpolation
     of \f$ g\f$, \f$ f_{(1)}\f$ respectively \f$ f_{(2)}\f$. 

     Full discretization (Elliott+Styles discretization)
     --------------------------------------------------

     We perform three steps.

     1. Given \f$\nodalValue{X}^n,\, \nodalValue{u}^n,\, \nodalValue{w}^n\f$
        solve for \f$\nodalValue{X}^{n+1}\f$
        \f{equation*}{
	  \parentheses[\big]{M_3^n + (\alpha + \varepsilon\tau) A_3^n}
	  \nodalValue{X}^{n+1} 
	  =  (M_3^n + \alpha A_3^n) \nodalValue{X}^n
	  + \tau \delta M_3^n(\nodalValue{u}^n, 
	  \nodalValue{\surfaceNormal}) + \tau G^n,
	\f}
       where \f$\nodalValue{\surfaceNormal}^n\f$ is elementwise normal.  
     2. Given \f$\nodalValue{X}^{n+1},\, \nodalValue{u}^n,\, \nodalValue{w}^n\f$
        solve for \f$\nodalValue{u}^{n+1}\f$
        \f{equation*}{
	  (1+ \tau \gamma) (M\nodalValue{u})^{n+1} + \tau (A \nodalValue{u})^{n+1}
	  - \tau \gamma M^{n+1}(\nodalValue{u}^n, \nodalValue{w}^n)
	  \nodalValue{u}^{n+1}
	  = (M\nodalValue{u})^n + \tau \gamma a M^{n+1} \nodalValue{1}
	  + \tau F^{n+1}_{(1)},
	\f} 
       where \f$M(a,b)\f$ a \f$4\f$ tensor is, namely 
       \f$M_{ijkl} = \int \chi_i \chi_j \chi_k \chi_l\f$ 
       and \f$\nodalValue{1}\f$ means the 
       constant \f$1\f$ finite element function.
     3. Given \f$\nodalValue{X}^{n+1},\, \nodalValue{u}^{n+1},\,
        \nodalValue{w}^n\f$
        solve for \f$\nodalValue{w}^{n+1}\f$
	\f{equation*}{
	  (M\nodalValue{w})^{n+1} + \tau D_c (A \nodalValue{w})^{n+1} 
	  + \tau \gamma M^{n+1}(\nodalValue{u}^{n+1}, \nodalValue{u}^{n+1})
	  \nodalValue{w}^{n+1}
	  = (M \nodalValue{w})^n + \tau \gamma b M^{n+1} \nodalValue{1}
	  + \tau F^{n+1}_{(2)}.
	\f}

    \author Christian Power
    \date 28. April 2016
    \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef BRUSSELATOR_ALGO_H
#define BRUSSELATOR_ALGO_H 

#include <string>
#include "esfem.h"

#ifndef FEF_PATH
#error Give full path to folder in FEF_PATH
#endif 

namespace Esfem{
  //! ESFEM algorithm.  Only this should be invoked by main.
  /*! \param argc `argc` from `main`
    \param[in] argv `argv` from `main`

    As exact solution for the scalar valued surface equation I chose
    \f{equation*}{
     u(x,y,z,t) = x y e^{-6t} \quad \text{and}\quad
     w(x,y,z,t) = y z e^{-6t}.
    \f}
    For the surface evolution I chose
    \f{equation*}{
     \Phi(x,t) := r(t) x, \quad
     r(t) := \frac{r_{end} r_0}{r_{end} e^{-kt} + r_0 (1-e^{-kt})},
    \f}
    where \f$ r(t)\f$ is the logistic growth function, which satisfies the
    following ODE
    \f{equation*}{
     \dot{r} = k \Bigl( 1 - \frac{r}{r_{end}}\Bigr) r,\quad r(0) = r_0.
    \f}
    From this it follows easily that the velocity is given by
    \f{equation*}{
     v(x,t) = k \Bigl(1 - \frac{r}{r_{end}}\Bigr) x.
    \f}
   */
  void brusselator_algo(int argc, char** argv);

  //! Implementation of the Elliott and Styles full discretization of the tumor problem
  class Brusselator_scheme{
  public:
    Brusselator_scheme(int argc, char** argv,
		       const std::string& parameter_fname);
    /*!< \brief The constructor that also performs the
                first part before the loop enters
      \param argc `argc` from `main`
      \param[in] argv `argv` from `main`
      \param parameter_fname Preferable absolute path to parameter file.
      \warning Absolute path differs on different operating systems.
     */

    /*! \name Prescribed ESFEM */
    //@{
    //! Standard Dziuk Elliott evolving surface finite element method 
    void standard_esfem();
    //@}

    /*! \name Loop action */
    //@{
    void prePattern_loop();
    /*!< \brief In some sense calculates the inital data for
                the solution driven problem.  Starts from \f$t_0\f$.
     */
    void intermediate_action();
    /*!< \brief To be used between prePattern_loop()
                and pattern_loop().

      The inital data has been created.  It will be saved
      and the right hand side for the surface PDE will be created. 
      \warning Do not use next_timeStep() afterwards.
     */
    void pattern_loop();
    /*!< \brief To be used in the second for-loop.
      If prePattern_loop() is off, then \f$t_0\f$ is here.

      At this stage the tumor is growing.
      \warning Do not use next_timeStep() afterwards.
     */
    void final_action();
    /*!< \brief To be used after the second for-loop to save some data. */
    //@}
  private:    
    struct Io{
      //! Functor for saving the current grid. 
      SecOrd_op::Identity identity {};
      
      const Esfem::Io::Dgf::Handler dgf_handler;
      /*!< \brief Converts finite element function into dgf file. */
      
      Esfem::Io::Error_stream u;
      /*!< \brief File to record errors of u. */
      Esfem::Io::Error_stream w;
      /*!< \brief File to record errors of w. */
      //! File to plot errors in the surface
      Esfem::Io::Error_stream surface;
      //! File to plot errors in the velocity
      Esfem::Io::Error_stream velocity;
      
      //! Get file name for error streams.
      Io(const Esfem::Io::Parameter&);
    };
    /*!< \brief Shortens Brusselator_scheme().

      Members are used for input and output of
      the nodal values of the finite element functions
      in `fef` or to calculate \f$L^2\f$- or \f$H^1\f$-norms of the errors.
    */
    //! Initial data respectively exact solution for numerical experiments 
    struct Init_data{
      //! Initial data respectively exact solution functor for \f$u\f$ 
      SecOrd_op::Init_data u;
      //! Initial data respectively exact solution functor for \f$w\f$
      SecOrd_op::Init_data w;
      //! Interpolation functor for the exact velocity
      SecOrd_op::Exact_velocity v;

      //! Provides analytically given initial data.
      explicit Init_data(const Grid::Grid_and_time&);
      //! Provides uniform distributed random inital values. 
      explicit Init_data(const Grid::Grid_and_time&, const Esfem::Io::Parameter&);
    };
    struct Fef{
      Grid::Scal_FEfun_set u;
      /*!< \brief Container for growth promoting numerical solution */
      Grid::Scal_FEfun_set w;
      /*!< \brief Container for growth inhibiting numerical solution */
      Grid::Vec_FEfun_set surface;
      /*!< \brief Container for the numerical solution of the surface
	\warning Only this container is allowed to write into and read from
	  a dgf file.
       */
      //! Analytically given exact velocity.
      Grid::Vec_FEfun_set velocity;
      const std::string tmpFile_path {FEF_PATH};
      /*!< \brief Directory which I have read and write access.
	\warning `FEF_PATH` is a macro variable which has be set by
	  the makefile. 
       */
      Fef(const Grid::Grid_and_time&);
    };
    /*!< \brief Shortens Brusselator_scheme().

      Collects all finite element functions and serves as a
      backup container. 
      \todo Add error checking for `tmpFile_path`.
    */

    Esfem::Io::Parameter data;
    /*!< \brief Contains parameter from `tumor_parameter.txt`. */
    //! Error streams and identity functor for the surface 
    Io io;
    //! Analytically given grid with absolute time provider.
    Esfem::Grid::Grid_and_time fix_grid;
    //! Norms on the analytically given grid
    Esfem::Io::L2H1_calculator norm;
    //! Finite element functions 
    Fef fef;
    //! Exact solution for \f$u\f$, \f$w\f$ and \f$v\f$
    Init_data exact;
    
    // ------------------------------------------------------------
    // Helper member functions
    
    /*! \name Helper classes for the for-loops */
    //@{
    friend class PreLoop_helper;
    /*!< \brief Used in pre_loop_action(). */
    friend class PrePattern_helper;
    /*!< \brief Used in prePattern_loop(). */
    friend class Pattern_helper;
    /*!< \brief Used in pattern_loop(). */
    friend class RhsAndSolve_helper;
    /*!< \brief Used in rhs_and_solve_SPDE(). */
    //@}

    //! Assign a new value to `fef.surface.exact` and `fef.surface.app`
    void update_surface();
    //! Assign a new value to `fef.velocity.exact` and `fef.velocity.fun`
    /*! \pre The `fef.surface.fun` represents the new surface and
      `fef.surface.app` represents the old surface.
    */
    void update_velocity();
    //! Assign a new value to `fef.u.exact` and `fef.w.exact`
    void update_scalar_solution();

    //! Plot error of `fef` on the interpolated surface
    /*! \pre The exact solution has been updated.
      \sa update_exact_surface(), update_exact_velocity(),
       update_scalar_solution() 
    */
    void error_on_intSurface();
    
    //! Constructor helper
    /*! \pre Should only be invoked by Brusselator_scheme().*/
    void pre_loop_action();
    //! Helper for pattern_loop().
    /*! Generates first part of the right-hand side for the scalar
      SPDE.  Then solves the vector SPDE and prints out a dgf file.
     */
    void rhs_and_solve_SPDE();
    
    /*! \name Flow control */
    //@{ 
    void next_timeStep(); 
    /*!< \brief Increments the next time step in `fix_grid`. */
    long prePattern_timeSteps() const; 
    /*!< \brief Maximum number of time steps for prePattern_loop(). */
    long pattern_timeSteps() const; 
    /*!< \brief Maximum number of time steps for pattern_loop(). */
    //! Absolute number of time steps 
    long time_steps() const;
    //@}
  };

  // ----------------------------------------------------------------------
  // Inline implementation

  inline void Brusselator_scheme::next_timeStep(){
    fix_grid.next_timeStep(data.global_timeStep());
  }
  inline long Brusselator_scheme::prePattern_timeSteps() const{
    return data.prePattern_timeSteps();
  }
  inline long Brusselator_scheme::pattern_timeSteps() const{
    return data.pattern_timeSteps();
  }
  inline long Brusselator_scheme::time_steps() const{
    return data.max_timeSteps();
  }
}

#endif // BRUSSELATOR_ALGO_H
