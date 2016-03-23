/*! \file brusselator_algo.h
    \author Christian Power
    \date 23. March 2016

    \brief Numerical experiment for the solution driven paper

     Revision history
     --------------------------------------------------

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Idea
     --------------------------------------------------

     This header provides model classes and operator classes to solve 
     a tumor growth model proposed by Elliott and Styles via the ESFEM.  


     Partial differential equation
     ==================================================

     Parameter
     --------------------------------------------------

     - \f$a,b, D_{c} \in \R\f$ for the scalar equation with  
       \f$\gamma \sim \vol(\surface)\f$.
     - \f$\varepsilon, \delta, a \in \R\f$ for the surface equation, where 
       \f$\varepsilon\f$ and \f$\alpha\f$ are small regularization parameter.

     Smooth problem
     --------------------------------------------------

     Find \f$u\colon \surface \to \R\f$ (growth-promoting), 
     \f$w\colon \surface \to \R\f$ (growth-inhibiting) and
     \f$X\colon \surface_{0} \times [0,T] \to \R^{m+1}\f$ such that
     \f{gather*}{
       \matd u + u \diver(v) - \laplaceBeltrami u = f_1(u,w), \\
       \matd w + w \diver(v) - D_c \laplaceBeltrami w = f_2(u,w),
     \f}
     where for \f$f_1,\, f_2\f$ we use the Brusselator model
     \f{equation*}{
       f_1(u,w) = \gamma (a - u + u^2 w) \quad \land \quad f_2(u,w)
       = \gamma (b - u^2 w),
     \f}
     and for the surface 
     \f{align*}{
       v - \alpha \laplaceBeltrami v = {} &  
       \parentheses[\big]{\varepsilon (-\meanCurvature) + \delta u} \surfaceNormal
       = \varepsilon \laplaceBeltrami X + \delta u \surfaceNormal, \\
       \dell_{t} X = {} & v(X).
     \f}

     Finite element discretization
     --------------------------------------------------

     Find \f$\nodalValue{u}\colon I \to \R^{N}\f$ (growth-promoting nodal values), 
     \f$\nodalValue{w}\colon I \to \R^{N}\f$ (growth-inhibiting 
     nodal values) and 
     \f$\nodalValue{X}\colon I \to \R^{3N}\f$ (surface nodal values) such that
     \f{gather*}{
       \parentheses[\big]{M(X) + \alpha A(X)} \dell_t X = 
       \varepsilon A(X) + \delta M(u,\nodalValue{\surfaceNormal}), \\
       \dell_t \parentheses[\big]{M(X) \nodalValue{u} } + A(X) \nodalValue{u}
       = \gamma \parentheses[\big]{b M(X) \nodalValue{1} 
       + M(X; \nodalValue{u},\nodalValue{w}) \nodalValue{u}}, \\
       \dell_t \parentheses[\big]{M(X) \nodalValue{w}} 
       + D_c A(X) \nodalValue{w}
       = \gamma \parentheses[\big]{b M(X) \nodalValue{1} 
       - M(X; \nodalValue{u},\nodalValue{u}) \nodalValue{w}}.
     \f}

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
	  \nodalValue{\surfaceNormal}),
	\f}
       where \f$\nodalValue{\surfaceNormal}^n\f$ is elementwise normal.  
     2. Given \f$\nodalValue{X}^{n+1},\, \nodalValue{u}^n,\, \nodalValue{w}^n\f$
        solve for \f$\nodalValue{u}^{n+1}\f$
        \f{equation*}{
	  (1+ \tau \gamma) (M\nodalValue{u})^{n+1} + \tau (A \nodalValue{u})^{n+1}
	  - \tau \gamma M^{n+1}(\nodalValue{u}^n, \nodalValue{w}^n)
	  \nodalValue{u}^{n+1}
	  = (M\nodalValue{u})^n + \tau \gamma a M^{n+1} \nodalValue{1},
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
	  = (M \nodalValue{w})^n + \tau \gamma b M^{n+1} \nodalValue{1}.
	\f}

         Created by Christian Power on 23.03.2016
         Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef BRUSSELATOR_ALGO_H
#define BRUSSELATOR_ALGO_H 

#include <string>
#include "esfem.h"

namespace Esfem{
  void brusselator_algo(int argc, char** argv);
  /*!< \brief ESFEM algorithm.  Only this should be invoked by main. */
  
  class Brusselator_scheme{
  public:
    using Scal_FEfun_set = Grid::FEfun_set<Grid::Scal_FEfun>;
    /*!< \brief Four functions of type \f$ f\colon \R^3 \to \R \f$ */
    using Vec_FEfun_set = Grid::FEfun_set<Grid::Vec_FEfun>;
    /*!< \brief Four functions of type \f$ f\colon \R^3 \to \R^3 \f$ */
    
    Brusselator_scheme(int argc, char** argv,
		       const std::string& parameter_fname);
    /*!< \brief The constructor that also performs the
                first part before the loop enters
      \param argc `argc` from `main`
      \param argv `argv` from `main`
      \param parameter_fname Preferable absolute path to parameter file.
      \warning Absolute path differs on different operating systems.
     */

    /*! \name Loop action */
    //@{
    void prePattern_loop();
    /*!< \brief In some sense calculates the inital data for
                the solution driven problem.
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

      At this stage the tumor is growing.
      \warning Do not forget to use next_timeStep() afterwards.
     */
    void final_action();
    /*!< \brief To be used after the second for-loop to save some data. */
    //@}
  private:    
    /*! \name Data members */
    //@{
    Io::Parameter data;
    /*!< \brief Contains parameter from `tumor_parameter.txt`. */
    struct Io{
      SecOrd_op::Identity identity {};
      const Esfem::Io::Dgf::Handler dgf_handler;
      Io::Error_stream u;
      Io::Error_stream w;
      Io(const Io::Parameter&);
    } io;
    /*!< \brief Members are used for input and output of
                the nodal values of the finite element functions
		in #fef.
    */
    /*! \struct Io
     \brief Shortens Brusselator_scheme().
     \sa `io`
    */
    Grid::Grid_and_time fix_grid; /*!< \brief Non evolving grid */
    struct Fef{
      Scal_FEfun_set u;
      Scal_FEfun_set w;
      Vec_FEfun_set surface;
      Fef(Grid::Grid_and_time&);
    } fef;
    /*!< \brief Collects all finite element functions */
    /*! \struct Fef
     \brief Shortens Brusselator_scheme().
     \sa `fef`
    */
    //@}

    // ------------------------------------------------------------
    // Helper member functions
    
    /*! \name Helper classes for the for-loops */
    //@{
    friend class PrePattern_helper;
    friend class Pattern_helper;    
    //@}

    void pre_loop_action();
    /*!< \warning Should only be invoked by Brusselator_scheme(). */
    void intermediate_surface_rhs();
    /*!< \warning Should only be invoked by Brusselator_scheme(). */
    void solve_surfacePDE();
    /*!< \brief Solves the surface PDE and prints out a dgf file. */
    
    /*! \name Flow control */
    //@{
    void next_timeStep(); 
    /*!< \brief Increments the next time step in #fix_grid. */
    long prePattern_timeSteps() const; 
    /*!< \brief Maximum number of time steps for prePattern_loop(). */
    long pattern_timeSteps() const; 
    /*!< \brief Maximum number of time steps for pattern_loop(). */
    //@}
  };
  /*!< \brief Implementation of the Elliott and Styles
              full discretization of the tumor problem
  */

  // ----------------------------------------------------------------------
  // Inline implementation

  inline void Brusselator_scheme::next_timeStep(){
    fix_grid.next_timeStep(data.global_timeStep());
  }
}

#endif // BRUSSELATOR_ALGO_H
