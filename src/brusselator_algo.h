/*! \file brusselator_algo.h
    \brief Numerical experiment for the solution driven paper

     Revision history
     --------------------------------------------------

          Revised by Christian Power June 2016
          Revised by Christian Power May 2016
          Revised by Christian Power April 2016
          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Idea
     --------------------------------------------------

     This header provides model classes and operator classes to solve 
     a tumor growth model proposed by Elliott and Styles via the ESFEM.
     Enter path to writable directory in the macro variable FEF_PATH.

    \author Christian Power
    \date 7. June 2016
    \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef BRUSSELATOR_ALGO_H
#define BRUSSELATOR_ALGO_H 

#include <memory>
#include <string>
#include "esfem.h"

#ifndef FEF_PATH
#error Give full path to folder in FEF_PATH
#endif 

namespace Esfem{
  //! ESFEM algorithm.  Only this should be invoked by main.
  /*! \param argc `argc` from `main`
    \param[in] argv `argv` from `main`
   */
  void brusselator_algo(int argc, char** argv);

  //! Implementation of the Elliott and Styles full discretization of the tumor problem
  /*!
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
       + g
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
       -\varepsilon A(X)X + \delta M(u,\nodalValue{\surfaceNormal})
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
   */
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
    //! Experiment for the maximum norm paper
    void maxnorm_esfem();
    //! Elliott and Styles ALE example
    /*! Standard linear parabolic equation 
      \f{equation*}{
      \matd u + u \diver(v) - \laplaceBeltrami u 
      = 
      0
      \f}
      with initial data \f$u(0)= 1\f$.  We consider no ALE movement.
      A levelset for our evolution reads as
      \f{equation*}
      d(x,t) := x_1^2 + x_2^2 + K(t)^2 G\Big( \frac{x_3^2}{L(t)^2} \Big)-K(t)^2,
      \f}
      with the functions \f$G\f$, \f$L\f$ and \f$K\f$ given as
      \f{align*}
      G(s)={}& 200s\Big( s - \frac{199}{200}\Big),\\
      L(t)={}& 1 + 0.2\sin(4\pi t),\\
      K(t)={}& 0.1 + 0.05\sin(2\pi t).
      \f}
      We define the velocity \f$v\f$ as the normal velocity 
      of the surface defined by
      the differential equation (formulated for the nodes): 
      \f{equation*}
      \dot{a}_j = V_j \nu_j, \quad \textnormal{where} 
      \quad V_j=\frac{-\partial_t d(a_j,t)}{|\nabla d(a_j,t)|}, 
      \quad \nu_j=\frac{\nabla d(a_j,t)}{|\nabla d(a_j,t)|}.
      \f}
     */
    void ale_aleMovement();
    //! Elliott and Styles ALE example
    /*! Standard linear parabolic equation 
      \f{equation*}{
      \matd u + u \diver(v) - \laplaceBeltrami u 
      = 
      0
      \f}
      with initial data \f$u(0)= 1\f$.  
      We consider the levelset function 
      \f{equation*}
      d(x,t) := x_1^2 + x_2^2 + K(t)^2 G\Big( \frac{x_3^2}{L(t)^2} \Big)-K(t)^2,
      \f}
      with the functions \f$G\f$, \f$L\f$ and \f$K\f$ given as
      \f{align*}
      G(s)={}& 200s\Big( s - \frac{199}{200}\Big),\\
      L(t)={}& 1 + 0.2\sin(4\pi t),\\
      K(t)={}& 0.1 + 0.05\sin(2\pi t).
      \f}
      The velocity \f$v\f$ is the normal velocity of the surface defined by
      the differential equation (formulated for the nodes): 
      \f{equation*}
      \dot{a}_j = V_j \nu_j, \quad \textnormal{where} 
      \quad V_j=\frac{-\partial_t d(a_j,t)}{|\nabla d(a_j,t)|}, 
      \quad \nu_j=\frac{\nabla d(a_j,t)}{|\nabla d(a_j,t)|}.
      \f}
      ALE version of the PDE reads as
      \f{equation*}{
      \matd_A u + u \diver(v_A) - \laplaceBeltrami u
      + \diver( (v - v_A) u )
      = 
      0.
      \f}
      As an ALE movement we define 
      \f{equation*}
      (a_i(t))_1= (a_0(t))_1 \frac{K(t)}{K(0)}, 
      \quad (a_i(t))_2= (a_0(t))_2 \frac{K(t)}{K(0)}, 
      \quad (a_i(t))_3= (a_0(t))_3 \frac{L(t)}{L(0)}, 
      \f}
      Note that \f$d(a_i(t) ,t)=0\f$ for every 
      \f$t\in[0,T]\f$, for \f$i=1,2,\dotsc,N\f$.
     */
    void ale_normalMovement();
    //@}

    //! Experiment, where right-hand side is calculated for a known solution
    /*! As exact solution for the scalar valued surface equation I chose
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
    void eoc_logisticSphere();

    //! Dziuk mean curvature flow ESFEM
    /*! (This is not any more true) Experiment reads as follows: 
      Stationary surface \f$\surface = S^2\f$
      with exact solution \f$X=e^{-t}(xy, yz, xz)\f$,  with the formula
      \f[
      \Delta f = \Delta_{\R^{n+1}} f - H D_{\R^{n+1}} f \cdot n 
      - D^2_{\R^{n+1}} f(n,n),
      \f]
      and 
      \f[
      H = div(n) = \frac{n}{|x|}
      \f]
      one easily sees
      \f[
      \Delta(xy) = -\frac{2n+2}{|x|^2} xy, 
      \f]
      which implies that the right-hand side of the surface PDE vanishes. 
    */
    void eoc_mcf();

    //! Sphere Dalquist
    /*! Initial surface is unit sphere in \f$\R^3\f$.  Exact solution is
     \f[
     \Phi(x,t) := r(t)x,\quad r(t):= e^t
     \f]
     The PDE is 
     \f[
     v - \alpha \Delta v - \varepsilon \Delta x 
     = f(x,t) = (1 + (\alpha + \varepsilon) 2 e^{-2t})x.
     \f]
    */
    void sd();
    //! Surface logistic sphere
    /*! As exact solution for the surface evolution I choose
    \f{equation*}{
     \Phi(x,t) := r(t) x, \quad
     r(t) := \frac{r_{end} r_0}{r_{end} e^{-kt} + r_0 (1-e^{-kt})}.
    \f}
    \f$r(t)\f$ satiesfies the ODE
    \f[
    \dot{r}(t) = k \left( 1 - \frac{r(t)}{r_{end}}\right) r(t),
    \f]
    which implies for the velocity
    \f[
    v(x,t) = k \left( 1 - \frac{r(t)}{r_{end}} \right)x 
    =: \tilde{a} x
    \f]
    I do not consider coupling.  For the surface PDE I choose
    \f[
    v - \alpha \Delta v - \varepsilon \Delta x - \delta u \vec{n} = f,
    \f]
    where \f$f\f$ must be
    \f[
    f = \left( \tilde{a} + 
    \frac{(\alpha \tilde{a} + \varepsilon)H - \delta u}{|x|}\right) x,
    \f]
    where we have used \f$ \Delta x = - H \vec{n}\f$, where \f$H\f$ 
    is the mean curvature  
    (without aritmetic mean) and \f$\vec{n}\f$ is the outwards pointing normal, and 
    \f$ \vec{n} = \frac{x}{|x|}\f$.  It holds \f$H = \frac{n}{|x|}\f$, 
    where \f$n\f$ is the dimension of the sphere.  Note that \f$|x| = r(t)\f$ 
    on the exact surface.
    \sa Esfem::SecOrd_op::vRhs::new_sls(), Esfem::SecOrd_op::vIdata::new_sls()
    */
    void eoc_sls();

    //! EOC experiment for solution driven paper
    /*! As exact solution for the surface evolution I choose
    \f{equation*}{
     \Phi(x,t) := r(t) x, \quad
     r(t) := \frac{r_{end} r_0}{r_{end} e^{-kt} + r_0 (1-e^{-kt})}.
    \f}
    \f$r(t)\f$ satiesfies the ODE
    \f[
    \dot{r}(t) = k \left( 1 - \frac{r(t)}{r_{end}}\right) r(t),
    \f]
    which implies for the velocity
    \f[
    v(x,t) = k \left( 1 - \frac{r(t)}{r_{end}} \right)x 
    =: \tilde{a} x
    \f]
    As exact solution for the scalar diffusion equation I choose
    \f[
    u(x,t) := x_1 x_2 e^{-6t}.
    \f]
    The coupled PDE equations reads as
    \f{gather*}{
    \partial^{\bullet} u + u \nabla\cdot v - \Delta u = f \\
    v - \Delta v - \Delta X - \delta u = g,
    \f}
    where \f$g\f$ is given via
    \f[
    g = \left( \tilde{a} + 
    \frac{(\alpha \tilde{a} + \varepsilon)H - \delta u}{|x|}\right) x,    
    \f]
    and the complicated \f$f\f$ is computed via Sage.
    \pre Parameter \f$\gamma\f$ has to be zero.
    \sa eoc_sls(), Esfem::SecOrd_op::vRhs::new_sls(), 
    Esfem::SecOrd_op::sRhs::new_sdp_u(), Esfem::SecOrd_op::vIdata::new_sls(),
    Esfem::SecOrd_op::sIdata::new_1ssef()
    */
    void eoc_sdp();

    //! Run a simple test
    void test();
    
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

      //! File to capture PDE parameter
      Esfem::Io::Error_stream para;
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
      //! Exact solution of the surface
      std::unique_ptr<SecOrd_op::vIdata> X_ptr;
      
      //! Provides analytically given initial data.
      explicit Init_data(const Grid::Grid_and_time&);
      //! Provides uniform distributed random inital values. 
      Init_data(const Grid::Grid_and_time&, const Esfem::Io::Parameter&);
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
    //! Assign a new value to `fef.u.exact` and `fef.w.exact`
    void update_scalar_solution();

    //! Lifted printing 
    template<class F>
    void print(Esfem::Io::Error_stream& os, const F& fem){
      os << fix_grid.time_provider().deltaT() << ' '
	 << norm.l2_err(fem.fun, fem.exact) << ' '
	 << norm.h1_err(fem.fun, fem.exact) << std::endl;
    }

    //! Calculate the velocity via simple differential quotient
    /*! \param[in] xn_first Iterator to first point of the new surface
      \param[in] xn_last Iterator to last point of the new surface
      \param[in] xo_first Iterator to first point of the previous surface
      \param[out] v_first Iterator to the first value of the velocity
      \pre Both `xo_first` and `v_first` point to as many elements as 
      `xn_first` does. 
    */
    template<class It1, class It2>
    void calculate_velocity(It1 xn_first, It1 xn_last, It1 xo_first, It2 v_first){
      const double dT = fix_grid.time_provider().deltaT();
      for(; xn_first != xn_last; ++xn_first, ++xo_first, ++v_first)
	*v_first = (*xn_first - *xo_first)/dT;
    }

    
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
