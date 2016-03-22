/*! \file brusselator_algo.h
    \author Christian Power
    \date 17. March 2016

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

         Created by Christian Power on 17.03.2016
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
    explicit Brusselator_scheme(int argc, char** argv,
				const std::string& parameter_fname);
    /*!< \brief The constructor that also performs the
                first part before the loop enters
      \param argc `argc` from `main`
      \param argv `argv` from `main`
      \param parameter_fname Preferable absolute path to parameter file.
      \warning Absolute path differs on different operating systems.
     */
    ~Brusselator_scheme();

    /*! \name Flow control */
    //@{
    void next_timeStep(); 
    /*!< \brief Increments the next time step for `time_provider`.
                Use this in every for-loop.
     */
    long prePattern_timeSteps() const; 
    /*!< \brief Maximum number of time steps for the first for-loop. */
    long pattern_timeSteps() const; 
    /*!< \brief Maximum number of time steps for the second for-loop. */
    //@}

    /*! \name Loop action */
    //@{
    void pre_pattern_action();
    /*!< \brief To be used in the first for-loop.  

      At this stage the tumor is not growing, but rather the
      inital data is created.
      \warning Do not forget to use `Brusselator_scheme::next_timeStep` afterwards.
     */
    void intermediate_action();
    /*!< \brief To be used between first and second for-loop.

      The inital data has been created.  It will be saved
      and the right hand side for the surface PDE will be created. 
      \warning Do not use `Brusselator_scheme::next_timeStep` afterwards.
     */
    void pattern_action();
    /*!< \brief To be used in the second for-loop.

      At this stage the tumor is growing.
      \warning Do not forget to use `Brusselator_scheme::next_timeStep` afterwards.
     */
    void final_action();
    /*!< \brief To be used after the second for-loop to save some data. */
    //@}
  private:
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
      FEfun_set& operator>>(const Esfem::Io::Dgf::Handler&);
      /*! \brief Specialized `FEfun_set::write` */
      FEfun_set& operator<<(const Esfem::Io::Dgf::Handler&);
      /*! \brief Specialized `FEfun-set::read` */
      void write(const Esfem::Io::Dgf::Handler&,
		 const std::string& dir = "/tmp/");
      /*! \brief Save nodal values in a dgf file. */
      void read(const Esfem::Io::Dgf::Handler&,
		const std::string& dir = "/tmp/");
      /*! \brief Read nodal values from a dgf file. */
    };
    /*!< \brief Minimal amount of finite element functions to perform
                higher order BDF-ESFEM
    */
    using Scal_FEfun_set = FEfun_set<Grid::Scal_FEfun>;
    using Vec_FEfun_set = FEfun_set<Grid::Vec_FEfun>;

    Io::Parameter data;
    Grid::Grid_and_time fix_grid;
    Scal_FEfun_set u;
    Scal_FEfun_set w;
    Vec_FEfun_set surface;

    const Esfem::Io::Dgf::Handler dgf_handler;
    friend class PrePattern_helper;
    friend class Pattern_helper;
    void pre_loop_action(); // To be invoked only in the constructor
  };
  /*!< \brief Implementation of the Elliott and Styles
              full discretization of the tumor problem
  */

  // ----------------------------------------------------------------------
  // Template implementation

  template<typename FEfun>
  FEfun_set<FEfun>::FEfun_set(const std::string& name,
			      const Esfem::Grid::Grid_and_time& gt)
    : fun {name, gt}, app {name + "_app", gt},
    nexact {name + "_exact", gt}, rhs_les {name + "_rhs_les", gt}
  {}
  template<typename FEfun>
  FEfun_set<FEfun>::FEfun_set(const FEfun_set& other, const Esfem::Grid::Grid_and_time& gt)
    : fun {other.fun, gt}, app {other.app, gt}, 
    exact {other.exact, gt}, rhs_les {other.rhs_les, gt}
  {}

  template<typename FEfun>
  FEfun_set<FEfun>& FEfun_set<FEfun>::operator>>(const Esfem::Io::Dgf::Handler& h){
    write(h);
    return *this;
  }
  template<typename FEfun>
  FEfun_set<FEfun>& FEfun_set<FEfun>::operator<<(const Esfem::Io::Dgf::Handler& h){
    read(h);
    return *this;
  }
  template<typename FEfun>
  void FEfun_set<FEfun>::write(const Esfem::Io::Dgf::Handler& h,
			       const std::string& dir){
    h.write(compose_dgfName(fun.name(), dir), fun);
    h.write(compose_dgfName(app.name(), dir), app);
    h.write(compose_dgfName(exact.name(), dir), exact);
    h.write(compose_dgfName(rhs_les.name(), dir), rhs_les);
  }
  template<typename FEfun>
  void FEfun_set<FEfun>::read(const Esfem::Io::Dgf::Handler& h, const std::string& dir){
    h.read(compose_dgfName(fun.name(), dir), fun);
    h.read(compose_dgfName(app.name(), dir), app);
    h.read(compose_dgfName(exact.name(), dir), exact);
    h.read(compose_dgfName(rhs_les.name(), dir), rhs_les);
  }

}

#endif // BRUSSELATOR_ALGO_H
