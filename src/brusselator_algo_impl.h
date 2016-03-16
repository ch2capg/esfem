/*! \file brusselator_algo.h

    \brief <Program Name>

     Revision history:

          Revised by Christian Power March 2016
          Originally written by Christian Power
               (power22c@gmail.com) February 2016

     Pseudocode:
       (1+τγ)(Mu)ⁿ⁺¹ + τ (Au)ⁿ⁺¹ - τγ Mⁿ⁺¹(uⁿ,wⁿ)uⁿ⁺¹ = (Mu)ⁿ + τγa Mⁿ⁺¹1
       (Mw)ⁿ⁺¹ + τDᶜ (Aw)ⁿ⁺¹ + τγ Mⁿ⁺¹(uⁿ⁺¹,uⁿ⁺¹)wⁿ⁺¹ = (Mw)ⁿ + τγb Mⁿ⁺¹1

     Created by Christian Power on 11.03.2016
     Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#ifndef BRUSSELATOR_ALGO_IMPL_H
#define BRUSSELATOR_ALGO_IMPL_H 

#include <string>
#include "esfem.h"

namespace Esfem{
  template<typename FEfun>
  struct FEfun_set{
    FEfun fun;
    FEfun app;
    FEfun exact;
    FEfun rhs_les;
    explicit FEfun_set(const std::string& name, const Esfem::Grid::Grid_and_time&);
    explicit FEfun_set(const FEfun_set&, const Esfem::Grid::Grid_and_time&);
    FEfun_set& operator>>(const Esfem::Io::Dgf::Handler&);
    FEfun_set& operator<<(const Esfem::Io::Dgf::Handler&);
    void write(const Esfem::Io::Dgf::Handler&,
	       const std::string& dir = "/tmp/");
    void read(const Esfem::Io::Dgf::Handler&,
	      const std::string& dir = "/tmp/");  
  };
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

  // ----------------------------------------------------------------------
  // helper functions

  void first_interpolate(const Init_data&,
			 FEfun_set<Esfem::Grid::Scal_FEfun>& u,
			 FEfun_set<Esfem::Grid::Scal_FEfun>& w);
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
  void clog_uw(const FEfun_set<Esfem::Grid::Scal_FEfun>& u, 
	       const FEfun_set<Esfem::Grid::Scal_FEfun>& w);

  // ----------------------------------------------------------------------
  // Template implementation

  template<typename FEfun>
  FEfun_set<FEfun>::FEfun_set(const std::string& name,
			      const Esfem::Grid::Grid_and_time& gt)
    : fun {name, gt}, app {name + "_app", gt},
			 exact {name + "_exact", gt}, rhs_les {name + "_rhs_les", gt}
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
  void FEfun_set<FEfun>::write(const Esfem::Io::Dgf::Handler& h, const std::string& dir){
    constexpr auto suffix = ".dgf";
    h.write(dir + fun.name() + suffix, fun);
    h.write(dir + app.name() + suffix, app);
    h.write(dir + exact.name() + suffix, exact);
    h.write(dir + rhs_les.name() + suffix, rhs_les);
  }
  template<typename FEfun>
  void FEfun_set<FEfun>::read(const Esfem::Io::Dgf::Handler& h, const std::string& dir){
    constexpr auto suffix = ".dgf";
    h.read(dir + fun.name() + suffix, fun);
    h.read(dir + app.name() + suffix, app);
    h.read(dir + exact.name() + suffix, exact);
    h.read(dir + rhs_les.name() + suffix, rhs_les);
  }

} // namespace Esfem

#endif // BRUSSELATOR_ALGO_IMPL_H
