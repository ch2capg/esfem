/*! \file tumor_growth_operator.h

    \brief Tumor growth model (2012-Elliott+Styles) + our regularizing term

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 22.11.2015

     This header provides model classes and operator classes to solve 
     a tumor growth model proposed by Elliott and Styles via the ESFEM.  
     Parameters: 
     - a,b,Dᶜ,γ ∈ R for the scalar equation with  γ ∼ vol(Γ).
     - ε,δ,α ∈ R for the surface equation, where ε and α are small 
       regularization parameter.
     PDE: 
     Find u:Γ→R (growth-promoting), w:Γ→R (growth-inhibiting) and
     X: Γ⁰×[0,T]→Rᵐ⁺¹ such that
  		∂•u + u div(v) - Δu   = f₁(u,w),
  		∂•w + w div(v) - DᶜΔw = f₂(u,w),
     where for f₁,f₂ we use the Brusselator model
  		f₁(u,w) = γ(a-u+u²w) ∧ f₂(u,w) = γ(b-u²w),
     and for the surface 
  		v - αΔv = (ε (-H) + δu) n = εΔX + δun,
  		    ∂ᵗX = v(X).
     FEM discretization:
     Find u:I→Rᴺ (growth-promoting nodal values), w:I→Rᴺ (growth-inhibiting 
     nodal values) and X:I→R³ᴺ (surface nodal values) such that
		(M(X) + α A(X))∂ᵗX = ε A(X) X + δ M(u,n)
		∂ᵗ(M(X)u) + A(X)u = γ(b M(X)1 + M(X;u,w)u)
		∂ᵗ(M(X)w) + Dᶜ A(X)w = γ(b M(X)1 - M(X;u,u)w)

     Pseudocode (Elliott+Styles discretization):
     - Given Xⁿ,uⁿ,wⁿ, solve for Xⁿ⁺¹
  
       (M₃ⁿ + (α+ετ)A₃ⁿ)Xⁿ⁺¹ = (M₃ⁿ+αA₃ⁿ)Xⁿ + τδ M₃ⁿ(uⁿ,nⁿ),
  
       where nⁿ is elementwise normal.  For the moment we choose an alternative
       normal namely nⁿ = A₃ⁿ id / ‖A₃ⁿ id‖.
     - Given Xⁿ⁺¹,uⁿ,wⁿ, solve for uⁿ⁺¹
  
       (1+τγ)(Mu)ⁿ⁺¹ + τ (Au)ⁿ⁺¹ - τγ Mⁿ⁺¹(uⁿ,wⁿ)uⁿ⁺¹ = (Mu)ⁿ + τγa Mⁿ⁺¹1
  
       where M(a,b) a 4 tensor is, namely Mⁱʲᵏˡ = ∫χⁱ̧χʲχᵏχˡ and 1 means the 
       constant 1 finite element function.
     - Given Xⁿ⁺¹,uⁿ⁺¹,wⁿ, solve for wⁿ⁺¹
  
       (Mw)ⁿ⁺¹ + τDᶜ (Aw)ⁿ⁺¹ + τγ Mⁿ⁺¹(uⁿ⁺¹,uⁿ⁺¹)wⁿ⁺¹ = (Mw)ⁿ + τγb Mⁿ⁺¹1


     Created by Christian Power on 22.11.2015
     Copyright (c) 2015 Christian Power.  All rights reserved.
 */

#ifndef TUMOR_GROWTH_OPERATOR_H
#define TUMOR_GROWTH_OPERATOR_H

#include <cmath>
#include <random>
#include <functional>

#include <dune/fem/function/common/function.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/quadrature/quadrature.hh>

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>

#include <heat.hh>

namespace Tumor_growth{
  // ----------------------------------------------------------------------
  // Initial data for the Brusselator model

  template<typename Function_space>
  struct Initial_data_u :
    Dune::Fem::Function<Function_space, Initial_data_u<Function_space> >{

    explicit Initial_data_u(const double mean_value, const double variance)
      : random_fun (std::bind(Random_dist {mean_value, variance},
			  Random_engine {}))
    {}

    using Domain = typename Function_space::DomainType;
    using Range = typename Function_space::RangeType;

    void evaluate(const Domain& p, Range& q) const{
      q = random_fun();
    }
  private:
    using Random_dist = std::normal_distribution<>;
    using Random_engine = std::default_random_engine;
    std::function<double()> random_fun;
  };

  template<typename Function_space>
  struct Initial_data_w :
    Dune::Fem::Function<Function_space, Initial_data_w<Function_space> >{

    explicit Initial_data_w(const double mean_value, const double variance)
      : base {mean_value}, pertubation {variance}
    {}
    
    using Domain = typename Function_space::DomainType;
    using Range = typename Function_space::RangeType;

    void evaluate(const Domain& p, Range& q) const{
      // const double  x = p[0], y = p[1], z = p[2];
      q = Range {base};
    }
  private:
    double base;
    double pertubation;
  };
  
  // ----------------------------------------------------------------------
  // Surface and Brusselator Models

  template <typename Vector_function_space, typename Scalar_function_space>
  struct Surface_evolution_model{
    using Vector = Vector_function_space;
    using Scalar = Scalar_function_space;

    template<typename T>
    using Domain = typename T::DomainType;
    template<typename T>
    using Range = typename T::RangeType;
    template<typename T>
    using JacobianRange = typename T::JacobianRangeType;
    template<typename T>
    using DomainField = typename T::DomainFieldType;
    template<typename T>
    using RangeField = typename T::RangeFieldType;
    static const int dim_vec_domain = Domain<Vector>::dimension;
    static const int dim_vec_range = Range<Vector>::dimension;
    static const int dim_scalar_domain = Domain<Scalar>::dimension;
    static const int dim_scalar_range = Range<Scalar>::dimension;

    explicit
    Surface_evolution_model (const Dune::Fem::TimeProviderBase& time_provider,
			     const double alpha_regularization_parameter,
			     const double epsilon_regularization_parameter,
			     const double delta_parameter) :
      tp {time_provider}, reg_alpha {alpha_regularization_parameter},
      reg_epsilon {epsilon_regularization_parameter}, delta {delta_parameter}
    {}

    // Pseudocode:
    // (M₃ⁿ + (α+ετ)A₃ⁿ)Xⁿ⁺¹ = (M₃ⁿ+αA₃ⁿ)Xⁿ + τδ M₃ⁿ(uⁿ,nⁿ)
    // alpha_regularization_parameter == α,
    // epsilon_regularization_parameter == ε,
    // delta_parameter == δ,
    // time_provider.deltaT() == τ.

    template <typename Entity, typename Point>
    void jac_ran_to_jac_ran_lhs(const Entity& entity, const Point& xx, 
				const JacobianRange<Vector>& dX_p,
				JacobianRange<Vector>& a_dX_p) const;
    // a_dX_p = (α+ετ)A₃ⁿdX_p
    template <typename Entity, typename Point>
    void jac_ran_to_jac_ran_rhs(const Entity& entity, const Point& xx, 
				const JacobianRange<Vector>& dX_p,
				JacobianRange<Vector>& a_dX_p) const;
    // a_dX_p = αA₃ⁿdX_p
    template <typename Entity, typename Point>
    void ran_to_ran(const Entity& entity, const Point& xx, 
		    const Range<Vector>& X_p, Range<Vector>& X_p_chi) const;
    // X_p_chi = M₃ⁿX_p
    template <typename Entity, typename Point>
    void ran_u_normal_rhs(const Entity& entity, const Point& xx, 
			  const Range<Scalar>& u_p,
			  Range<Vector>& n_p_chi) const;
    // n_p_chi = τδM₃ⁿ(u_p,nⁿ)
  private:
    const Dune::Fem::TimeProviderBase& tp;
    const double reg_alpha;
    const double reg_epsilon;
    const double delta;
  };

  enum class Growth {promoting, inhibiting};

  template <typename Scalar_function_space>
  class Brusselator_model{
  public:
    using Scalar = Scalar_function_space;
    using Range = typename Scalar::RangeType;
    using JacobianRange = typename Scalar::JacobianRangeType;

    explicit
    Brusselator_model(const double uniform_time_step,
		      const double gamma, 
		      const double a_or_b, 
		      const Growth scalar_type,
		      const double diffusion = 1.);
    // Pseudocode:
    // growth-promoting
    // (1+τγ)(Mu)ⁿ⁺¹ + τ (Au)ⁿ⁺¹ - τγ Mⁿ⁺¹(uⁿ,wⁿ)uⁿ⁺¹ = (Mu)ⁿ + τγa Mⁿ⁺¹1
    // growth-inhibiting
    // (Mw)ⁿ⁺¹ + τDᶜ (Aw)ⁿ⁺¹ + τγ Mⁿ⁺¹(uⁿ⁺¹,uⁿ⁺¹)wⁿ⁺¹ = (Mw)ⁿ + τγb Mⁿ⁺¹1

    template <typename Entity, typename Point>
    void ran_to_ran(const Entity& entity, const Point& xx, 
		    const Range& u_p, Range& u_p_chi) const;
    // u_p_chi = (Mu)ⁿ || u_p_chi = (Mw)ⁿ 

    template <typename Entity, typename Point>
    void ran_to_ran_ones(const Entity& entity, const Point& xx, 
			 Range& ones_chi) const;
    // ones_chi = τγa Mⁿ⁺¹1 || ones_chi = τγb Mⁿ⁺¹1 

    template <typename Entity, typename Point>
    void ran_to_ran_lhs(const Entity& entity, const Point& xx, 
			 const Range& u_p, Range& u_p_chi) const;
    // u_p_chi = (1+τγ)(Mu)ⁿ⁺¹ || u_p_chi = (Mw)ⁿ⁺¹ 

    template <typename Entity, typename Point>
    void jac_ran_to_jac_ran_lhs(const Entity& entity, const Point& xx, 
				const JacobianRange& du_p,
				JacobianRange& du_p_dchi) const;
    // du_p_dchi = τ(Au)ⁿ⁺¹ || du_p_dchi = τDᶜ(Aw)ⁿ⁺¹ 

    template <typename Entity, typename Point>
    void tri_ran_to_ran_lhs(const Entity& entity, const Point& xx,
			    const Range& first, const Range& second,
			    const Range& u_p, Range& uwu_p_chi) const;
    // uwu_p_chi = - τγ Mⁿ⁺¹(uⁿ,wⁿ)uⁿ⁺¹ || uwu_p_chi = τγ Mⁿ⁺¹(uⁿ⁺¹,uⁿ⁺¹)wⁿ⁺¹

  private:
    using Domain = typename Scalar::DomainType;

    // ----------
    // data members
    double mass_matrix_lhs {0.};
    double stiffness_matrix_lhs {0.};
    double quadrilinear_mass_matrix_lhs {0.};
    double mass_matrix_rhs {0.};
  };

  // ----------------------------------------------------------------------
  // Parabolic operators

  template <typename Vec_discrete_function, typename Scalar_discrete_function, 
	    typename Model>
  /* requires Model:
     template <typename Entity, typename Point>
     void jac_ran_to_jac_ran_lhs(const Entity& entity, const Point& xx, 
     const JacobianRange<Vector>& dX_p,
     JacobianRange<Vector>& a_dX_p) const;
     template <typename Entity, typename Point>
     void jac_ran_to_jac_ran_rhs(const Entity& entity, const Point& xx, 
     const JacobianRange<Vector>& dX_p,
     JacobianRange<Vector>& a_dX_p) const;
     template <typename Entity, typename Point>
     void ran_to_ran(const Entity& entity, const Point& xx, 
     const Range<Vector>& X_p, Range<Vector>& X_p_chi) const;
     template <typename Entity, typename Point>
     void ran_u_normal_rhs(const Entity& entity, const Point& xx, 
     const Range<Scalar>& u_p,
     Range<Vector>& n_p_chi) const;
  */
  struct Surface_parabolic_operator
    : Dune::Fem::Operator<Vec_discrete_function>{
    using Vector = Vec_discrete_function;
    using Scalar = Scalar_discrete_function;
  protected:
    using Discrete_function_space = typename Vector::DiscreteFunctionSpaceType;
    using Iterator = typename Discrete_function_space::IteratorType;
    using Entity = typename Iterator::Entity;
    using Geometry = typename Entity::Geometry;
    using GridPart = typename Discrete_function_space::GridPartType;
    using QuadratureType = Dune::Fem::CachingQuadrature<GridPart, 0>;

    template<typename T>
    using Local_function = typename T::LocalFunctionType;
    template<typename T>
    using Range = typename Local_function<T>::RangeType;
    template<typename T>
    using Jacobian_range = typename Local_function<T>::JacobianRangeType;

  public:
    explicit
    Surface_parabolic_operator(const Model& model_input,
			       const Scalar_discrete_function& u_previous,
			       const double bdf_alpha_lead_coeff = 1.)
      : model (model_input), u_approx (u_previous), bdf_factor {bdf_alpha_lead_coeff}
    {}
  
    virtual void operator() (const Vec_discrete_function& u, 
			     Vec_discrete_function& w) const;
    // Pseudocode:
    // w = (M₃ⁿ + (α+ετ)A₃ⁿ)u

    void generate_rhs(const Vec_discrete_function& u, Vec_discrete_function& w) const;
    // Pseudocode:
    // w = (M₃ⁿ+αA₃ⁿ)u + τδ M₃ⁿ(u_previous,nⁿ)

  private:
    const Model model;
    const Scalar_discrete_function& u_approx;
    const double bdf_factor;
  };

  template <typename Scalar_discrete_function, typename Model>
  /* requires Model:
     template <typename Entity, typename Point>
     void ran_to_ran(const Entity& entity, const Point& xx, 
		     const Range& u_p, Range& u_p_chi) const;
     template <typename Entity, typename Point>
     void ran_to_ran_ones(const Entity& entity, const Point& xx, 
                          Range& ones_chi) const;
     template <typename Entity, typename Point>
     void ran_to_ran_lhs(const Entity& entity, const Point& xx, 
			 const Range& u_p, Range& u_p_chi) const;
     template <typename Entity, typename Point>
     void jac_ran_to_jac_ran_lhs(const Entity& entity, const Point& xx, 
				 const JacobianRange& du_p,
				 JacobianRange& du_p_dchi) const;
     template <typename Entity, typename Point>
     void tri_ran_to_ran_lhs(const Entity& entity, const Point& xx,
 			     const Range& first, const Range& second,
			     const Range& u_p, Range& uwu_p_chi) con st;
  */
  struct Scalar_parabolic_operator
    : Dune::Fem::Operator<Scalar_discrete_function>{
  public:
    using FE_function = Scalar_discrete_function;

    explicit
    Scalar_parabolic_operator(const Model& model_input,
			      const FE_function& u_previous,
			      const FE_function& u_or_w_previous,
			      const double bdf_alpha_lead_coeff = 1.)
      : model (model_input), u_approx (u_previous), u_or_w_approx (u_or_w_previous),
	bdf_factor {bdf_alpha_lead_coeff}
    {}

    virtual void operator() (const FE_function& u, FE_function& w) const;
    // Pseudocode:
    // w = (1+τγ)(Mu)ⁿ⁺¹ + τ (Au)ⁿ⁺¹ - τγ Mⁿ⁺¹(uⁿ,wⁿ)uⁿ⁺¹  ||
    // w = (Mw)ⁿ⁺¹ + τDᶜ (Aw)ⁿ⁺¹ + τγ Mⁿ⁺¹(uⁿ⁺¹,uⁿ⁺¹)wⁿ⁺¹

    void generate_rhs_old_surface(const FE_function& u, FE_function& w) const;
    // Pseudocode: w = (Mu)ⁿ || w = (Mw)ⁿ

    void generate_rhs_new_surface(FE_function& w) const;
    // Pseudocode: w = τγa Mⁿ⁺¹1 || w = τγb Mⁿ⁺¹1
  private:
    using Discrete_function_space = typename FE_function::DiscreteFunctionSpaceType;
    using Iterator = typename Discrete_function_space::IteratorType;
    using Entity = typename Iterator::Entity;
    using Geometry = typename Entity::Geometry;
    using GridPart = typename Discrete_function_space::GridPartType;
    using QuadratureType = Dune::Fem::CachingQuadrature<GridPart, 0>;

    using Local_function = typename FE_function::LocalFunctionType;
    using Range = typename Local_function::RangeType;
    using Jacobian_range = typename Local_function::JacobianRangeType;
    
    // ----------
    // Data members
    const Model model;
    const FE_function& u_approx;
    const FE_function& u_or_w_approx;
    const double bdf_factor;
  };

  // ======================================================================
  // All implementation details

  // ----------------------------------------------------------------------
  // Implementation details for 'Surface_evolution_model'
  template <typename Vector_function_space, typename Scalar_function_space>
  template <typename Entity, typename Point>
  inline void Surface_evolution_model<Vector_function_space, Scalar_function_space>::
  jac_ran_to_jac_ran_lhs(const Entity& entity, const Point& xx, 
			 const JacobianRange<Vector>& dX_p,
			 JacobianRange<Vector>& a_dX_p) const{
    // DomainType xGlobal = entity.geometry().global(Dune::coordinate(xx));
    // xGlobal.size() == dim_domain == 3;
    // DomainFieldType x = xGlobal[0], y = xGlobal[1], z = xGlobal[2];
    // double t = tp.time(), dT = tp.deltaT();
    // 'dX_p' is a R²ˣ³ Matrix, so it has two indices.		
    a_dX_p = dX_p;
    a_dX_p *= (reg_alpha + reg_epsilon * tp.deltaT());
  }
  template <typename Vector_function_space, typename Scalar_function_space>
  template <typename Entity, typename Point>
  inline void Surface_evolution_model<Vector_function_space, Scalar_function_space>::
  jac_ran_to_jac_ran_rhs(const Entity& entity, const Point& xx, 
			 const JacobianRange<Vector>& dX_p,
			 JacobianRange<Vector>& a_dX_p) const{
    // DomainType xGlobal = entity.geometry().global(Dune::coordinate(xx));
    // xGlobal.size() == dim_domain == 3;
    // DomainFieldType x = xGlobal[0], y = xGlobal[1], z = xGlobal[2];
    // double t = tp.time(), dT = tp.deltaT();
    // 'dX_p' is a R²ˣ³ Matrix, so it has two indices.		
    a_dX_p = dX_p;
    a_dX_p *= reg_alpha;
  }
  template <typename Vector_function_space, typename Scalar_function_space>
  template <typename Entity, typename Point>
  inline void Surface_evolution_model<Vector_function_space, Scalar_function_space>::
  ran_to_ran(const Entity& entity, const Point& xx, 
	     const Range<Vector>& X_p, Range<Vector>& X_p_chi) const{
    // DomainType xGlobal = entity.geometry().global(Dune::coordinate(xx));
    // xGlobal.size() == dim_domain == 3;
    // DomainFieldType x = xGlobal[0], y = xGlobal[1], z = xGlobal[2];
    // double t = tp.time(), dT = tp.deltaT();
    X_p_chi = X_p;
  }
  template <typename Vector_function_space, typename Scalar_function_space>
  template <typename Entity, typename Point>
  void Surface_evolution_model<Vector_function_space, Scalar_function_space>::
  ran_u_normal_rhs(const Entity& entity, const Point& xx, const Range<Scalar>& u_p,
		   Range<Vector>& n_p_chi) const{
    // Domain<Vector> xGlobal = entity.geometry().global(Dune::coordinate(xx));
    // xGlobal.size() == dim_domain == 3;
    // DomainFieldType x = xGlobal[0], y = xGlobal[1], z = xGlobal[2];
    // double t = tp.time(), dT = tp.deltaT();

    // Consistency check.  Comment this later out.
    // if(entity.geometry().corners() != 3) throw std::runtime_error
    //   {"Error in Surface_evolution_model::ran_u_normal_rhs().\n"
    // 	  "Number of corners unequal to 3."};
    Domain<Vector> p0 = entity.geometry().corner(0);
    Domain<Vector> p1 = entity.geometry().corner(1);
    Domain<Vector> p2 = entity.geometry().corner(2);
    
    // Calculating normal vector

    // Tangent vector components
    Domain<Vector> v = p2 - p0;
    Domain<Vector> w = p1 - p0;

    auto eucl_norm = [](const Range<Vector>& vv){
      return std::sqrt(vv[0] * vv[0] + vv[1] * vv[1] + vv[2] * vv[2]);
    };

    // Formula for cross product
    // n = v × w / ‖v × w‖
    n_p_chi[0] = v[1] * w[2] - v[2] * w[1];
    n_p_chi[1] = - v[0] * w[2] + v[2] * w[0];
    n_p_chi[2] = v[0] * w[1] - v[1] * w[0];
    n_p_chi /= eucl_norm(n_p_chi);

    // The final result
    n_p_chi *= u_p * tp.deltaT() * delta; 
  }

  // ----------------------------------------------------------------------
  // Implementation details for 'Brusselator_model'

  template <typename Scalar_function_space>
  Brusselator_model<Scalar_function_space>::
  Brusselator_model(const double uniform_time_step,
		    const double gamma, const double a_or_b, const Growth scalar_type,
		    const double diffusion)
  {
    switch(scalar_type){
    case Growth::promoting:
      mass_matrix_lhs = 1. + uniform_time_step * gamma;
      stiffness_matrix_lhs = uniform_time_step;
      quadrilinear_mass_matrix_lhs = (-1.) * uniform_time_step * gamma;
      break;
    case Growth::inhibiting:
      mass_matrix_lhs = 1.;
      stiffness_matrix_lhs = uniform_time_step * diffusion;
      quadrilinear_mass_matrix_lhs = uniform_time_step * gamma;
      break;
    default:
      throw std::runtime_error {"Error in constructor of Brusselator_model.\n"
	  "Your scalar_type is not implemented."};
      break;
    };
    mass_matrix_rhs = uniform_time_step * gamma * a_or_b;
  }

  template <typename Scalar_function_space>
  template <typename Entity, typename Point>
  inline void Brusselator_model<Scalar_function_space>::
  ran_to_ran(const Entity& entity, const Point& xx, 
	     const Range& u_p, Range& u_p_chi) const{
    u_p_chi = u_p;
  }
  template <typename Scalar_function_space>
  template <typename Entity, typename Point>
  inline void Brusselator_model<Scalar_function_space>::
  ran_to_ran_ones(const Entity& entity, const Point& xx, 
		  Range& ones_chi) const{
    // ones_chi = τγa Mⁿ⁺¹1 || ones_chi = τγb Mⁿ⁺¹1     
    ones_chi = Range {mass_matrix_rhs};
  }
  template <typename Scalar_function_space>
  template <typename Entity, typename Point>
  inline void Brusselator_model<Scalar_function_space>::
  ran_to_ran_lhs(const Entity& entity, const Point& xx, 
		 const Range& u_p, Range& u_p_chi) const{
    // u_p_chi = (1+τγ)(Mu)ⁿ⁺¹ || u_p_chi = (Mw)ⁿ⁺¹
    u_p_chi = u_p;
    u_p_chi *= mass_matrix_lhs;
  }
  template <typename Scalar_function_space>
  template <typename Entity, typename Point>
  inline void Brusselator_model<Scalar_function_space>::
  jac_ran_to_jac_ran_lhs(const Entity& entity, const Point& xx, 
			 const JacobianRange& du_p,
			 JacobianRange& du_p_dchi) const{
    // du_p_dchi = τ(Au)ⁿ⁺¹ || du_p_dchi = τDᶜ(Aw)ⁿ⁺¹
    du_p_dchi = du_p;
    du_p_dchi *= stiffness_matrix_lhs;
  }
  template <typename Scalar_function_space>
  template <typename Entity, typename Point>
  inline void Brusselator_model<Scalar_function_space>::
  tri_ran_to_ran_lhs(const Entity& entity, const Point& xx,
		     const Range& first, const Range& second,
		     const Range& u_p, Range& uwu_p_chi) const{
    // uwu_p_chi = - τγ Mⁿ⁺¹(uⁿ,wⁿ)uⁿ⁺¹ || uwu_p_chi = τγ Mⁿ⁺¹(uⁿ⁺¹,uⁿ⁺¹)wⁿ⁺¹
    uwu_p_chi = u_p;
    uwu_p_chi *= first * second * quadrilinear_mass_matrix_lhs;
  }
  
  // ----------------------------------------------------------------------
  // Implementation details for 'Surface_parabolic_operator'

  template <typename Vec_discrete_function, typename Scalar_discrete_function, 
  	    typename Model>
  void Surface_parabolic_operator<Vec_discrete_function, 
  				  Scalar_discrete_function, Model>::
  operator() (const Vec_discrete_function& u, Vec_discrete_function& w) const{
    w.clear();

    const Discrete_function_space& dfSpace = w.space();
    for(const Entity& entity : dfSpace){
      // auto == tmp class Iterator<Vector>
      // auto& == tmp class Entity<Vector> convertible to Entity<Vector>
      // The constructor of 'Entity<Vector>' is protected, hence auto&
      // Remark: Entity<Vector> == Entity<Scalar>!
  
      const Geometry& geometry = entity.geometry();
      
      const Local_function<Vector> uLocal = u.localFunction(entity);
      Local_function<Vector> wLocal = w.localFunction(entity);
    
      QuadratureType quadrature(entity, uLocal.order() + wLocal.order() + 1);
    
      for(size_t pt = 0; pt < quadrature.nop(); ++pt){
    
  	Range<Vector> u_x;	// u_x == X(x);
  	uLocal.evaluate(quadrature[pt], u_x);	// initialize
    
  	Jacobian_range<Vector> du_x;	// du_x == ∇X(x) (jacobian matrix)
  	uLocal.jacobian(quadrature[pt], du_x);	// initialize
    
  	Range<Vector> u_x_psi;
  	// u_x_psi ≈ ∫ X Ψ
  	model.ran_to_ran(entity, quadrature[pt], u_x, u_x_psi);	
  	// initialize
  	// u_x_chi *= bdf_factor;
  
  	Jacobian_range<Vector> du_x_dpsi;
  	// du_x_dpsi ≈ (α+ετ)∫ ∇X · ∇Ψ  (∇ stands for gradient matrix)
  	model.jac_ran_to_jac_ran_lhs(entity, quadrature[pt], du_x, du_x_dpsi);
  	// initialize 
    
  	const typename QuadratureType::CoordinateType& x = quadrature.point(pt);
  	const double weight = quadrature.weight(pt) * geometry.integrationElement(x);
  	// weigth = ω_pt (det(dψᵀ·dψ))^½ 
    		
  	u_x_psi *= weight;
  	du_x_dpsi *= weight;
  	wLocal.axpy(quadrature[pt], u_x_psi, du_x_dpsi);
  	// calculating matrix * vector without assembling the matrices
      }
    }
    w.communicate();
  }
  template <typename Vec_discrete_function, typename Scalar_discrete_function, 
	    typename Model>
  void Surface_parabolic_operator<Vec_discrete_function, 
				  Scalar_discrete_function, Model>::
  generate_rhs(const Vec_discrete_function& u, Vec_discrete_function& w) const{
    w.clear();

    const Discrete_function_space& dfSpace = w.space();
    for(const Entity& entity : dfSpace){ 
      // auto& == tmp class Entity<Vector> convertible to Entity<vector>&
      // Entity<Vector> == Entity<Scalar>!

      const Geometry& geometry = entity.geometry();

      const Local_function<Vector> uLocal = u.localFunction(entity);
      const Local_function<Scalar> u_approx_local 
	= u_approx.localFunction(entity);
      Local_function<Vector> wLocal = w.localFunction(entity);

      QuadratureType quadrature(entity, uLocal.order() + wLocal.order() + 1);

      for(size_t pt = 0; pt < quadrature.nop(); ++pt){
	Range<Vector> u_x;	// u_x == X(x);
	uLocal.evaluate(quadrature[pt], u_x);	// initialize

	Range<Scalar> u_approx_x;	// u_approx_x = u(x);
	u_approx_local.evaluate(quadrature[pt], u_approx_x);

	Jacobian_range<Vector> du_x;	// du_x == ∇X(x) (jacobian matrix)
	uLocal.jacobian(quadrature[pt], du_x);	// initialize

	Range<Vector> u_x_psi;
	// u_x_psi ≈ ∫ X Ψ
	model.ran_to_ran(entity, quadrature[pt], u_x, u_x_psi);	// initialize
	// u_x_chi *= bdf_factor;

	Jacobian_range<Vector> du_x_dpsi;
	// du_x_dpsi ≈ α ∫ ∇X · ∇Ψ  (∇ stands for gradient matrix)
	model.jac_ran_to_jac_ran_rhs(entity, quadrature[pt], du_x, du_x_dpsi);
	// initialize 

	Range<Vector> un_x_psi;
	// un_x_psi ≈ τδ ∫ u_approx n Ψ
	model.ran_u_normal_rhs(entity, quadrature[pt], u_approx_x, un_x_psi);

	const typename QuadratureType::CoordinateType& x = quadrature.point(pt);
	const double weight = quadrature.weight(pt) * geometry.integrationElement(x);
	// weigth = ω_pt (det(dψᵀ·dψ))^½ 
		
	u_x_psi *= weight;
	un_x_psi *= weight;
	du_x_dpsi *= weight;

	wLocal.axpy(quadrature[pt], u_x_psi, du_x_dpsi);
	wLocal.axpy(quadrature[pt], un_x_psi, 0);
	// calculating matrix * vector without assembling the matrices
      }
    }
    w.communicate();
  }

  // ----------------------------------------------------------------------
  // Implementation details for 'Scalar_parabolic_operator'

  template <typename Scalar_discrete_function, typename Model>
  void Scalar_parabolic_operator<Scalar_discrete_function, Model>::
  operator() (const FE_function& u, FE_function& w) const{
    w.clear();

    const Discrete_function_space& dfSpace = w.space();
    for(const Entity& entity : dfSpace){  
      const Geometry& geometry = entity.geometry();

      const Local_function uLocal = u.localFunction(entity);
      const auto u_approx_local = u_approx.localFunction(entity);
      const auto u_or_w_approx_local = u_or_w_approx.localFunction(entity);
      Local_function wLocal = w.localFunction(entity);

      QuadratureType quadrature(entity, uLocal.order() + wLocal.order() + 1);
      for(size_t pt = 0; pt < quadrature.nop(); ++pt){
  	Range u_x;	// u_x == u(x);
  	uLocal.evaluate(quadrature[pt], u_x);	// initialize
    
  	Jacobian_range du_x;	// du_x == ∇u(x) (jacobian matrix)
  	uLocal.jacobian(quadrature[pt], du_x);	// initialize
    
  	Range u_x_chi;
  	// u_x_chi ≈ (1+τγ) ∫ u χ || u_x_chi ≈ ∫ u χ
  	model.ran_to_ran_lhs(entity, quadrature[pt], u_x, u_x_chi);	
  	// initialize
  	// u_x_chi *= bdf_factor;
  
  	Jacobian_range du_x_dchi;
  	// du_x_dchi ≈ τ ∫ ∇u · ∇χ  (∇ stands for gradient matrix)
	// || du_x_dchi ≈ τDᶜ ∫ ∇u · ∇χ
  	model.jac_ran_to_jac_ran_lhs(entity, quadrature[pt], du_x, du_x_dchi);
  	// initialize 
    
  	const typename QuadratureType::CoordinateType& x = quadrature.point(pt);
  	const double weight = quadrature.weight(pt) * geometry.integrationElement(x);
  	// weigth = ω_pt (det(dψᵀ·dψ))^½ 
    		
  	u_x_chi *= weight;
  	du_x_dchi *= weight;
  	wLocal.axpy(quadrature[pt], u_x_chi, du_x_dchi);
  	// calculating matrix * vector without assembling the matrices
      }

      QuadratureType 
	bigger_quadrature(entity, uLocal.order() + u_approx_local.order() 
			  + u_or_w_approx_local.order() + wLocal.order() + 1);
      for(size_t pt = 0; pt < bigger_quadrature.nop(); ++pt){
  	Range u_x;	// u_x == u(x);
  	uLocal.evaluate(bigger_quadrature[pt], u_x);	// initialize
  	Range u_app;	// u_app == u_approx(x);
	u_approx_local.evaluate(bigger_quadrature[pt], u_app);	// initialize
  	Range u_or_w_app;	// u_or_w_app == u_or_w_approx(x);
	u_or_w_approx_local.evaluate(bigger_quadrature[pt], u_or_w_app);
	// initialize

	Range uwu_x_chi;	
	// uwu_x_chi ≈ -τγ ∫ uⁿ wⁿ uⁿ⁺¹ χ || uwu_x_chi ≈ τγ ∫ uⁿ⁺¹uⁿ⁺¹wⁿ⁺¹χ
	model.tri_ran_to_ran_lhs(entity, bigger_quadrature[pt],
				 u_app, u_or_w_app, u_x, uwu_x_chi);

  	const typename QuadratureType::CoordinateType& x = bigger_quadrature.point(pt);
  	const double weight
	  = bigger_quadrature.weight(pt) * geometry.integrationElement(x);
  	// weigth = ω_pt (det(dψᵀ·dψ))^½ 

  	wLocal.axpy(bigger_quadrature[pt], uwu_x_chi, 0);
      }
    }
    w.communicate();
  }  
  template <typename Scalar_discrete_function, typename Model>
  void Scalar_parabolic_operator<Scalar_discrete_function, Model>::
  generate_rhs_old_surface (const FE_function& u, FE_function& w) const{
    w.clear();
  
    const Discrete_function_space& dfSpace = w.space();
    for(const Entity& entity : dfSpace){
  
      const Geometry& geometry = entity.geometry();
      
      const Local_function uLocal = u.localFunction(entity);
      Local_function wLocal = w.localFunction(entity);
    
      QuadratureType quadrature(entity, uLocal.order() + wLocal.order() + 1);
    
      for(size_t pt = 0; pt < quadrature.nop(); ++pt){
    
  	Range u_x;	// u_x == X(x);
  	uLocal.evaluate(quadrature[pt], u_x);	// initialize
    
  	Range u_x_chi;
  	// u_x_chi ≈ ∫ X χ
  	model.ran_to_ran(entity, quadrature[pt], u_x, u_x_chi);	
  	// initialize
  
  	const typename QuadratureType::CoordinateType& x = quadrature.point(pt);
  	const double weight = quadrature.weight(pt) * geometry.integrationElement(x);
  	// weigth = ω_pt (det(dψᵀ·dψ))^½ 
    		
  	u_x_chi *= weight;
  	wLocal.axpy(quadrature[pt], u_x_chi, 0);
  	// calculating matrix * vector without assembling the matrices
      }
    }
    w.communicate();
  }
  template <typename Scalar_discrete_function, typename Model>
  void Scalar_parabolic_operator<Scalar_discrete_function, Model>::
  generate_rhs_new_surface (FE_function& w) const{
    w.clear();
  
    const Discrete_function_space& dfSpace = w.space();
    for(const Entity& entity : dfSpace){
      const Geometry& geometry = entity.geometry();      

      Local_function wLocal = w.localFunction(entity);
   
      QuadratureType quadrature(entity, wLocal.order() + 1);
    
      for(size_t pt = 0; pt < quadrature.nop(); ++pt){

  	Range ones_chi;
  	// ones_chi ≈ τγa ∫ 1 χ || ones_chi ≈ τγb ∫ 1 χ
  	model.ran_to_ran_ones(entity, quadrature[pt], ones_chi);	
  	// initialize

 	const typename QuadratureType::CoordinateType& x = quadrature.point(pt);
  	const double weight = quadrature.weight(pt) * geometry.integrationElement(x);
  	// weigth = ω_pt (det(dψᵀ·dψ))^½ 
    		
  	ones_chi *= weight;
  	wLocal.axpy(quadrature[pt], ones_chi, 0);
  	// calculating matrix * vector without assembling the matrices
      }
    }
    w.communicate();
  }
}	// namespace Tumor_growth

#endif // TUMOR_GROWTH_OPERATOR_H

/*! Log:
    2015.11.27: Save a promising tumor_growth_code()

                I have run the experiment with u_base_value = 1, w_base_value = 0, 
		no perturbation values

    2015.11.26: Crush a bug in the Scalar_parabolic_operator::operator() function
    
                I received a assertion error from an dune array class.  The reason was 
		that a bigger_quadrature object had to replace a quadrature object
		in the second for-loop.  

    2015.11.25: Add Brusselator_operator class without testing it
    2015.11.25: Add Brusselator_model class without testing it
    2015.11.25: Corrected mistake in typedef/using declaration of Surface_operator
    2015.11.24: Give operator() and generate_rhs()
    2015.11.21: Save tested version of the discrete normal vector
 */
