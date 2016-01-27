#ifndef ELLIPTIC_HH
#define ELLIPTIC_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>

template< class DiscreteFunction, class Model >
struct EllipticOperator
: public Dune::Fem::Operator< DiscreteFunction >
{
  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model ModelType;

protected:
  typedef typename DiscreteFunctionType
  ::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  //typedef Dune::CachingQuadrature< GridPartType, 0 > QuadratureType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

public:
	explicit EllipticOperator(const ModelType& model_input = Model {} )
		: model{model_input} 
		{}
      
	virtual void
	operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const;
	/**
	 *	EllipticOperator::operator() takes FE-function u and calculates the 
	 *	FE-function w = L(u), which is meant in the weak sense.  We don't assemble 
	 *	the matrix, but we calculate the result 
	 *		"dof_vector(w) = matrix * dof_vector(u)" 
	 *	locally for each element.  The matrix is determined by the member functions 
	 *	'model.massFlux()' and 'model.diffusiveFlux()'.  
	 */

	template< class JacobianOperator >
	void jacobian ( const DiscreteFunction &u, JacobianOperator &jOp ) const;

	void mass_matrix(const DiscreteFunctionType& u, 
					 DiscreteFunctionType& w, const double cdbl = 1) const;
	// calculates: dof(w) = cdbl * M * dof(u), where 'M' is the mass matrix.
	void stiffness_matrix(const DiscreteFunctionType& u, 
						  DiscreteFunctionType& w, const double cdbl = 1) const;
	// calculates: dof(w) = cdbl * A * dof(u), where 'A' is the stiffness matrix.
	void load_vector(DiscreteFunctionType& w) const;
	// calculates: dof(w) = M * dof( Interpolation of 'f'), where 'f' is the RHS
	// of the PDE.

	void set_timeStepFactor(const double d) {model.set_timeStepFactor(d);}
	double read_timeStepFactor() const {return model.read_timeStepFactor();} 
private:
	ModelType model;
};

template < class DiscreteFunction, class Model >
// requires 
// Model:
//  template <typename Entity, typename Point>
//  void jac_ran_to_jac_ran(const Entity entity, const Point xx, 
//			    const JacobianRangeType du_x, const RangeType xi_x, 
//			    JacobianRangeType a_xi_du_x) const;
//  template <typename Entity, typename Point>
//  void ran_to_jac_ran(const Entity entity, const Point xx, const RangeType u_x, 
//    		        JacobianRangeType vec_u_x) const;
//  template <typename Entity, typename Point>
//  void jac_ran_to_ran(const Entity entity, const Point xx, 
//			const JacobianRangeType du_x, RangeType du_x_vec) const;
//  template <typename Entity, typename Point>
//  void ran_to_ran(const Entity entity, const Point xx, const RangeType u_x, 
//		    RangeType u_x_chi) const;
struct NonlinearHeatOperator
  : Dune::Fem::Operator<DiscreteFunction>{
  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model ModelType;
protected:
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType 
  DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
public:
  NonlinearHeatOperator(const DiscreteFunctionType& xi_input,
			const ModelType& model_input,
			const double bdf_factor = 1.)
    : xi {xi_input}, model {model_input}, bdf_alpha_lead_coeff {bdf_factor}
  {}
  
  virtual void
  operator() ( const DiscreteFunctionType& u, DiscreteFunctionType& w ) const;
   // EllipticOperator::operator() takes FE-function u and calculates the 
   // FE-function w = L(ξ)u, which is meant in the weak sense.  ξ is a priori 
   // given FE-function.  We don't assemble the matrix, but we calculate the result 
   //		"dof_vector(w) = matrix(dof-vector(ξ)) * dof_vector(u)" 
   // locally for each element.  The matrix is determined by the member functions 
   // 'model.massFlux()' and 'model.diffusiveFlux()'.  
	
  void mass_matrix(const DiscreteFunctionType& u, 
		   DiscreteFunctionType& w) const;
   // calculates: dof(w) = cdbl * M * dof(u), where 'M' is the mass matrix.

private:
  const DiscreteFunctionType& xi;
  const ModelType model;
  const double bdf_alpha_lead_coeff;
};

// All implementation details

// Implementation details for 'EllipticOperator'

template< class DiscreteFunction, class Model >
void EllipticOperator< DiscreteFunction, Model >
  ::operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const 
/*
  Pseudo code:
  for_each entity in mesh{
  E = entity;		// E = (xᴱ_n)_n; n iterates throw all knots of E
  vu = local(u,E);	// vu == (u(xᴱ_n))_n; part of dof_vector(u)
  avu = massFlux(vu);	// avu == (m(xᴱ_n) φⱼᴱ(xᴱ_n))_n; 
			// m is a reell valued function;
  du = local_jac(u,E);	// du == (∇u(xᴱ_n))_n
  adu = diffusiveFlux(u,E); // adu == (D(xᴱ_n) ∇φⱼᴱ(xᴱ_n))_n; 
			    // D is a diffusion tensor
    for_each point in entity{ 
      pt = point;				// pt = xᴱ_n
      weight = E.weight(pt) * E.integrationElement(pt);	
		// weight == ωᴱ_n = ω (xᴱ_n) * (det(dψᵀ·dψ))^½(xᴱ_n)
      avu *= weight; adu *= weight;
      wLocal.axpy( pt, avu, adu);
      }	
      // after this for loop 'axpy' has implicitly calculated the matrix element:
      // aᴱ(φⱼᴱ,φᵢᴱ):= ∑_n ωᴱ_n [D∇φⱼᴱ(xᴱ_n)·∇φᵢᴱ(xᴱ_n) + m φⱼᴱ(xᴱ_n) φᵢᴱ (xᴱ_n)]
      // it multipied the a part of dof_vector(u) with this matrix element and
      // assigned the result to the same part of dof_vector(w).
      }
      The matrtix element reads	
      aᴱ(φⱼᴱ,φᵢᴱ):= ∑_n ωᴱ_n [D∇φⱼᴱ(xᴱ_n)·∇φᵢᴱ(xᴱ_n) + m φⱼᴱ(xᴱ_n) φᵢᴱ (xᴱ_n)],
      with quadratur[ pt ] == xᴱ_n, avu == m φⱼᴱ(xᴱ_n), adu == D∇φⱼᴱ(xᴱ_n),
      ωᴱ_n == ω_pt (det(dψᵀ·dψ))^½, D = id_2x2;
*/
{
  w.clear();

  const DiscreteFunctionSpaceType &dfSpace = w.space();
  
  // original code
  //const IteratorType end = dfSpace.end();
  //for( IteratorType it = dfSpace.begin(); it != end; ++it )

  // equivalent but clearer; note: for -O1 or higher this is equivalent optimal
  //for( IteratorType it = dfSpace.begin(); it != dfSpace.end(); ++it )

  // c++11 style
  for(auto& it : dfSpace) // the constructor of 'it' is protected, hence auto&
  {
    // const EntityType &entity = *it;
    const EntityType& entity = it;
    const GeometryType& geometry = entity.geometry();

    const LocalFunctionType uLocal = u.localFunction( entity );
    LocalFunctionType wLocal = w.localFunction( entity );

    QuadratureType quadrature( entity, uLocal.order() + wLocal.order() );

    
    for( size_t pt = 0; pt < quadrature.nop(); ++pt ){
      typename LocalFunctionType::RangeType vu, avu;
      uLocal.evaluate( quadrature[ pt ], vu );
      // assign 'vu': vu = uLocal(pt);
      model.massFlux( entity, quadrature[ pt ], vu, avu );
       // assign 'avu':
       // avu = vu, for explicit and implicit; cf. heat.hh
       // CAPG: massFlux because this should be m_k φᵢᴱ(x_nᴱ)

       // typename LocalFunctionType::RangeType new_avu;
       // model.massFlux_2( entity, quadrature[ pt ], vu, new_avu );
       // avu += new_avu;
	  
      typename LocalFunctionType::JacobianRangeType du, adu, Bu;
      uLocal.jacobian( quadrature[ pt ], du );
       // assign 'du': du = grad(uLocal)(pt);
      model.diffusiveFlux( entity, quadrature[ pt ], du, adu );
       // adu *= timeStepFactor;
       // assign 'adu': adu == Δt·du (resp. adu == deltaT * du)
      model.aleFlux( entity, quadrature[pt], vu, Bu);
      adu += Bu;
	  
      const typename QuadratureType::CoordinateType& x = quadrature.point( pt );
      const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );
       // weigth = ω_pt (det(dψᵀ·dψ))^½ 
       // (resp. weight = \omega_pt * \sqrt{det((d\psi)^T (d\psi))});
      avu *= weight; // weight = ω_pt (det(dψᵀ·dψ))^½ 
      adu *= weight; 
      wLocal.axpy( quadrature[ pt ], avu, adu );
      /***********************************************************************
       * CAPG explanation:
       * elliptic operator takes u and calculates w == L(u);
       * we don't assemble the matrix, but we calculate  the result
       * matrix * vector; the matrix element reads
       * aᴱ(φⱼᴱ,φᵢᴱ):= ωᴱ_n D∇φⱼᴱ(xᴱ_n)·∇φᵢᴱ(xᴱ_n) + m φⱼᴱ(xᴱ_n) φᵢᴱ (xᴱ_n);
       * quadratur[ pt ] == xᴱ_n, avu == m φⱼᴱ(xᴱ_n), adu == D∇φⱼᴱ(xᴱ_n),
       * ωᴱ_n == ω_pt (det(dψᵀ·dψ))^½, D = id_2x2;
       ***********************************************************************/
       // 1. variable means which coordinate,
       // 2. variable must be a RangeFieldType (scalar)
       // 3. variable must be a DiscreteFunctionInterfaceType (vector)
    }
  }
  w.communicate();
}

template< class DiscreteFunction, class Model >
template< class JacobianOperator >
void EllipticOperator< DiscreteFunction, Model >
::jacobian ( const DiscreteFunction &u, JacobianOperator &jOp ) const
{
	typedef typename JacobianOperator::LocalMatrixType LocalMatrixType;
	typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
	
	const DiscreteFunctionSpaceType& dfSpace = u.space();
	
	jOp.reserve();
	jOp.clear();
	
	std::vector<typename LocalFunctionType::RangeType> 
		phi( dfSpace.mapper().maxNumDofs() );
	std::vector<typename LocalFunctionType::JacobianRangeType> 
		dphi( dfSpace.mapper().maxNumDofs() );

	/*
	const IteratorType end = dfSpace.end();
	for( IteratorType it = dfSpace.begin(); it != end; ++it )
	{
	*/
	//for(IteratorType it = dfSpace.begin(); it != dfSpace.end(); ++it ){
	for(auto& it : dfSpace){
		const EntityType& entity = *it;
		const GeometryType& geometry = entity.geometry();

		LocalMatrixType jLocal = jOp.localMatrix( entity, entity );

		const BaseFunctionSetType& baseSet = jLocal.domainBaseFunctionSet();
		const unsigned int numBaseFunctions = baseSet.numBaseFunctions();
          
		QuadratureType quadrature( entity, 2*dfSpace.order() );

		/*
		const size_t numQuadraturePoints = quadrature.nop();
		for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
		{
		*/
		for(size_t pt = 0; pt < quadrature.nop(); ++pt ){
			const typename QuadratureType::CoordinateType&
				x = quadrature.point( pt );
			const double 
				weight = quadrature.weight( pt ) * geometry.integrationElement( x );

			const typename GeometryType::Jacobian& 
				gjit = geometry.jacobianInverseTransposed( x );
			baseSet.evaluateAll( quadrature[ pt ], phi );
			baseSet.jacobianAll( quadrature[ pt ], gjit, dphi );

			for( unsigned int i = 0; i < numBaseFunctions; ++i ){
				typename LocalFunctionType::RangeType aphi;
				model.massFlux( entity, quadrature[ pt ], phi[ i ], aphi );

				typename LocalFunctionType::JacobianRangeType adphi;
				model.diffusiveFlux( entity, quadrature[ pt ], dphi[ i ], adphi );
				//dphi *= timeStepFactor;
				//cout << "Bäm!" << endl;

				for( unsigned int j = 0; j < numBaseFunctions; ++j ){
					double value = aphi * phi[ j ];
					for( int k = 0; k < adphi.rows; ++k )
						value += adphi[ k ] * dphi[ j ][ k ];
					jLocal.add( j, i, weight * value );
				}
			}
		}
	}
}

template< class DiscreteFunction, class Model >
void EllipticOperator< DiscreteFunction, Model >
::mass_matrix(const DiscreteFunctionType& u, 
			  DiscreteFunctionType& w, const double cdbl) const
{
	w.clear();

	const DiscreteFunctionSpaceType &dfSpace = w.space();
  
	for(auto& it : dfSpace){ // the constructor of 'it' is protected, hence auto&
		const EntityType& entity = it;
		const GeometryType& geometry = entity.geometry();

		const LocalFunctionType uLocal = u.localFunction( entity );
		LocalFunctionType wLocal = w.localFunction( entity );

		QuadratureType quadrature( entity, uLocal.order() + wLocal.order() );

		for( size_t pt = 0; pt < quadrature.nop(); ++pt ){
			typename LocalFunctionType::RangeType vu, avu;
			uLocal.evaluate( quadrature[ pt ], vu );
			// assign 'vu': vu = uLocal(pt);
			model.massFlux( entity, quadrature[ pt ], vu, avu );
			// assign 'avu': avu = vu;
			// CAPG: massFlux because this should be m_k φᵢᴱ(x_nᴱ)

			const typename QuadratureType::CoordinateType& x = quadrature.point( pt );
			const double weight = 
				quadrature.weight( pt ) * geometry.integrationElement( x );
			// weigth = ω_pt (det(dψᵀ·dψ))^½ 
			// (resp. weight = \omega_pt * \sqrt{det((d\psi)^T (d\psi))});
			avu *= weight * cdbl; // weight = ω_pt (det(dψᵀ·dψ))^½ 
			wLocal.axpy( quadrature[ pt ], avu, 0 );
			/***********************************************************************
			 * CAPG explanation:
			 * elliptic operator takes u and calculates w == L(u);
			 * we don't assemble the matrix, but we calculate  the result
			 * matrix * vector; the matrix element reads
			 * aᴱ(φⱼᴱ,φᵢᴱ):= ωᴱ_n D∇φⱼᴱ(xᴱ_n)·∇φᵢᴱ(xᴱ_n) + m φⱼᴱ(xᴱ_n) φᵢᴱ (xᴱ_n);
			 * quadratur[ pt ] == xᴱ_n, avu == m φⱼᴱ(xᴱ_n), adu == D∇φⱼᴱ(xᴱ_n),
			 * ωᴱ_n == ω_pt (det(dψᵀ·dψ))^½, D = id_2x2;
			 ***********************************************************************/
			// 1. variable means which coordinate,
			// 2. variable must be a RangeFieldType (scalar)
			// 3. variable must be a DiscreteFunctionInterfaceType (vector)
		}
	}
	w.communicate();
}

template< class DiscreteFunction, class Model >
void EllipticOperator<DiscreteFunction,Model>
::stiffness_matrix(const DiscreteFunctionType& u, 
				   DiscreteFunctionType& w, const double cdbl) const{
	w.clear();
	const DiscreteFunctionSpaceType &dfSpace = w.space();

	for(auto& it : dfSpace){ // the constructor of 'it' is protected, hence auto&
		const EntityType& entity = it;
		const GeometryType& geometry = entity.geometry();
		const LocalFunctionType uLocal = u.localFunction( entity );
		LocalFunctionType wLocal = w.localFunction( entity );
		QuadratureType quadrature( entity, uLocal.order() + wLocal.order() );

		for( size_t pt = 0; pt < quadrature.nop(); ++pt ){
			typename LocalFunctionType::JacobianRangeType du, adu;
			uLocal.jacobian( quadrature[ pt ], du );
			// assign 'du': du = grad(uLocal)(pt);
			model.diffusiveFlux( entity, quadrature[ pt ], du, adu );
			// assign 'adu':
			// adu == Δt·du (resp. adu == deltaT * du) for implicit elliptic
			// adu == 0, for explicit elliptic

			const typename QuadratureType::CoordinateType& x = quadrature.point( pt );
			const double 
				weight = quadrature.weight( pt ) * geometry.integrationElement( x );
			// weigth = ω_pt (det(dψᵀ·dψ))^½ 
			// (resp. weight = \omega_pt * \sqrt{det((d\psi)^T (d\psi))});
			adu *= weight; 
			wLocal.axpy( quadrature[ pt ], 0, adu );
		}
	}
	w.communicate();
}

template< class DiscreteFunction, class Model >
void EllipticOperator<DiscreteFunction,Model>
::load_vector(DiscreteFunctionType& w) const{
	w.clear();
	const DiscreteFunctionSpaceType &dfSpace = w.space();

	for(auto& it : dfSpace){ // the constructor of 'it' is protected, hence auto&
		const EntityType& entity = it;
		const GeometryType& geometry = entity.geometry();
		LocalFunctionType wLocal = w.localFunction( entity );
		QuadratureType quadrature( entity, 2 * wLocal.order() + 1);

		for( size_t pt = 0; pt < quadrature.nop(); ++pt ){
			typename LocalFunctionType::RangeType vec_w;
			model.rhs( entity, quadrature[ pt ], vec_w );
			// assign 'avu': vec_w = f( quadratur[pt] );
			
			const typename QuadratureType::CoordinateType& x = quadrature.point( pt );
			const double weight = 
				quadrature.weight( pt ) * geometry.integrationElement( x );
			// weigth = ω_pt (det(dψᵀ·dψ))^½ 
			// (resp. weight = \omega_pt * \sqrt{det((d\psi)^T (d\psi))});
			vec_w *= weight; // weight = ω_pt (det(dψᵀ·dψ))^½ 
			wLocal.axpy( quadrature[ pt ], vec_w, 0 );
		}
	}
	w.communicate();
}

// Implementation details for 'NonlinearHeatOperator'

template < class DiscreteFunction, class Model >
void NonlinearHeatOperator<DiscreteFunction, Model>
  ::operator() (const DiscreteFunctionType& u, DiscreteFunctionType& w) const 
/*
  Pseudo code:
  for_each entity in mesh{
  E = entity;		// E = (xᴱ_n)_n; n iterates throw all knots of E
  vu = local(u,E);	// vu == (u(xᴱ_n))_n; part of dof_vector(u)
  avu = massFlux(vu);	// avu == (m(xᴱ_n) φⱼᴱ(xᴱ_n))_n; 
			// m is a reell valued function
  du = local_jac(u,E);	// du == (∇u(xᴱ_n))_n
  adu = diffusiveFlux(u,E);	// adu == (D(xᴱ_n) ∇φⱼᴱ(xᴱ_n))_n; 
				// D is a diffusion tensor
  for_each point in entity{ 
    pt = point;		// pt = xᴱ_n
    weight = E.weight(pt) * E.integrationElement(pt); 
     // weight == ωᴱ_n = ω (xᴱ_n) * (det(dψᵀ·dψ))^½(xᴱ_n)
    avu *= weight; adu *= weight;
    wLocal.axpy( pt, avu, adu);
   }	// after this for loop 'axpy' has implicitly calculated the matrix element:
   // aᴱ(φⱼᴱ,φᵢᴱ):= ∑_n ωᴱ_n [D∇φⱼᴱ(xᴱ_n)·∇φᵢᴱ(xᴱ_n) + m φⱼᴱ(xᴱ_n) φᵢᴱ (xᴱ_n)]
   // it multipied the a part of dof_vector(u) with this matrix element and
   // assigned the result to the same part of dof_vector(w).
  }
  The matrtix element reads	
  aᴱ(φⱼᴱ,φᵢᴱ):= ∑_n ωᴱ_n [D∇φⱼᴱ(xᴱ_n)·∇φᵢᴱ(xᴱ_n) + m φⱼᴱ(xᴱ_n) φᵢᴱ (xᴱ_n)],
  with quadratur[ pt ] == xᴱ_n, avu == m φⱼᴱ(xᴱ_n), adu == D∇φⱼᴱ(xᴱ_n),
  ωᴱ_n == ω_pt (det(dψᵀ·dψ))^½, D = id_2x2;
*/
{
  w.clear();

  const DiscreteFunctionSpaceType &dfSpace = w.space();
  
  for(auto& it : dfSpace) // the constructor of 'it' is protected, hence auto&
  {
    const EntityType& entity = it;
    const GeometryType& geometry = entity.geometry();

    const LocalFunctionType uLocal = u.localFunction(entity);
    const LocalFunctionType xiLocal = xi.localFunction(entity);
    LocalFunctionType wLocal = w.localFunction(entity);

    QuadratureType quadrature(entity, uLocal.order() + wLocal.order() + 1);

    for(size_t pt = 0; pt < quadrature.nop(); ++pt){

      typename LocalFunctionType::RangeType u_x, xi_x;
      // u_x == u(x), xi_x == ξ(x)
      uLocal.evaluate(quadrature[pt], u_x);	 // initialize u_x
      // xiLocal.evaluate(quadrature[pt], xi_x); // initialize xi_x

      typename LocalFunctionType::JacobianRangeType du_x;
      // du_x == u'(x) (jacobian matrix)
      uLocal.jacobian(quadrature[pt], du_x);	 // initialize du_x

      typename LocalFunctionType::JacobianRangeType a_xi_du_x;
      // a_du_x ≈ ∫ A(ξ) ∇u · ∇χ  (∇ stands for gradient matrix)
      model.jac_ran_to_jac_ran(entity, quadrature[pt], du_x, xi_x, a_xi_du_x);
      // initialize a_xi_du_x

      typename LocalFunctionType::JacobianRangeType vec_u_x;
      // vec_u_x ≈ ∫ (u V) · ∇χ  (∇ stands for gradient matrix, 
      // V is an appropriate vectorfield)
      model.ran_to_jac_ran(entity, quadrature[pt], u_x, vec_u_x);
      // initialize vec_u_x

      typename LocalFunctionType::RangeType du_x_vec;
      // vec_du_x ≈ ∫ (V · ∇u) χ  (∇ stands for gradient matrix, 
      // V is an appropriate vectorfield)
      model.jac_ran_to_ran(entity, quadrature[pt], du_x, du_x_vec);
      // initialize du_x_vec

      typename LocalFunctionType::RangeType u_x_chi;
      // chi_u_x ≈ ∫ u χ
      model.ran_to_ran(entity, quadrature[pt], u_x, u_x_chi);
      u_x_chi *= bdf_alpha_lead_coeff;
      // initialize u_x_chi

      const typename QuadratureType::CoordinateType& x = quadrature.point(pt);
      const double weight = quadrature.weight(pt) * geometry.integrationElement(x);
      // weigth = ω_pt (det(dψᵀ·dψ))^½ 
		
      du_x_vec *= weight;
      u_x_chi *= weight;
      a_xi_du_x *= weight;
      vec_u_x *= weight;
      wLocal.axpy(quadrature[pt], du_x_vec , a_xi_du_x);
      wLocal.axpy(quadrature[pt], u_x_chi , vec_u_x);
      // calculating matrix * vector without assembling the matrices
    }
  }
  w.communicate();
}

template <class DiscreteFunction, class Model>
void NonlinearHeatOperator<DiscreteFunction, Model>
::mass_matrix(const DiscreteFunctionType& u, 
			  DiscreteFunctionType& w) const{
  w.clear();

  const DiscreteFunctionSpaceType &dfSpace = w.space();
  
  for(auto& it : dfSpace) // the constructor of 'it' is protected, hence auto&
    {
      const EntityType& entity = it;
      const GeometryType& geometry = entity.geometry();

      const LocalFunctionType uLocal = u.localFunction(entity);
      LocalFunctionType wLocal = w.localFunction(entity);

      QuadratureType quadrature(entity, uLocal.order() + wLocal.order() + 1);

      for(size_t pt = 0; pt < quadrature.nop(); ++pt){
	typename LocalFunctionType::RangeType u_x;
	// u_x == u(x)
	uLocal.evaluate(quadrature[pt], u_x);		// initialize u_x

	typename LocalFunctionType::RangeType u_x_chi;
	// chi_u_x ≈ ∫ u χ
	model.ran_to_ran(entity, quadrature[pt], u_x, u_x_chi);
	// initialize u_x_chi

	const typename QuadratureType::CoordinateType& x = quadrature.point(pt);
	const double weight = quadrature.weight(pt) * geometry.integrationElement(x);
	// weigth = ω_pt (det(dψᵀ·dψ))^½ 
		
	wLocal.axpy(quadrature[pt], u_x_chi * weight, 0);
	// calculating matrix * vector * cdbl without assembling the matrices
      }
    }
  w.communicate();
}

#endif // #ifndef ELLIPTIC_HH
