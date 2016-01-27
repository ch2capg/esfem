/*********************************************************************
 *  dune_typedef_heat.hpp                                            *
 *                                                                   *
 *  This header just contains lots of typedef by the original author *
 *  of heat.cc, heat.hh, etc.                                        *
 *                                                                   *
 *  Revision history:                                                *
 *  none                                                             *
 *                                                                   *
 *                                                                   *
 *  Created by Christian Power on 19.06.14.                          *
 *  Copyright (c) 2014 Christian Power. All rights reserved.         *
 *                                                                   *
 *********************************************************************/

#ifndef DUNE_TYPEDEF_HEAT_HPP
#define DUNE_TYPEDEF_HEAT_HPP

#include "deformation.hh"
#include "tumor_growth.h"

typedef Dune::GridSelector::GridType HostGridType;
typedef Dune::GeometryGrid< HostGridType, DeformationCoordFunction > GridType;

// CAPG: Perhaps this is usefull
// create host grid part consisting of leaf level elements
// typedef Dune::Fem::LeafGridPart< HGridType > HostGridPartType;
// HostGridPartType hostGridPart( grid );

typedef Dune::Fem::AdaptiveLeafGridPart< GridType, Dune::InteriorBorder_Partition > 
GridPartType; 
typedef Dune::Fem::FunctionSpace< double, double, GridType::dimensionworld, 1 > 
FunctionSpaceType;

typedef Dune::Fem::L2Norm< GridPartType > L2NormType;
typedef Dune::Fem::H1Norm< GridPartType > H1NormType;

typedef Problem::RHSFunction< FunctionSpaceType > RHSFunctionType;
typedef Problem::InitialData< FunctionSpaceType > InitialDataType;

// choose type of discrete function space
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, 
						  GridPartType, POLORDER > 
DiscreteFunctionSpaceType;

// choose type of discrete function
typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > 
DiscreteFunctionType;

// choose inverse operator
typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType > 
LinearInverseOperatorType;

// select Lagrange interpolation
typedef Dune::LagrangeInterpolation< DiscreteFunctionType > InterpolationType;

// define Laplace operator
typedef Problem::HeatModel< FunctionSpaceType > HeatModelType;
// typedef EllipticOperator< DiscreteFunctionType, HeatModelType > EllipticOperatorType;

// define nonlinear operator
// typedef Problem::NonlinearHeatModel<FunctionSpaceType> NonlinearModel;
// typedef NonlinearHeatOperator<DiscreteFunctionType, 
// 			      NonlinearModel>		   NonlinearOperator;

// type of input/output
typedef Dune::tuple<DiscreteFunctionType*, DiscreteFunctionType*> IOTupleType;
// type of data writer
typedef Dune::Fem::DataOutput< GridType, IOTupleType > DataOutputType;

// ----------------------------------------------------------------------
// Vector valued stuff
using Vec_Fun_Space = Dune::Fem::
  FunctionSpace<double, double, GridType::dimensionworld, 3>;
// using Vec_Init_Data = Problem::InitialData<Vec_Fun_Space>;  
// Instead of defining a costume InitialData_Vec we take identity.  
using Vec_FE_Space = Dune::Fem::
  LagrangeDiscreteFunctionSpace<Vec_Fun_Space, GridPartType, POLORDER>; 
using Vec_FE_Fun = Dune::Fem::AdaptiveDiscreteFunction<Vec_FE_Space>;
using Vec_CG_Solver = Dune::Fem::CGInverseOperator<Vec_FE_Fun>;
using R3Interpolation = Dune::LagrangeInterpolation<Vec_FE_Fun>;
using R3Identity = Tumor_growth::Identity<Vec_Fun_Space>;
// ----------------------------------------------------------------------

// MCF model
// using Vec_Model = MCF_Model<Vec_Fun_Space>;
// using Vec_Operator = Vec_EllipticOperator<Vec_FE_Fun, Vec_Model>

// reg MCF model
// using Vec_Model = Reg_MCF_Model<Vec_Fun_Space>;
// using Vec_Operator = Reg_Vec_EllipticOperator<Vec_FE_Fun, Vec_Model>;

// Tumor growth model
using Tumor_surface_model 
= Tumor_growth::Surface_evolution_model<Vec_Fun_Space, FunctionSpaceType>;
using Tumor_surface_operator 
= Tumor_growth::Surface_parabolic_operator<Vec_FE_Fun, DiscreteFunctionType,
					   Tumor_surface_model>;
using Tumor_Brusselator_model
= Tumor_growth::Brusselator_model<FunctionSpaceType>;
using Tumor_Brusselator_operator
= Tumor_growth::Scalar_parabolic_operator<DiscreteFunctionType, Tumor_Brusselator_model>;

using Initial_data_gauss = Tumor_growth::Initial_data_gauss<FunctionSpaceType>;
using Initial_data_uniform = Tumor_growth::Initial_data_uniform<FunctionSpaceType>;

#endif // #ifndef DUNE_TYPEDEF_HEAT_HPP
