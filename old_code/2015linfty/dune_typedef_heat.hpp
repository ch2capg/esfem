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

typedef Dune::GridSelector::GridType HostGridType;
typedef Dune::GeometryGrid< HostGridType, DeformationCoordFunction > GridType;

// CAPG: Perhaps this is usefull
// create host grid part consisting of leaf level elements
// typedef Dune::Fem::LeafGridPart< HGridType > HostGridPartType;
// HostGridPartType hostGridPart( grid );

typedef Dune::Fem::AdaptiveLeafGridPart< GridType, Dune::InteriorBorder_Partition > GridPartType;
typedef Dune::Fem::FunctionSpace< double, double, GridType::dimensionworld, 1 > FunctionSpaceType;

typedef Dune::Fem::L2Norm< GridPartType > L2NormType;
typedef Dune::Fem::H1Norm< GridPartType > H1NormType;

typedef Problem::RHSFunction< FunctionSpaceType > RHSFunctionType;
typedef Problem::InitialData< FunctionSpaceType > InitialDataType;

// choose type of discrete function space
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;

// choose type of discrete function
typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

// choose inverse operator
typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType > LinearInverseOperatorType;

// select Lagrange interpolation
typedef Dune::LagrangeInterpolation< DiscreteFunctionType > InterpolationType;

// define Laplace operator
typedef Problem::HeatModel< FunctionSpaceType > HeatModelType;
typedef EllipticOperator< DiscreteFunctionType, HeatModelType > EllipticOperatorType;

// define nonlinear operator
typedef Problem::NonlinearHeatModel<FunctionSpaceType> NonlinearModel;
typedef NonlinearHeatOperator<DiscreteFunctionType, 
			      NonlinearModel>		   NonlinearOperator;

// type of input/output
typedef Dune::tuple< DiscreteFunctionType* > IOTupleType;
// type of data writer
typedef Dune::Fem::DataOutput< GridType, IOTupleType > DataOutputType;

#endif // #ifndef DUNE_TYPEDEF_HEAT_HPP
