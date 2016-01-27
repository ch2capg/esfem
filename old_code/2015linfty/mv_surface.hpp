/*********************************************************************
 *  mv_surface.hpp                                                   *
 *                                                                   *
 *  Defines auxiliar classes like a string vector with the correct   *
 *  names.                                                           *
 *                                                                   *
 *  Created by Christian Power on 19.06.14.                          *
 *  Copyright (c) 2014 Christian Power. All rights reserved.         *
 *                                                                   *
 *********************************************************************/

#ifndef MV_SURFACE_HPP
#define MV_SURFACE_HPP

namespace BDF
{
  struct  TimeProviderType
  {
    const double beginTime;
    const double endTime;
    const double femStepSize;
    const int femCounter;
    double time;
    //CAPG: do we really need this variable???
        
    TimeProviderType(double,double,double);
    // initialize beginTime, endTime, femStepSize;
    // calculates femCounter (== (endTime - beginTime)/ femStepSize)
    // sets time == beginTime
  };
}

//#include "mv_surface.hxx"  //include template related definitions

#endif // #ifndef MV_SURFACE_HPP
