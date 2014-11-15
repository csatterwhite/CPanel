//
//  GeomTest.h
//  CPanel
//
//  Created by Chris Satterwhite on 3/27/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__GeomTest__
#define __CPanel__GeomTest__

#include <iostream>
#include "cpptest/cpptest.h"
#include "Eigen/Dense"
#include <fstream>
#include "geometry.h"
#include "surface.h"
#include "liftingSurf.h"
#include "wake.h"

class GeomTest : public Test::Suite
{
    void test_readGeom();
    int countTEpanels(surface* surf);
    
public:
    GeomTest();
};

#endif /* defined(__CPanel__GeomTest__) */
