//
//  InfluenceTests.h
//  CPanel
//
//  Created by Chris Satterwhite on 10/7/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__InfluenceTests__
#define __CPanel__InfluenceTests__

#include <iostream>
#include <Eigen/Dense>
#include "cpptest/cpptest.h"
#include "panel.h"
#include "convexHull.h"

class influenceTests : public Test::Suite
{
    void testPntClose();
    void testPntFar();
    // void testPntOnBoundary();
    
public:
    influenceTests();
};

#endif /* defined(__CPanel__InfluenceTests__) */
