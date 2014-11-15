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
#include <fstream>
#include "bodyPanel.h"
#include "convexHull.h"

class influenceTests : public Test::Suite
{
//    void testPntClose();
//    void testPntFar();
//    void testCollocationPnt();
//    void testDoublet();
//    void testSource();
    bool isEqual(double var1, double var2, int decimalPrecision);
    
public:
    influenceTests();
};

#endif /* defined(__CPanel__InfluenceTests__) */
