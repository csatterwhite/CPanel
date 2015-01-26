//
//  influenceTests.h
//  CPanel
//
//  Created by Chris Satterwhite on 1/24/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__influenceTests__
#define __CPanel__influenceTests__

#include <stdio.h>
#include <Eigen/Dense>
#include "cpptest/cpptest.h"
#include "bodyPanel.h"

class influenceTests : public Test::Suite
{
    Eigen::MatrixXd testNodes;
    bodyPanel* testPan;
    
    bool isEqual(double var1, double var2, int decimalPrecision);
    
    void createTestPanel();
    void testPntClose();
    void testPntFar();
//    void testCollocationPnt();
public:
    influenceTests();
    
    ~influenceTests() {delete testPan;}
};

#endif /* defined(__CPanel__influenceTests__) */
