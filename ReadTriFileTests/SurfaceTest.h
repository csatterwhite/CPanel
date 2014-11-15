//
//  SurfaceTest.h
//  CPanel
//
//  Created by Chris Satterwhite on 3/27/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__SurfaceTest__
#define __CPanel__SurfaceTest__

#include <iostream>
#include "cpptest/cpptest.h"
#include "Eigen/Dense"
#include "surface.h"

class SurfaceTest : public Test::Suite
{
    void test_addPanel();
    
public:
    SurfaceTest();
};

#endif /* defined(__CPanel__SurfaceTest__) */
