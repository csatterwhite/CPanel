//
//  convexHullTests.h
//  CPanel
//
//  Created by Chris Satterwhite on 9/24/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__convexHullTests__
#define __CPanel__convexHullTests__

#include <iostream>
#include "cpptest/cpptest.h"
#include "convexHull.h"

class convexHullTests : public Test::Suite
{
    void testPntInside();
    void testPntOutside();
    void testPntOnBoundary();
    
public:
    convexHullTests();
};

#endif /* defined(__CPanel__convexHullTests__) */
