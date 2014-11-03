//
//  main.cpp
//  InfluenceCoefficientTests
//
//  Created by Chris Satterwhite on 5/19/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include <iostream>
#include "cpptest/cpptest.h"
#include "CoordTransformTest.h"
#include "convexHullTests.h"
#include "InfluenceTests.h"

int main()
{
    Test::Suite tests;
//    tests.add(std::auto_ptr<Test::Suite>(new CoordTransformTest));
    tests.add(std::auto_ptr<Test::Suite>(new convexHullTests));
    tests.add(std::auto_ptr<Test::Suite>(new influenceTests));
    
    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
}

