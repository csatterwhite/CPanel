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

int main()
{
    Test::Suite tests;
    tests.add(std::auto_ptr<Test::Suite>(new CoordTransformTest));
    
    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
}

