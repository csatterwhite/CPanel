//
//  main.cpp
//  OctreeUnitTests
//
//  Created by Chris Satterwhite on 4/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include <iostream>
#include "cpptest/cpptest.h"
#include "OctreeTests.h"
#include "NodeTests.h"

int main()
{
    Test::Suite tests;
    tests.add(std::auto_ptr<Test::Suite>(new OctreeTests));
    tests.add(std::auto_ptr<Test::Suite>(new NodeTests));
    
    Test::TextOutput output(Test::TextOutput::Terse);
    return tests.run(output);
}

