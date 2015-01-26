//
//  main.cpp
//  CPanelTests
//
//  Created by Chris Satterwhite on 1/22/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include <iostream>
#include "geometryTests.h"
#include "influenceTests.h"

int main(int argc, const char * argv[])
{
    Test::Suite tests;
    tests.add(std::auto_ptr<Test::Suite>(new GeomTests));
    tests.add(std::auto_ptr<Test::Suite>(new influenceTests));
    
    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
}
