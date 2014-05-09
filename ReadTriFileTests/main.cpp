//
//  main.cpp
//  ReadTriFileTests
//
//  Created by Chris Satterwhite on 4/28/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include <iostream>
#include "PanelTest.h"
#include "GeomTest.h"
#include "SurfaceTest.h"

int main()
{
    Test::Suite tests;
    tests.add(std::auto_ptr<Test::Suite>(new GeomTest));
    tests.add(std::auto_ptr<Test::Suite>(new SurfaceTest));
    tests.add(std::auto_ptr<Test::Suite>(new PanelTest));
    
    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
}

