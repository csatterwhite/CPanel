//
//  main.cpp
//  Solver
//
//  Created by Chris Satterwhite on 4/30/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include <iostream>
#include "geometry.h"
#include "runCase.h"

int main(int argc, const char * argv[])
{
    time_t ts;
    time(&ts);
    time_t tf;
    geometry temp("ellipse.tri");
    
    runCase case1(&temp,1,0,0);
    
    time(&tf);
    double t = difftime(tf,ts);
    std::cout << "Elapsed time for program execution : " << t << " seconds" << std::endl;
}
