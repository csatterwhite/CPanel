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
    std::string path = "/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files/";
    std::string inFile = "sphere_coarse.tri";
    std::string outFile = "sphere_coarse_testVTUfile.vtu";
    time_t ts;
    time(&ts);
    time_t tf;
    geometry temp(path+inFile);
    
    runCase case1(&temp,1,0,0,path+outFile);
    
    time(&tf);
    double t = difftime(tf,ts);
    std::cout << "Elapsed time for program execution : " << t << " seconds" << std::endl;
}
