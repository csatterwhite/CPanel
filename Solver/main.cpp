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
    std::string inFile = "NACA4412_fineTE.tri";
    std::string outFile = "NACA4412_alpha5_CGA_v100.vtu";
    time_t ts;
    time(&ts);
    time_t tf;
    geometry temp(path+inFile);
    
    runCase case1(&temp,100,5,0,path+outFile);
    
    time(&tf);
    std::cout << "Elapsed time for program execution : " << difftime(tf,ts) << " seconds" << std::endl;
}
