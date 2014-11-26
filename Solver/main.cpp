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
<<<<<<< HEAD
<<<<<<< HEAD
    std::string inFile = "NACA4412_fineTE.tri";
    std::string outFile = "NACA4412_alpha5_CGA_v100.vtu";
=======
    std::string inFile = "sphere_coarse.tri";
    std::string outFile = "sphere_coarse_testVTUfile.vtu";
>>>>>>> OutputFileWriting
=======
    std::string inFile = "genericAC.tri";
    std::string outFile = "generic_AC_neighborTest.vtu";
>>>>>>> CHTSLS
    time_t ts;
    time(&ts);
    time_t tf;
    geometry temp(path+inFile);
    
    std::string neighborFile = "neighborCheck.txt";
    std::ofstream fid;
    fid.open(path+neighborFile);
    std::vector<panel*> pans = temp.getPanels();
    for (int i=0; i<pans.size(); i++)
    {
        fid << pans[i]->getCenter()(0) << "\t" << pans[i]->getCenter()(1) << "\t" << pans[i]->getCenter()(2) << "\t" << pans[i]->getNeighbors().size() << "\n";
    }
    fid.close();
    
    runCase case1(&temp,1,0,0,path+outFile);
    
    time(&tf);
    std::cout << "Elapsed time for program execution : " << difftime(tf,ts) << " seconds" << std::endl;
}
