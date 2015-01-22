//
//  main.cpp
//  Solver
//
//  Created by Chris Satterwhite on 4/30/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include "CPanelMgr.h"
#include "inputParams.h"
#include "geometry.h"

void usage()
{
    printf("\n");
    printf("CPanel is an unstructured panel code developed by students at California Polytechnic State University - San Luis Obispo.\n");
    printf("\n");
    printf("Usage: CPanel -i<infile> \n");
    printf("   -i<infile>   Required parameter to specify settings input file.\n");
    printf("\n");
    printf("The input file contains the geometry filename, reference values, and freestream conditions.\n\n");
    printf("Currently supported geometry file formats are:\n");
    printf("\t.tri (Cart3D format)\n");
    printf("\t.tricp (Modified .tri file including normal vectors from underlying bezier surfaces)\n\n");
    printf("This example below may be used as a template for the input file.\n");
    printf("\n");
    printf("--CPanel Input File--\n");
    printf("GeomFile = wing.tri\n");
    printf("S_ref =  6.0\n");
    printf("b_ref =  6.0\n");
    printf("c_ref =  1.0\n");
    printf("X_cg =  0.25\n");
    printf("Y_cg =  0.0\n");
    printf("Z_cg =  0.0\n");
    printf("Velocity (ft/s)\n");
    printf("3\n");
    printf("120\n");
    printf("160\n");
    printf("200\n");
    printf("Angle of Attack (degrees)\n");
    printf("6\n");
    printf("0\n");
    printf("2\n");
    printf("4\n");
    printf("6\n");
    printf("8\n");
    printf("10\n");
    printf("Angle of Sideslip (degrees)\n");
    printf("1\n");
    printf("0\n");
    printf("Mach Number\n");
    printf("1\n");
    printf("1.4\n");
    printf("\n");
    printf("The number following Velocity, Angle of Attack, Angle of Sideslip, and Mach Number indicates the number of different values for those variables, followed by the actual values\n");
    printf("\n");
    printf("For questions, contact :\n");
    printf("\tDavid D. Marshall\n");
    printf("\tCal Poly, San Luis Obispo\n");
    printf("\tddmarsha@calpoly.edu\n\n");
    printf("v. 1.0 - 1/24/2014\n");
    printf("\n");
    exit(1);
}

int main(int argc, const char * argv[])
{
    // Start Timer
    time_t ts,tf;
    time(&ts);
    
    // Check arguments
    if (argc != 2)
    {
        printf("ERROR : Incorrect Usage");
        usage();
        exit(EXIT_FAILURE);
    }
    
    // Check for file existence
    const char* inFile = argv[1];
    std::ifstream fid;
    fid.open(inFile);
    
    if (!fid)
    {
        std::cout << "ERROR : File could not be opened" << std::endl;
        usage();
        exit(EXIT_FAILURE);
    }
    
    // Check that file is right type
    std::string line1;
    std::getline(fid,line1);
    
    if (line1 != "%%% CPanel Input File %%%")
    {
        std::cout << "ERROR : Incorrect Format - Use Template Below" << std::endl;
        usage();
        exit(EXIT_FAILURE);
    }
    
    inputParams params(inFile);
    params.set(fid);
    std::cout << "Running CPanel with the following inputs..." << std::endl;
    params.print(std::cout);
    geometry geom(params.geomFile,params.normFlag,params.Sref,params.bref,params.cref,params.cg);
    caseMgr cm(params,&geom);
    
    
    time(&tf);
    std::cout << "Elapsed time for program execution : " << difftime(tf,ts) << " seconds" << std::endl;
}
