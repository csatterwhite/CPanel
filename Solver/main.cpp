//
//  main.cpp
//  Solver
//
//  Created by Chris Satterwhite on 4/30/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include <iostream>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include "geometry.h"
#include "runCase.h"

int main(int argc, const char * argv[])
{
    
    // Start Timer
    
    time_t ts;
    time(&ts);
    
    
    // Set Input File Path and Name
    
    std::string path;
    std::string geomName;
    std::string inFile = "/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files (Horseshoe)/coarseNormCompare.tri";
    
    double Vinf = 100;
    double alpha = 5;
    double beta = 0;
    
    // Create Geometry Name for Subdirectory Creation
    
    std::size_t pathEnd = inFile.find_last_of("/")+1;
    std::size_t nameEnd = inFile.find_last_of(".");
    if (pathEnd != 0)
    {
        path = inFile.substr(0,pathEnd);
        geomName = inFile.substr(pathEnd,nameEnd-pathEnd);
    }
    else
    {
        boost::filesystem::path p = boost::filesystem::current_path();
        path = p.string();
        geomName = inFile.substr(0,nameEnd);
    }
    
    
    // Check file type
    bool normFlag; // Are panel normals included in file?
    std::string ext = inFile.substr(nameEnd,inFile.size()-nameEnd);
    if (ext == ".tri")
    {
        normFlag = false;
    }
    else if (ext == ".tricp")
    {
        normFlag = true;
    }
    else
    {
        std::cout << "ERROR : Unsupported File Type \nAccepted Filetypes : '.tri','.tricp'" << std::endl;
        exit(EXIT_FAILURE);
    }

    
    time_t tf;
    geometry temp(inFile,normFlag);
    runCase case1(&temp,Vinf,alpha,beta,path,geomName);
    
    time(&tf);
    std::cout << "Elapsed time for program execution : " << difftime(tf,ts) << " seconds" << std::endl;
}
