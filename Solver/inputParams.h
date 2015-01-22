//
//  inputParams.h
//  CPanel
//
//  Created by Chris Satterwhite on 1/16/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__inputParams__
#define __CPanel__inputParams__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <iomanip>
#include <limits>

struct inputParams
{
    std::string inputFile;
    std::string inputPath;
    std::string inputName;
    std::string geomFile;
    std::string geomPath;
    std::string geomName;
    bool normFlag; //Normals included in input file
    double Sref;
    double bref;
    double cref;
    Eigen::Vector3d cg;
    Eigen::VectorXd velocities,alphas,betas,machs;
    
    inputParams(const char* inFile)
    {
        inputFile = inFile;
    }
    
    void set(std::ifstream &fid);
    
    void checkGeomFile();
    
    void print(std::ostream &stream);
    
    void printVec(Eigen::VectorXd &vec,std::ostream &stream);
};


#endif /* defined(__CPanel__inputParams__) */
