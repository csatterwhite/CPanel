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
#include "cpFile.h"

struct inputParams
{
    cpFile* inputFile;
    cpFile* geomFile;
    bool normFlag; //Normals included in input file
    double Sref;
    double bref;
    double cref;
    double streamSpacing;
    Eigen::Vector3d cg;
    Eigen::VectorXd velocities,alphas,betas,machs;
    
    inputParams(cpFile* inFile) : inputFile(inFile) {}
    
    ~inputParams()
    {
        delete geomFile;
    }
    
    bool set();
    void print(std::ostream &stream);

private:
    bool checkGeomFile();
    void makeWorkingDir();
    void writeInputFile();
    void printVec(Eigen::VectorXd &vec,std::ostream &stream);
};


#endif /* defined(__CPanel__inputParams__) */
