//
//  runCase.h
//  CPanel
//
//  Created by Chris Satterwhite on 10/13/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__runCase__
#define __CPanel__runCase__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include "geometry.h"

class runCase
{
    geometry *geom;
    std::string outFile;
    Eigen::Vector3d Vinf;
    Eigen::MatrixXd velocities;
    Eigen::VectorXd sigmas;
    std::vector<bodyPanel*> bPanels;
    std::vector<wakePanel*> wPanels;
    std::vector<panel*> allPanels;
    Eigen::MatrixXd Ab; // Potential Influence Coefficient Matrix for doublets on just body panels
    Eigen::MatrixXd Bb; // Potential Influence Coefficient Matrix for sources on just body panels
    
    Eigen::Vector3d windToBody(double V,double alpha,double beta);
    void setSourceStrengths();
    void solveMatrixEq();
    Eigen::Vector4i getIndices(std::vector<bodyPanel*> interpPans);
    void writeVTU(std::string filename);
    
public:
    runCase(geometry *geom,double V,double alpha,double beta,std::string outFile);
    
    
};
#endif /* defined(__CPanel__runCase__) */
