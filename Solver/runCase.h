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
#include <boost/filesystem/operations.hpp>
#include "geometry.h"
#include "VTUfile.h"

class runCase
{
    geometry *geom;
    std::string path;
    std::string geomName;
    double Vmag;
    double alpha;
    double beta;
    Eigen::Vector3d Vinf;
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
    void writeFiles();
    void writeBodyData(boost::filesystem::path path);
    void writeWakeData(boost::filesystem::path path);
    
public:
    runCase(geometry *geom,double V,double alpha,double beta,std::string path,std::string geomName) : geom(geom), Vmag(V), alpha(alpha), beta(beta), path(path), geomName(geomName)
    {
        Vinf = windToBody(V,alpha,beta);
        setSourceStrengths();
        solveMatrixEq();
    }
    
};
#endif /* defined(__CPanel__runCase__) */
