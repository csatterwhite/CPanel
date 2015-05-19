//
//  streamline.h
//  CPanel - Unstructured Panel Code
//
//  Created by Chris Satterwhite on 1/26/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel___Unstructured_Panel_Code__streamline__
#define __CPanel___Unstructured_Panel_Code__streamline__

#include <stdio.h>
#include <vector>
#include <Eigen/Dense>
#include "bodyPanel.h"
#include "wakePanel.h"
#include "geometry.h"

class streamline
{
    std::vector<Eigen::Vector3d> pnts;
    std::vector<Eigen::Vector3d> velocities;
    Eigen::Vector3d Vinf;
    double PG;
    geometry* geom;
    
    Eigen::Matrix<double,1,6> coeff5,coeff4;
    
    
    Eigen::Vector3d rkf(const Eigen::Vector3d &x0,double dt,double &error);
    
public:
    streamline(const Eigen::Vector3d &startPnt, double xMax, double tol, const Eigen::Vector3d &Vinf, double PG, geometry* geom);
    
    std::vector<Eigen::Vector3d> getPnts() {return pnts;}
    std::vector<Eigen::Vector3d> getVelocities() {return velocities;}
};

#endif /* defined(__CPanel___Unstructured_Panel_Code__streamline__) */
