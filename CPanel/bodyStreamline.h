//
//  bodyStreamline.h
//  CPanel - Unstructured Panel Code
//
//  Created by Chris Satterwhite on 1/27/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel___Unstructured_Panel_Code__bodyStreamline__
#define __CPanel___Unstructured_Panel_Code__bodyStreamline__

#include <stdio.h>
#include <vector>
#include <Eigen/Dense>
#include "bodyPanel.h"
#include "wakePanel.h"
#include "edge.h"
#include "cpNode.h"

class geometry;

class bodyStreamline
{
    std::vector<Eigen::Vector3d> pnts;
    std::vector<Eigen::Vector3d> velocities;
    Eigen::Vector3d Vinf;
    double Vmag;
    geometry* geom;
    int pntsPerPanel;
    
//    Eigen::Vector3d trailingEdgePnt(bodyPanel* p);
    bool edgeIntersection(edge* e,const Eigen::Vector3d &pnt, const Eigen::Vector3d &vel, double &dt, Eigen::Vector3d &pntOnEdge);
    
    bool stagnationPnt(const Eigen::Vector3d vel, const Eigen::Vector3d &velOld, double maxAngle);
    
    double safeInvCos(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2);
public:
    bodyStreamline(Eigen::Vector3d startPnt, bodyPanel* startPan, const Eigen::Vector3d &Vinf, geometry* geom, int pntsPerPanel, bool marchFwd);
    
    std::vector<Eigen::Vector3d> getPnts() {return pnts;}
    std::vector<Eigen::Vector3d> getVelocities() {return velocities;}
    
};

#endif /* defined(__CPanel___Unstructured_Panel_Code__bodyStreamline__) */
