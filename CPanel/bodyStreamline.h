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

class bodyStreamline
{
    std::vector<Eigen::Vector3d> pnts;
    std::vector<Eigen::Vector3d> velocities;
    std::vector<double> potentials;
    std::vector<bodyPanel*>* bPanels;
    std::vector<wakePanel*>* wPanels;
    int pntsPerPanel;
    
    Eigen::Vector3d trailingEdgePnt(bodyPanel* p);
    double pntPotential(const Eigen::Vector3d &pnt, const Eigen::Vector3d Vinf, bodyPanel* currentPan);
    bool edgeIntersection(edge* e,const Eigen::Vector3d &pnt, const Eigen::Vector3d &vel, double &dt, Eigen::Vector3d &pntOnEdge);
    
    bool stagnationPnt(const Eigen::Vector3d vel, const Eigen::Vector3d &velOld);
public:
    bodyStreamline(bodyPanel* startPan,std::vector<bodyPanel*>* bPanels,std::vector<wakePanel*>* wPanels, const Eigen::Vector3d &Vinf, int pntsPerPanel);
    
    std::vector<Eigen::Vector3d> getPnts() {return pnts;}
    
};

#endif /* defined(__CPanel___Unstructured_Panel_Code__bodyStreamline__) */
